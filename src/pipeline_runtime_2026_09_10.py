"""Portable, bounded model transport and isolated execution for CodingAgent.

No local credential service is required. Scientific acceptance is handled by the
audit, separately from successful software execution. Linux/bubblewrap is the
supported execution environment; missing isolation never falls back to host exec.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import signal
import subprocess
import sys
import time

from audit_system.agent_workflow_2026_09_10 import sha, write_new


class PipelineIntegrityError(ValueError):
    """A frozen input, implementation or receipt was changed."""


_OMITTED_SERVICE_TIER = object()


def service_tier_options(settings):
    """Keep omission unchanged; an explicit transport tier must be default."""
    if 'service_tier' not in settings:
        return {}
    value = settings['service_tier']
    if type(value) is not str or value != 'default':
        raise ValueError('Explicit service_tier must be default; omit it for legacy transport')
    return {'service_tier': value}


def openai_request(model, messages, output_tokens, *, json_response=False,
                   service_tier=_OMITTED_SERVICE_TIER):
    """Build the same wire object for transport, reservations and replay hashes."""
    request = {'model': model, 'messages': messages, 'max_completion_tokens': output_tokens,
               'reasoning_effort': 'low', 'store': False}
    if json_response:
        request['response_format'] = {'type': 'json_object'}
    if service_tier is not _OMITTED_SERVICE_TIER:
        request.update(service_tier_options({'service_tier': service_tier}))
    return request


def verify_service_tier(request, response):
    """Reject a bound tier mismatch before content use or accounting settlement."""
    try:
        required = service_tier_options(request)
    except ValueError as error:
        raise PipelineIntegrityError('Bound request has an invalid service tier') from error
    if required and (not isinstance(response, dict) or
                     response.get('service_tier') != required['service_tier']):
        raise PipelineIntegrityError('Returned service tier differs from the explicit request; no fallback')


def identity(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _redact_provider_error_field(value, *, key, request_text):
    """Return a bounded diagnostic string without credentials or request text."""
    if value is None:
        return None
    if not isinstance(value, (str, int, float, bool)):
        return None
    text = str(value)[:4096]
    if key:
        text = text.replace(key, '[REDACTED_CREDENTIAL]')
    text = re.sub(r'(?i)authorization\s*:\s*bearer\s+[^\s,;]+',
                  'Authorization: Bearer [REDACTED_CREDENTIAL]', text)
    text = re.sub(r'(?i)\b(?:sk|sess|gh[opsu]|github_pat)-[A-Za-z0-9_-]{8,}\b',
                  '[REDACTED_CREDENTIAL]', text)
    text = re.sub(r'(?i)\b(?:api[_ -]?key|access[_ -]?token|secret)\s*[:=]\s*[^\s,;]+',
                  '[REDACTED_CREDENTIAL]', text)
    text = re.sub(r"(?<![A-Za-z0-9_])/(?:Entropy|home|Users|tmp|var/tmp)/[^\s\"'<>]*",
                  '[REDACTED_LOCAL_PATH]', text)
    # Provider errors occasionally quote an offending request value. Keep the
    # diagnostic useful, but never persist content copied from the prompt.
    for delimiter, quoted in re.findall(r'([\"\'`])(.{8,}?)\1', text):
        if quoted in request_text:
            text = text.replace(delimiter + quoted + delimiter, '[REDACTED_REQUEST_TEXT]')
    if len(text) >= 24 and any(text[index:index + 24] in request_text
                               for index in range(len(text) - 23)):
        text = '[REDACTED_REQUEST_TEXT]'
    text = ''.join(character if character >= ' ' else ' ' for character in text)
    return text[:2048]


def _provider_error_record(body, http_code, *, key, request):
    """Extract only the documented OpenAI error fields from a failed reply."""
    request_text = '\n'.join(str(message.get('content', '')) for message in request['messages'])
    try:
        parsed = json.loads(body) if len(body) <= 1024 * 1024 else None
        candidate = parsed.get('error') if isinstance(parsed, dict) else None
        if not isinstance(candidate, dict):
            raise ValueError('no structured provider error')
        error = {field: _redact_provider_error_field(
                     candidate.get(field), key=key, request_text=request_text)
                 for field in ('type', 'code', 'param', 'message')}
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError):
        error = {'type': None, 'code': None, 'param': None,
                 'message': 'Provider returned an unstructured error body; content was not retained'}
    return {'schema': 'clockbase-openai-error-v1', 'provider': 'openai',
            'http_status': http_code, 'billing_status': 'uncertain',
            'automatic_retry': False, 'raw_response_retained': False, 'error': error}


class StageStore:
    def __init__(self, root, config):
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        path = self.root / 'config.json'
        if path.exists():
            if json.loads(path.read_text()) != config:
                raise PipelineIntegrityError('Configuration, source or implementation changed; use a new run directory')
        else:
            write_new(path, config)

    def run(self, name, function):
        directory = self.root / name
        directory.mkdir(exist_ok=True)
        complete = directory / 'complete.json'
        if complete.exists():
            if not (directory / 'complete.sha256').is_file() or sha(complete) != (directory / 'complete.sha256').read_text().strip():
                raise PipelineIntegrityError('Completed stage receipt changed: ' + name)
            record = json.loads(complete.read_text())
            for path, checksum in record['outputs'].items():
                if sha(directory / path) != checksum:
                    raise PipelineIntegrityError('Completed stage output changed: ' + name + '/' + path)
            print('Resume verified: ' + name, flush=True)
            return record['result']
        attempt = directory / ('attempt-%03d' % len(list(directory.glob('attempt-*'))))
        attempt.mkdir()
        started = time.monotonic()
        print('Running: ' + name, flush=True)
        try:
            result = function(attempt)
            outputs = {str(p.relative_to(directory)): sha(p) for p in attempt.rglob('*') if p.is_file()}
            write_new(complete, {'result': result, 'outputs': outputs,
                                 'elapsed_seconds': time.monotonic() - started})
            with (directory / 'complete.sha256').open('x') as handle:
                handle.write(sha(complete) + '\n')
            return result
        except Exception as error:
            write_new(attempt / 'failure.json', {'status': 'CANNOT VERIFY',
                      'error': f'{type(error).__name__}: {error}',
                      'elapsed_seconds': time.monotonic() - started})
            raise


class BoundedOpenAI:
    """LangChain-compatible invoke with a shared persistent budget and no retry.

    A request is reserved before transport. Ambiguous/failed billing remains
    reserved. Successful replies can be replayed only for identical requests.
    The cost estimate is not a provider invoice. Rates are explicit run inputs.
    """
    def __init__(self, directory, actor, *, model='gpt-5.6-sol', max_calls=8,
                  max_cost=2.0, output_tokens=2500, prompt_bytes=48000,
                  input_per_million=5.0, output_per_million=20.0,
                  require_model_match=False, json_response=False,
                  service_tier=_OMITTED_SERVICE_TIER):
        if type(prompt_bytes) is not int or prompt_bytes <= 0:
            raise ValueError('prompt_bytes must be a positive integer')
        if type(json_response) is not bool:
            raise ValueError('json_response must be a boolean')
        self.tier_options = ({} if service_tier is _OMITTED_SERVICE_TIER else
                             service_tier_options({'service_tier': service_tier}))
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.actor, self.model, self.index = actor, model, 0
        self.max_calls, self.max_cost = max_calls, max_cost
        self.output_tokens, self.prompt_bytes = output_tokens, prompt_bytes
        self.input_rate, self.output_rate = input_per_million, output_per_million
        self.require_model_match = require_model_match
        self.json_response = json_response
        if min(max_calls, max_cost, output_tokens, prompt_bytes, input_per_million, output_per_million) <= 0:
            raise ValueError('All resource limits and tariff rates must be positive')

    def _verify_returned_tier(self, slot, request, response):
        diagnostic = slot / 'service-tier-error.json'
        if diagnostic.exists():
            raise PipelineIntegrityError('Prior service-tier rejection retained; use a new run directory')
        try:
            verify_service_tier(request, response)
        except PipelineIntegrityError:
            write_new(diagnostic, {
                'schema': 'clockbase-openai-service-tier-error-v1', 'provider': 'openai',
                'category': 'OPENAI_SERVICE_TIER_MISMATCH',
                'requested_service_tier': request.get('service_tier'),
                'returned_service_tier': response.get('service_tier') if isinstance(response, dict) else None,
                'request_sha256': sha(slot / 'request.json'),
                'response_sha256': sha(slot / 'response.json'),
                'billing_status': 'uncertain', 'automatic_retry': False,
                'raw_response_retained': True})
            raise

    def invoke(self, messages):
        from langchain_core.messages import AIMessage
        roles = {'human': 'user', 'ai': 'assistant', 'system': 'system'}
        wire = [{'role': roles[m.type], 'content': m.content} for m in messages]
        request = openai_request(self.model, wire, self.output_tokens,
                                 json_response=self.json_response, **self.tier_options)
        raw = json.dumps(request, ensure_ascii=False).encode()
        if len(raw) > self.prompt_bytes:
            raise ValueError(f'Prompt exceeds the frozen run byte limit: {len(raw)} bytes > '
                             f'permitted {self.prompt_bytes}; nothing truncated or sent')
        slot = self.directory / (self.actor + '-%03d' % self.index)
        self.index += 1
        if slot.exists():
            if json.loads((slot / 'request.json').read_text()) != request:
                raise PipelineIntegrityError('Request changed on resume; use a new run directory')
            if (slot / 'service-tier-error.json').exists():
                raise PipelineIntegrityError('Prior service-tier rejection retained; use a new run directory')
            if not (slot / 'receipt.json').exists():
                if self.tier_options and (slot / 'response.json').is_file():
                    self._verify_returned_tier(slot, request, json.loads((slot / 'response.json').read_text()))
                raise RuntimeError('Prior request has uncertain billing; no automatic repeat')
            receipt = json.loads((slot / 'receipt.json').read_text())
            if sha(slot / 'response.json') != receipt['response_sha256']:
                raise PipelineIntegrityError('Model response changed on resume')
            response = json.loads((slot / 'response.json').read_text())
        else:
            # A generic key may belong to a different provider. This transport
            # only targets OpenAI and must never forward another provider's key.
            key = os.environ.get('OPENAI_API_KEY')
            if not key:
                raise RuntimeError('Set OPENAI_API_KEY before starting a model stage')
            if key.encode() in raw:
                raise ValueError('Credential detected in prompt; request rejected')
            slots = list(self.directory.glob('*/reservation.json'))
            charged = 0.0
            for p in slots:
                receipt = p.parent / 'receipt.json'
                settled = receipt.exists() and not (p.parent / 'service-tier-error.json').exists()
                charged += json.loads(receipt.read_text())['estimated_cost_usd'] if settled else json.loads(p.read_text())['reserved_usd']
            # UTF-8 bytes plus message framing conservatively bound input tokens.
            reserve = ((len(raw) + 100 * len(wire)) * self.input_rate + self.output_tokens * self.output_rate) / 1e6
            if len(slots) >= self.max_calls or charged + reserve > self.max_cost:
                raise RuntimeError('Run API call/cost limit reached before another request')
            slot.mkdir()
            write_new(slot / 'request.json', request)
            write_new(slot / 'reservation.json', {'reserved_usd': reserve, 'started_unix': time.time(),
                      'model': self.model, 'automatic_retries': 0, 'concurrency': 1})
            started = time.monotonic()
            command = ['curl', '--silent', '--show-error', '--max-time', '180', '--connect-timeout', '30',
                       '--request', 'POST', '--header', '@-', '--data-binary', '@' + str(slot / 'request.json'),
                       '--write-out', '\n%{http_code}', 'https://api.openai.com/v1/chat/completions']
            timeout_seconds = 190
            try:
                transport = subprocess.run(command, input=('Authorization: Bearer ' + key + '\nContent-Type: application/json\n').encode(),
                                           capture_output=True, timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                # The provider may have received the request. Keep its reservation
                # and exclude the exception's command and captured output.
                write_new(slot / 'transport-error.json', {
                    'schema': 'clockbase-openai-transport-error-v1', 'provider': 'openai',
                    'category': 'OPENAI_TRANSPORT_TIMEOUT', 'error_type': 'TimeoutExpired',
                    'timeout_seconds': timeout_seconds,
                    'elapsed_seconds': max(0.0, time.monotonic() - started),
                    'request_sha256': sha(slot / 'request.json'),
                    'billing_status': 'uncertain', 'automatic_retry': False,
                    'raw_output_retained': False})
                raise RuntimeError('OPENAI_TRANSPORT_TIMEOUT: OpenAI transport exceeded its deadline; '
                                   'billing uncertain; no automatic retry') from None
            if transport.returncode:
                raise RuntimeError('OpenAI transport failed; billing uncertain; response suppressed')
            body, code = transport.stdout.rsplit(b'\n', 1)
            if code != b'200':
                code_text = code.decode(errors='replace')
                http_code = int(code_text) if re.fullmatch(r'[0-9]{3}', code_text) else None
                write_new(slot / 'provider-error.json',
                          _provider_error_record(body, http_code, key=key, request=request))
                raise RuntimeError('OpenAI HTTP ' + code_text + '; response suppressed')
            if key.encode() in body:
                raise RuntimeError('Response echoed a credential; discarded')
            response = json.loads(body)
            write_new(slot / 'response.json', response)
            self._verify_returned_tier(slot, request, response)
            usage = response.get('usage', {})
            if any(type(usage.get(k)) is not int or usage[k] < 0 for k in ('prompt_tokens', 'completion_tokens')):
                raise ValueError('Missing token usage; billing reservation retained')
            estimate = (usage['prompt_tokens'] * self.input_rate + usage['completion_tokens'] * self.output_rate) / 1e6
            write_new(slot / 'receipt.json', {'usage': usage, 'estimated_cost_usd': estimate,
                      'billed_cost_usd': None, 'response_sha256': sha(slot / 'response.json'),
                       'requested_model': self.model, 'returned_model': response.get('model'),
                       'elapsed_seconds': time.monotonic() - started})
        self._verify_returned_tier(slot, request, response)
        if self.require_model_match and response.get('model') != self.model:
            raise PipelineIntegrityError('Returned model differs from the requested model; no fallback')
        choices = response.get('choices', [])
        if len(choices) != 1 or choices[0].get('finish_reason') != 'stop':
            raise ValueError('Incomplete model response; no invented closing tags or code execution')
        return AIMessage(content=choices[0]['message']['content'])


class IsolatedExecutor:
    def __init__(self, workspace, *, timeout=120):
        self.workspace, self.timeout = Path(workspace).resolve(), timeout
        self.index = 0

    def __call__(self, code):
        import shutil
        if not shutil.which('bwrap'):
            raise RuntimeError('bubblewrap is required; host code execution is disabled')
        if code.lstrip().startswith(('#!R', '#!BASH', '#!CLI', '# R ', '# Bash')):
            return 'This bounded run supports Python execution only. Use pandas, numpy and scipy.'
        script = self.workspace / ('agent-step-%03d.py' % self.index)
        self.index += 1
        script.write_text(code, encoding='utf-8')
        command = ['bwrap', '--unshare-all', '--die-with-parent', '--new-session', '--clearenv']
        for folder in ('/usr', '/lib', '/lib64', '/bin'):
            if Path(folder).exists(): command += ['--ro-bind', folder, folder]
        command += ['--ro-bind', sys.prefix, '/runtime', '--bind', str(self.workspace), '/work',
                    '--tmpfs', '/tmp', '--proc', '/proc', '--dev', '/dev', '--chdir', '/work',
                    '--setenv', 'PATH', '/runtime/bin:/usr/bin:/bin', '--setenv', 'OPENBLAS_NUM_THREADS', '1',
                    '--setenv', 'OMP_NUM_THREADS', '1', '--setenv', 'MPLBACKEND', 'Agg']
        for name in ('metadata.json', 'clock_predictions.csv'):
            if (self.workspace / name).exists():
                command += ['--ro-bind', str(self.workspace / name), '/work/' + name]
        command += ['/runtime/bin/python', '/work/' + script.name]
        def limits():
            resource.setrlimit(resource.RLIMIT_CPU, (self.timeout, self.timeout))
            resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))
            resource.setrlimit(resource.RLIMIT_FSIZE, (32 * 1024**2, 32 * 1024**2))
        log_path = self.workspace / ('execution-%03d.log' % (self.index - 1))
        with log_path.open('xb') as log:
            process = subprocess.Popen(command, stdout=log, stderr=log,
                                       env={'PATH': '/usr/bin:/bin'}, start_new_session=True, preexec_fn=limits)
            try:
                process.wait(timeout=self.timeout + 5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                return 'ERROR: isolated code exceeded its execution time limit'
        with log_path.open('rb') as log:
            result = log.read(10000).decode(errors='replace')
        write_new(self.workspace / ('execution-%03d.json' % (self.index - 1)),
                  {'code_sha256': sha(script), 'returncode': process.returncode, 'output': result,
                   'network_access': False, 'host_home_access': False})
        return ('ERROR: ' if process.returncode else '') + result
