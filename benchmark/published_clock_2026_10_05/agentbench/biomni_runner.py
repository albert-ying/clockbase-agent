"""Run one Biomni A1 session on a prepared run directory (PROTOCOL.md of groupbench, item 24).

  BIOMNI_PYTHON -I -B biomni_runner.py MODEL        (working directory: the run directory)

run_biomni.sh starts it inside bubblewrap with OPENAI_API_KEY from the credential launcher. The task is the run's
TASK.md. RUN/events.jsonl receives one JSON object per line, in order:
  {"type": "config", ...}                                   Biomni version, model and options
  {"type": "system_prompt", "content": ...}                 A1's prompt after its resource retrieval
  {"type": "message", "role": "human" | "ai", ...}          every message of A1's graph, as it is produced
  {"type": "usage", ...}                                    token usage of each model call
  {"type": "end", "messages": n} or {"type": "error", ...}
Adaptations of item 24 (b), made here and not in Biomni's files: biomni.llm's handling of gpt-5 models (Responses API,
no stop or temperature parameter) applied to gpt-6 with medium reasoning effort and a usage callback; no data-lake
download (A1's expected_data_lake_files=[]); the loop of A1.go driven here so that each message is written as it is
produced; the API key removed from the environment once the client is built, since executed code shares the process.
Item 25: the Responses API does not apply A1's stop sequences, so the client cuts each reply at the first of them, as
an API that accepts stop sequences does.
"""
import json
import os
import sys
import time
import traceback
from pathlib import Path

RUN = Path.cwd()
EVENTS = RUN / 'events.jsonl'
STOPS = ('</execute>', '</solution>')  # A1's stop sequences (biomni/agent/a1.py)


def truncate_at_stop(text, stops=STOPS):
    """The text before the first stop sequence, as an API that applies stop sequences returns it."""
    cuts = [text.find(s) for s in stops if s in text]
    return text[:min(cuts)] if cuts else text


def emit(event):
    with EVENTS.open('a') as sink:
        sink.write(json.dumps({**event, 'time': time.time()}, default=str) + '\n')


def model_client(model):
    """biomni.llm's gpt-5 client (Responses API without stop or temperature), with reasoning effort and usage."""
    from langchain_core.callbacks import BaseCallbackHandler
    from langchain_openai import ChatOpenAI

    class Usage(BaseCallbackHandler):
        def on_llm_end(self, response, **kwargs):
            for generations in response.generations:
                for generation in generations:
                    meta = getattr(getattr(generation, 'message', None), 'usage_metadata', None) or {}
                    if meta:
                        details = meta.get('input_token_details') or {}
                        emit({'type': 'usage', 'input_tokens': meta.get('input_tokens', 0),
                              'cached_input_tokens': details.get('cache_read', 0) or 0,
                              'output_tokens': meta.get('output_tokens', 0)})

    class ResponsesWithoutStop(ChatOpenAI):
        def _get_request_payload(self, input_, *, stop=None, **kwargs):
            payload = super()._get_request_payload(input_, stop=stop, **kwargs)
            payload.pop('stop', None)
            payload.pop('temperature', None)
            return payload

        def invoke(self, input, config=None, **kwargs):  # noqa: A002 (langchain's signature)
            message = super().invoke(input, config=config, **kwargs)
            content = message.content
            if isinstance(content, list):  # Responses API blocks, joined as A1's generate joins them
                content = ''.join(b.get('text') or '' for b in content
                                  if isinstance(b, dict) and b.get('type') in ('text', 'output_text', 'redacted_text'))
            return message.model_copy(update={'content': truncate_at_stop(str(content))})

    return ResponsesWithoutStop(model=model, temperature=1, stop_sequences=['</execute>', '</solution>'],
                                use_responses_api=True, output_version='v0', reasoning={'effort': 'medium'},
                                callbacks=[Usage()])


def main():
    model = sys.argv[1]
    task = (RUN / 'TASK.md').read_text()
    import biomni.agent.a1 as a1_module
    from biomni.version import __version__
    from langchain_core.messages import HumanMessage

    original = a1_module.get_llm

    def get_llm(name=None, *args, **kwargs):
        return model_client(name) if name and name.startswith('gpt-6') else original(name, *args, **kwargs)

    a1_module.get_llm = get_llm
    emit({'type': 'config', 'biomni_version': __version__, 'model': model, 'reasoning_effort': 'medium',
          'data_lake': 'not downloaded (expected_data_lake_files=[])', 'run': str(RUN)})
    try:
        agent = a1_module.A1(path=str(RUN / 'biomni_home'), llm=model, source='OpenAI', expected_data_lake_files=[])
        os.environ.pop('OPENAI_API_KEY', None)
        # A1.go: retrieval, then the compiled graph with the same inputs and configuration.
        agent.critic_count = 0
        agent.user_task = task
        if agent.use_tool_retriever:
            agent.update_system_prompt_with_selected_resources(agent._prepare_resources_for_retrieval(task))
        emit({'type': 'system_prompt', 'content': agent.system_prompt})
        inputs = {'messages': [HumanMessage(content=task)], 'next_step': None}
        config = {'recursion_limit': 500, 'configurable': {'thread_id': 42}}
        seen = 0
        for state in agent.app.stream(inputs, stream_mode='values', config=config):
            messages = state['messages']
            for message in messages[seen:]:
                emit({'type': 'message', 'role': 'human' if isinstance(message, HumanMessage) else 'ai',
                      'content': message.content})
            seen = len(messages)
        emit({'type': 'end', 'messages': seen})
    except Exception as error:  # recorded as the run's outcome
        emit({'type': 'error', 'error': repr(error), 'traceback': traceback.format_exc()})
        sys.exit(1)


if __name__ == '__main__':
    main()
