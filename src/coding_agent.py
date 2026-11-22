"""
Self-contained coding agent distilled from BioAgentOS codeact_agent.
Provides a minimal LangGraph-based agent with Python/R/Bash code execution.
"""

import os, re, sys, subprocess, tempfile, threading, queue, ctypes
from uuid import uuid4
from io import StringIO
from typing import List, Optional, Literal, TypedDict

from langchain_core.messages import BaseMessage, HumanMessage, AIMessage, SystemMessage
from langgraph.graph import StateGraph, START, END
from langgraph.checkpoint.memory import MemorySaver


# ---------------------------------------------------------------------------
# LLM
# ---------------------------------------------------------------------------

def get_llm(model='claude-3-7-sonnet-20250219', temperature=0.7, stop_sequences=None):
    """Create a ChatAnthropic LLM instance."""
    # Keep the historical default; injected OpenAI clients need no Anthropic package.
    from langchain_anthropic import ChatAnthropic
    return ChatAnthropic(
        model=model,
        temperature=temperature,
        max_tokens=8192,
        stop_sequences=stop_sequences,
    )


# ---------------------------------------------------------------------------
# Code execution tools
# ---------------------------------------------------------------------------

_persistent_namespace = {}


def run_python_repl(command: str) -> str:
    """Execute Python in a persistent namespace (variables survive across calls)."""
    old_stdout = sys.stdout
    sys.stdout = buf = StringIO()
    global _persistent_namespace
    try:
        exec(command.strip("```").strip(), _persistent_namespace)
        output = buf.getvalue()
    except Exception as e:
        output = f"Error: {e}"
    finally:
        sys.stdout = old_stdout
    return output


def run_r_code(code: str) -> str:
    """Run R code via Rscript subprocess."""
    try:
        with tempfile.NamedTemporaryFile(suffix='.R', mode='w', delete=False) as f:
            f.write(code)
            tmp = f.name
        result = subprocess.run(['Rscript', tmp], capture_output=True, text=True, check=False)
        os.unlink(tmp)
        return result.stderr if result.returncode != 0 else result.stdout
    except Exception as e:
        return f"Error running R code: {e}"


def run_bash_script(script: str) -> str:
    """Run a Bash script via subprocess."""
    try:
        script = script.strip()
        if not script:
            return "Error: Empty script"
        with tempfile.NamedTemporaryFile(suffix='.sh', mode='w', delete=False) as f:
            if not script.startswith("#!/"):
                f.write("#!/bin/bash\n")
            if "set -e" not in script:
                f.write("set -e\n")
            f.write(script)
            tmp = f.name
        os.chmod(tmp, 0o755)
        result = subprocess.run(
            [tmp], shell=True, capture_output=True, text=True,
            check=False, env=os.environ.copy(), cwd=os.getcwd(),
        )
        os.unlink(tmp)
        if result.returncode != 0:
            return f"Error running Bash script (exit code {result.returncode}):\n{result.stderr}"
        return result.stdout
    except Exception as e:
        return f"Error running Bash script: {e}"


# ---------------------------------------------------------------------------
# Timeout wrapper
# ---------------------------------------------------------------------------

def run_with_timeout(func, args=None, kwargs=None, timeout=600):
    """Run *func* with a timeout (seconds) using a daemon thread."""
    args = args or []
    kwargs = kwargs or {}
    result_queue = queue.Queue()

    def _worker():
        try:
            result_queue.put(("ok", func(*args, **kwargs)))
        except Exception as e:
            result_queue.put(("err", str(e)))

    t = threading.Thread(target=_worker, daemon=True)
    t.start()
    t.join(timeout)

    if t.is_alive():
        print(f"TIMEOUT: Code execution timed out after {timeout} seconds")
        try:
            tid = t.ident
            if tid:
                res = ctypes.pythonapi.PyThreadState_SetAsyncExc(
                    ctypes.c_long(tid), ctypes.py_object(SystemExit))
                if res > 1:
                    ctypes.pythonapi.PyThreadState_SetAsyncExc(ctypes.c_long(tid), None)
        except Exception:
            pass
        return f"ERROR: Code execution timed out after {timeout} seconds. Please try with simpler inputs or break your task into smaller steps."

    try:
        status, result = result_queue.get(block=False)
        return result if status == "ok" else f"Error in execution: {result}"
    except queue.Empty:
        return "Error: Execution completed but no result was returned"


# ---------------------------------------------------------------------------
# Pretty-print helper
# ---------------------------------------------------------------------------

def pretty_print(message, printout=True):
    """Format and optionally print a LangChain message."""
    if isinstance(message, tuple):
        return message
    header = f"=== {message.type.title()} Message ==="
    name_line = f"\nName: {message.name}" if getattr(message, 'name', None) else ""
    if isinstance(message.content, list):
        parts = []
        for item in message.content:
            if item.get('type') == 'text':
                parts.append(item['text'])
            elif item.get('type') == 'tool_use':
                parts.append(f"Tool: {item['name']}\nInput: {item['input']}")
        body = "\n".join(parts)
    else:
        body = str(message.content)
    title = f"{header}{name_line}\n\n{body}"
    if printout:
        print(title)
    return title


# ---------------------------------------------------------------------------
# Default system prompt
# ---------------------------------------------------------------------------

DEFAULT_SYSTEM_PROMPT = """\
You are a helpful assistant assigned with the task of problem-solving.
To achieve this, you will be using an interactive coding environment equipped with Python, R, and Bash.

Given a task, make a plan first. The plan should be a numbered list of steps that you will take to solve the task. Be specific and detailed.
Format your plan as a checklist with empty checkboxes like this:
1. [ ] First step
2. [ ] Second step
3. [ ] Third step

Follow the plan step by step. After completing each step, update the checklist by replacing the empty checkbox with a checkmark:
1. [✓] First step (completed)
2. [ ] Second step
3. [ ] Third step

If a step fails or needs modification, mark it with an X and explain why:
1. [✓] First step (completed)
2. [✗] Second step (failed because...)
3. [ ] Modified second step
4. [ ] Third step

Always show the updated plan after each step so the user can track progress.

At each turn, you should first provide your thinking and reasoning given the conversation history.
After that, you have two options:

1) Interact with a programming environment and receive the corresponding output within <observe></observe>. Your code should be enclosed using "<execute>" tag, for example: <execute> print("Hello World!") </execute>. IMPORTANT: You must end the code block with </execute> tag.
   - For Python code (default): <execute> print("Hello World!") </execute>
   - For R code: <execute> #!R
library(ggplot2)
print("Hello from R") </execute>
   - For Bash scripts and commands: <execute> #!BASH
echo "Hello from Bash"
ls -la </execute>
   - For CLI softwares, use Bash scripts.

2) When you think it is ready, directly provide a solution that adheres to the required format for the given task to the user. Your solution should be enclosed using "<solution>" tag, for example: The answer is <solution> A </solution>. IMPORTANT: You must end the solution block with </solution> tag.

You have many chances to interact with the environment to receive the observation. So you can decompose your code into multiple steps.
Don't overcomplicate the code. Keep it simple and easy to understand.
When writing the code, please print out the steps and results in a clear and concise manner, like a research log.
YOU MUST SAVE THE OUTPUT and PRINT OUT the result of each step.
For example, result = some_function(XXX); print(result)
Otherwise the system will not be able to know what has been done.

For R code, use the #!R marker at the beginning of your code block to indicate it's R code.
For Bash scripts and commands, use the #!BASH marker at the beginning of your code block. This allows for both simple commands and multi-line scripts with variables, conditionals, loops, and other Bash features.

In each response, you must include EITHER <execute> or <solution> tag. Not both at the same time. Do not respond with messages without any tags. No empty messages.
"""


# ---------------------------------------------------------------------------
# LangGraph state
# ---------------------------------------------------------------------------

class AgentState(TypedDict):
    messages: List[BaseMessage]
    next_step: Optional[str]


# ---------------------------------------------------------------------------
# CodingAgent
# ---------------------------------------------------------------------------

class CodingAgent:
    """Minimal CodeAct-style agent: generate → execute → observe → loop."""

    def __init__(self, llm='claude-3-7-sonnet-20250219', system_prompt=None, timeout=600,
                 executor=None, strict_tags=False, recursion_limit=500, printout=True):
        self.llm = get_llm(llm, stop_sequences=['</execute>', '</solution>']) if isinstance(llm, str) else llm
        self.system_prompt = system_prompt or DEFAULT_SYSTEM_PROMPT
        self.timeout_seconds = timeout
        self.executor = executor
        self.strict_tags = strict_tags
        self.recursion_limit = recursion_limit
        self.printout = printout
        self.app = None
        self._run_lock = threading.Lock()

    def configure(self):
        """Build the LangGraph workflow."""
        agent = self  # capture for closures

        def generate(state: AgentState) -> AgentState:
            messages = [SystemMessage(content=agent.system_prompt)] + state['messages']
            response = agent.llm.invoke(messages)
            msg = str(response.content)

            # Fix incomplete tags
            if not agent.strict_tags and '<execute>' in msg and '</execute>' not in msg:
                msg += '</execute>'
            if not agent.strict_tags and '<solution>' in msg and '</solution>' not in msg:
                msg += '</solution>'
            if not agent.strict_tags and '<think>' in msg and '</think>' not in msg:
                msg += '</think>'

            think_match = re.search(r'<think>(.*?)</think>', msg, re.DOTALL)
            execute_match = re.search(r'<execute>(.*?)</execute>', msg, re.DOTALL)
            answer_match = re.search(r'<solution>(.*?)</solution>', msg, re.DOTALL)
            if agent.strict_tags:
                if bool(execute_match) == bool(answer_match):
                    raise ValueError('CodingAgent requires exactly one complete execute or solution block')
                if len(re.findall(r'<(?:execute|solution)>', msg)) != 1:
                    raise ValueError('CodingAgent returned multiple action blocks')

            state['messages'].append(AIMessage(content=msg.strip()))

            if answer_match:
                state['next_step'] = 'end'
            elif execute_match:
                state['next_step'] = 'execute'
            elif think_match:
                state['next_step'] = 'generate'
            else:
                error_count = sum(
                    1 for m in state['messages']
                    if isinstance(m, AIMessage) and "There are no tags" in m.content
                )
                if error_count >= 2:
                    print('Detected repeated parsing errors, ending conversation')
                    state['next_step'] = 'end'
                    state['messages'].append(AIMessage(
                        content="Execution terminated due to repeated parsing errors. Please check your input and try again."
                    ))
                else:
                    state['messages'].append(AIMessage(
                        content="There are no tags (e.g. <execute><solution>). Please follow the instruction, fix and update."
                    ))
                    state['next_step'] = 'generate'
            return state

        def execute(state: AgentState) -> AgentState:
            last_message = state['messages'][-1].content
            if '<execute>' in last_message and '</execute>' not in last_message:
                last_message += '</execute>'

            match = re.search(r'<execute>(.*?)</execute>', last_message, re.DOTALL)
            if match:
                code = match.group(1)
                timeout = agent.timeout_seconds

                if agent.executor is not None:
                    result = agent.executor(code)
                elif code.strip().startswith(("#!R", "# R code", "# R script")):
                    r_code = re.sub(r'^#!R|^# R code|^# R script', '', code, 1).strip()
                    result = run_with_timeout(run_r_code, [r_code], timeout=timeout)
                elif code.strip().startswith(("#!BASH", "# Bash script", "#!CLI")):
                    if code.strip().startswith("#!CLI"):
                        cli_cmd = re.sub(r'^#!CLI', '', code, 1).strip().replace('\n', ' ')
                        result = run_with_timeout(run_bash_script, [cli_cmd], timeout=timeout)
                    else:
                        bash_script = re.sub(r'^#!BASH|^# Bash script', '', code, 1).strip()
                        result = run_with_timeout(run_bash_script, [bash_script], timeout=timeout)
                else:
                    result = run_with_timeout(run_python_repl, [code], timeout=timeout)

                if len(result) > 10000:
                    result = 'The output is too long to be added to context. Here are the first 10K characters...\n' + result[:10000]

                observation = f"\n<observation>{result}</observation>"
                state['messages'].append(AIMessage(content=observation.strip()))
            return state

        def routing_function(state: AgentState) -> Literal['execute', 'generate', 'end']:
            ns = state.get('next_step')
            if ns in ('execute', 'generate', 'end'):
                return ns
            raise ValueError(f"Unexpected next_step: {ns}")

        workflow = StateGraph(AgentState)
        workflow.add_node('generate', generate)
        workflow.add_node('execute', execute)
        workflow.add_conditional_edges(
            'generate', routing_function,
            path_map={'execute': 'execute', 'generate': 'generate', 'end': END},
        )
        workflow.add_edge('execute', 'generate')
        workflow.add_edge(START, 'generate')

        self.app = workflow.compile()
        checkpointer = MemorySaver()
        self.app.checkpointer = checkpointer

    def go(self, prompt):
        """Return a per-call trace, serializing reuse of this mutable agent.

        Checkpoint and log isolation do not certify source identity or isolate
        the legacy default code executor's process-wide Python namespace.
        """
        with self._run_lock:
            return self._go_locked(prompt)

    def _go_locked(self, prompt):
        if self.app is None:
            raise RuntimeError("Call .configure() before .go()")
        inputs = {'messages': [HumanMessage(content=prompt)], 'next_step': None}
        config = {"recursion_limit": self.recursion_limit, "configurable": {"thread_id": str(uuid4())}}
        log = []
        message = None
        for s in self.app.stream(inputs, stream_mode="values", config=config):
            message = s["messages"][-1]
            out = pretty_print(message, printout=self.printout)
            log.append(out)
        final = message.content if message else ""
        # Retain the last-run compatibility attribute without sharing the
        # returned list with the next invocation or its caller.
        self.log = list(log)
        return log, final
