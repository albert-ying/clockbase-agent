import httpx
import time
import re
import json
import ast
from openai import AsyncOpenAI
from dataclasses import dataclass, field
from pydantic import BaseModel, Field
from typing import AsyncGenerator, Union, Optional

# Local imports
from configs.global_config import SYS_LLM_CONFIG
from configs.config_cls import LLMConfig


class ChatResponse(BaseModel):
    content: str
    prompt_tokens: int = Field(0)
    completion_tokens: int = Field(0)
    cost: float = Field(0.0)
    time_cost: float = Field(0.0)


class OpenAIWrapper:
    def __init__(self, config: LLMConfig = SYS_LLM_CONFIG):
        self.config = config
        self.client = AsyncOpenAI(
            base_url=config.base_url,
            api_key=config.api_key.get_secret_value(),
            http_client=httpx.AsyncClient(
                proxy=config.proxy,
                timeout=config.timeout
            ) if self.config.proxy else None
        )

    async def chat(self, messages: list, stop: Optional[list] = None) -> ChatResponse:
        st = time.time()
        response = await self.client.chat.completions.create(
            model=self.config.model,
            messages=messages,
            temperature=self.config.temperature,
            max_tokens=self.config.max_tokens,
            top_p=self.config.top_p,
            frequency_penalty=self.config.frequency_penalty,
            presence_penalty=self.config.presence_penalty,
            timeout=self.config.timeout,
            stop=stop
        )
        et = time.time()
        return ChatResponse(
            content=response.choices[0].message.content or '',
            prompt_tokens=response.usage.prompt_tokens or 0,
            completion_tokens=response.usage.completion_tokens or 0,
            cost=response.usage.prompt_tokens * self.config.prompt_cost + response.usage.completion_tokens * self.config.completion_cost,
            time_cost=et - st
        )

    async def stream_chat(self, messages: list, stop: Optional[list] = None):
        response = await self.client.chat.completions.create(
            model=self.config.model,
            messages=messages,
            temperature=self.config.temperature,
            max_tokens=self.config.max_tokens,
            top_p=self.config.top_p,
            frequency_penalty=self.config.frequency_penalty,
            presence_penalty=self.config.presence_penalty,
            timeout=self.config.timeout,
            stream=True,
            stop=stop
        )
        async for chunk in response:
            if chunk.choices[0].finish_reason == "stop":
                break
            yield chunk

    async def chat_once(self, prompt: str, sys_prompt: Optional[str] = None, stop: Optional[list] = None) -> ChatResponse:
        messages = []
        if sys_prompt:
            messages.append({'role': 'system', 'content': sys_prompt})
        messages.append({'role': 'user', 'content': prompt})
        response: ChatResponse = await self.chat(messages, stop)
        return response

    @staticmethod
    def parse_llm_response(content: str) -> Optional[Union[dict, list]]:
        """Parse LLM response and extract JSON content."""
        json_match = re.search(r'```json\n(.*)\n```', content, re.DOTALL)
        if json_match:
            return json.loads(json_match.group(1))
        return None


def parse_llm_json_output(llm_output: str, verbose: bool = False) -> Optional[Union[dict, list]]:
    """
    Attempt to parse JSON object or array from potentially irregular LLM output.

    This function tries the following strategies:
    1. Directly parse the entire input string.
    2. Look for JSON in Markdown code blocks (e.g., ```json ... ```).
    3. Use general regex to find the first or largest JSON structure.
    4. Process the extracted candidate string:
        a. Try parsing with `ast.literal_eval` (handles Python-like dict/list,
           including single quotes, None/True/False, etc.).
        b. If `ast.literal_eval` fails, apply string cleaning:
           - Convert Python's None/True/False to JSON's null/true/false.
           - Remove trailing commas.
           - Remove C-style (//) and Python-style (#) comments.
           - Try parsing with `json.loads`.
        c. As a last resort for `json.loads`, if specific errors indicate quote issues,
           try replacing single quotes with double quotes and parse again.

    Args:
        llm_output: String output from LLM.
        verbose: If True, print diagnostic information.

    Returns:
        If parsing succeeds, returns a Python dict or list; otherwise returns None.
    """
    if not llm_output or not isinstance(llm_output, str):
        if verbose:
            print("Input is not a valid string or is empty.")
        return None

    # Strategy 1: Try parsing the entire input directly
    try:
        return json.loads(llm_output)
    except json.JSONDecodeError:
        if verbose:
            print("Direct parsing of entire input failed, attempting extraction and cleaning...")
        pass

    json_str_candidate = None

    # Strategy 2: Look for JSON in Markdown code blocks
    # Match ```json ... ``` or ``` ... ``` (if containing valid JSON object/array)
    # re.DOTALL makes '.' match newlines
    markdown_match = re.search(r'```(?:json)?\s*(\{[\s\S]*?\}|\[[\s\S]*?\])\s*```', llm_output, re.DOTALL)
    if markdown_match:
        json_str_candidate = markdown_match.group(1).strip()
        if verbose:
            print(f"Extracted from Markdown code block: {json_str_candidate[:200]}...")
    else:
        # Strategy 3: Use general regex to find first/largest JSON structure
        # This regex tries to find the most complete object {} or array []
        # It's greedy and will try to match the largest valid structure.
        # [\s\S] is used to match any character including newlines.
        general_match = re.search(r'(\{[\s\S]*\}|\[[\s\S]*\])', llm_output)
        if general_match:
            json_str_candidate = general_match.group(0).strip()
            if verbose:
                print(f"Extracted via general regex: {json_str_candidate[:200]}...")

    if not json_str_candidate:
        if verbose:
            print("Failed to find potential JSON structure via regex.")
        return None

    # Strategy 4: Process the extracted candidate string

    # 4.a. Try parsing with ast.literal_eval
    # ast.literal_eval is effective for Python-like dict/list (single quotes, None/True/False)
    try:
        if verbose:
            print(f"Attempting to parse with ast.literal_eval: {json_str_candidate[:200]}...")
        # ast.literal_eval needs Python's None, True, False
        # If candidate string is already in JSON style null, true, false, convert them first
        temp_for_ast = json_str_candidate
        temp_for_ast = re.sub(r'\bnull\b', 'None', temp_for_ast)
        temp_for_ast = re.sub(r'\btrue\b', 'True', temp_for_ast)
        temp_for_ast = re.sub(r'\bfalse\b', 'False', temp_for_ast)

        # ast.literal_eval is sensitive to trailing commas, try to remove them
        temp_for_ast = re.sub(r',\s*([\}\]])', r'\1', temp_for_ast)
        
        # Remove comments, ast.literal_eval doesn't support comments
        temp_for_ast = re.sub(r"//.*", "", temp_for_ast)
        temp_for_ast = re.sub(r"#.*", "", temp_for_ast)
        temp_for_ast = temp_for_ast.strip()

        parsed_object = ast.literal_eval(temp_for_ast)
        if verbose:
            print("ast.literal_eval parsing successful.")
        # Ensure result is dict or list, as expected for JSON
        if isinstance(parsed_object, (dict, list)):
            return parsed_object
        else:
            if verbose:
                print(f"ast.literal_eval produced non-dict/list type: {type(parsed_object)}")
            # If not a collection type, continue trying as it might just be a simple string parsed by ast.literal_eval
    except (ValueError, SyntaxError, TypeError) as e_ast:
        if verbose:
            print(f"ast.literal_eval failed: {e_ast}. Continuing with json.loads attempt.")
        pass # If ast.literal_eval fails, continue with subsequent string cleaning and json.loads attempt

    # 4.b. String cleaning for json.loads
    cleaned_str = json_str_candidate

    # Convert Python's None/True/False to JSON's null/true/false (if not already converted)
    # Use \b to ensure we're replacing whole words
    replacements = {
        r'\bNone\b': 'null',
        r'\bTrue\b': 'true',
        r'\bFalse\b': 'false'
    }
    for py_val, json_val in replacements.items():
        cleaned_str = re.sub(py_val, json_val, cleaned_str)

    # Remove C-style (//) and Python-style (#) comments
    cleaned_str = re.sub(r"//.*", "", cleaned_str)  # Remove // to end of line
    cleaned_str = re.sub(r"#.*", "", cleaned_str)   # Remove # to end of line

    # Remove trailing commas (before '}' or ']')
    cleaned_str = re.sub(r',\s*([\}\]])', r'\1', cleaned_str)
    cleaned_str = cleaned_str.strip()  # Remove leading/trailing whitespace

    # Try parsing cleaned string with json.loads
    try:
        if verbose:
            print(f"Attempting to parse cleaned string with json.loads: {cleaned_str[:200]}...")
        return json.loads(cleaned_str)
    except json.JSONDecodeError as e_json:
        if verbose:
            print(f"json.loads failed on cleaned string: {e_json}")

        # 4.c. Last resort: If error indicates quote issues, try global single quote replacement
        # Similar to "Expecting property name enclosed in double quotes" error message
        # Or "Unterminated string starting at" (if single quotes exist in candidate string)
        if ("double quotes" in str(e_json) or \
            ("Unterminated string" in str(e_json) and "'" in cleaned_str)) and \
           "'" in json_str_candidate:  # Check if single quotes exist in original candidate string

            if verbose:
                print("Error suggests quote issues. Attempting global single to double quote replacement...")

            # Operate on original extracted json_str_candidate to avoid accumulated errors
            desperate_str = json_str_candidate.replace("'", '"')

            # Reapply other cleaning steps
            for py_val, json_val in replacements.items():  # None, True, False
                desperate_str = re.sub(py_val, json_val, desperate_str)
            desperate_str = re.sub(r"//.*", "", desperate_str)
            desperate_str = re.sub(r"#.*", "", desperate_str)
            desperate_str = re.sub(r',\s*([\}\]])', r'\1', desperate_str)
            desperate_str = desperate_str.strip()

            try:
                if verbose:
                    print(f"Attempting to parse last-resort cleaned string with json.loads: {desperate_str[:200]}...")
                return json.loads(desperate_str)
            except json.JSONDecodeError as e_desperate:
                if verbose:
                    print(f"Last-resort cleaning attempt also failed: {e_desperate}")
                pass

        if verbose:
            print("All parsing attempts failed.")
        return None



if __name__ == "__main__":
    import asyncio

    llm = OpenAIWrapper()
    res = asyncio.run(llm.chat([{"role": "user", "content": "Hello!"}]))
    print(res)
    ...
    