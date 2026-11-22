import os
from dotenv import load_dotenv

load_dotenv()

# Local imports
from configs.config_cls import LLMConfig


LLM_CONFIG = LLMConfig(
    base_url=os.getenv("LLM_BASE_URL", ""),
    api_key=os.getenv("LLM_API_KEY", ""),
    model="google/gemini-2.5-flash-preview",
    temperature=0.7,
    max_tokens=8192,
    top_p=1.0,
    frequency_penalty=0.0,
    presence_penalty=0.0,
    proxy=os.getenv("LLM_PROXY"),
    completion_cost=1e-5,
    prompt_cost=2.5e-6,
    timeout=180.0
)