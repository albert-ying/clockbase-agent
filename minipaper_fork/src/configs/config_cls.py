from typing import Optional, Literal, Dict, List
from pathlib import Path
from pydantic import BaseModel
from pydantic_settings import BaseSettings, SettingsConfigDict
from pydantic import SecretStr


BASE_PATH = Path(__file__).parent.parent


# Init cache directory
CACHE = BASE_PATH.joinpath('.cache')


# Logger config
class LoggerConfig(BaseSettings):
    name: str
    file_name: str = 'task_graph.log'
    level: Literal['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'] = 'INFO'
    when: str = 'midnight'    # S-Seconds, M-Minutes, H-Hours, D-Days, midnight, W{0-6}-certain day
    interval: int = 1
    backup_count: int = 30
    delay: bool = False
    utc: bool = True


class LLMConfig(BaseSettings):
    model_config = SettingsConfigDict(env_file=".env", env_prefix="LLM_", extra="ignore")

    base_url: Optional[str] = None
    model: str = "gpt-4o-mini"
    temperature: float = 0.7
    max_tokens: int = 4096
    top_p: float = 1.0
    frequency_penalty: float = 0.0
    presence_penalty: float = 0.0
    api_key: SecretStr
    completion_cost: float = 1e-5    # The cost of completion tokens
    prompt_cost: float = 2.5e-6        # The cost of prompt tokens
    proxy: Optional[str] = None
    timeout: float = 300.0
    system_prompt: Optional[str] = None


class QueryConfig(BaseSettings):
    endpoint: str
    threshold: float = 0.4
    topk: int = 5
    output_fields: list[str] = ["text", "file_name", "agg_index", "url"]


class ToKnowledgeConfig(BaseModel):
    semaphore: int = 4
    domain: str    # Knowledge base name
    rag_insert_api: str
    ocr_apis: List[str]
    emb_apis: List[str]


class Neo4jConfig(BaseSettings):
    model_config = SettingsConfigDict(env_file=".env", env_prefix="NEO4J_", extra="ignore")

    uri: str = "bolt://localhost:7687"
    user: str = "neo4j"
    password: SecretStr
    database: str = "neo4j"
    
    
    