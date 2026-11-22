import os
from pathlib import Path
from configs.config_cls import (
    LoggerConfig, 
    LLMConfig, 
    QueryConfig, 
    ToKnowledgeConfig, 
    Neo4jConfig
)
from dotenv import load_dotenv


load_dotenv()   

# Logger config
LOGGER_CONFIG = LoggerConfig(
    name='clockbase',
    file_name='clockbase.log',
    level='INFO',
    when='midnight',
    interval=1,
    backup_count=30,
    delay=False,
    utc=True
)

# LLM config
SYS_LLM_CONFIG = LLMConfig(
    base_url=os.getenv("LLM_BASE_URL", ""),
    api_key=os.getenv("LLM_API_KEY", ""),
    # model="openai/o4-mini-2025-04-16",
    # model="google/gemini-2.5-pro-preview-05-06",
    # model="deepseek-reasoner",
    # model="google/gemini-2.5-flash-preview",
    model="claude-sonnet-4-20250514",
    temperature=0.7,
    max_tokens=8192,
    top_p=1.0,
    frequency_penalty=0.0,
    presence_penalty=0.0,
    proxy=os.getenv("LLM_PROXY"),
    completion_cost=1.5e-5,
    prompt_cost=3e-6,
    timeout=180.0
)




QUERY_CONFIG = QueryConfig(
    endpoint=os.getenv("RAG_ENDPOINT", ""),
    threshold=0.4,
    topk=5,
    output_fields=["text", "file_name", "agg_index", "url"]
)


TO_KNOWLEDGE_CONFIG = ToKnowledgeConfig(
    semaphore=16,
    domain='clockbase_paper_202506',
    rag_insert_api=f'{os.getenv("RAG_ENDPOINT")}/api/v1/insert/insert_pdf',
    ocr_apis=[url for url in os.getenv('OCR_ENDPOINTS', '').split(',') if url],
    emb_apis=[url for url in os.getenv('EMBEDDING_ENDPOINTS', '').split(',') if url]
)

NEO4J_CONFIG = Neo4jConfig()


# Which CSV modality: "mm_rna" (mouse RNA, original), "hs_rna" (human RNA),
# or "hs_dnam" (human DNA methylation, multi-clock panel).
# run.py picks the IDENTIFY_RELEVANCE / GENERATE_PAPER_2 variant accordingly.
CSV_MODALITY = os.getenv("CSV_MODALITY", "mm_rna")

CSV_PATH = Path(__file__).parent.parent.joinpath("data", "sample.csv")
PAPER_DOWNLOAD_DIR = Path(__file__).parent.parent.joinpath("data", "papers")
PAPER_DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)