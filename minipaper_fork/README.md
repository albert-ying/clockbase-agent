# ClockBase

An LLM-powered pipeline for analyzing biological aging clock data from GEO datasets. ClockBase retrieves experimental metadata, searches a curated knowledge base of scientific literature via RAG, and uses large language models to assess the relevance of interventions to aging — ultimately producing structured scientific reports.

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                      Batch Pipeline (run.py)                    │
│                                                                 │
│  CSV Data ──► GEO Metadata ──► LLM Query Gen ──► RAG Retrieval │
│                                       │                         │
│                                       ▼                         │
│                              LLM Relevance Analysis             │
│                                       │                         │
│                                       ▼                         │
│                              Report Generation                  │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                      QA API (api.py)                            │
│                                                                 │
│  User Question ──► Intent + Query Gen ──► Parallel RAG Retrieve │
│                                                  │              │
│                                                  ▼              │
│                                       Evidence Deduction        │
│                                                  │              │
│                                                  ▼              │
│                                          Final Answer           │
└─────────────────────────────────────────────────────────────────┘
```

## Project Structure

```
clockbase/
├── src/
│   ├── orchestrate.py              # E2E pipeline orchestrator (6 stages)
│   ├── run.py                      # Batch processing pipeline
│   ├── api.py                      # FastAPI QA service
│   ├── get_data.py                 # GEO / PubMed metadata
│   ├── configs/
│   │   ├── config_cls.py           # Pydantic config models
│   │   └── global_config.py        # Runtime config instances
│   ├── nlp_process/
│   │   ├── llm.py                  # OpenAI-compatible LLM wrapper
│   │   ├── query.py                # RAG query engine
│   │   ├── prompt_template.py      # Prompts for batch pipeline
│   │   └── to_knowledge.py         # PDF ingestion into knowledge base
│   ├── paper_qa/
│   │   ├── pipe.py                 # Paper QA multi-agent pipeline
│   │   ├── prompt.py               # Prompts for QA pipeline
│   │   └── config.py               # QA-specific LLM config
│   ├── result_parse/
│   │   ├── result_process.py       # Scoring, ranking, FDR filtering
│   │   ├── clean.py                # Markdown reference cleanup
│   │   ├── convert.py              # Markdown → HTML conversion
│   │   └── run.py                  # Batch HTML conversion
│   ├── store/
│   │   └── neo4j_store.py          # Neo4j graph database interface
│   └── utils/
│       ├── logger.py               # Rotating file logger
│       ├── request.py              # Semantic Scholar / PDF download
│       └── helpers.py              # Retry decorator and utilities
├── tests/
│   └── test_api.py                 # API smoke test
├── data/                           # Runtime data (gitignored)
├── pyproject.toml
└── .env                            # Secrets & endpoints (gitignored)
```

## Prerequisites

- Python >= 3.13
- [uv](https://docs.astral.sh/uv/) (recommended package manager)
- Access to an OpenAI-compatible LLM API
- A running RAG service endpoint (for document retrieval)
- Neo4j database (optional, for paper metadata storage)
- Network proxy (for GEO/PubMed API access, if needed)

## Setup

1. **Clone and install dependencies:**

```bash
git clone <repo-url> && cd clockbase
uv sync
```

> This creates a `.venv` automatically and installs all dependencies from `pyproject.toml`.

2. **Configure environment variables:**

Create a `.env` file in the project root:

```env
# LLM
LLM_BASE_URL=https://api.openai.com/v1
LLM_API_KEY=sk-...

# RAG knowledge base
RAG_ENDPOINT=http://your-rag-host:port

# Network proxy (for GEO/PubMed access)
PROXY_HOST=proxy.example.com
PROXY_PORT=7777
PROXY_USER=your_user
PROXY_PWD=your_password

# Neo4j
NEO4J_PASSWORD=your_neo4j_password

# API auth token
LOCAL_TOKEN=your_secret_token
```

3. **Prepare input data:**

Place your main dataset CSV at `data/main.csv` (or any path you specify). The CSV should contain columns: `gse_id`, `comparison_name`, `condition_name`, `condition_category`, `effect_size`, `direction`, `pvalue`, `outlier_status`, `fdr`.

## Usage

### End-to-End Pipeline

The orchestrator chains all 6 processing stages into a single command.

**Quick start with demo data** (single sample, GSE196908 — Lanifibranor anti-aging study):

```bash
uv run --directory src orchestrate.py --df-path ../demo_data/demo_sample.csv --workdir ../demo_output
```

**Full run with your own data:**

```bash
uv run --directory src orchestrate.py --df-path ../data/main.csv
```

**Stages:**

| Stage | Description |
|-------|-------------|
| 1 | **Batch Processing** — Fetch GEO metadata, generate queries, retrieve documents (with CrossRef reference resolution), analyze relevance |
| 2 | **Scoring & FDR Filter** — Rank results by `final_score`, filter by FDR threshold |
| 3 | **Case Selection** — Copy top-N cases to report directory |
| 4 | **Report Generation** — Generate scientific reports via LLM |
| 5 | **Reference Cleanup** — Fix citation numbering and remove invalid references |
| 6 | **HTML Conversion** — Convert Markdown reports to styled HTML |

**Resume from a specific stage** (e.g., stages 1-2 already done):

```bash
uv run --directory src orchestrate.py --df-path ../demo_data/demo_sample.csv --workdir ../demo_output --start-stage 3
```

**Custom parameters:**

```bash
uv run --directory src orchestrate.py \
  --df-path ../data/main.csv \
  --workdir ../my_run \
  --top-n 500 \
  --fdr-threshold 0.01 \
  --num-semaphore 16
```

> By default, `--top-n` is unlimited (processes all cases) and `--fdr-threshold` is `1.0` (no filtering), so the pipeline works correctly even with a single data entry.

### Running Individual Stages

Each stage can also be run independently via its own script:

```bash
# Stage 1: Batch processing
uv run --directory src run.py

# Stage 2: Scoring & filtering
uv run --directory src result_parse/result_process.py

# Stage 5-6: Clean & convert
uv run --directory src result_parse/clean.py
uv run --directory src result_parse/run.py
```

### QA API Server

Start the FastAPI server for interactive paper-based question answering:

```bash
uv run --directory src api.py
# or
uv run --directory src uvicorn api:app --host 0.0.0.0 --port 8000
```

**Endpoints:**

| Method | Path | Description |
|--------|------|-------------|
| `POST` | `/qa` | Ask a question (body: `{"token": "...", "user_input": "..."}`) |
| `GET` | `/health` | Health check |

**Example request:**

```bash
curl -X POST http://localhost:8000/qa \
  -H "Content-Type: application/json" \
  -d '{"token": "your_token", "user_input": "What biomarkers are used for brain aging?"}'
```

## External Dependencies

| Service | Purpose | Config Key |
|---------|---------|------------|
| LLM API | Query generation, relevance analysis, report writing | `LLM_BASE_URL`, `LLM_API_KEY` |
| RAG Endpoint | Document retrieval from curated paper knowledge base | `RAG_ENDPOINT` |
| Neo4j | Paper metadata storage (DOI, PDF links) — **optional**, only used by the standalone `add_doi_url_to_top()` helper, not required for the E2E pipeline | `NEO4J_*` |
| CrossRef | Reference resolution — resolves RAG file names to canonical paper metadata (title, DOI, authors, journal) using a 3-level fallback: DOI extraction → title search → text search | no key needed |
| NCBI GEO | Experimental metadata for GSE datasets | via proxy |
| PubMed | Citation details (title, DOI, authors) | via proxy |

## Status

> **This project is under active development and not yet feature-complete.**

### Roadmap

- [ ] **Dynamic literature knowledge base updates** — Integrate automated paper discovery, download, and ingestion so the RAG knowledge base stays current without manual intervention (leveraging `nlp_process/to_knowledge.py`)
- [ ] Automated test suite (unit + integration)
- [ ] Pre-flight health checks for external services (LLM, RAG, GEO API)

## License

All rights reserved.
