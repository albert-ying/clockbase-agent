# Minipaper source snapshot

This folder holds the source of the biological interpreter, scorer and manuscript
agent. It is taken from https://github.com/Ard-Skelling/minipaper at base commit
`b4570d446f7a1aa0a0131010ffb217846f744a93`, with the changes made for ClockBase to
seven tracked source files, plus `src/utils/http_client.py` and
`run_both_modalities.sh`.

Two files differ further from that source. `src/get_data.py` keeps the GEO and
PubMed metadata functions that the pipeline uses; the standalone paper-download
utility is not part of this release. `src/configs/global_config.py` reads the
addresses of the OCR and embedding services from the environment variables
`OCR_ENDPOINTS` and `EMBEDDING_ENDPOINTS` (comma-separated).

Credentials, local environments, logs, status files and generated outputs are
not included. The dependency manifest and lockfile are included. This copy has
not been validated as a complete RAG or writing-agent run.
