import requests
from collections import defaultdict
import hashlib
import re
from pathlib import Path
from typing import List, Optional

# Local modules
from configs.global_config import QUERY_CONFIG
from configs.config_cls import QueryConfig
from utils.crossref import resolve_reference, PaperMetadata
from utils.logger import logger
from utils.http_client import safe_get, safe_post  # timeout + retry + proxy fix


generate_md5 = lambda text: hashlib.md5(text.encode('utf-8')).hexdigest()

TEMPLATE_PLAIN = '---\n**File Name**: {file_name}\n\n**Reference**:\n\n{text}\n\n---\n\n'
TEMPLATE_ENRICHED = '---\n{header}\n\n**Reference**:\n\n{text}\n\n---\n\n'


class QueryEngine:
    def __init__(self, config: QueryConfig = None):
        self.config = config or QUERY_CONFIG
        self._crossref_cache: dict[str, PaperMetadata] = {}

    def test_query(self, query, domain, threshold=0.4, topk=5, output_fields=["text", "file_name", "agg_index", "url"]):
        data = {
            "csid": "test",
            "query": query,
            "domain": domain,
            "topk": topk,
            "output_fields": output_fields,
            "threshold": threshold
        }
        response = safe_post(f"{self.config.endpoint}/api/v1/query/search", json=data)
        return response.json()

    def parse(self, documents: List[dict], domain: str, endpoint: str, enrich: bool = False) -> tuple[list[str], dict[str, str]]:
        result = defaultdict(list)
        url_map = dict()
        index = 0
        for record in documents:
            file_name = record['file_name']
            agg_index = record['agg_index']
            text = record['text']
            index = 0
            file_md5 = generate_md5(Path(file_name).stem)
            raw_urls = re.findall(r"images/[a-fA-F0-9]{64}\.jpg", text)
            img_md5s = re.findall(r"images/([a-fA-F0-9]{64})\.jpg", text)
            new_urls = [f'{endpoint}/api/v1/storage/image?domain={domain}&file={file_md5}&image={image_md5}' for image_md5 in img_md5s]
            for i, url in enumerate(new_urls):
                url_sub = f'<image>{file_md5}-{agg_index}-{index}</image>'
                text = re.sub(raw_urls[i], url_sub, text)
                url_map[f"{file_md5}-{agg_index}-{index}"] = url
                index += 1
            result[file_name].append((agg_index, text))

        for k, v in result.items():
            result[k] = sorted(v)
            result[k] = '\n\n'.join([ele[1] for ele in v])

        if enrich:
            formatted = []
            for file_name, text in result.items():
                meta = self._resolve_cached(file_name, text)
                formatted.append(TEMPLATE_ENRICHED.format(header=meta.header_block, text=text))
            return formatted, url_map

        formatted = [TEMPLATE_PLAIN.format(file_name=fn, text=t) for fn, t in result.items()]
        return formatted, url_map

    def _resolve_cached(self, file_name: str, document_text: str) -> PaperMetadata:
        if file_name in self._crossref_cache:
            return self._crossref_cache[file_name]
        try:
            meta = resolve_reference(file_name, document_text)
        except Exception as e:
            logger.warning(f'CrossRef lookup failed for {file_name}: {e}')
            meta = PaperMetadata(title=file_name)
        self._crossref_cache[file_name] = meta
        return meta

    def pipeline(self, query, domain, threshold=0.4, topk=5, output_fields=["text", "file_name", "agg_index", "url"], enrich: bool = False):
        result = self.test_query(query, domain, threshold, topk, output_fields)
        result, url_map = self.parse(result['result'], domain, self.config.endpoint, enrich=enrich)
        return result, url_map

    
    def count_docs(self, domain):
        url = f"{self.config.endpoint}/api/v1/maintain/list_paper"
        body = {
            "domain": domain,
            "csid": "test",
        }
        response = safe_post(url, json=body)
        return len(response.json()['result'])


if __name__ == "__main__":
    
    domain = "clockbase_paper_202506"
    query_engine = QueryEngine()

    query = 'wearable device data'
    result, url_map = query_engine.pipeline(query, "clockbase_paper_202506")
    print(result)
    print(url_map)
    ...
    # print(query_engine.count_docs(domain))
