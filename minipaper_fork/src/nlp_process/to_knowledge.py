import uuid
import asyncio
import aiohttp
import sys
import requests
from pathlib import Path
from typing import Tuple, List, Optional
from concurrent.futures.thread import ThreadPoolExecutor


# Local modules
from utils.helpers import generate_md5, bytes_to_b64, RoundRobin
from utils.logger import logger
from configs.config_cls import ToKnowledgeConfig
from configs.global_config import TO_KNOWLEDGE_CONFIG
from utils.http_client import safe_get, safe_post  # timeout + retry + proxy fix


mac = uuid.UUID(int=uuid.getnode()).hex[-12:].upper()
MAC = "-".join([mac[e:e+2] for e in range(0, 11, 2)])


class ToKnowledge:
    def __init__(self, config: Optional[ToKnowledgeConfig] = None):
        self.config = config or TO_KNOWLEDGE_CONFIG
        self.ocr_apis = RoundRobin(self.config.ocr_apis)
        self.emb_apis = RoundRobin(self.config.emb_apis)

    def _insert(self, domain:str, file_name:str, file_bytes:bytes):
        body = {
            'file_bs64': bytes_to_b64(file_bytes),
            'domain': domain,
            'file_name': file_name,
            'csid': MAC,
            'ocr_api': self.ocr_apis.get_next(),
            'emb_api': self.emb_apis.get_next()
        }
        try:
            resp = safe_post(self.config.rag_insert_api, json=body)
            # Request succeed
            if resp.status_code == 200:
                result = resp.json()
                # Task succeed
                if result['status'] == 200:
                    logger.info(f'Insert successfully: {file_name}')
                # Task fail
                else:
                    logger.error(f"Insert paper failed: {file_name}\nMessage: {result['message']}")
            # Request fail
            else:
                logger.error(f'Insert paper failed - {resp.status_code}: {file_name}')
        # Other error
        except Exception as err:
            logger.error(f'Insert paper Error: {repr(err)}')


    async def a_insert(self, domain:str, file_name:str, file_bytes:bytes):
        body = {
        'file_bs64': bytes_to_b64(file_bytes),
        'domain': domain,
        'file_name': file_name,
        'csid': MAC,
        'ocr_api': self.ocr_apis.get_next(),
        'emb_api': self.emb_apis.get_next()
    }
        try:
            async with aiohttp.ClientSession() as sess:
                async with sess.post(self.config.rag_insert_api, json=body) as resp:
                    # Request succeed
                    if resp.status == 200:
                        result = await resp.json()
                        # Task succeed
                        if result['status'] == 200:
                            # Remove successful task from Redis
                            logger.info(f'Insert successfully: {file_name}')
                        # Task fail
                        else:
                            logger.error(f"Insert paper failed: {file_name}\nMessage: {result['message']}")
                    # Request fail
                    else:
                        logger.error(f'Insert paper failed - {resp.status}: {file_name}')
        # Other error
        except Exception as err:
            logger.error(f'Insert paper Error: {repr(err)}')
                

    def insert(self, pdf_info: Tuple[Path, Optional[str]]):
        pdf_path, file_name = pdf_info
        with open(str(pdf_path), "rb") as f:
            pdf_bytes = f.read()
        file_name = file_name or pdf_path.stem
        # Skip large pdf more than 20M
        size = sys.getsizeof(pdf_bytes)
        if size <= 20000000:
            self._insert(self.config.domain, file_name, pdf_bytes)
        else:
            logger.warning(f'Skip large file: {file_name}')


    def pipe(self, pdfs: List[Tuple[Path, Optional[str]]]):
        logger.info(f'Tasks count: {len(pdfs)}')
        with ThreadPoolExecutor(max_workers=self.config.semaphore) as executor:
            executor.map(self.insert, pdfs)
        # for doi in dois:
        #     self.insert(doi)


if __name__ == '__main__':
    import json

    to_knowledge = ToKnowledge()

    types = ['significant_df_Drug', 'significant_df_Disease']
    data_dir = Path(__file__).parent.parent.joinpath('data')
    paper_dir = Path(__file__).parent.parent.joinpath('data', 'papers')
    for t in types:
        processed_file = data_dir.joinpath(t, 'processed.json')
        if processed_file.exists():
            with open(str(processed_file), 'r') as f:
                processed_data = json.load(f)
        else:
            processed_data = dict()

        success_file = data_dir.joinpath(t, 'success.json')
        with open(success_file, 'r') as f:
            success_data = json.load(f)
        for doi, info in success_data.items():
            if doi in processed_data:
                continue
            title = info['title']
            pdf_path = paper_dir.joinpath(f'{title}.pdf')
            if pdf_path.exists():
                try:
                    to_knowledge.insert((pdf_path, title))
                    processed_data[doi] = True
                    with open(str(processed_file), 'w') as f:
                        json.dump(processed_data, f, indent=4)
                except Exception as err:
                    logger.error(f'Insert paper failed: {doi}\nError: {repr(err)}')
        