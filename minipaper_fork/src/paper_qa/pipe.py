import asyncio
import json
import pandas as pd
from copy import deepcopy
from pathlib import Path


# Local imports
from get_data import get_geo_metadata
from nlp_process.query import QueryEngine
from nlp_process.llm import OpenAIWrapper, parse_llm_json_output, ChatResponse
from utils.logger import logger
from paper_qa.prompt import GENERATE_QUERY, DEDUCT_RETRIEVAL, ANSWER
from paper_qa.config import LLM_CONFIG


class PaperQAPipe:
    def __init__(self):
        self.llm = OpenAIWrapper(config=LLM_CONFIG)
        self.domain = "clockbase_paper_202506"
        self.query_engine = QueryEngine()
        self.history = []


    async def generate_query(self, user_input: str):
        messages = deepcopy(self.history)
        prompt = GENERATE_QUERY.format(user_input=user_input)
        messages.append({'role': 'user', 'content': prompt})
        resp = await self.llm.chat(messages)
        resp = parse_llm_json_output(resp.content)
        return resp

    async def retrieve(self, query: str):
        result, url_map = self.query_engine.pipeline(query, self.domain)
        return result, url_map

    async def parallel_retrieve(self, queries: list):
        results = await asyncio.gather(*[self.retrieve(query) for query in queries])
        final_url_map = dict()
        final_result = []
        for result, url_map in results:
            final_url_map.update(url_map)
            final_result.extend(result)
        return final_result, final_url_map

    async def deduct_evidence(self, question: str, retrieved: str):
        """
        Extract evidence from retrieved text.
        """
        prompt = DEDUCT_RETRIEVAL.format(question=question, retrieved=retrieved)
        resp = await self.llm.chat_once(prompt)
        if 'NO RELATED' in resp.content.upper():
            return ''
        return resp.content

    async def parallel_deduct_evidence(self, question: str, retrieved: list):
        results = await asyncio.gather(*[self.deduct_evidence(question, r) for r in retrieved])
        return results


    async def answer(self, question: str, evidence: str):
        prompt = ANSWER.format(question=question, evidence=evidence)
        resp = await self.llm.chat_once(prompt)
        return resp.content


    async def pipeline(self, user_input: str):
        resp = await self.generate_query(user_input)
        if resp is None:
            return 'I am sorry, I cannot answer your question.'
        if resp['retrieval'] == 'YES':
            queries = resp['queries']
            question = resp['question']
            results, url_map = await self.parallel_retrieve(queries)
            evidences = await self.parallel_deduct_evidence(question, results)
            # for e in evidences:
            #     if '<image>' in e:
            #         print('****************')
            #         print(e)
            #         print('****************')
            valid_evidences ='\n-----\n'.join([e for e in evidences if e])
            answer = await self.answer(question, valid_evidences)
            for key, url in url_map.items():
                answer = answer.replace(f'<image>{key}</image>', f'![]({url})')
                answer = answer.replace(f'{key}', f'![]({url})')
                answer = answer.replace(f'<table>', '')
                answer = answer.replace(f'</table>', '')
            print(answer)
            return answer
        else:
            resp = await self.llm.chat_once(user_input)
            return resp.content



if __name__ == "__main__":
    pipe = PaperQAPipe()
    question = 'Frequently mentioned potential “routine laboratory” biomarkers'
    # question = "How to evaluate the brain's age?"
    result = asyncio.run(pipe.pipeline(question))
    print(result)




