"""
ClockBase data processing and analysis main module

This module handles circadian clock related data, including:
- Fetching metadata from the GEO database
- Using LLMs to generate queries and analyze relevance
- Retrieving related literature and generating reports
- Managing paper information in the Neo4j database
"""

import asyncio
import json
import re
import shutil
from pathlib import Path
import pandas as pd

# Local module imports
from get_data import get_geo_metadata  # GEO database metadata retrieval
from nlp_process.query import QueryEngine  # Query engine
from nlp_process.prompt_template import (
    IDENTIFY_RELEVANCE, GENERATE_QUERY, GENERATE_PAPER_2,
    IDENTIFY_RELEVANCE_HUMAN_RNA, IDENTIFY_RELEVANCE_HUMAN_DNAM,
    GENERATE_PAPER_2_HUMAN_RNA, GENERATE_PAPER_2_HUMAN_DNAM,
)  # Prompt templates
from nlp_process.llm import OpenAIWrapper, parse_llm_json_output, ChatResponse  # LLM wrapper
from configs.global_config import SYS_LLM_CONFIG, CSV_PATH, CSV_MODALITY  # Global configuration

# Modality → (identify template, paper template) lookup
_IR_TEMPLATES = {
    "mm_rna":  IDENTIFY_RELEVANCE,
    "hs_rna":  IDENTIFY_RELEVANCE_HUMAN_RNA,
    "hs_dnam": IDENTIFY_RELEVANCE_HUMAN_DNAM,
}
_PAPER_TEMPLATES = {
    "mm_rna":  GENERATE_PAPER_2,
    "hs_rna":  GENERATE_PAPER_2_HUMAN_RNA,
    "hs_dnam": GENERATE_PAPER_2_HUMAN_DNAM,
}
# The 9 extra fields the DNAm prompt needs from the CSV row
_DNAM_EXTRA_FIELDS = [
    "PhenoAge_effect_size", "PhenoAge_pvalue",
    "Horvathv1_effect_size", "Horvathv1_pvalue",
    "YingDamAge_effect_size", "YingDamAge_pvalue",
    "DunedinPACE_effect_size", "DunedinPACE_pvalue",
    "n_support",
]
from utils.logger import logger  # Logger
from store.neo4j_store import Neo4jStore  # Neo4j data store


def load_all_data():
    """
    Load data from a CSV file
    
    Returns:
        list: A list of dictionaries containing all records
    """
    df = pd.read_csv(CSV_PATH)
    records = df.to_dict(orient="records")
    return records

async def generate_query(title, summary, study_design, condition, comparison, clock, condition_cate):
    """
    Generate query statements using an LLM
    
    Args:
        title (str): Study title
        summary (str): Study summary
        study_design (str): Study design
        condition (str): Experimental condition
        comparison (str): Control condition
        clock (str): Biological clock type
        condition_cate (str): Condition category
        
    Returns:
        list: List of generated query statements
    """
    llm = OpenAIWrapper(config=SYS_LLM_CONFIG)
    # Prepare keyword arguments for query generation
    gq_kw = {
        'title': title,
        'summary': summary,
        'study_design': study_design,
        'condition': condition,
        'comparison': comparison,
        'clock': clock,
        'condition_cate': condition_cate
    }
    # Format prompt template and send to LLM
    generate_query_prompt = GENERATE_QUERY.format(**gq_kw)
    queries: ChatResponse = await llm.chat_once(generate_query_prompt)
    # Parse the JSON-formatted query results returned by LLM
    queries = parse_llm_json_output(queries.content)
    return queries

async def identify_relevance(gse_id, title, summary, study_design, condition, comparison, condition_cate, effect_size, effect_direction, effect_pvalue, extreme_value, retrieved_documents, intervention_source_literature, modality: str = "mm_rna", extra_cross_clock: dict | None = None):
    """
    Use LLM to identify the relevance of a study to biological clocks
    
    Args:
        gse_id (str): GEO dataset ID
        title (str): Study title
        summary (str): Study summary
        study_design (str): Study design
        condition (str): Experimental condition
        comparison (str): Control condition
        condition_cate (str): Condition category
        effect_size (float): Effect size
        effect_direction (str): Effect direction
        effect_pvalue (float): P-value
        extreme_value (str): Outlier status
        retrieved_documents (str): Retrieved documents
        intervention_source_literature (str): Intervention source literature
        
    Returns:
        tuple: (relevance result, prompt content, response content, cost)
    """
    llm = OpenAIWrapper(config=SYS_LLM_CONFIG)
    # Prepare arguments for relevance identification
    ir_kw = {
        'gse_id': gse_id,
        'title': title,
        'summary': summary,
        'study_design': study_design,
        'condition': condition,
        'comparison': comparison,
        'condition_cate': condition_cate,
        'effect_size': effect_size,
        'effect_direction': effect_direction,
        'effect_pvalue': effect_pvalue,
        'extreme_value': extreme_value,
        'retrieved_documents': retrieved_documents,
        'intervention_source_literature': intervention_source_literature
    }
    # Inject the 9 cross-clock fields the DNAm prompt expects
    if modality == "hs_dnam" and extra_cross_clock:
        ir_kw.update({
            'phenoage_effect_size':    extra_cross_clock.get('PhenoAge_effect_size'),
            'phenoage_pvalue':         extra_cross_clock.get('PhenoAge_pvalue'),
            'horvathv1_effect_size':   extra_cross_clock.get('Horvathv1_effect_size'),
            'horvathv1_pvalue':        extra_cross_clock.get('Horvathv1_pvalue'),
            'yingdamage_effect_size':  extra_cross_clock.get('YingDamAge_effect_size'),
            'yingdamage_pvalue':       extra_cross_clock.get('YingDamAge_pvalue'),
            'dunedinpace_effect_size': extra_cross_clock.get('DunedinPACE_effect_size'),
            'dunedinpace_pvalue':      extra_cross_clock.get('DunedinPACE_pvalue'),
            'n_support':               extra_cross_clock.get('n_support'),
        })
    # Pick the right template for this modality
    ir_template = _IR_TEMPLATES.get(modality, IDENTIFY_RELEVANCE)
    # Format the relevance identification prompt and send to LLM
    identify_relevance_prompt = ir_template.format(**ir_kw)
    response: ChatResponse = await llm.chat_once(identify_relevance_prompt)
    cost = response.cost
    # Parse the relevance analysis result
    relevance = parse_llm_json_output(response.content)
    return relevance, identify_relevance_prompt, response.content, cost

async def generate_report(identify_relevance_prompt, response_content, report_prompt):
    """
    Generate a report based on relevance analysis results
    
    Args:
        identify_relevance_prompt (str): Relevance identification prompt
        response_content (str): LLM response content
        report_prompt (str): Report generation prompt
        
    Returns:
        tuple: (report content, cost)
    """
    llm = OpenAIWrapper(config=SYS_LLM_CONFIG)
    # Build the conversation message sequence to maintain context coherence
    messages = [
        {'role': 'user', 'content': identify_relevance_prompt},
        {'role': 'assistant', 'content': response_content},
        {'role': 'user', 'content': report_prompt}
    ]
    response: ChatResponse = await llm.chat(messages)
    return response.content, response.cost


async def generate_paper(identify_relevance_prompt, response_content, paper_prompt):
    """
    Generate an academic paper based on relevance analysis results
    
    Args:
        identify_relevance_prompt (str): Relevance identification prompt
        response_content (str): LLM response content
        paper_prompt (str): Paper generation prompt
        
    Returns:
        tuple: (paper content, cost)
    """
    llm = OpenAIWrapper(config=SYS_LLM_CONFIG)
    # Build the conversation message sequence to maintain context coherence
    messages = [
        {'role': 'user', 'content': identify_relevance_prompt},
        {'role': 'assistant', 'content': response_content},
        {'role': 'user', 'content': paper_prompt}
    ]
    response: ChatResponse = await llm.chat(messages)
    return response.content, response.cost


# Commented-out relevance judgment function - replaced by identify_relevance
# async def judge_relevance(intervention_name, comparison_context, biomarker_clock, scientific_report_content):
#     """
#     Judge the relevance of an intervention to biological clocks (deprecated)
#     This function has been replaced by identify_relevance, which provides more comprehensive relevance analysis
#     """
#     llm = OpenAIWrapper(config=SYS_LLM_CONFIG)
#     jr_kw = {
#         'intervention_name': intervention_name,
#         'comparison_context': comparison_context,
#         'biomarker_clock': biomarker_clock,
#         # 'aging_relevance_score': aging_relevance_score,
#         'scientific_report_content': scientific_report_content
#     }
#     judge_relevance_prompt = JUDGE_RELEVANCE.format(**jr_kw)
#     response: ChatResponse = await llm.chat_once(judge_relevance_prompt)
#     return response.content


async def pipeline(sample_index: int, sample: dict, output_base: Path = None):
    """
    Complete data processing pipeline
    
    Full processing flow for a single sample:
    1. Data preprocessing and field mapping
    2. Create output directory
    3. Fetch GEO metadata
    4. Generate query statements
    5. Retrieve related documents
    6. Analyze relevance
    7. Save results and costs
    
    Args:
        sample_index (int): Sample index
        sample (dict): Sample data dictionary
        output_base (Path): Output base directory, defaults to <project_root>/output
        
    Returns:
        float: Processing cost
    """
    if output_base is None:
        output_base = Path(__file__).parent.parent / 'output'
    cost = 0.0
    # Map sample data to database field names.
    # For DNAm, also keep the 9 cross-clock fields the prompt needs.
    fields = ['gse_id', 'comparison_name', 'condition_name', 'condition_category', 'effect_size', 'direction', 'pvalue', 'outlier_status']
    if CSV_MODALITY == "hs_dnam":
        fields = fields + _DNAM_EXTRA_FIELDS
    chosen = dict()
    for k in fields:
        if k in sample:
            chosen[k] = sample[k]
    sample = chosen
    del chosen

    # Create output directory, format: GSE_ID-index
    output_dir = output_base / f'{sample["gse_id"]}-{sample_index}'
    output_dir.mkdir(parents=True, exist_ok=True)
    # If the output directory already contains retrieved documents, the pipeline can be skipped (commented out)
    # if output_dir.joinpath('retrieved_documents.md').exists():
    #     logger.info(f'Output directory {output_dir} already exists')
    #     return None, None

    # Fetch GEO metadata for the sample
    metadata = get_geo_metadata(sample['gse_id'])

    # Generate query statements
    llm = OpenAIWrapper(config=SYS_LLM_CONFIG)
    # Prepare keyword arguments for query generation
    gq_kw = {
        'title': metadata['title'],
        'summary': metadata['summary'],
        'study_design': metadata['overall_design'],
        'condition': sample['condition_name'],
        'comparison': sample['comparison_name'],
        'condition_cate': sample['condition_category'],
        'effect_size': sample['effect_size'],
        'effect_direction': 'pro-aging' if sample['direction'] > 0 else 'anti-aging',
        'effect_pvalue': sample['pvalue'],
    }

    # Save sample information including metadata
    sample.update(metadata)
    with open(output_dir.joinpath('sample.json'), 'w') as f:
        json.dump(sample, f, indent=4)

    # Use LLM to generate query statements
    generate_query_prompt = GENERATE_QUERY.format(**gq_kw)
    queries: ChatResponse = await llm.chat_once(generate_query_prompt)
    cost += queries.cost  # Accumulate cost
    queries = parse_llm_json_output(queries.content)
    # Save the generated query statements
    with open(output_dir.joinpath('queries.json'), 'w') as f:
        json.dump(queries, f, indent=4)

    # Retrieve related documents
    query = QueryEngine()
    all_retrieved_documents = ''
    for q in queries:
        try:
            # Use the query engine to retrieve documents
            retrieved_documents, url_map = query.pipeline(q, 'clockbase_paper_202506', enrich=True)
            # Join retrieved documents with separators
            q_document = '\n\n-----\n\n'.join(retrieved_documents)
            # Replace image tags with markdown image links
            for k, v in url_map.items():
                q_document = q_document.replace(f'<image>{k}</image>', f'![]({v})')
        except Exception as e:
            logger.error(f'Error retrieving documents for query {q}: {e}')
            q_document = ''
        all_retrieved_documents += q_document
    # Save all retrieved documents
    with open(output_dir.joinpath('retrieved_documents.md'), 'w') as f:
        f.write(all_retrieved_documents)

    # If no documents were retrieved, skip subsequent processing
    if not all_retrieved_documents:
        logger.info(f'No retrieved documents for {sample["gse_id"]}')
        return None, None

    # Identify the relevance of the study to biological clocks
    # Extract intervention source literature from metadata
    intervention_source_literature = '\n\n'.join([f'{c["title"]}: {c["abstract"]}' for c in metadata['citations']]) if metadata['citations'] else ''

    # Truncate oversized retrieval blocks so we stay under Anthropic's 200K-token
    # context ceiling. A conservative char budget of 450000 ≈ 112K tokens for the
    # RAG block, leaving ~80K tokens for the rest of the prompt. Some rows end up
    # with huge RAG enrichment (hundreds of CrossRef-augmented snippets), which
    # otherwise triggers a 400 error that cascades into the event loop.
    MAX_RAG_CHARS = 450_000
    if len(all_retrieved_documents) > MAX_RAG_CHARS:
        trim = all_retrieved_documents[:MAX_RAG_CHARS]
        all_retrieved_documents = trim + "\n\n[... truncated: retrieval exceeded context budget ...]"
        logger.warning(
            f"{sample['gse_id']}-{sample_index}: truncated retrieved_documents from "
            f"{len(trim) + 1000}+ chars to {MAX_RAG_CHARS} chars to fit context window"
        )
    MAX_LIT_CHARS = 80_000
    if len(intervention_source_literature) > MAX_LIT_CHARS:
        intervention_source_literature = (
            intervention_source_literature[:MAX_LIT_CHARS]
            + "\n\n[... truncated: intervention source literature too long ...]"
        )
    # Pull the 9 cross-clock fields out of sample for DNAm (ignored for others)
    extra_cross_clock = {k: sample[k] for k in _DNAM_EXTRA_FIELDS if k in sample}
    # Human-friendly direction label for the prompt (positive → pro-aging)
    effect_direction_label = 'pro-aging' if float(sample['direction']) > 0 else 'anti-aging'
    # Call the relevance identification function
    relevance, identify_relevance_prompt, response_content, ir_cost = await identify_relevance(
        gse_id=sample['gse_id'],
        title=metadata['title'],
        summary=metadata['summary'],
        study_design=metadata['overall_design'],
        condition=sample['condition_name'],
        comparison=sample['comparison_name'],
        condition_cate=sample['condition_category'],
        effect_size=sample['effect_size'],
        effect_direction=effect_direction_label,
        effect_pvalue=sample['pvalue'],
        extreme_value=sample['outlier_status'],
        retrieved_documents=all_retrieved_documents,
        intervention_source_literature=intervention_source_literature,
        modality=CSV_MODALITY,
        extra_cross_clock=extra_cross_clock,
    )
    cost += ir_cost  # Accumulate relevance identification cost
    # If deemed irrelevant, skip subsequent processing
    if not relevance:
        logger.info(f'Irrelevant: {sample["gse_id"]}')
        return None, None
    # Save relevance identification results
    with open(output_dir.joinpath('identified_relevance.json'), 'w') as f:
        json.dump(relevance, f, indent=4)
    with open(output_dir.joinpath('identify_relevance_prompt.txt'), 'w') as f:
        f.write(identify_relevance_prompt)
    with open(output_dir.joinpath('identify_relevance_response.txt'), 'w') as f:
        f.write(response_content)

    # This section generates the academic report
    # paper_prompt = GENERATE_PAPER_2.format(
    #     gse_id=sample['gse_id'], 
    #     title=metadata['title'], 
    #     condition=sample['condition_name'], 
    #     comparison=sample['comparison_name'], 
    #     effect_size=sample['effect_size'], 
    #     effect_direction=sample['direction'], 
    #     effect_pvalue=sample['pvalue']
    # )
    # report, report_cost = await generate_paper(identify_relevance_prompt, response_content, paper_prompt)
    # cost += report_cost
    # with open(output_dir.joinpath('report.md'), 'w') as f:
    #     f.write(report)

    # Output total cost and save to file
    print(f'Cost: {sample["gse_id"]}: {cost}')
    with open(output_dir.joinpath('cost.txt'), 'w') as f:
        f.write(f'Cost: {cost}')

    return cost


def find_papers(ref_prompt: str, store: Neo4jStore):
    """
    Find paper information in the Neo4j database and add DOI and PDF links
    
    Args:
        ref_prompt (str): Prompt text containing file names
        store (Neo4jStore): Neo4j data store instance
        
    Returns:
        str: Prompt text with DOI and PDF links added
    """
    # Extract file names using regex
    file_name_pattern = r'\*\*File Name\*\*: (.*)\n\n'
    file_names = re.findall(file_name_pattern, ref_prompt)
    # Look up paper information for each file name
    for file_name in file_names:
        paper = store.find_nodes("Paper", {"title": file_name})
        if paper:
            doi = paper[0]['doi']
            pdf_link = paper[0]['pdf_link']
            # Add DOI and PDF link information to the original text
            ref_prompt = ref_prompt.replace(f'**File Name**: {file_name}\n\n', f'**File Name**: {file_name}\n**DOI**: {doi}\n**PDF URL**: {pdf_link}\n\n')
    return ref_prompt


def add_doi_url_to_top(output_dir: Path, csv_path: str, top_n: int = 1000):
    """
    Add DOI and PDF link information for the top-N ranked cases
    
    Args:
        output_dir (Path): Output directory path
        csv_path (str): CSV file path
        top_n (int): Number of cases to process, default 1000
    """
    ir_prompt_file = 'identify_relevance_prompt.txt'
    # Read and sort data
    df = pd.read_csv(csv_path)
    df = df.sort_values(by='final_score', ascending=False)
    df = df[:top_n]
    # Use Neo4j store to find paper information
    with Neo4jStore() as store:
        for i, row in df.iterrows():
            case_id = row['case_id']
            # Read the relevance identification prompt file
            ref_prompt = output_dir.joinpath(case_id, ir_prompt_file).read_text()
            # Add paper DOI and PDF links
            ref_prompt = find_papers(ref_prompt, store)
            # Save the updated prompt file
            with open(output_dir.joinpath(case_id, ir_prompt_file), 'w') as f:
                f.write(ref_prompt)

async def process_one(output_dir: Path, case_id: str):
    """
    Process a single case to generate an academic report
    
    This function reads existing relevance analysis results and generates
    a complete academic report based on those results.
    
    Args:
        output_dir (Path): Output directory path
        case_id (str): Case ID
    """
    # Define file name constants
    ir_prompt_file = 'identify_relevance_prompt.txt'
    ir_response_file = 'identify_relevance_response.txt'
    sample_file = 'sample.json'
    cost_file = 'cost.txt'
    
    # Read sample data
    with open(output_dir.joinpath(case_id, sample_file), 'r') as f:
        sample = json.load(f)
    # Read relevance identification prompt
    with open(output_dir.joinpath(case_id, ir_prompt_file), 'r') as f:
        identify_relevance_prompt = f.read()
    # Read relevance identification response
    with open(output_dir.joinpath(case_id, ir_response_file), 'r') as f:
        response_content = f.read()
    # Read existing cost
    with open(output_dir.joinpath(case_id, cost_file), 'r') as f:
        cost = float(f.read().split(':')[-1].strip())
    
    # Generate paper prompt template
    paper_prompt = GENERATE_PAPER_2.format(
        gse_id=sample['gse_id'], 
        title=sample['title'], 
        condition=sample['condition_name'], 
        comparison=sample['comparison_name'], 
        effect_size=sample['effect_size'], 
        effect_direction='pro-aging' if sample['direction'] > 0 else 'anti-aging', 
        effect_pvalue=sample['pvalue']
    )
    # Generate the academic report
    report, report_cost = await generate_paper(identify_relevance_prompt, response_content, paper_prompt)
    cost += report_cost
    # Save the report
    with open(output_dir.joinpath(case_id, 'report.md'), 'w') as f:
        f.write(report)
    # Update the total cost
    with open(output_dir.joinpath(case_id, cost_file), 'w') as f:
        f.write(f'Cost: {cost + report_cost}')

async def process_top(output_dir: Path, csv_path: str, top_n: int = 1000, workdir: Path = None):
    """
    Batch process the top-N ranked cases to generate academic reports
    
    This function sorts by final_score, processes the top-N cases,
    and generates an academic report for each. Supports resumption
    from interruption by skipping already successfully processed cases.
    
    Args:
        output_dir (Path): Output directory path
        csv_path (str): CSV file path
        top_n (int): Number of cases to process, default 1000
        workdir (Path): Working directory for status files, defaults to project root
    """
    workdir = workdir or Path(__file__).parent.parent
    # Read and sort data, select the top N cases
    df = pd.read_csv(csv_path)
    df = df.sort_values(by='final_score', ascending=False)
    df = df[:top_n]
    
    # Load or create status file for resumption support
    status_file = workdir / 'report_status.json'
    if status_file.exists():
        with open(status_file, 'r') as f:
            status = json.load(f)
    else:
        status = {}
    
    # Process cases one by one
    for i, row in df.iterrows():
        case_id = row['case_id']
        # Skip already successfully processed cases
        if status.get(case_id, {}).get('status') == 'success':
            continue
        try:
            # Process a single case
            await process_one(output_dir, case_id)
            # Update status to success
            status[case_id] = {'status': 'success'}
            with open(status_file, 'w') as f:
                json.dump(status, f, indent=4)
            logger.info(f'Processed: {case_id} -- {i}/{top_n}')
        except Exception as e:
            # Record error status on failure
            logger.error(f'Error processing {case_id}: {e}')
            status[case_id] = {'status': 'failed'}
            with open(status_file, 'w') as f:
                json.dump(status, f, indent=4)

def load_data(df_path: str):
    """
    Load and preprocess the main dataset
    
    This function loads data from the specified path, filters records by
    specific condition categories, separately processes pro-aging and
    anti-aging effects, and deduplicates.
    
    Returns:
        list: List of processed records
        
    Note: This is the second load_data function, with different functionality from the first
    """
    df = pd.read_csv(df_path)
    # Define condition categories of interest
    category_set = ['Drug', 'Disease', 'Genetic', 'Environment']
    
    # Separate pro-aging and anti-aging effect data
    df_pro = df.query('direction > 0.5')  # Pro-aging effects
    df_anti = df.query('direction < -0.5')  # Anti-aging effects
    
    # Filter pro-aging data containing specified condition categories
    max_pro = df_pro[df_pro['condition_category'].str.contains('|'.join(category_set), case=False, na=False)]
    # Sort by GSE_ID and p-value, keep the record with the smallest p-value per GSE_ID
    sorted_pro = max_pro.sort_values(by=['gse_id', 'pvalue'], ascending=[True, True])
    max_pro = sorted_pro.drop_duplicates(subset=['gse_id'], keep='first')
    
    # Filter anti-aging data containing specified condition categories
    max_anti = df_anti[df_anti['condition_category'].str.contains('|'.join(category_set), case=False, na=False)]
    # Sort by GSE_ID and p-value, keep the record with the smallest p-value per GSE_ID
    sorted_anti = max_anti.sort_values(by=['gse_id', 'pvalue'], ascending=[True, True])
    max_anti = sorted_anti.drop_duplicates(subset=['gse_id'], keep='first')
    
    # Merge pro-aging and anti-aging data
    records = pd.concat([max_pro, max_anti])
    records = records.to_dict(orient='records')
    return records


async def producer(records: list, queue: asyncio.Queue):
    """
    Producer coroutine that puts records into the queue
    
    Args:
        records (list): List of records to process
        queue (asyncio.Queue): Async queue
    """
    for i, record in enumerate(records):
        await queue.put((i, record))


async def consumer(queue: asyncio.Queue, status: dict, status_file: Path, output_base: Path = None):
    """
    Consumer coroutine that takes records from the queue and processes them
    
    This function continuously retrieves records from the queue, calls the
    pipeline function to process them, and updates processing status.
    Supports resumption from interruption and error handling.
    
    Args:
        queue (asyncio.Queue): Async queue
        status (dict): Status dictionary
        status_file (Path): Status file path
        output_base (Path): Output base directory
    """
    while True:
        # Get a record from the queue
        i, record = await queue.get()
        await asyncio.sleep(0.05)  # Brief delay to avoid excessive concurrency
        
        # Skip already successfully processed records
        if status.get(f'{record["gse_id"]}-{i}', {}).get('status') == 'success':
            queue.task_done()
            continue
            
        try:
            # Call the processing pipeline
            cost = await pipeline(i, record, output_base=output_base)
            # Update success status
            status[f'{record["gse_id"]}-{i}'] = {'cost': cost, 'status': 'success'}
            with open(status_file, 'w') as f:
                json.dump(status, f, indent=4)
        except Exception as e:
            # Record processing failure
            logger.error(f'Error processing {i}th sample: {e}')
            status[f'{record["gse_id"]}-{i}'] = {'cost': None, 'status': 'failed'}
            with open(status_file, 'w') as f:
                json.dump(status, f, indent=4)
        queue.task_done()


async def main(num_semaphore: int = 8, df_path: str = None, workdir: Path = None):
    """
    Main async function that coordinates the entire processing flow
    
    This function sets up a producer-consumer pattern, uses semaphores
    to control concurrency, enabling efficient batch data processing.
    
    Args:
        num_semaphore (int): Concurrency control semaphore count, default 8
        df_path (str): Data file path
        workdir (Path): Working directory for intermediate artifacts, defaults to project root
    """
    workdir = workdir or Path(__file__).parent.parent
    # Load or create status file
    status_file = workdir / 'status.json'
    if status_file.exists():
        with open(status_file, 'r') as f:
            status = json.load(f)
    else:
        status = {}
    
    # Load records to process
    records = load_data(df_path)

    output_base = workdir / 'output'
    queue = asyncio.Queue()
    producer_task = asyncio.create_task(producer(records, queue))
    # Spawn `num_semaphore` consumers so they process the queue concurrently.
    # The original code created only 1 consumer, serialising everything.
    consumer_tasks = [
        asyncio.create_task(consumer(queue, status, status_file, output_base=output_base))
        for _ in range(max(1, num_semaphore))
    ]
    logger.info(f"Spawned {len(consumer_tasks)} consumer tasks (num_semaphore={num_semaphore})")
    await producer_task
    await queue.join()
    for ct in consumer_tasks:
        ct.cancel()

def select_cases(csv_path: str, from_dir: str, to_dir: str, top_n: int = 1000):
    """
    Select the top-N ranked cases and copy them to the target directory
    
    Sorts by final_score, selects the top-N cases, and copies them from
    the source directory to the target directory. Supports directory overwriting.
    
    Args:
        csv_path (str): CSV file path
        from_dir (str): Source directory path
        to_dir (str): Target directory path
        top_n (int): Number of cases to select, default 1000
    """
    # Read and sort data
    df = pd.read_csv(csv_path)
    df = df.sort_values(by='final_score', ascending=False)
    df = df[:top_n]
    
    # Ensure target directory exists
    Path(to_dir).mkdir(parents=True, exist_ok=True)
    
    # Copy case directories one by one
    for i, row in df.iterrows():
        case_id = row['case_id']
        src_path = Path(from_dir).joinpath(case_id)
        dst_path = Path(to_dir).joinpath(case_id)
        
        # If the target directory already exists, remove it first
        if dst_path.exists():
            shutil.rmtree(dst_path)
        
        # Copy the entire directory
        shutil.copytree(src_path, dst_path)


async def whole_pipe_one(case_josn_path: str, multi_turns: int = 1, num_semaphore: int = 8):
    """
    Complete data processing pipeline
    
    Full processing flow for a single case:
    1. Data preprocessing and field mapping
    2. Create output directory
    3. Fetch GEO metadata
    4. Generate query statements
    5. Retrieve related documents
    6. Analyze relevance
    7. Save results and costs
    """

    async def process(case: dict, i: int):
        output_base = Path(__file__).parent.parent / 'output'
        await pipeline(i, case, output_base=output_base)
        case_id = f'{case["gse_id"]}-{i}'
        await process_one(output_base, case_id)

    with open(case_josn_path, 'r') as f:
        case = json.load(f)
    async with asyncio.Semaphore(num_semaphore):
        tasks = []
        for i in range(multi_turns):
            tasks.append(process(case, i))
        await asyncio.gather(*tasks)


if __name__ == "__main__":
    """
    Main program entry point
    
    Provides three running modes:
    1. main(16) - Async batch processing mode with 16 concurrent workers
    2. select_cases() - Case selection and copy mode, copies top-N ranked cases to target directory
    3. process_top() - Report generation mode for top-N cases (currently active), generates reports for top-N cases and saves to target directory
    """

    workdir = Path(__file__).parent.parent

    # --- Full scoring run (scoring only; no paper generation) -------------
    import os
    modality = os.environ.get("CSV_MODALITY", "hs_dnam")
    csv_name = f"main_{modality}.csv"
    asyncio.run(main(
        num_semaphore=int(os.environ.get("NUM_WORKERS", "4")),
        df_path=str(Path(__file__).parent / 'data' / csv_name),
        workdir=workdir,
    ))
    # ----------------------------------------------------------------------

    # Mode 1: Async batch processing (full run)
    # asyncio.run(main(16, df_path='/path/to/for_scoring_hs_rna.csv'))

    # Mode 2: Case selection and copy
    # select_cases(
    #     csv_path=str(workdir / 'output_df_fdr_1000.csv'),
    #     from_dir=str(workdir / 'output'),
    #     to_dir=str(workdir / 'report_output'),
    #     top_n=1000
    # )

    # Mode 3: Report generation for top-N cases
    # asyncio.run(process_top(workdir / 'report_output', str(workdir / 'output_df_fdr_1000.csv'), top_n=1000))

    # Mode 4: Single case report generation, generates report for a single case and saves to target directory
    # asyncio.run(whole_pipe_one('data/GSE122080.json', multi_turns=10, num_semaphore=8))
