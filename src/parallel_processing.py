import os, sys, json, time, signal, pickle, warnings, requests, datetime, pandas as pd, concurrent.futures
from typing import Dict, List, Any, Optional
warnings.filterwarnings("ignore")

# Global variables for interruption handling
interrupted = False
all_results = {}

def signal_handler(sig, frame):
    """Handle interruption signals like Ctrl+C"""
    global interrupted
    print("\nInterruption detected. Saving current results before exiting...")
    interrupted = True
    # Note: save_results will be called by the main script
    print("Partial results saved. Exiting gracefully.")
    sys.exit(0)

def setup_signal_handlers(save_func, results_dict, timestamp, results_dir):
    """Register signal handlers for graceful interruption"""
    def custom_signal_handler(sig, frame):
        global interrupted
        print("\nInterruption detected. Saving current results before exiting...")
        interrupted = True
        save_func(results_dict, timestamp, results_dir, final=False)
        print("Partial results saved. Exiting gracefully.")
        sys.exit(0)
    
    signal.signal(signal.SIGINT, custom_signal_handler)
    signal.signal(signal.SIGTERM, custom_signal_handler)

def save_results(results: Dict[str, Any], timestamp: str, results_dir: str, final: bool = True):
    """
    Save the current results to disk
    
    Args:
        results (dict): The results to save
        timestamp (str): Timestamp for filename
        results_dir (str): Directory to save results in
        final (bool): Whether this is the final save or an intermediate one
    """
    try:
        # Create results directory if it doesn't exist
        os.makedirs(results_dir, exist_ok=True)
        
        # Determine prefix based on whether this is a final or partial save
        prefix = "final" if final else "partial"
        
        # Save summary results as JSON
        summary_file = os.path.join(results_dir, f"{prefix}_analysis_summary_{timestamp}.json")
        with open(summary_file, 'w') as f:
            # Create a simplified version of results for JSON serialization
            serializable_results = {}
            for study_id, result in results.items():
                serializable_results[study_id] = {
                    "status": result.get("status", "unknown"),
                    "message": result.get("message", "") if result.get("status") != "success" else "",
                    "metadata": result.get("metadata", {}) if result.get("status") == "success" else {},
                    "execution_time": result.get("execution_time", 0)
                }
            
            json.dump(serializable_results, f, indent=2)
        
        if not final:
            print(f"\nUpdated partial results in: {summary_file}")
        else:
            print(f"\nSaved final results to: {summary_file}")
        
        # Save full results (including agent outputs)
        full_results_file = os.path.join(results_dir, f"{prefix}_analysis_full_{timestamp}.pkl")
        with open(full_results_file, 'wb') as f:
            pickle.dump(results, f)
        
        if final:
            print(f"Saved full results to: {full_results_file}")
            
    except Exception as e:
        print(f"Error saving results: {str(e)}")

def get_geo_metadata(gse_id: str, max_retries: int = 3, retry_delay: int = 60) -> Dict[str, str]:
    """
    Fetch metadata for a GEO dataset given its GSE ID.

    Args:
        gse_id (str): GEO Series accession number (e.g. 'GSE122080')
        max_retries (int): Maximum number of retry attempts
        retry_delay (int): Delay in seconds between retries

    Returns:
        dict: Dictionary containing:
            - title: Study title
            - summary: Study description/summary
            - overall_design: Study design information

    Raises:
        ValueError: If GSE ID is invalid or data cannot be retrieved
        ConnectionError: If there are network issues after all retries
    """
    from bs4 import BeautifulSoup

    # Validate GSE ID format
    if not gse_id.startswith("GSE") or not gse_id[3:].isdigit():
        raise ValueError("Invalid GSE ID format. Must start with 'GSE' followed by numbers")
    
    # Use GEO's direct API
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}&targ=self&form=xml&view=quick"

    retry_count = 0
    last_exception = None

    while retry_count < max_retries:
        try:
            # Fetch the record
            print(f"Fetching metadata for {gse_id} (attempt {retry_count + 1}/{max_retries})...")
            response = requests.get(url)
            response.raise_for_status()
            soup = BeautifulSoup(response.content, "xml")

            # Find the Series element
            series = soup.find("Series")
            if not series:
                raise ValueError(f"No metadata found for GSE ID: {gse_id}")

            # Extract metadata with error handling
            try:
                title = series.find("Title").text.strip()  # type: ignore
            except (AttributeError, TypeError):
                title = ""

            try:
                summary = series.find("Summary").text.strip()  # type: ignore
            except (AttributeError, TypeError):
                summary = ""

            try:
                design = series.find("Overall-Design").text.strip()  # type: ignore
            except (AttributeError, TypeError):
                design = ""

            # Validate that we got at least some data
            if not any([title, summary, design]):
                raise ValueError(f"No metadata content found for GSE ID: {gse_id}")

            return {
                "title": title,
                "summary": summary,
                "overall_design": design,
            }

        except requests.exceptions.RequestException as e:
            last_exception = e
            retry_count += 1
            if retry_count < max_retries:
                print(f"  Request failed: {str(e)}. Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
            
        except Exception as e:
            last_exception = e
            retry_count += 1
            if retry_count < max_retries:
                print(f"  Error parsing response: {str(e)}. Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
    
    # If we've exhausted all retries
    if isinstance(last_exception, requests.exceptions.RequestException):
        raise ConnectionError(f"Failed to fetch data from GEO after {max_retries} attempts: {str(last_exception)}")
    else:
        raise ValueError(f"Failed to process GEO response after {max_retries} attempts: {str(last_exception)}")

def should_skip_study(study_id: str, results_dir: str) -> bool:
    """
    Determine if a study should be skipped based on existing results
    
    Args:
        study_id (str): The GSE ID to check
        results_dir (str): Directory where results are stored
        
    Returns:
        bool: True if study should be skipped, False otherwise
    """
    results_file = f"{results_dir}/{study_id}/aging_analysis.csv"
    
    # Check if the results file already exists
    if os.path.exists(results_file):
        print(f"Skipping {study_id}: Results file already exists at {results_file}")
        return True
    
    return False

def process_study(study_id: str, agent, generate_task_intention_func, results_dir: str) -> Dict[str, Any]:
    """
    Process a single study using the agent
    
    Args:
        study_id (str): The GSE ID to process
        agent: The configured agent instance
        generate_task_intention_func: Function to generate task intention
        results_dir (str): Directory where results are stored
        
    Returns:
        dict: Results of the study processing
    """
    start_time = datetime.datetime.now()
    try:
        # Skip if results already exist
        if should_skip_study(study_id, results_dir):
            return {
                "study_id": study_id,
                "status": "skipped",
                "message": "Results file already exists",
                "execution_time": 0
            }
        
        print(f"Processing study {study_id}...")
        study_design = get_geo_metadata(study_id)
        task_intention = generate_task_intention_func(study_id, study_design)
        
        # Run the agent
        result = agent.go(task_intention)
        print(f"  ✓ Agent completed analysis for {study_id}")
        
        end_time = datetime.datetime.now()
        execution_time = (end_time - start_time).total_seconds()
        
        return {
            "study_id": study_id,
            "status": "success",
            "metadata": study_design,
            "result": result,
            "execution_time": execution_time
        }
    except Exception as e:
        end_time = datetime.datetime.now()
        execution_time = (end_time - start_time).total_seconds()
        print(f"  ✗ Error processing {study_id}: {str(e)}")
        return {
            "study_id": study_id,
            "status": "error",
            "message": str(e),
            "execution_time": execution_time
        }

def parallel_process_studies(study_ids: List[str], agent, generate_task_intention_func, results_dir: str, timestamp: str, max_workers: int = 3) -> Dict[str, Any]:
    """
    Process multiple studies in parallel using the agent
    
    Args:
        study_ids (list): List of GSE IDs to analyze
        agent: The configured agent instance
        generate_task_intention_func: Function to generate task intention
        results_dir (str): Directory where results are stored
        timestamp (str): Timestamp for this run
        max_workers (int): Maximum number of parallel workers
        
    Returns:
        dict: Mapping of study IDs to their analysis results
    """
    global all_results, interrupted
    
    # Track the number of completed studies for periodic saving
    completed_count = 0
    save_interval = 5  # Save every 5 completed studies
    last_save_time = time.time()
    save_time_interval = 1800  # 30 minutes in seconds
    
    # Use ThreadPoolExecutor for parallel processing
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all processing tasks
        future_to_study = {
            executor.submit(process_study, study_id, agent, generate_task_intention_func, results_dir): study_id 
            for study_id in study_ids
        }
        
        # Process results as they complete
        for future in concurrent.futures.as_completed(future_to_study):
            study_id = future_to_study[future]
            try:
                result = future.result()
                all_results[study_id] = result
                status = result['status']
                if status == "skipped":
                    print(f"Skipped analysis for {study_id}: {result.get('message', '')}")
                else:
                    print(f"Completed analysis for {study_id}: {status}")
                
                # Increment completed count
                completed_count += 1
                
                # Check if we should save intermediate results (by count or time)
                current_time = time.time()
                if (completed_count % save_interval == 0) or (current_time - last_save_time >= save_time_interval):
                    print(f"\nSaving intermediate results after {completed_count} studies...")
                    save_results(all_results, timestamp, results_dir, final=False)
                    last_save_time = current_time
                
                # Check if we've been interrupted
                if interrupted:
                    print("\nInterruption detected. Stopping further processing.")
                    break
                    
            except Exception as e:
                all_results[study_id] = {
                    "study_id": study_id,
                    "status": "error",
                    "message": f"Exception during analysis: {str(e)}",
                    "execution_time": 0
                }
                print(f"Error processing {study_id}: {str(e)}")
                
                # Still count as completed for saving purposes
                completed_count += 1
    
    return all_results



def print_summary(results: Dict[str, Any]):
    """Print summary of analysis results"""
    print("\nSummary of analysis results:")
    success_count = 0
    skipped_count = 0
    error_count = 0
    
    for study_id, result in results.items():
        status = result.get("status", "unknown")
        if status == "success":
            print(f"✓ {study_id}: Success")
            success_count += 1
        elif status == "skipped":
            print(f"⏭ {study_id}: Skipped - {result.get('message', 'Results already exist')}")
            skipped_count += 1
        else:
            print(f"✗ {study_id}: {status.capitalize()} - {result.get('message', 'Unknown error')}")
            error_count += 1
    
    print(f"\nProcessed {len(results)} studies: {success_count} successful, {skipped_count} skipped, {error_count} failed")
    return success_count, skipped_count, error_count 