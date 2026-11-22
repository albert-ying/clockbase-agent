import sys
from pathlib import Path

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_system import paths

import json
import pickle
import random
from pprint import pprint

import pandas as pd

# Load the traces
with open(
    paths.TRACES / "successful_traces_mouse.pkl",
    "rb",
) as f:
    traces = pickle.load(f)

print(f"Total number of traces: {len(traces)}")

# Set random seed for reproducibility
random.seed(42)

# Randomly select 2 traces
all_keys = list(traces.keys())
selected_keys = random.sample(all_keys, 2)

print(f"\n{'='*80}")
print(f"RANDOMLY SELECTED TRACES FOR VALIDATION")
print(f"{'='*80}")

for idx, key in enumerate(selected_keys, 1):
    print(f"\n\n{'#'*80}")
    print(f"TRACE {idx}: {key}")
    print(f"{'#'*80}\n")

    trace = traces[key]

    # Basic info
    print(f"Study ID: {trace['study_id']}")
    print(f"Status: {trace['status']}")
    print(f"Execution Time: {trace['execution_time']:.2f} seconds")

    # Metadata
    print(f"\n{'='*70}")
    print("STUDY METADATA")
    print(f"{'='*70}")
    print(f"Title: {trace['metadata']['title']}")
    print(f"\nSummary: {trace['metadata']['summary'][:500]}...")
    print(f"\nOverall Design: {trace['metadata']['overall_design'][:500]}...")

    # Result analysis
    print(f"\n{'='*70}")
    print("ANALYSIS TRACE")
    print(f"{'='*70}")

    result_tuple = trace["result"]
    messages = result_tuple[0] if isinstance(result_tuple, tuple) else result_tuple

    # Parse the messages/conversation
    if isinstance(messages, list):
        print(f"\nTotal messages in trace: {len(messages)}")

        # Show first few and last few messages
        print(f"\n{'-'*70}")
        print("INITIAL PROMPT:")
        print(f"{'-'*70}")
        if len(messages) > 0:
            print(messages[0][:1000])
            if len(messages[0]) > 1000:
                print("...[truncated]...")

        # Look for tool calls and analysis steps
        print(f"\n{'-'*70}")
        print("ANALYSIS STEPS (Tool Calls and Key Actions):")
        print(f"{'-'*70}")

        tool_calls = []
        ai_responses = []

        for i, msg in enumerate(messages):
            msg_str = str(msg)

            # Look for code execution or tool calls
            if "Tool Message" in msg_str or "tool_call" in msg_str.lower():
                # Extract relevant parts
                snippet = msg_str[:500]
                tool_calls.append((i, snippet))

            # Look for AI responses with analysis
            if "AI Message" in msg_str or "assistant" in msg_str.lower():
                if len(msg_str) > 100:  # Only substantial responses
                    ai_responses.append((i, msg_str))

        print(f"\nFound {len(tool_calls)} tool-related messages")
        print(f"Found {len(ai_responses)} AI response messages")

        # Show some key tool calls
        if tool_calls:
            print(f"\n{'-'*70}")
            print("Sample Tool Calls:")
            for i, (msg_idx, snippet) in enumerate(tool_calls[:5]):
                print(f"\nMessage {msg_idx}:")
                print(snippet)
                print("...")

        # Show final analysis/conclusion
        print(f"\n{'-'*70}")
        print("FINAL MESSAGES (Last 3):")
        print(f"{'-'*70}")
        for msg in messages[-3:]:
            msg_str = str(msg)
            print(f"\n{msg_str[:800]}")
            if len(msg_str) > 800:
                print("...[truncated]...")

    # If result has additional structure
    if isinstance(result_tuple, tuple) and len(result_tuple) > 1:
        print(f"\n{'-'*70}")
        print("ADDITIONAL RESULT DATA:")
        print(f"{'-'*70}")
        print(f"Result tuple length: {len(result_tuple)}")
        for i, item in enumerate(result_tuple[1:], 1):
            print(f"\nItem {i} type: {type(item)}")
            if isinstance(item, dict):
                print(f"Keys: {item.keys()}")
            elif isinstance(item, pd.DataFrame):
                print(f"DataFrame shape: {item.shape}")
                print(f"Columns: {item.columns.tolist()}")

    print(f"\n{'='*80}\n")

print("\n" + "=" * 80)
print("VALIDATION SUMMARY")
print("=" * 80)
print("\nTo validate these traces, check:")
print("1. Does the analysis follow the stated workflow?")
print("2. Are appropriate data files loaded?")
print("3. Are statistical methods correctly applied?")
print("4. Are the conclusions supported by the data?")
print("5. Are age-related comparisons made correctly?")
