import re
import traceback
from pathlib import Path


def correct_markdown_references(markdown_text: str) -> str:
    """
    Corrects formatting errors in a Markdown scientific report.

    This function cleans reference list entries, removes invalid references,
    and updates the corresponding in-text citations to be sequential and accurate.

    Args:
        markdown_text: A string containing the full Markdown report.

    Returns:
        A string with the corrected Markdown text.
    """
    # Define patterns for unwanted prefixes and entirely invalid lines.
    # These are case-insensitive.
    PREFIXES_TO_CLEAN = r'^(OPEN ACCESS|Review|REVIEW|Article|ARTICLE|BASIC SCIENCES|Symposium\s+\w+:)\s*[-.]?\s*'
    INVALID_LINE_PATTERNS = [
        r'^\s*OPEN ACCESS\s*$',
        r'^\s*REVIEW\.\s*DOI:.*$',
        r'^\s*Article\s*$',
        r'^\s*ARTICLE\s*$',
        r'^\s*BASIC SCIENCES\s*$',
        r'^\s*Symposium\s+\w+:\s*$',
        r'^\s*DOI:\s*Not\s+provided\.?\s*$',
        r'^\s*DOI:\s*N/A\.?\s*$'
    ]

    # Split the document into the main body and the references section.
    # We use a regex split to handle variations in whitespace and capitalization.
    parts = re.split(r'\n##\s+References\s*\n', markdown_text, flags=re.IGNORECASE)
    if len(parts) < 2:
        # No reference section found, return the original text.
        return markdown_text

    body, references_section = parts[0], parts[1]
    
    original_lines = references_section.strip().split('\n')
    
    valid_references = []
    invalid_indices = set()
    seen_texts = {}  # Track unique references and their original indices

    # 1. First pass: Parse references, clean them, and identify invalid ones.
    for line in original_lines:
        line = line.strip()
        if not line:
            continue

        # Match reference lines like "1. Some text" or "[2] Some text"
        match = re.match(r'^\s*\[?(\d+)\]?\.?\s*(.*)', line)
        if not match:
            continue

        index = int(match.group(1))
        text = match.group(2).strip()

        # Check if the line matches any of the invalid patterns.
        is_invalid = any(re.match(pattern, text, re.IGNORECASE) for pattern in INVALID_LINE_PATTERNS)

        if is_invalid:
            invalid_indices.add(index)
        else:
            # If valid, clean any unwanted prefixes.
            cleaned_text = re.sub(PREFIXES_TO_CLEAN, '', text, flags=re.IGNORECASE).strip()
            
            # Strip "DOI: Not provided" or similar patterns
            cleaned_text = re.sub(r'\.\s*DOI:\s*Not\s+provided\.?$', '', cleaned_text, flags=re.IGNORECASE).strip()
            cleaned_text = re.sub(r'\.\s*DOI:\s*N/A\.?$', '', cleaned_text, flags=re.IGNORECASE).strip()
            
            # Check if the reference is too short to be meaningful (less than 10 characters after cleaning)
            # This filters out entries that are just "DOI: Not provided" or similar
            if len(cleaned_text) < 20:
                invalid_indices.add(index)
                continue
            
            # Check for duplicates
            if cleaned_text in seen_texts:
                # This is a duplicate - map this index to the first occurrence
                seen_texts[cleaned_text].append(index)
            else:
                # First occurrence of this reference
                seen_texts[cleaned_text] = [index]
                valid_references.append({
                    'original_indices': [index],  # Store all indices that map to this reference
                    'text': cleaned_text
                })
    
    # Update valid_references with all duplicate indices
    for ref in valid_references:
        ref_text = ref['text']
        ref['original_indices'] = seen_texts[ref_text]

    # 2. Process the main body of the text.
    cleaned_body = body
    
    # Find all citations in the body text
    all_citations_in_body = set()
    # Pattern to match both single citations [1] and multiple citations [1,2,3]
    citation_pattern = r'\[\s*(\d+(?:\s*,\s*\d+)*)\s*\]'
    for match in re.finditer(citation_pattern, cleaned_body):
        # Split by comma to handle multiple citations in one bracket
        citation_nums = match.group(1).replace(' ', '').split(',')
        for num in citation_nums:
            all_citations_in_body.add(int(num))
    
    # Determine which citations are dangling (not in valid references)
    valid_reference_indices = set()
    for ref in valid_references:
        valid_reference_indices.update(ref['original_indices'])
    
    dangling_citations = all_citations_in_body - valid_reference_indices - invalid_indices

    # Remove citations pointing to invalid references AND dangling citations.
    citations_to_remove = invalid_indices.union(dangling_citations)
    
    # Handle multiple citations in brackets like [1,2,3]
    def replace_citations(match):
        full_match = match.group(0)
        citations = match.group(1).replace(' ', '').split(',')
        remaining_citations = []
        for cit in citations:
            if int(cit) not in citations_to_remove:
                remaining_citations.append(cit)
        
        if not remaining_citations:
            # If no citations remain, return empty string
            return ''
        elif remaining_citations != citations:
            # If some citations remain, return updated bracket
            return '[' + ','.join(remaining_citations) + ']'
        else:
            # If all citations remain, return original
            return full_match
    
    # Replace citations, handling both single and multiple citations
    cleaned_body = re.sub(r'\[\s*(\d+(?:\s*,\s*\d+)*)\s*\]', replace_citations, cleaned_body)
    
    # Clean up any leftover spaces before punctuation
    cleaned_body = re.sub(r'\s+([.,;!?])', r'\1', cleaned_body)

    # Create a mapping from the old reference indices to the new sequential index.
    old_to_new_map = {}
    for i, ref in enumerate(valid_references):
        new_index = i + 1
        for old_index in ref['original_indices']:
            old_to_new_map[old_index] = new_index

    # Update the remaining valid citations to their new numbers.
    def renumber_citations(match):
        full_match = match.group(0)
        citations = match.group(1).replace(' ', '').split(',')
        renumbered = []
        for cit in citations:
            old_num = int(cit)
            if old_num in old_to_new_map:
                renumbered.append(str(old_to_new_map[old_num]))
            else:
                renumbered.append(cit)
        return '[' + ','.join(renumbered) + ']'
    
    # Renumber citations if needed
    if old_to_new_map:
        cleaned_body = re.sub(r'\[\s*(\d+(?:\s*,\s*\d+)*)\s*\]', renumber_citations, cleaned_body)
            
    # 3. Reconstruct the final document.
    final_references_list = [f"{i + 1}. {ref['text']}" for i, ref in enumerate(valid_references)]
    final_references_section = "\n".join(final_references_list)

    # Combine the cleaned body with the new, corrected reference list.
    final_text = cleaned_body.strip() + "\n\n## References\n\n" + final_references_section

    return final_text

# --- Example Usage ---
# Here is the sample Markdown text provided in the prompt.
def example():
    sample_markdown = """
# Smchd1 Knockout Accelerates Epigenetic Aging Through Disruption of Chromatin Organization

## Introduction

Structural maintenance of chromosomes flexible hinge domain containing 1 (Smchd1) is an epigenetic regulator that plays a critical role in maintaining heterochromatin organization and gene silencing. This chromatin-associated protein belongs to the SMC family and functions as a key component in X-chromosome inactivation and genomic imprinting processes. The biological significance of Smchd1 extends beyond developmental regulation, as emerging evidence suggests its involvement in cellular aging mechanisms through epigenetic control pathways. This report synthesizes existing scientific evidence regarding Smchd1's potential impact on aging-related processes and longevity, with particular focus on its role in chromatin stability and epigenetic maintenance.

## Results

Analysis of the GSE44665 dataset revealed that Smchd1 knockout (Smchd1MommeD1/MommeD1 vs Smchd1+/+) in E9.5 mouse embryos significantly accelerated the mouse mortality aging clock with an effect size of 3.10 (p = 0.021), indicating a pro-aging effect. The RNA-seq analysis demonstrated substantial transcriptomic alterations in homozygous mutant embryos compared to wild-type controls.

The mechanistic basis for Smchd1's role in aging centers on its function in maintaining heterochromatin organization. Heterochromatin domains are structurally inaccessible and transcriptionally inactive regions that become progressively lost during aging, resulting in de-repression of normally silenced genes [1]. This loss of heterochromatin represents a fundamental aging hallmark, with specific proteins like heterochromatin protein 1 (HP1) and histone modifications such as H3K9me3 decreasing alongside heterochromatin deterioration [1].

Epigenetic alterations constitute a primary mechanism through which aging manifests at the cellular level. DNA methylation patterns, histone modifications, and chromatin remodeling undergo systematic changes during the aging process [2]. These modifications profoundly influence transcriptional accessibility and cellular state determination, with methylation at CpG islands typically leading to gene repression while gene body methylation associates with active transcription [3].

The relationship between chromatin organization and aging extends to multiple tissue types and cellular contexts. Age-related epigenetic changes affect posttranslational modifications of histones, including lysine methylations and cytosine methylations at CpG-rich DNA regions [2]. Repressive histone marks involving lysine trimethylations, particularly H3K9me3 and H3K27me3, show general reduction during physiological aging, contributing to heterochromatin loss [4].

Experimental evidence demonstrates that maintenance of heterochromatin organization can effectively counter aging processes. Studies using human mesenchymal stem cells revealed that proteins involved in heterochromatin stability, including nuclear envelope protein Lamin B1 and heterochromatin-associated proteins KAP1 and HP1γ, can reverse premature senescent phenotypes when overexpressed [1]. This suggests that chromatin organizational integrity directly impacts cellular aging trajectories.

The broader implications of epigenetic dysregulation in aging encompass multiple organ systems and disease processes. Age-related changes in DNA methylation and histone modifications contribute to the development of neurodegenerative diseases, cardiovascular disorders, and metabolic dysfunction [2]. The progressive nature of these epigenetic alterations correlates with functional decline across tissues, establishing chromatin organization as a fundamental determinant of organismal aging.

## Discussion

The collective evidence indicates that Smchd1 knockout produces a clear pro-aging effect through disruption of fundamental chromatin organizational mechanisms. The significant acceleration of the epigenetic aging clock observed in Smchd1-deficient embryos aligns with established understanding of heterochromatin loss as a primary aging hallmark. This finding demonstrates remarkable consistency with broader literature documenting the critical role of chromatin stability in maintaining cellular youth and preventing age-related dysfunction.

The plausibility of Smchd1's anti-aging function stems from its established role in heterochromatin maintenance and gene silencing. Since heterochromatin loss represents one of the most well-characterized features of cellular aging, proteins that preserve chromatin organization would logically exhibit anti-aging properties. The observed pro-aging effect of Smchd1 knockout therefore represents the expected consequence of disrupting this protective mechanism.

Several limitations warrant consideration in interpreting these findings. The analysis was conducted in embryonic tissue, which may not fully reflect aging dynamics in mature organisms. Additionally, the specific mechanisms through which Smchd1 loss accelerates epigenetic aging require further elucidation, particularly regarding tissue-specific effects and temporal progression of chromatin deterioration.

Future research should focus on investigating Smchd1's role in adult tissue aging and determining whether therapeutic interventions targeting chromatin organization can reverse age-related epigenetic changes. Understanding the precise molecular pathways through which Smchd1 maintains chromatin stability could reveal novel targets for anti-aging interventions based on epigenetic restoration.

## References

[1] Stabilizing heterochromatin by DGCR8 alleviates senescence and osteoarthritis

[2] Review Epigenetic Mechanisms of Aging and Aging-Associated Diseases

[3] DNA Methylation in Epidermal Differentiation, Aging, and Cancer

[4] OPEN ACCESS
"""

    # Run the correction function
    corrected_text = correct_markdown_references(sample_markdown)

    # Print the result to see the changes
    print("----------- ORIGINAL TEXT -----------")
    print(sample_markdown.strip())
    print("\n" + "="*40 + "\n")
    print("----------- CORRECTED TEXT -----------")
    print(corrected_text)


def main(raw_dir: str, parsed_dir: str):
    for gse in Path(raw_dir).glob('*'):
        gse_id = gse.name
        report_path = gse.joinpath('report.md')
        if report_path.exists():
            with open(str(report_path), 'r') as f:
                report = f.read()
            try:
                parsed_report = correct_markdown_references(report)
                save_path = parsed_dir.joinpath(f'{gse_id}/parsed_report.md')
                save_path.parent.mkdir(parents=True, exist_ok=True)
                with open(str(save_path), 'w') as f:
                    f.write(parsed_report)
            except Exception as e:
                traceback.print_exc()
                raise e
                message = f'Error parsing {gse_id}: {e}\n\n'
                print(message)
                with open('error_parsed.log', 'a') as f:
                    f.write(message)


if __name__ == '__main__':
    ...

