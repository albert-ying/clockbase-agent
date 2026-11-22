import re
import time
from typing import Dict, List, Optional

import requests
from bs4 import BeautifulSoup


def get_pubmed_details(pmid: str) -> Dict[str, str]:
    """
    Fetch paper details from PubMed using PMID.

    Args:
        pmid (str): PubMed ID

    Returns:
        dict: Dictionary containing title, doi, authors, journal, etc.
    """
    try:
        # Use PubMed eutils API
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pmid}&retmode=xml"

        response = requests.get(url)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, "xml")

        article = soup.find("PubmedArticle")
        if not article:
            return {"title": "", "doi": "", "journal": "", "authors": ""}

        # Extract title
        title_elem = article.find("ArticleTitle")
        title = title_elem.text.strip() if title_elem else ""

        # Extract DOI
        doi = ""
        article_ids = article.find_all("ArticleId")
        for article_id in article_ids:
            if article_id.get("IdType") == "doi":
                doi = article_id.text.strip()
                break

        # Extract journal
        journal_elem = article.find("Title")
        journal = journal_elem.text.strip() if journal_elem else ""

        # Extract authors (first few)
        authors = []
        author_list = article.find("AuthorList")
        if author_list:
            for author in author_list.find_all("Author")[:3]:  # First 3 authors
                last_name = author.find("LastName")
                first_name = author.find("ForeName")
                if last_name:
                    author_name = last_name.text.strip()
                    if first_name:
                        author_name += f" {first_name.text.strip()}"
                    authors.append(author_name)

        authors_str = ", ".join(authors)
        if len(author_list.find_all("Author")) > 3:
            authors_str += " et al."

        return {"title": title, "doi": doi, "journal": journal, "authors": authors_str}

    except Exception as e:
        print(f"Error fetching PubMed details for PMID {pmid}: {e}")
        return {"title": "", "doi": "", "journal": "", "authors": ""}


def get_geo_metadata(gse_id: str) -> Dict[str, str]:
    """
    Fetch metadata for a GEO dataset given its GSE ID.

    Args:
        gse_id (str): GEO Series accession number (e.g. 'GSE122080')

    Returns:
        dict: Dictionary containing:
            - title: Study title
            - summary: Study description/summary
            - overall_design: Study design information
            - citations: List of citations with DOI and PMID information

    Raises:
        ValueError: If GSE ID is invalid or data cannot be retrieved
        ConnectionError: If there are network issues
    """

    # Validate GSE ID format
    if not gse_id.startswith("GSE") or not gse_id[3:].isdigit():
        raise ValueError(
            "Invalid GSE ID format. Must start with 'GSE' followed by numbers"
        )

    # Use GEO's direct API
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}&targ=self&form=xml&view=quick"

    try:
        # Fetch the record
        response = requests.get(url)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, "xml")

        # Find the Series element
        series = soup.find("Series")
        if not series:
            raise ValueError(f"No metadata found for GSE ID: {gse_id}")

        # Extract metadata with error handling
        try:
            title = series.find("Title").text.strip()
        except (AttributeError, TypeError):
            title = ""

        try:
            summary = series.find("Summary").text.strip()
        except (AttributeError, TypeError):
            summary = ""

        try:
            design = series.find("Overall-Design").text.strip()
        except (AttributeError, TypeError):
            design = ""

            # Extract citations from Pubmed-ID elements
        citations = []
        try:
            pubmed_elements = series.find_all("Pubmed-ID")

            for pubmed_element in pubmed_elements:
                if pubmed_element:
                    pmid = pubmed_element.text.strip()

                    # For now, just store the PMID
                    # In a more complete implementation, you could fetch paper details from PubMed API
                    citations.append(
                        {
                            "pmid": pmid,
                            "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                            "doi": "",  # Would need to fetch from PubMed API
                            "title": "",  # Would need to fetch from PubMed API
                            "full_citation": f"PMID: {pmid}",
                        }
                    )
        except (AttributeError, TypeError):
            citations = []

        # Validate that we got at least some data
        if not any([title, summary, design]):
            raise ValueError(f"No metadata content found for GSE ID: {gse_id}")

        return {
            "title": title,
            "summary": summary,
            "overall_design": design,
            "citations": citations,
        }

    except requests.exceptions.RequestException as e:
        raise ConnectionError(f"Failed to fetch data from GEO: {str(e)}")
    except Exception as e:
        raise ValueError(f"Failed to parse GEO response: {str(e)}")


# Example usage:
if __name__ == "__main__":
    try:
        gse_id = "GSE122088"
        metadata = get_geo_metadata(gse_id)
        print("GSE ID:", gse_id)
        print("URL:", f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}")
        print("Title:", metadata["title"])
        print("\nSummary:", metadata["summary"])
        print("\nStudy Design:", metadata["overall_design"])

        print("\nCitations:")
        if metadata["citations"]:
            for i, citation in enumerate(metadata["citations"], 1):
                print(f"\n  Paper {i}:")
                print(f"    PMID: {citation['pmid']}")
                print(f"    PubMed URL: {citation['pubmed_url']}")
                if citation["doi"]:
                    print(f"    DOI: {citation['doi']}")
                if citation["title"]:
                    print(f"    Title: {citation['title']}")
        else:
            print("  No citations available")
    except Exception as e:
        print(f"Error: {str(e)}")
