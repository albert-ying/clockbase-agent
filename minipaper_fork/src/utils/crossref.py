"""
CrossRef API client for resolving paper references.

Three-level fallback:
  1. Extract DOI from document text → CrossRef DOI lookup (100% accurate)
  2. Use file_name as title → CrossRef bibliographic search
  3. Use document text fragment → CrossRef bibliographic search
  Fallback: return original file_name as title
"""

import re
import time
import requests
from dataclasses import dataclass
from typing import Optional

from utils.logger import logger
from utils.helpers import retry


DOI_PATTERN = re.compile(r'\b(10\.\d{4,}/[^\s,;\]})\"\']+)')
CROSSREF_API = 'https://api.crossref.org/works'
REQUEST_HEADERS = {
    'User-Agent': 'ClockBase/0.1 (https://github.com/clockbase; mailto:clockbase@avinasi.com)',
}
MIN_TITLE_LEN = 15
MIN_CROSSREF_SCORE = 40
RATE_LIMIT_INTERVAL = 0.2


@dataclass
class PaperMetadata:
    title: str = ''
    doi: str = ''
    authors: str = ''
    journal: str = ''
    year: str = ''
    url: str = ''

    @property
    def citation_line(self) -> str:
        """Format as a single-line citation for Markdown references."""
        parts = []
        if self.authors:
            parts.append(self.authors)
        if self.title:
            parts.append(self.title)
        if self.journal:
            j = f'*{self.journal}*'
            if self.year:
                j += f', {self.year}'
            parts.append(j)
        if self.doi:
            parts.append(f'DOI: {self.doi}')
        return '. '.join(parts) if parts else self.title

    @property
    def header_block(self) -> str:
        """Format as a multi-line header block for retrieved documents."""
        lines = [f'**Title**: {self.title}']
        if self.authors:
            lines.append(f'**Authors**: {self.authors}')
        if self.journal:
            j = self.journal
            if self.year:
                j += f', {self.year}'
            lines.append(f'**Journal**: {j}')
        if self.doi:
            lines.append(f'**DOI**: {self.doi}')
        return '\n'.join(lines)


def extract_dois(text: str) -> list[str]:
    """Extract unique DOIs from text via regex."""
    matches = DOI_PATTERN.findall(text)
    seen = set()
    result = []
    for doi in matches:
        doi = doi.rstrip('.,;:')
        if doi not in seen:
            seen.add(doi)
            result.append(doi)
    return result


def _parse_item(item: dict) -> PaperMetadata:
    """Parse a CrossRef API response item into PaperMetadata."""
    title = item.get('title', [''])[0] if item.get('title') else ''

    authors_raw = item.get('author', [])
    if authors_raw:
        names = []
        for a in authors_raw[:3]:
            family = a.get('family', '')
            given = a.get('given', '')
            names.append(f'{family} {given}'.strip() if given else family)
        authors = ', '.join(names)
        if len(authors_raw) > 3:
            authors += ' et al.'
    else:
        authors = ''

    journal = ''
    for key in ('container-title', 'short-container-title'):
        if item.get(key):
            journal = item[key][0]
            break

    year = ''
    for key in ('published-print', 'published-online', 'created'):
        date_parts = item.get(key, {}).get('date-parts')
        if date_parts and date_parts[0]:
            year = str(date_parts[0][0])
            break

    doi = item.get('DOI', '')
    url = f'https://doi.org/{doi}' if doi else ''

    return PaperMetadata(
        title=title, doi=doi, authors=authors,
        journal=journal, year=year, url=url,
    )


@retry(times=2, delay=3, exceptions=(Exception,))
def _query_by_doi(doi: str) -> Optional[PaperMetadata]:
    resp = requests.get(
        f'{CROSSREF_API}/{doi}',
        headers=REQUEST_HEADERS,
        timeout=15,
    )
    if resp.status_code != 200:
        return None
    item = resp.json().get('message', {})
    return _parse_item(item) if item.get('title') else None


@retry(times=2, delay=3, exceptions=(Exception,))
def _query_by_text(text: str) -> Optional[PaperMetadata]:
    text = text[:500].strip()
    if len(text) < 10:
        return None
    resp = requests.get(
        CROSSREF_API,
        params={'query.bibliographic': text, 'rows': 1},
        headers=REQUEST_HEADERS,
        timeout=15,
    )
    if resp.status_code != 200:
        return None
    items = resp.json().get('message', {}).get('items', [])
    if not items:
        return None
    item = items[0]
    if item.get('score', 0) < MIN_CROSSREF_SCORE:
        return None
    return _parse_item(item)


def resolve_reference(file_name: str, document_text: str) -> PaperMetadata:
    """
    Resolve a RAG file_name + document text to structured paper metadata.

    Three-level fallback:
      1. Extract DOI from document text → CrossRef DOI lookup
      2. Use file_name as title → CrossRef bibliographic search
      3. Use document text fragment → CrossRef bibliographic search
    """
    # Level 1: DOI in document text
    dois = extract_dois(document_text)
    for doi in dois:
        try:
            meta = _query_by_doi(doi)
            if meta and meta.title:
                logger.info(f'CrossRef [DOI]: {doi} → {meta.title[:80]}')
                time.sleep(RATE_LIMIT_INTERVAL)
                return meta
        except Exception:
            continue

    # Level 2: file_name as title
    if file_name and len(file_name) >= MIN_TITLE_LEN:
        try:
            meta = _query_by_text(file_name)
            if meta and meta.title:
                logger.info(f'CrossRef [title]: {file_name[:50]} → {meta.title[:80]}')
                time.sleep(RATE_LIMIT_INTERVAL)
                return meta
        except Exception:
            pass

    # Level 3: document text fragment
    try:
        meta = _query_by_text(document_text[:500])
        if meta and meta.title:
            logger.info(f'CrossRef [text]: → {meta.title[:80]}')
            time.sleep(RATE_LIMIT_INTERVAL)
            return meta
    except Exception:
        pass

    logger.warning(f'CrossRef unresolved: {file_name}')
    return PaperMetadata(title=file_name)
