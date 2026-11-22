"""
Robust HTTP helpers for minipaper.

The original pipeline used bare `requests.get(...)` / `requests.post(...)` with
no timeout, which caused the entire asyncio event loop to freeze when NCBI
stalled mid-SSL-handshake (read() with no timeout blocks forever). This module
provides `safe_get` / `safe_post` wrappers that:

  - apply a sensible default timeout (connect 10s, read 30s)
  - retry transient network errors with exponential backoff (3 attempts)
  - fix the Oxylabs proxy dict (the original used a 'socks5' key, which
    `requests` ignores — it expects 'http' / 'https' keys).

Import pattern (drop-in replacement):
    from utils.http_client import safe_get, safe_post
    resp = safe_get(url)
    resp = safe_post(url, json=body)
"""

import logging
import os
import time
from typing import Optional

import requests
from requests.exceptions import (
    ConnectionError,
    ReadTimeout,
    ConnectTimeout,
    Timeout,
    ChunkedEncodingError,
)

logger = logging.getLogger(__name__)

# (connect, read) timeouts in seconds. Read is generous so we don't kill
# legitimate slow-but-progressing responses, but no longer "infinite".
DEFAULT_TIMEOUT = (10, 30)

# How many times to retry transient failures. Total wall clock on a fully
# timing-out call: ~ 30 + 30 + 30 + backoff = ~100s before giving up.
DEFAULT_RETRIES = 3

# Errors worth retrying. Everything else propagates.
_RETRIABLE = (
    ConnectionError,
    ReadTimeout,
    ConnectTimeout,
    Timeout,
    ChunkedEncodingError,
)


def get_proxies() -> Optional[dict]:
    """Build a correct `requests` proxies dict from PROXY_* env vars, or None.

    The original `get_data.py` used a dict keyed on 'socks5' which `requests`
    never honours. This returns properly-keyed http/https entries so proxies
    actually engage when PROXY_HOST is set.
    """
    host = os.getenv("PROXY_HOST")
    port = os.getenv("PROXY_PORT")
    user = os.getenv("PROXY_USER")
    pwd = os.getenv("PROXY_PWD")
    if not host or not port:
        return None
    auth = f"{user}:{pwd}@" if user and pwd else ""
    url = f"http://{auth}{host}:{port}"
    return {"http": url, "https": url}


def _retry_request(method: str, url: str, **kwargs):
    """Core retry loop shared by safe_get / safe_post.

    - Sets timeout if caller didn't.
    - Injects proxies from env if caller didn't and USE_PROXY env var is truthy.
    - Retries on the _RETRIABLE exception tuple with backoff 1s, 2s, 4s.
    """
    kwargs.setdefault("timeout", DEFAULT_TIMEOUT)
    if "proxies" not in kwargs and os.getenv("USE_PROXY", "").lower() in ("1", "true", "yes"):
        kwargs["proxies"] = get_proxies()

    retries = kwargs.pop("_retries", DEFAULT_RETRIES)
    last_exc = None
    for attempt in range(1, retries + 1):
        try:
            return requests.request(method, url, **kwargs)
        except _RETRIABLE as e:
            last_exc = e
            if attempt >= retries:
                logger.warning(f"HTTP {method} {url} failed after {attempt} attempts: {e}")
                raise
            backoff = 2 ** (attempt - 1)
            logger.info(f"HTTP {method} {url} attempt {attempt} failed ({type(e).__name__}), retrying in {backoff}s")
            time.sleep(backoff)
    # unreachable, but make the type checker happy
    raise last_exc if last_exc else RuntimeError("retry loop exited without exception")


def safe_get(url: str, **kwargs):
    """`requests.get` with timeout + retry + correct proxy handling."""
    return _retry_request("GET", url, **kwargs)


def safe_post(url: str, **kwargs):
    """`requests.post` with timeout + retry + correct proxy handling."""
    return _retry_request("POST", url, **kwargs)
