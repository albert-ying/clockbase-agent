"""Cache public curated-source GETs with byte hashes and retrieval receipts."""
import argparse
import hashlib
import json
import logging
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path('/Entropy/kying/clockbase-benchmark-2026-10-05/curated')


def retrieve(url, name, root=ROOT):
    """Return a cached response, recording successful and failed GETs."""
    if Path(name).name != name or name in {'', '.', '..'}:
        raise ValueError('Cache name must be one filename')
    cache = root / 'downloads'
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / name
    receipt = cache / (name + '.receipt.json')
    if path.exists() and receipt.exists():
        saved = json.loads(receipt.read_text())
        if saved['url'] != url or saved.get('sha256') != hashlib.sha256(path.read_bytes()).hexdigest():
            raise ValueError(f'Cache identity mismatch: {path}')
        return path
    if path.exists() or receipt.exists():
        raise FileExistsError(f'Incomplete cache entry: {path}')
    time.sleep(0.5)
    entry = {'url': url, 'retrieved_at': datetime.now(timezone.utc).isoformat(),
             'path': str(path), 'method': 'GET'}
    try:
        request = urllib.request.Request(url, headers={'User-Agent': 'ClockBase-curated-literature-audit/1.0'})
        with urllib.request.urlopen(request, timeout=90) as response:
            payload = response.read()
            entry.update(status=response.status, final_url=response.url,
                         content_type=response.headers.get('Content-Type'))
        path.write_bytes(payload)
        entry.update(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
        receipt.write_text(json.dumps(entry, indent=2) + '\n')
        logging.info('Cached %s (%s bytes)', name, len(payload))
        return path
    except (OSError, ValueError) as exc:
        entry.update(error=str(exc))
        if isinstance(exc, urllib.error.HTTPError):
            payload = exc.read()
            error_path = cache / (name + '.error')
            error_path.write_bytes(payload)
            entry.update(status=exc.code, error_path=str(error_path),
                         error_sha256=hashlib.sha256(payload).hexdigest())
        receipt.write_text(json.dumps(entry, indent=2) + '\n')
        logging.error('GET failed for %s: %s', url, exc)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('url')
    parser.add_argument('name')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    print(retrieve(args.url, args.name))


if __name__ == '__main__':
    main()
