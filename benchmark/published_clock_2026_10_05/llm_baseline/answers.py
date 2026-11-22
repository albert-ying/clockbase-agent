"""Read S4 JSON answers without inventing missing predictions."""
import json
import math
import re


def _boolean(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {'yes', 'true', 'significant'}:
            return True
        if text in {'no', 'false', 'not significant', 'non-significant'}:
            return False
    return None


def _unique_fields(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'Duplicate answer field: {key}')
        result[key] = value
    return result


def parse_answer(text):
    """Accept JSON or fenced JSON; reject incomplete or ambiguous answers."""
    if not isinstance(text, str):
        return None
    text = text.strip()
    fenced = re.fullmatch(r'```(?:json)?\s*(.*?)\s*```', text, re.S | re.I)
    if fenced:
        text = fenced.group(1)
    try:
        answer = json.loads(text, object_pairs_hook=_unique_fields)
    except (ValueError, TypeError):
        return None
    if not isinstance(answer, dict):
        return None
    direction = answer.get('prediction', answer.get('direction'))
    if not isinstance(direction, str):
        return None
    direction = direction.strip().lower().replace('_', ' ')
    direction = {'no difference': 'none', 'not different': 'none', 'unchanged': 'none'}.get(direction, direction)
    if direction not in {'higher', 'lower', 'none'}:
        return None
    significant = _boolean(answer.get('significance', answer.get('significant')))
    recognized = _boolean(answer.get('recognizes_publication'))
    confidence = answer.get('confidence')
    if isinstance(confidence, bool):
        return None
    try:
        confidence = float(confidence)
    except (TypeError, ValueError):
        return None
    if significant is None or recognized is None or not math.isfinite(confidence) or not 0 <= confidence <= 1:
        return None
    publication = answer.get('publication', answer.get('publication_reference'))
    if recognized:
        if not isinstance(publication, str) or not publication.strip():
            return None
        publication = publication.strip()
        if publication.lower() in {'none', 'null', 'unknown', 'n/a', 'no'}:
            return None
    else:
        if publication is not None and (not isinstance(publication, str) or
                                        publication.strip().lower() not in {'', 'none', 'null', 'n/a', 'no'}):
            return None
        publication = None
    if direction == 'none' and significant:
        return None
    return {'prediction': direction, 'significance': significant, 'confidence': confidence,
            'recognizes_publication': recognized, 'publication': publication}
