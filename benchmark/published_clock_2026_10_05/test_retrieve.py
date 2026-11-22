"""Retrieval checks using synthetic prose/API fixtures, never benchmark measurements."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

from benchmark.published_clock_2026_10_05.retrieve import (
    CachedGET, OUTPUT, accession_ids, fetch_annotations, fetch_ncbi_links, load_atlas, ncbi_accessions, search, sentences, text_blocks, triage,
)


class RetrievalTest(unittest.TestCase):
    def setUp(self):
        folder = OUTPUT / 'checks'
        folder.mkdir(parents=True, exist_ok=True)
        self.workspace = tempfile.TemporaryDirectory(dir=folder)
        self.root = Path(self.workspace.name)
        self.addCleanup(self.workspace.cleanup)

    def test_http_failure_is_cached_without_offline_refetch(self):
        response = Mock(status_code=404, content=b'not available', headers={})
        with patch('requests.get', return_value=response) as get:
            body, meta = CachedGET(self.root).get('https://example.invalid/fullTextXML')
            replay, replay_meta = CachedGET(self.root, offline=True).get('https://example.invalid/fullTextXML')
        self.assertEqual(get.call_count, 1)
        self.assertEqual(body, replay)
        self.assertEqual(meta, replay_meta)
        self.assertEqual(meta['status'], 404)
        self.assertEqual(meta['sha256'], hashlib.sha256(body).hexdigest())

    def test_cache_corruption_stops_without_network(self):
        response = Mock(status_code=200, content=b'original', headers={})
        with patch('requests.get', return_value=response):
            _, meta = CachedGET(self.root).get('https://example.invalid/paper')
        Path(meta['cache_body']).write_bytes(b'changed')
        with patch('requests.get') as get, self.assertRaisesRegex(RuntimeError, 'checksum mismatch'):
            CachedGET(self.root).get('https://example.invalid/paper')
        get.assert_not_called()

    def test_transient_responses_are_preserved_before_retry(self):
        responses = [Mock(status_code=503, content=b'retry', headers={}),
                     Mock(status_code=200, content=b'ok', headers={})]
        with patch('requests.get', side_effect=responses), patch('time.sleep'):
            body, meta = CachedGET(self.root).get('https://example.invalid/retry')
        self.assertEqual(body, b'ok')
        self.assertEqual(meta['attempt'], 2)
        attempts = list((self.root / 'cache').glob('*.attempt-1.body'))
        self.assertEqual(len(attempts), 1)
        self.assertEqual(attempts[0].read_bytes(), b'retry')

    def test_inline_markup_abstracts_tables_and_reference_exclusion(self):
        blocks = text_blocks(b'<article><front><article-meta><abstract><p>Clock abstract.</p>'
                            b'</abstract></article-meta></front><body><sec><title>Results</title>'
                            b'<p><italic>Dunedin</italic>PACE in GSE123.</p><table-wrap id="t1">'
                            b'<label>Table 1</label><table><tr><td>GSE456</td><td>control</td></tr>'
                            b'</table></table-wrap></sec></body><back><ref-list><ref><p>GSE789'
                            b'</p></ref></ref-list></back></article>')
        self.assertEqual(blocks[0]['section'], 'Abstract')
        self.assertIn('DunedinPACE', blocks[1]['text'])
        self.assertEqual(blocks[2]['label'], 'Table 1')
        self.assertEqual(accession_ids(blocks[2]['text']), ['GSE456'])
        self.assertFalse(any('GSE789' in b['text'] for b in blocks))

    def classify(self, paragraphs, publication_types=()):
        blocks = [{'text': text, 'section': section, 'label': None, 'element': 'p',
                   'block_index': i} for i, (section, text) in enumerate(paragraphs)]
        paper = {'paper_id': 'TEST:fixture', 'ids': {}, 'title': 'Synthetic test prose', 'year': None,
                 'publication_types': publication_types, 'fulltext_status': 'ok', 'text_path': 'fixture',
                 'accession_evidence': []}
        atlas = {'modalities': ['human DNAm'], 'atlas_rows': 1, 'atlas_clocks': ['PhenoAge'],
                 'atlas_clock_count': 1, 'atlas_by_modality': {}}
        return triage(paper, 'GSE123', blocks, atlas)

    def test_accession_local_reported_clock_comparison(self):
        candidate = self.classify([('Results', 'In GSE123, we calculated PhenoAge; patients had '
                                               'significantly higher epigenetic age than controls.')])
        self.assertEqual(candidate['triage_class'], 'LIKELY')
        self.assertTrue(candidate['passages'])

    def test_unrelated_result_does_not_bind_to_accession(self):
        candidate = self.classify([('Data availability', 'RNA data are in GSE123.'),
                                   ('Results', 'We measured PhenoAge in a separate cohort: patients '
                                               'had significantly higher age than controls.')])
        self.assertEqual(candidate['triage_class'], 'POSSIBLE')
        self.assertFalse(candidate['rule_features']['direct_blocks'])

    def test_training_and_reviews_are_not_automatically_likely(self):
        text = 'GSE123 was the training set for PhenoAge; patients had significantly higher age than controls.'
        self.assertNotEqual(self.classify([('Methods', text)])['triage_class'], 'LIKELY')
        self.assertEqual(self.classify([('Results', text)], ['Review'])['triage_class'], 'UNLIKELY')

    def test_cpg_group_result_does_not_become_a_clock_result(self):
        candidate = self.classify([('Results', 'In GSE123, CpG levels were significantly higher in patients '
                                               'than controls. These observations inform our epigenetic age predictions.')])
        self.assertNotEqual(candidate['triage_class'], 'LIKELY')

    def test_figures_and_tables_inside_paragraphs_keep_their_locations(self):
        blocks = text_blocks(b'<article><body><sec><title>Results</title><p>First sentence.'
                            b'<table-wrap><label>Table 2</label><table><tr><td>GSE123</td><td>PhenoAge</td>'
                            b'</tr></table></table-wrap> Last sentence.<fig><label>Figure 1</label>'
                            b'<caption><p>Clock caption.</p></caption></fig></p></sec></body></article>')
        self.assertEqual(blocks[0]['text'], 'First sentence. Last sentence.')
        self.assertEqual(blocks[1]['label'], 'Table 2')
        self.assertEqual(blocks[2]['label'], 'Figure 1')
        self.assertNotIn('GSE123', blocks[0]['text'])

    def test_citation_numbers_do_not_extend_a_geo_accession(self):
        blocks = text_blocks(b'<article><body><p>Data are GSE123<xref ref-type="bibr">24</xref>.'
                            b'</p></body></article>')
        self.assertEqual(accession_ids(blocks[0]['text']), ['GSE123'])

    def test_geo_dataset_summary_retains_its_originating_series(self):
        self.assertEqual(ncbi_accessions({'accession': 'GDS123', 'gse': '456', 'gpl': '0'}), ['GSE456'])
        self.assertEqual(ncbi_accessions({'accession': 'GSE789', 'gse': '789'}), ['GSE789'])
        self.assertEqual(ncbi_accessions({'accession': 'GDS123', 'gse': '0'}), [])

    def test_sentences_keep_comparison_and_figure_abbreviations(self):
        text = 'Clock age was higher in OIS vs. EP (Fig. S2). Next result.'
        self.assertEqual(sentences(text), ['Clock age was higher in OIS vs. EP (Fig. S2).', 'Next result.'])

    def test_ncbi_links_keep_each_pmid_separate(self):
        client = Mock()
        meta = {'sha256': 'synthetic test hash', 'url': 'https://example.invalid/fixture'}
        link_data = {'linksets': [{'ids': ['1'], 'linksetdbs': [{'dbto': 'gds', 'links': ['101']}]},
                                 {'ids': ['2'], 'linksetdbs': [{'dbto': 'gds', 'links': ['202']}]}]}
        summary_data = {'result': {'101': {'accession': 'GSE123'}, '202': {'accession': 'GSE456'}}}
        client.json.side_effect = [(link_data, meta), (summary_data, meta)]
        papers = {'MED:1': {'record': {'pmid': '1'}}, 'MED:2': {'record': {'pmid': '2'}}}
        links, errors = fetch_ncbi_links(client, papers, self.root)
        self.assertFalse(errors)
        self.assertEqual([e['accession'] for e in links['MED:1']['accession_evidence']], ['GSE123'])
        self.assertEqual([e['accession'] for e in links['MED:2']['accession_evidence']], ['GSE456'])
        self.assertEqual(client.json.call_args_list[0].args[1]['id'], ['1', '2'])

    def test_ncbi_merged_pmid_links_are_rejected(self):
        client = Mock()
        client.json.return_value = {'linksets': [{'ids': ['1', '2'], 'linksetdbs': []}]}, {}
        papers = {'MED:1': {'record': {'pmid': '1'}}, 'MED:2': {'record': {'pmid': '2'}}}
        with self.assertLogs(level='WARNING'):
            links, errors = fetch_ncbi_links(client, papers, self.root)
        self.assertEqual(len(errors), 1)
        self.assertTrue(all(entry['status'] == 'failed' for entry in links.values()))
        self.assertTrue(all(not entry['accession_evidence'] for entry in links.values()))

    def test_missing_empirical_atlas_fails_explicitly(self):
        with self.assertRaisesRegex(RuntimeError, 'Cannot proceed - required data file not found:'):
            load_atlas(self.root / 'absent.csv')

    def test_cursor_paging_and_cross_query_deduplication(self):
        client = Mock()
        first = {'hitCount': 2, 'nextCursorMark': 'next',
                 'resultList': {'result': [{'id': '1', 'pmid': '1', 'source': 'MED'}]}}
        second = {'hitCount': 2, 'nextCursorMark': 'next',
                  'resultList': {'result': [{'id': '2', 'pmid': '2', 'source': 'MED'}]}}
        meta = {'fetched_at': 'synthetic test timestamp'}
        client.json.side_effect = [(first, meta), (second, meta)] * 3
        papers, queries = search(client, self.root)
        self.assertEqual(set(papers), {'MED:1', 'MED:2'})
        self.assertEqual(len(queries), 3)
        self.assertTrue(all(q['retrieved_records'] == 2 for q in queries))
        cursors = [call.args[1]['cursorMark'] for call in client.json.call_args_list]
        self.assertEqual(cursors, ['*', 'next'] * 3)

    def test_stalled_cursor_fails_before_partial_outputs(self):
        client = Mock()
        data = {'hitCount': 2, 'nextCursorMark': '*',
                'resultList': {'result': [{'id': '1', 'pmid': '1', 'source': 'MED'}]}}
        client.json.return_value = data, {'fetched_at': 'synthetic test timestamp'}
        with self.assertRaisesRegex(RuntimeError, 'Paging stopped before hit count'):
            search(client, self.root)

    def test_empty_page_cannot_silently_truncate_search(self):
        client = Mock()
        client.json.return_value = {'hitCount': 2, 'resultList': {'result': []}}, {'fetched_at': 'synthetic test timestamp'}
        with self.assertRaisesRegex(RuntimeError, 'Incomplete search'):
            search(client, self.root)

    def test_duplicate_article_annotation_records_are_combined(self):
        annotation = {'exact': 'GSE123', 'subType': 'GEO', 'tags': []}
        client = Mock()
        client.json.return_value = [
            {'pmcid': 'PMC1', 'source': 'MED', 'extId': '1', 'annotations': [annotation]},
            {'pmcid': 'PMC1', 'source': 'PMC', 'extId': 'PMC1', 'annotations': []},
            {'pmcid': 'PMC1', 'source': 'MED', 'extId': '1', 'annotations': [annotation]},
        ], {'sha256': 'synthetic test response'}
        papers = {'MED:1': {'record': {'source': 'MED', 'id': '1', 'pmid': '1', 'pmcid': 'PMC1'}}}
        annotations, errors = fetch_annotations(client, papers, self.root)
        self.assertEqual(annotations['MED:1']['annotations'], [annotation])
        self.assertEqual(errors, [])


if __name__ == '__main__':
    unittest.main()
