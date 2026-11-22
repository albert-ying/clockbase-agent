"""Explicit SOFTWARE FIXTURES, not empirical studies or audit performance estimates.

The schema is the unchanged sealed contract. All SOFT values, GSM memberships,
reference interpretations and release declarations below are software fixtures.
No empirical source payload is needed or read by this test module.
"""
import copy
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from audit_system import evaluation_2026_09_10 as subject
from audit_system import agent_workflow_2026_09_10 as workflow
from src import pipeline_runtime_2026_09_10 as runtime


SCHEMA = Path(__file__).parent / 'fixtures/source_readiness_judgment_2026_09_24.schema.json'


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')
    return {'path': str(path), 'sha256': workflow.sha(path)}


class SourceReadinessTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='RNA80_SOFTWARE_FIXTURE_')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.schema = json.loads(SCHEMA.read_text())
        self.allowed_ids = self.schema['properties']['study_id']['enum']
        ids = self.allowed_ids[:]
        ids.insert(60, 'GSE295890'); ids.insert(68, 'GSE298220')
        self.enrollment = {'schema': 'clockbase-rna80-source-readiness-enrollment-v1',
            'state': 'SOFTWARE_FIXTURE_NOT_EMPIRICAL', 'split': 'development_exposed',
            'source_manifest_sha256': hashlib.sha256(b'SOFTWARE FIXTURE enrollment').hexdigest(),
            'enrolled_studies': 80, 'enrolled_source_records': 320,
            'source_access_allowed_studies': 78, 'source_access_allowed_records': 312,
            'administrative_holds': 2, 'administrative_hold_source_records': 8, 'studies': []}
        self.allowed = {'schema': 'clockbase-rna80-source-readiness-allowed-sources-v1',
            'state': 'SOFTWARE_FIXTURE_NOT_EMPIRICAL', 'split': 'development_exposed',
            'source_manifest_sha256': self.enrollment['source_manifest_sha256'], 'studies': []}
        for index, sid in enumerate(ids, 1):
            hold = sid in subject.SOURCE_READINESS_HOLDS
            self.enrollment['studies'].append({'study_id': sid, 'manifest_ordinal': index,
                'source_record_count_from_fixed_manifest': 4, 'reference_judgment': None,
                'new_source_access_allowed': not hold,
                'administrative_state': subject.SOURCE_READINESS_HOLDS[sid] if hold else 'ENROLLED_NOT_ASSESSED',
                'reason': 'SOFTWARE FIXTURE administrative state, not an empirical judgment'})
            if hold: continue
            samples = ['GSM' + str(index * 10 + i) for i in range(1, 5)]
            sources = {}
            for role in ('series', 'samples'):
                path = self.root / 'SOFTWARE_FIXTURE_RAW' / sid / (role + '.soft')
                if len(self.allowed['studies']) < 5:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(self.soft(sid, samples, role).encode())
                    sources[role] = {**{'path': str(path), 'sha256': workflow.sha(path)}, 'bytes': path.stat().st_size}
                else:
                    # No unselected raw files exist; the projection must never scan them.
                    sources[role] = {'path': str(path), 'sha256': '0' * 64, 'bytes': 1}
            self.allowed['studies'].append({'study_id': sid, 'manifest_ordinal': index,
                                           'sample_ids': samples, 'sources': sources})
        self.schedule = {'schema': 'clockbase-rna80-source-readiness-schedule-v1', 'batch_size': 5, 'batches': []}
        for offset in range(0, 78, 5):
            members = self.allowed['studies'][offset:offset+5]
            self.schedule['batches'].append({'batch_id': 'source-reference-batch-%03d' % (offset // 5 + 1),
                'study_ids': [r['study_id'] for r in members],
                'manifest_ordinals': [r['manifest_ordinal'] for r in members], 'source_record_count': 4 * len(members)})
        self.batch = self.schedule['batches'][0]['batch_id']
        self.sid = self.allowed_ids[0]
        self.refresh_release()
        self.projection = None

    @staticmethod
    def soft(sid, samples, role):
        if role == 'series':
            return ('^SERIES = ' + sid + '\n!Series_geo_accession = ' + sid + '\n'
                '!Series_title = SOFTWARE FIXTURE study, not empirical data\n'
                '!Series_summary = EXCLUDED SOFTWARE FIXTURE SOURCE SUMMARY\n'
                '!Series_overall_design = SOFTWARE FIXTURE paired units\n'
                '  continuation with two leading spaces  \n\n'
                '\tlast continuation\n'
                '!Series_overall_design = Repeated design declaration\n'
                '!Series_relation = BioProject: PRJ_SOFTWARE_FIXTURE1\n'
                '!Series_contact_email = excluded@example.invalid\n'
                '!Series_contact_zip/postal_code = SOFTWARE_FIXTURE_POSTAL\n'
                '  SOFTWARE_FIXTURE_POSTAL continuation\n' +
                ''.join('!Series_sample_id = ' + s + '\n' for s in samples))
        records = []
        for i, sample in enumerate(samples):
            records.append('^SAMPLE = ' + sample + '\n!Sample_geo_accession = ' + sample + '\n'
                '!Sample_title = SOFTWARE FIXTURE ' + ('control' if i < 2 else 'case') + str(i % 2 + 1) + '\n'
                '!Sample_characteristics_ch1 = condition: ' + ('control' if i < 2 else 'case') + '\n'
                '!Sample_characteristics_ch1 = unit: unit-' + str(i % 2 + 1) + '\n'
                '!Sample_characteristics_ch1 = age: ' + str(20 + i % 2) + '\n'
                '!Sample_description = SOFTWARE FIXTURE first line  \n'
                '  second line\n\n\tlast description line\n'
                '!Sample_contact_zip/postal_code = SOFTWARE_FIXTURE_POSTAL\n'
                '  SOFTWARE_FIXTURE_POSTAL continuation\n'
                '!Sample_data_processing = SOFTWARE FIXTURE processing declaration\n'
                '!Sample_library_strategy = RNA-Seq\n!Sample_organism_ch1 = Homo sapiens\n'
                '!Sample_relation = BioProject: PRJ_SOFTWARE_FIXTURE1\n'
                '!Sample_series_id = ' + sid + '\n!Sample_contact_email = excluded@example.invalid\n')
        return ''.join(records)

    def refresh_release(self, *, phase='SOURCE_ONLY_FROZEN'):
        docs = self.root / 'SOFTWARE_FIXTURE_PROTOCOL'; docs.mkdir(exist_ok=True)
        protocol = docs / 'PROTOCOL.md'
        protocol.write_text('SOFTWARE FIXTURE protocol; no empirical/source release is authorized by these tests.\n')
        bindings = {'protocol': {'path': str(protocol), 'sha256': workflow.sha(protocol)}}
        for key, value in [('enrollment', self.enrollment), ('allowed_sources', self.allowed),
                           ('schedule', self.schedule), ('judgment_schema', self.schema)]:
            bindings[key] = write_json(docs / (key + '.json'), value)
        self.release_data = {'schema': 'clockbase-rna80-source-readiness-release-v1', 'phase': phase,
            'split': 'development_exposed', 'protocol_commit': 'f' * 40,
            'reference_runtime': {'runtime': 'SOFTWARE FIXTURE no model runtime', 'model': 'SOFTWARE_FIXTURE_NOT_A_MODEL'},
            'jsonschema_version': version('jsonschema'), 'bindings': bindings,
            'implementation_sha256': {'audit_system/evaluation_2026_09_10.py': workflow.sha(subject.__file__),
                'audit_system/agent_workflow_2026_09_10.py': workflow.sha(workflow.__file__),
                'src/pipeline_runtime_2026_09_10.py': workflow.sha(runtime.__file__)}}
        self.release = write_json(docs / 'release.json', self.release_data)

    def context(self):
        return subject.source_readiness_context(self.release['path'], self.release['sha256'])

    def project(self):
        if self.projection is None:
            self.projection = subject.prepare_source_readiness_batch(self.release['path'], self.release['sha256'],
                                                                     self.batch, self.root / 'projection')
        binding = self.projection['packets'][0]
        return {'path': binding['path'], 'sha256': binding['sha256']}, json.loads(Path(binding['path']).read_text())

    def submission(self, *, comparison=False):
        binding, packet = self.project()
        data = json.loads(Path(self.projection['packets'][0]['form']).read_text())
        data.update(completed=True, comparison_inventory_status='UNRESOLVED',
                    unassigned_sample_ids=packet['sample_ids'], unit_unresolved_sample_ids=packet['sample_ids'])
        data['reviewer'].update(name='SOFTWARE FIXTURE source assessor', independent_from_methods=True,
            runtime_model_as_reported='SOFTWARE_FIXTURE_NOT_A_MODEL', known_prior_exposure='SOFTWARE FIXTURE only',
            source_narrative_exposure_note='SOFTWARE FIXTURE narratives, no model or empirical review')
        for domain in subject.SOURCE_READINESS_DOMAINS:
            key = 'u-' + domain
            data['unresolved'].append({'unresolved_id': key, 'domain': domain,
                'question': 'SOFTWARE FIXTURE unresolved scope', 'reason': 'SOFTWARE FIXTURE missing prerequisite',
                'searched_source_roles': ['series'], 'searched_fields': ['!Series_overall_design'],
                'sample_ids': packet['sample_ids'], 'comparison_ids': [], 'related_evidence_ids': [],
                'needed_evidence': 'SOFTWARE FIXTURE additional source evidence'})
            section = 'domains' if domain in subject.REFERENCE_DESIGN_DOMAINS else 'readiness_domains'
            data[section][domain] = {'status': 'CANNOT VERIFY', 'material': None,
                'rationale': 'SOFTWARE FIXTURE unresolved judgment; not an empirical label',
                'evidence_ids': [], 'fact_ids': [], 'unresolved_ids': [key]}
        if not comparison: return data
        data['evidence'] = {k: v for k, v in packet['literal_witnesses'].items() if v['quote'].strip()}
        samples = packet['sample_ids']
        def quotes(ids, field, prefix=None):
            return [k for k, v in data['evidence'].items() if v['entity'] in ids and v['field'] == field
                    and (prefix is None or v['quote'].startswith(prefix))]
        def fact(fid, domain, ids, witnesses):
            data['facts'].append({'fact_id': fid, 'domain': domain, 'sample_ids': ids,
                'comparison_ids': [], 'qualification': 'LITERAL_SOURCE_DECLARATION',
                'statement': 'SOFTWARE FIXTURE declared source fact', 'evidence_ids': witnesses})
        comparisons_id = workflow.comparison_identity(self.sid, samples[2:], samples[:2], 'SOFTWARE FIXTURE contrast')
        core_facts = {
            'sample_identity': ('f-identity', quotes(samples, '!Sample_geo_accession')),
            'comparison_definition': ('f-arms', quotes(samples, '!Sample_characteristics_ch1', 'condition:')),
            'confounding': ('f-age', quotes(samples, '!Sample_characteristics_ch1', 'age:'))}
        for domain, (fid, witnesses) in core_facts.items():
            fact(fid, domain, samples, witnesses)
            data['domains'][domain] = {'status': 'PASS', 'material': False,
                'rationale': 'SOFTWARE FIXTURE scoped source declaration',
                'evidence_ids': witnesses, 'fact_ids': [fid], 'unresolved_ids': []}
        for i in range(2):
            members = [samples[i], samples[i+2]]
            witnesses = quotes(members, '!Sample_characteristics_ch1', 'unit:')
            fact('f-unit-' + str(i), 'biological_units', members, witnesses)
            fact('f-pair-' + str(i), 'pairing', members, witnesses)
            data['unit_map'].append({'unit_id': 'unit-' + str(i+1), 'unit_type': 'PARTICIPANT',
                'sample_ids': members, 'fact_ids': ['f-unit-' + str(i)], 'evidence_ids': witnesses})
            data['pair_map'].append({'comparison_id': comparisons_id, 'unit_id': 'unit-' + str(i+1),
                'case_sample_id': samples[i+2], 'control_sample_id': samples[i],
                'fact_ids': ['f-pair-' + str(i)], 'evidence_ids': witnesses})
        for domain, prefix in [('biological_units', 'f-unit-'), ('pairing', 'f-pair-')]:
            data['domains'][domain] = {'status': 'PASS', 'material': False,
                'rationale': 'SOFTWARE FIXTURE explicit source unit key',
                'evidence_ids': quotes(samples, '!Sample_characteristics_ch1', 'unit:'),
                'fact_ids': [prefix + '0', prefix + '1'], 'unresolved_ids': []}
        data['comparison_inventory_status'] = 'SOURCE_DEFINED_CONTRASTS_LISTED'
        data['comparisons'] = [{'comparison_id': comparisons_id, 'target': 'SOFTWARE FIXTURE contrast',
            'contrast_type': 'INTERVENTIONAL', 'case_ids': samples[2:], 'control_ids': samples[:2],
            'evidence_ids': core_facts['comparison_definition'][1], 'domains': copy.deepcopy(data['domains']),
            'conditional_test_family': 'PAIRED', 'analysis_requirements': [{'requirement': 'SOFTWARE FIXTURE outcome assumptions remain open',
                'fact_ids': ['f-pair-0', 'f-pair-1'], 'unresolved_ids': ['u-test_requirements']}],
            'excluded_sample_ids': [], 'exclusion_reason': 'SOFTWARE FIXTURE no exclusions'}]
        data['unassigned_sample_ids'] = []; data['unit_unresolved_sample_ids'] = []
        data['covariates'] = [{'name': 'SOFTWARE FIXTURE age', 'scope': 'STUDY', 'observed_sample_ids': samples,
            'missing_sample_ids': [], 'fact_ids': ['f-age'], 'evidence_ids': core_facts['confounding'][1], 'unresolved_ids': []}]
        data['family_links'] = [{'link_type': 'BIOPROJECT', 'target_identifier': 'PRJ_SOFTWARE_FIXTURE1',
            'evidence_ids': quotes(samples[:1], '!Sample_relation'),
            'interpretation': 'SOFTWARE FIXTURE literal link; does not close study families'}]
        return data

    def freeze(self, data, name='version-001', *, supersedes=None, output=None):
        packet, _ = self.project()
        path = self.root / 'submissions' / (name + '.json')
        write_json(path, data)
        return subject.freeze_source_readiness_submission(self.release['path'], self.release['sha256'], self.batch,
            self.sid, packet, path, name, output or self.root / name, supersedes=supersedes)

    @staticmethod
    def frozen_binding(result):
        return {'path': result['frozen_submission'], 'sha256': result['sha256']}

    def assert_rejected(self, data, match, name='version-001'):
        result = self.freeze(data, name)
        self.assertEqual(result['validation_status'], 'REJECTED')
        frozen = json.loads(Path(result['frozen_submission']).read_text())
        self.assertRegex(frozen['validation']['reason'], match)
        self.assertEqual(frozen['validation']['readiness_status'], 'NOT_ASSESSED')
        self.assertNotIn('domains', frozen['validation'])
        return result

    def test_sealed_schema_is_unchanged(self):
        self.assertEqual(workflow.sha(SCHEMA), '33fc5f4d3dc151499b865efec4a7825245f319809902183c4231cae61ab27f2a')

    def test_draft_or_changed_release_rejects_before_raw_access(self):
        self.refresh_release(phase='DRAFT_NOT_RELEASED')
        with patch.object(subject, '_readiness_raw', side_effect=AssertionError('must not read raw')):
            with self.assertRaisesRegex(ValueError, 'SOURCE_ONLY_FROZEN'):
                subject.prepare_source_readiness_batch(self.release['path'], self.release['sha256'], self.batch, self.root/'no-projection')
            with self.assertRaisesRegex(ValueError, 'hash changed'):
                subject.source_readiness_context(self.release['path'], '0'*64)

    def test_model_declaration_and_code_hash_are_bound_before_sources(self):
        for mutation in ('model', 'code'):
            self.refresh_release()
            if mutation == 'model': self.release_data['reference_runtime']['model'] = ''
            else: self.release_data['implementation_sha256']['audit_system/evaluation_2026_09_10.py'] = '0'*64
            self.release = write_json(Path(self.release['path']), self.release_data)
            with patch.object(subject, '_readiness_raw', side_effect=AssertionError('must not read raw')):
                with self.assertRaises(ValueError): self.context()

    def test_batch_and_raw_path_guards_precede_all_file_operations(self):
        context = self.context(); allowed = context['allowed_by_id'][self.sid]['sources']['series']
        cases = [('not-a-batch', self.sid, allowed), (self.batch, 'GSE295890', allowed),
                 (self.batch, 'GSE298220', allowed), (self.batch, self.allowed_ids[5], allowed),
                 (self.batch, self.sid, {**allowed, 'path': '/NEVER_OPEN_UNLISTED'}),
                 (self.batch, self.sid, {**allowed, 'sha256': '0'*64})]
        with patch('os.open', side_effect=AssertionError('opened denied source')), \
             patch.object(Path, 'is_symlink', side_effect=AssertionError('statted denied source')):
            for batch, sid, binding in cases:
                with self.subTest(batch=batch, sid=sid, binding=binding):
                    with self.assertRaisesRegex(ValueError, 'Unlisted|held'):
                        subject._readiness_raw(context, batch, sid, 'series', binding)

    def test_held_submission_rejected_before_packet_or_submission_open(self):
        for sid in subject.SOURCE_READINESS_HOLDS:
            with patch.object(subject, '_readiness_raw', side_effect=AssertionError('opened raw')):
                with self.assertRaisesRegex(ValueError, 'not released'):
                    subject.freeze_source_readiness_submission(self.release['path'], self.release['sha256'], self.batch,
                        sid, {'path': '/NEVER_OPEN_PACKET', 'sha256': '0'*64}, '/NEVER_OPEN_SUBMISSION',
                        'held-test', self.root/'no-hold-freeze')

    def test_enrollment_holds_and_fixed_batch_order_cannot_be_changed(self):
        self.enrollment['studies'][60]['new_source_access_allowed'] = True
        self.refresh_release()
        with self.assertRaisesRegex(ValueError, 'hold'): self.context()
        self.enrollment['studies'][60]['new_source_access_allowed'] = False
        self.schedule['batches'][0]['study_ids'].reverse()
        self.refresh_release()
        with self.assertRaisesRegex(ValueError, 'order'): self.context()

    def test_projection_keeps_repeated_multiline_values_and_no_default_statuses(self):
        binding, packet = self.project()
        self.assertEqual([r['study_id'] for r in self.projection['packets']], self.allowed_ids[:5])
        self.assertEqual(self.projection['source_files_read'], 10)
        sample = packet['source_packet']['samples'][0]
        self.assertEqual(sample['characteristics_ch1'], 'condition: control')
        self.assertEqual(sample['characteristics_ch1.1'], 'unit: unit-1')
        self.assertEqual(sample['characteristics_ch1.2'], 'age: 20')
        self.assertEqual(sample['description'], 'SOFTWARE FIXTURE first line  \n  second line\n\n\tlast description line')
        fields = packet['source_packet']['study_documents'][0]['fields']
        self.assertEqual(fields['Series_overall_design.1'], 'Repeated design declaration')
        self.assertIn('\n\n\tlast continuation', fields['Series_overall_design'])
        self.assertNotIn('Series_summary', fields)
        self.assertNotIn('Series_contact_email', fields); self.assertNotIn('contact_email', sample)
        self.assertIn('data_processing', sample)
        witnesses = [w for w in packet['literal_witnesses'].values() if w['field'] == '!Sample_description']
        self.assertTrue(all(w['line_end'] - w['line_start'] == 3 for w in witnesses))
        form = json.loads(Path(self.projection['packets'][0]['form']).read_text())
        self.assertFalse(form['completed']); self.assertTrue(all(v is None for v in form['domains'].values()))
        self.assertTrue(all(v is None for v in form['readiness_domains'].values()))

    def test_crlf_and_empty_continuations_are_literal(self):
        context = self.context(); row = context['allowed_by_id'][self.sid]
        raw = self.soft(self.sid, row['sample_ids'], 'samples').replace('\n', '\r\n').encode()
        records = subject._readiness_soft(raw, 'samples', self.sid, row['sample_ids'], row['sources']['samples'])
        descriptions = [r for r in records if r['field'] == '!Sample_description']
        self.assertTrue(all('\r\n\r\n\t' in r['quote'] for r in descriptions))

    def test_slash_named_excluded_fields_preserve_neighboring_scientific_values(self):
        row = self.context()['allowed_by_id'][self.sid]
        for role in ('series', 'samples'):
            raw = self.soft(self.sid, row['sample_ids'], role).encode()
            records = subject._readiness_soft(raw, role, self.sid, row['sample_ids'], row['sources'][role])
            self.assertTrue(records)
            self.assertTrue(all('/' not in r['field'] for r in records))
            self.assertTrue(all('SOFTWARE_FIXTURE_POSTAL' not in r['quote'] for r in records))
            if role == 'samples':
                descriptions = [r['quote'] for r in records if r['field'] == '!Sample_description']
                self.assertEqual(descriptions, ['SOFTWARE FIXTURE first line  \n  second line\n\n\tlast description line'] * len(row['sample_ids']))
                processing = [r['quote'] for r in records if r['field'] == '!Sample_data_processing']
                self.assertEqual(processing, ['SOFTWARE FIXTURE processing declaration'] * len(row['sample_ids']))
            else:
                members = [r['quote'] for r in records if r['field'] == '!Series_sample_id']
                self.assertEqual(members, row['sample_ids'])
                designs = [r['quote'] for r in records if r['field'] == '!Series_overall_design']
                self.assertEqual(designs[-1], 'Repeated design declaration')

    def test_slash_field_names_still_reject_malformed_segments_roles_and_delimiters(self):
        row = self.context()['allowed_by_id'][self.sid]
        for role in ('series', 'samples'):
            prefix = '!Series_' if role == 'series' else '!Sample_'
            other = '!Sample_' if role == 'series' else '!Series_'
            raw = self.soft(self.sid, row['sample_ids'], role)
            field = prefix + 'contact_zip/postal_code = '
            replacements = [prefix + name + ' = ' for name in
                            ('/contact_zip', 'contact_zip/', 'contact_zip//postal_code')]
            replacements += [other + 'contact_zip/postal_code = ', prefix + 'contact_zip/postal_code= ']
            for replacement in replacements:
                with self.subTest(role=role, replacement=replacement):
                    with self.assertRaisesRegex(ValueError, 'Malformed SOFT field boundary'):
                        subject._readiness_soft(raw.replace(field, replacement, 1).encode(),
                                                role, self.sid, row['sample_ids'], row['sources'][role])

    def test_malformed_and_unlisted_soft_boundaries_rejected(self):
        row = self.context()['allowed_by_id'][self.sid]
        raw = self.soft(self.sid, row['sample_ids'], 'samples')
        changes = [raw.replace('!Sample_title = ', '!Sample_title=', 1),
                   raw.replace('^SAMPLE = ' + row['sample_ids'][0], '^SAMPLE = GSM999999999', 1),
                   raw.replace('!Sample_geo_accession = ' + row['sample_ids'][0], '!Sample_geo_accession = GSM999999999', 1),
                   raw + '^SAMPLE = ' + row['sample_ids'][0] + '\n']
        for altered in changes:
            with self.assertRaisesRegex(ValueError, 'SOFT|accession'):
                subject._readiness_soft(altered.encode(), 'samples', self.sid, row['sample_ids'], row['sources']['samples'])

    def test_series_missing_extra_or_duplicate_membership_is_rejected(self):
        row = self.context()['allowed_by_id'][self.sid]
        raw = self.soft(self.sid, row['sample_ids'], 'series')
        for altered in (raw.replace('!Series_sample_id = ' + row['sample_ids'][0] + '\n', ''),
                        raw.replace('!Series_sample_id = ' + row['sample_ids'][0], '!Series_sample_id = GSM99999999'),
                        raw + '!Series_sample_id = ' + row['sample_ids'][0] + '\n'):
            with self.assertRaisesRegex(ValueError, 'membership'):
                subject._readiness_soft(altered.encode(), 'series', self.sid, row['sample_ids'], row['sources']['series'])

    def test_intact_raw_record_permutations_preserve_identity_fields_and_literal_lines(self):
        row = self.allowed['studies'][0]
        sample_path = Path(row['sources']['samples']['path'])
        original = sample_path.read_bytes()
        blocks = [b'^SAMPLE = ' + block for block in original.split(b'^SAMPLE = ')[1:]]
        permuted = b''.join(reversed(blocks))
        before = subject._readiness_soft(original, 'samples', self.sid, row['sample_ids'], row['sources']['samples'])
        after = subject._readiness_soft(permuted, 'samples', self.sid, row['sample_ids'], row['sources']['samples'])
        identity = lambda rows: {(r['entity'], r['projection_field']): r['quote'] for r in rows}
        self.assertEqual(identity(before), identity(after))
        self.assertEqual(after[0]['entity'], row['sample_ids'][-1])
        self.assertEqual(after[0]['line_start'], 2)
        sample_path.write_bytes(permuted)
        series_path = Path(row['sources']['series']['path'])
        lines = series_path.read_bytes().splitlines(keepends=True)
        member_lines = iter(reversed([line for line in lines if line.startswith(b'!Series_sample_id = ')]))
        series_path.write_bytes(b''.join(next(member_lines) if line.startswith(b'!Series_sample_id = ') else line for line in lines))
        for role in ('series', 'samples'):
            path = Path(row['sources'][role]['path'])
            row['sources'][role].update(sha256=workflow.sha(path), bytes=path.stat().st_size)
        self.refresh_release()
        _, packet = self.project()
        self.assertEqual([r['sample_id'] for r in packet['source_packet']['samples']], row['sample_ids'])
        title = next(r for r in before if r['entity'] == row['sample_ids'][0] and r['field'] == '!Sample_title')
        self.assertEqual(packet['source_packet']['samples'][0]['title'], title['quote'])
        self.assertEqual(self.freeze(self.submission(comparison=True))['validation_status'], 'ACCEPTED_SOURCE_ATTRIBUTION')

    def two_comparison_scope_fixture(self):
        data = self.submission(comparison=True)
        first = data['comparisons'][0]
        second = copy.deepcopy(first)
        second['case_ids'] = first['case_ids'][:1]
        second['excluded_sample_ids'] = first['case_ids'][1:]
        second['exclusion_reason'] = 'SOFTWARE FIXTURE only the first case belongs to this additional contrast'
        second['conditional_test_family'] = 'CANNOT_SPECIFY'
        second['comparison_id'] = workflow.comparison_identity(self.sid, second['case_ids'], second['control_ids'], second['target'])
        data['comparisons'].append(second)
        return data, first['comparison_id'], second['comparison_id']

    def test_comparison_domain_rejects_fact_scoped_only_to_other_comparison(self):
        data, first, second = self.two_comparison_scope_fixture()
        next(f for f in data['facts'] if f['fact_id'] == 'f-unit-0')['comparison_ids'] = [second]
        self.assert_rejected(data, 'scoped to a different comparison')

    def test_conditional_requirement_rejects_foreign_fact_and_unresolved_scope(self):
        for number, kind in enumerate(('fact', 'unresolved')):
            data, first, second = self.two_comparison_scope_fixture()
            requirement = data['comparisons'][0]['analysis_requirements'][0]
            if kind == 'fact':
                foreign = copy.deepcopy(next(f for f in data['facts'] if f['fact_id'] == 'f-unit-0'))
                foreign.update(fact_id='f-foreign-requirement', comparison_ids=[second])
                data['facts'].append(foreign)
                requirement.update(fact_ids=[foreign['fact_id']], unresolved_ids=[])
            else:
                foreign = next(u for u in data['unresolved'] if u['unresolved_id'] == 'u-test_requirements')
                foreign['comparison_ids'] = [second]
                requirement.update(fact_ids=[], unresolved_ids=[foreign['unresolved_id']])
            self.assert_rejected(data, 'scoped to a different comparison', 'scope-requirement-%03d' % number)

    def test_pair_binding_rejects_fact_scoped_only_to_other_comparison(self):
        data, first, second = self.two_comparison_scope_fixture()
        foreign = copy.deepcopy(next(f for f in data['facts'] if f['fact_id'] == 'f-pair-0'))
        foreign.update(fact_id='f-foreign-pair', comparison_ids=[second])
        data['facts'].append(foreign)
        data['pair_map'][0]['fact_ids'] = [foreign['fact_id']]
        self.assert_rejected(data, 'scoped to a different comparison')

    def test_explicit_shared_comparison_scope_remains_valid(self):
        data, first, second = self.two_comparison_scope_fixture()
        for fact in data['facts']:
            fact['comparison_ids'] = [first, second]
        for unresolved in data['unresolved']:
            unresolved['comparison_ids'] = [first, second]
        self.assertEqual(self.freeze(data)['validation_status'], 'ACCEPTED_SOURCE_ATTRIBUTION')

    def test_schema_template_and_unknown_property_are_not_assessments(self):
        data = self.submission(); data['completed'] = False
        with patch.object(subject, '_readiness_raw', side_effect=AssertionError('must reject before raw')):
            self.assert_rejected(data, 'schema')
        data = self.submission(); data['family_id'] = 'INVENTED_SINGLETON'
        self.assert_rejected(data, 'schema', 'version-002')

    def test_packet_binding_injection_rejected_before_raw_read(self):
        binding, packet = self.project()
        packet['sources']['series']['path'] = '/NEVER_OPEN_OTHER_STUDY'
        forged = write_json(self.root / 'forged-packet.json', packet)
        with patch.object(subject, '_readiness_raw', side_effect=AssertionError('must reject before raw')):
            with self.assertRaisesRegex(ValueError, 'Unlisted packet'):
                subject.freeze_source_readiness_submission(self.release['path'], self.release['sha256'], self.batch,
                    self.sid, forged, '/NEVER_OPEN_SUBMISSION', 'version-001', self.root/'no-freeze')

    def test_witness_unlisted_path_or_hash_rejected_before_raw_read(self):
        for field, value in [('source_path', '/NEVER_OPEN_OTHER_STUDY'), ('source_sha256', '0'*64)]:
            data = self.submission(comparison=True)
            next(iter(data['evidence'].values()))[field] = value
            with patch.object(subject, '_readiness_raw', side_effect=AssertionError('must reject before raw')):
                self.assert_rejected(data, 'Unlisted witness', 'version-' + field)

    def test_full_multiline_literal_valid_and_substrings_or_crossing_boundaries_rejected(self):
        data = self.submission(comparison=True)
        self.assertEqual(self.freeze(data)['validation_status'], 'ACCEPTED_SOURCE_ATTRIBUTION')
        for number, mutation in enumerate(('substring', 'trim', 'line', 'field'), 2):
            changed = copy.deepcopy(data)
            witness = next(w for w in changed['evidence'].values() if w['field'] == '!Sample_description')
            if mutation == 'substring': witness['quote'] = 'SOFTWARE FIXTURE first line'
            elif mutation == 'trim': witness['quote'] = witness['quote'].replace('  ', ' ')
            elif mutation == 'line': witness['line_end'] += 1
            else: witness['field'] = '!Sample_title'
            self.assert_rejected(changed, 'complete literal field', 'version-%03d' % number)

    def test_unknown_cross_reference_rejected_before_raw(self):
        data = self.submission()
        data['domains']['pairing']['unresolved_ids'] = ['missing-id']
        with patch.object(subject, '_readiness_raw', side_effect=AssertionError('must reject before raw')):
            self.assert_rejected(data, 'Unknown unresolved_ids')

    def test_duplicate_fact_and_unresolved_ids_rejected(self):
        data = self.submission(comparison=True)
        data['facts'].append(copy.deepcopy(data['facts'][0]))
        self.assert_rejected(data, 'Duplicate fact_id')
        data = self.submission(); data['unresolved'].append(copy.deepcopy(data['unresolved'][0]))
        self.assert_rejected(data, 'Duplicate unresolved_id', 'version-002')

    def test_domain_cross_reference_cannot_borrow_other_domain_fact(self):
        data = self.submission(comparison=True)
        data['domains']['pairing']['fact_ids'] = ['f-age']
        self.assert_rejected(data, 'different domain')

    def test_fail_requires_source_fact_and_explicit_materiality(self):
        data = self.submission()
        data['domains']['pairing']['status'] = 'FAIL'
        self.assert_rejected(data, 'schema')
        data = self.submission(comparison=True)
        data['domains']['pairing']['status'] = 'FAIL'; data['domains']['pairing']['material'] = None
        self.assert_rejected(data, 'materiality', 'version-002')

    def test_comparison_arm_witness_and_stable_identity_guards(self):
        for number, mutation in enumerate(('overlap', 'missing-witness', 'identity', 'excluded')):
            data = self.submission(comparison=True); comparison = data['comparisons'][0]
            if mutation == 'overlap': comparison['case_ids'].append(comparison['control_ids'][0])
            elif mutation == 'missing-witness': comparison['evidence_ids'] = comparison['evidence_ids'][:1]
            elif mutation == 'identity': comparison['comparison_id'] = 'invented-name'
            else: comparison['excluded_sample_ids'] = comparison['case_ids'][:1]
            self.assert_rejected(data, 'Comparison identity|Pair map', 'version-%03d' % number)

    def test_shared_control_membership_is_preserved(self):
        data = self.submission(comparison=True)
        original = data['comparisons'][0]
        second = copy.deepcopy(original)
        second['case_ids'] = original['case_ids'][:1]
        second['conditional_test_family'] = 'CANNOT_SPECIFY'
        second['comparison_id'] = workflow.comparison_identity(self.sid, second['case_ids'], second['control_ids'])
        data['comparisons'].append(second)
        result = self.freeze(data)
        self.assertEqual(result['validation_status'], 'ACCEPTED_SOURCE_ATTRIBUTION')
        frozen = json.loads(Path(result['frozen_submission']).read_text())
        self.assertEqual(frozen['validation']['comparisons'], 2)

    def test_unassigned_and_unknown_units_keep_explicit_dependencies(self):
        data = self.submission(); data['unassigned_sample_ids'] = []
        self.assert_rejected(data, 'unassigned sample accounting')
        data = self.submission()
        next(r for r in data['unresolved'] if r['domain'] == 'biological_units')['sample_ids'] = []
        self.assert_rejected(data, 'biological_units dependencies', 'version-002')

    def test_pair_wrong_arm_and_wrong_unit_fail(self):
        for number, mutation in enumerate(('arm', 'unit', 'reuse', 'missing')):
            data = self.submission(comparison=True)
            pair = data['pair_map'][0]
            if mutation == 'arm': pair['case_sample_id'] = pair['control_sample_id']
            elif mutation == 'unit': pair['unit_id'] = 'unit-2'
            elif mutation == 'reuse': data['pair_map'].append(copy.deepcopy(pair))
            else: data['pair_map'] = []
            self.assert_rejected(data, 'Pair map|paired conditional plan', 'version-%03d' % number)

    def test_title_suffix_cannot_supply_unit_or_pair_key_evidence(self):
        for number, target in enumerate(('unit_map', 'pair_map')):
            data = self.submission(comparison=True)
            titles = [key for key, value in data['evidence'].items() if value['field'] == '!Sample_title']
            data[target][0]['evidence_ids'] = titles
            self.assert_rejected(data, 'Unit map|Pair map|omits its source-fact', 'version-%03d' % number)

    def test_unit_and_pair_assertions_cannot_exchange_valid_but_unrelated_witnesses(self):
        for number, target in enumerate(('unit_map', 'pair_map')):
            data = self.submission(comparison=True)
            ages = [key for key, value in data['evidence'].items() if value['field'] == '!Sample_characteristics_ch1'
                    and value['quote'].startswith('age:')]
            data[target][0]['evidence_ids'] = ages
            self.assert_rejected(data, 'omits its source-fact', 'version-%03d' % number)

    def test_unknown_collaboration_model_is_recorded_as_null_not_an_api_model(self):
        self.release_data['reference_runtime'] = {'runtime': 'SOFTWARE FIXTURE collaboration-runtime', 'model': None}
        self.release = write_json(Path(self.release['path']), self.release_data)
        data = self.submission(); data['reviewer']['runtime_model_as_reported'] = None
        self.assertEqual(self.freeze(data)['validation_status'], 'ACCEPTED_SOURCE_ATTRIBUTION')

    def test_unsupported_covariate_and_method_requirements_fail(self):
        data = self.submission(comparison=True); data['covariates'][0]['fact_ids'] = []
        self.assert_rejected(data, 'source facts')
        data = self.submission(comparison=True)
        data['comparisons'][0]['analysis_requirements'][0].update(fact_ids=[], unresolved_ids=[])
        self.assert_rejected(data, 'source facts', 'version-002')

    def test_family_identifier_and_unresolved_closure_are_bound(self):
        data = self.submission(comparison=True); data['family_links'][0]['target_identifier'] = 'NOT_IN_SOURCE'
        self.assert_rejected(data, 'Family identifier')
        data = self.submission(); data['family_closure'] = 'ALL_SINGLETONS'
        self.assert_rejected(data, 'schema', 'version-002')

    def test_projection_resume_checks_raw_hash_without_overwriting(self):
        binding, _ = self.project(); before = workflow.sha(binding['path'])
        result = subject.prepare_source_readiness_batch(self.release['path'], self.release['sha256'], self.batch, self.root/'projection')
        self.assertEqual(result, self.projection)
        source = Path(self.allowed['studies'][0]['sources']['series']['path'])
        source.write_bytes(source.read_bytes() + b'SOFTWARE FIXTURE changed input\n')
        with self.assertRaisesRegex(ValueError, 'size/hash changed'):
            subject.prepare_source_readiness_batch(self.release['path'], self.release['sha256'], self.batch, self.root/'projection')
        self.assertEqual(workflow.sha(binding['path']), before)

    def test_source_unavailable_is_preserved_rejection_not_scientific_cv(self):
        data = self.submission()
        source = Path(self.allowed['studies'][0]['sources']['series']['path'])
        source.rename(source.with_suffix('.preserved-software-fixture'))
        self.assert_rejected(data, 'Cannot proceed - required data file not found:')

    def test_freeze_resume_identical_and_changed_submission_needs_new_version(self):
        data = self.submission(); result = self.freeze(data)
        self.assertEqual(self.freeze(data), result)
        checksum = workflow.sha(result['frozen_submission'])
        data['reviewer']['known_prior_exposure'] += ' changed software fixture declaration'
        with self.assertRaisesRegex(ValueError, 'changed; use a new run directory'):
            self.freeze(data)
        self.assertEqual(workflow.sha(result['frozen_submission']), checksum)
        newer = self.freeze(data, 'version-002', supersedes=self.frozen_binding(result))
        self.assertNotEqual(result['sha256'], newer['sha256'])

    def test_rejected_then_corrected_versions_and_all_denominators_survive(self):
        bad = self.submission(); bad['completed'] = False
        rejected = self.assert_rejected(bad, 'schema')
        correct = self.freeze(self.submission(), 'version-002', supersedes=self.frozen_binding(rejected))
        bindings = [self.frozen_binding(rejected), self.frozen_binding(correct)]
        with patch.object(subject, '_readiness_raw', side_effect=AssertionError('report must not read raw')):
            result = subject.report_source_readiness(self.release['path'], self.release['sha256'], bindings, self.root/'counts')
        table = json.loads(Path(result['table']).read_text())
        self.assertEqual(len(table['studies']), 80)
        self.assertEqual((table['allowed_domain_slots'], table['held_domain_slots'], table['enrolled_domain_slots']), (702,18,720))
        self.assertEqual(len(table['studies'][0]['versions']), 2)
        for count in table['domains'].values():
            self.assertEqual((count['CANNOT VERIFY'], count['NOT_ASSESSED'], count['administratively_held']), (1,77,2))
            self.assertEqual(count['source_validated_coverage_enrolled'], {'numerator':1,'denominator':80})
        self.assertEqual({r['administrative_state'] for r in table['studies'] if not r['source_access_allowed']}, set(subject.SOURCE_READINESS_HOLDS.values()))
        self.assertNotIn('accepted_result_error_rate', table)
        self.assertNotIn('scientific_reference_complete', json.dumps(table))

    def test_latest_invalid_submission_returns_to_not_assessed_without_erasing_history(self):
        good = self.freeze(self.submission())
        bad = self.submission(); bad['completed'] = False
        rejected = self.freeze(bad, 'version-002', supersedes=self.frozen_binding(good))
        result = subject.report_source_readiness(self.release['path'], self.release['sha256'],
            [self.frozen_binding(good), self.frozen_binding(rejected)], self.root/'counts')
        table = json.loads(Path(result['table']).read_text())
        self.assertEqual(table['domains']['sample_identity']['NOT_ASSESSED'], 78)
        self.assertEqual(len(table['studies'][0]['versions']), 2)
        self.assertIn('schema', table['studies'][0]['not_assessed_reason'])

    def test_readiness_requires_previous_versions_and_intact_outputs(self):
        result = self.freeze(self.submission())
        newer = self.freeze(self.submission(), 'version-002', supersedes=self.frozen_binding(result))
        with self.assertRaisesRegex(ValueError, 'complete ordered supersession chain'):
            subject.report_source_readiness(self.release['path'], self.release['sha256'], [self.frozen_binding(newer)], self.root/'bad-counts')
        source = Path(result['frozen_submission']).parent/'submission.json'
        source.write_text(source.read_text() + ' ')
        with self.assertRaisesRegex(ValueError, 'output changed'):
            subject.report_source_readiness(self.release['path'], self.release['sha256'], [self.frozen_binding(result)], self.root/'tampered-counts')

    def test_duplicate_json_keys_are_rejected_and_original_bytes_retained(self):
        binding, _ = self.project()
        raw = '{"completed": true, "completed": false}\n'
        path = self.root/'duplicate-keys.json'; path.write_text(raw)
        result = subject.freeze_source_readiness_submission(self.release['path'],self.release['sha256'],self.batch,self.sid,
            binding,path,'version-001',self.root/'version-001')
        frozen = Path(result['frozen_submission'])
        self.assertEqual((frozen.parent/'submission.json').read_text(), raw)
        self.assertIn('Duplicate JSON key', json.loads(frozen.read_text())['validation']['reason'])

    def test_metadata_only_initial_denominator_has_no_scientific_verdicts(self):
        with patch.object(subject, '_readiness_raw', side_effect=AssertionError('must not read raw')):
            table = subject.source_readiness_table(self.enrollment)
        self.assertEqual(table['enrolled_source_records'], 320)  # SOFTWARE FIXTURE only.
        self.assertEqual(len(table['studies']), 80)
        for counts in table['domains'].values():
            self.assertEqual((counts['PASS'], counts['FAIL'], counts['CANNOT VERIFY'], counts['NOT_ASSESSED']), (0,0,0,78))
            self.assertEqual(counts['resolved_coverage_allowed'], {'numerator':0,'denominator':78})


if __name__ == '__main__':
    unittest.main()
