"""Assemble source-bound literature records without estimating missing statistics."""
import argparse
import csv
import hashlib
import json
import logging
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path

import openpyxl

from benchmark.published_clock_2026_10_05.curated.retrieve_2026_10_05 import ROOT


ATLAS = Path('/Entropy/kying/clockbase-verifier-v3-2026-09-27/final-2026-10-05/'
             'atlas-v3/atlas_v3_combined.csv')
SOURCE_DEFINITIONS = {
    'tysh2023': ('10.1016/j.cell.2023.05.002', '37269831', 'PMC11192172'),
    'tysh2019': ('10.1016/j.cmet.2019.06.018', '31353263', 'PMC6907080'),
    'computage': ('10.1101/2024.06.06.597715', None, None),
    'biolearn': ('10.1038/s43587-025-00987-y', '41188602', 'PMC13476091'),
    'down': ('10.1111/acel.12325', '25678027', 'PMC4406678'),
    'hivhorvath': ('10.1093/infdis/jiv277', '25969563', 'PMC4621253'),
    'hivgross': ('10.1016/j.molcel.2016.03.019', '27105112', 'PMC4995115'),
    'obesity': ('10.1073/pnas.1412759111', '25313081', 'PMC4217403'),
    'huntington': ('10.18632/aging.101005', '27479945', 'PMC4993344'),
    'parkinson': ('10.18632/aging.100859', '26655927', 'PMC4712337'),
    'werner': ('10.18632/aging.101217', '28377537', 'PMC5425119'),
    'phenoage': ('10.18632/aging.101414', '29676998', 'PMC5940111'),
    'grimage': ('10.18632/aging.101684', '30669119', 'PMC6366976'),
    'pace': ('10.7554/elife.73420', '35029144', 'PMC8853656'),
}
CLOCK_MAP = {'Horvath 2013 (353-CpG epigenetic clock)': 'Horvathv1',
             'Horvath model': 'Horvathv1', 'Horvathv1': 'Horvathv1',
             'DNAm PhenoAge': 'PhenoAge', 'PhenoAge': 'PhenoAge',
             'DunedinPACE': 'DunedinPACE', 'YingDamAge': 'YingDamAge'}


def normalized(text):
    return re.sub(r'\s+', ' ', text).strip()


class DocumentParser(HTMLParser):
    """Keep paragraph and table structure in the cached PMC HTML."""
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.root = ET.Element('document')
        self.stack = [self.root]

    def handle_starttag(self, tag, attrs):
        node = ET.SubElement(self.stack[-1], tag, dict(attrs))
        if tag not in {'area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input',
                       'link', 'meta', 'param', 'source', 'track', 'wbr'}:
            self.stack.append(node)

    def handle_startendtag(self, tag, attrs):
        ET.SubElement(self.stack[-1], tag, dict(attrs))

    def handle_endtag(self, tag):
        for i in range(len(self.stack) - 1, 0, -1):
            if self.stack[i].tag == tag:
                self.stack = self.stack[:i]
                break

    def handle_data(self, data):
        node = self.stack[-1]
        if len(node):
            node[-1].tail = (node[-1].tail or '') + data
        else:
            node.text = (node.text or '') + data


def document(path):
    if path.suffix == '.xml':
        return ET.parse(path).getroot()
    parser = DocumentParser()
    parser.feed(path.read_text())
    return parser.root


def text_of(node):
    return normalized(''.join(node.itertext()))


def evidence(path, selector, text):
    return {'file': path.name, 'selector': selector, 'text': text}


def paragraph(path, containing, excerpt=None):
    tree = document(path)
    for i, node in enumerate(tree.findall('.//p'), 1):
        value = text_of(node)
        if containing in value:
            quote = excerpt if excerpt is not None else value
            if quote not in value:
                raise ValueError(f'Quote is absent from {path.name}: {quote}')
            return evidence(path, f'p[{i}] contains {containing!r}', quote)
    raise ValueError(f'Paragraph is absent from {path.name}: {containing}')


def table(path, table_id):
    tree = document(path)
    for node in tree.iter():
        if node.tag in {'table-wrap', 'table'} and node.get('id') == table_id:
            return [[text_of(cell) for cell in row if cell.tag in {'td', 'th'}]
                    for row in node.findall('.//tr')]
    raise ValueError(f'Table is absent from {path.name}: {table_id}')


def statistic(name, value, units, relation='=', method=None):
    return {'name': name, 'value': value, 'relation': relation,
            'units': units, 'method': method}


def record(source, gse, tissue, clock, case, control, kind, location, quotes,
           stats=None, direction=None, n_case=None, n_control=None, **extra):
    doi, pmid, pmcid = SOURCE_DEFINITIONS[source]
    row = {'source_id': source, 'doi': doi, 'pmid': pmid, 'pmcid': pmcid,
           'source_location': location, 'gse_id': gse, 'modality': 'human DNAm',
           'tissue': tissue, 'clock_source': clock, 'clock_clockbase': CLOCK_MAP.get(clock),
           'case_definition': case, 'control_definition': control,
           'n_case': n_case, 'n_control': n_control, 'n_total': None,
           'ground_truth_type': kind, 'reported_statistics': stats or [],
           'reported_direction': direction, 'expected_direction': None,
           'quote': '\n'.join(q['text'] for q in quotes), 'quote_evidence': quotes,
           'notes': []}
    row.update(extra)
    row['statistic_fields'] = {name: next((s for s in row['reported_statistics'] if s['name'] == name), None)
                              for name in ['effect', 'mean_difference', 'log_fold_change', 'beta', 'p_value', 'fdr']}
    return row


def classic_records(cache):
    result = []
    path = cache / 'down.xml'
    rows = table(path, 'tbl1')[1:]
    regression = table(path, 'tbl2')
    ds = next(r for r in regression if r[0] == 'DS')
    acceleration = next(r for r in regression if r[0] == 'Age accel. for DS')
    pvalues = [0.0048, 5.5e-7, 0.00056, 0.26]
    beta = [3.74593, 8.53045, 3.65564, 2.6495]
    beta_se = [1.37252, 1.23934, 1.08842, 1.7776]
    beta_p = [0.0086, 2.3e-9, 0.0012, 0.15]
    years = [3.9, 11.5, 4.6, 2.8]
    for i, cells in enumerate(rows):
        gse = next(c for c in cells if c.startswith('GSE'))
        tissue = ['leukocytes', 'brain', 'whole blood', 'buccal epithelium'][i]
        q = evidence(path, f'tbl1 row {i + 2}', ' | '.join(cells))
        counts = {'n_case': int(cells[4]), 'n_total': int(cells[2]),
                  'notes': ['Control count is not filled by subtracting cases from total.',
                            'Age acceleration is a residual relative to the control age regression.']}
        result.append(record('down', gse, tissue, 'Horvath 2013 (353-CpG epigenetic clock)',
                             'Down syndrome (trisomy 21)', 'Non-Down-syndrome samples in this dataset',
                             'reported_case_control_statistic', 'Table 1; Figure 1e-h', [q],
                             [statistic('p_value', pvalues[i], 'probability', method='Kruskal-Wallis')],
                             'case higher' if i < 3 else 'no significant difference', **counts))
        q2 = evidence(path, f'tbl2 DS and Age accel. for DS, dataset {i + 1}',
                      ' | '.join([ds[0], ds[1 + i * 2], ds[2 + i * 2],
                                  acceleration[0], acceleration[1 + i]]))
        result.append(record('down', gse, tissue, 'Horvath 2013 (353-CpG epigenetic clock)',
                             'Down syndrome (trisomy 21)', 'Non-Down-syndrome samples in this dataset',
                             'reported_case_control_statistic', 'Table 2, dataset ' + str(i + 1), [q, q2],
                             [statistic('beta', beta[i], 'DNAm years', method='DNAm age ~ age + DS'),
                              statistic('beta_standard_error', beta_se[i], 'DNAm years'),
                              statistic('p_value', beta_p[i], 'probability', method='Wald test for DS coefficient'),
                              statistic('age_acceleration', years[i], 'years', method='Reported DS/age coefficient ratio')],
                             'case higher' if i < 3 else 'no significant difference', **counts))

    path = cache / 'hivhorvath.xml'
    mapping = table(path, 'JIV277TB1')
    for dataset, gse, tissue, pvalue, effect, fragment, location, cases in [
        (1, 'GSE59457', 'multiple brain regions', 8e-8, None,
         'On average, brain samples from cases exhibited significant age acceleration effects', 'Figure 1E', 99),
        (3, 'GSE67748', 'cerebellum', 0.031, 7.4,
         'Using brain data set 3, we observed again', 'Figure 2C; Supplementary 1', 8),
        (4, 'GSE53841', 'blood (PBMCs in Methods; whole blood in Table 1)', 0.0048, 3.7,
         'Using 2 independent blood data sets', 'Figure 3C; Supplementary 1', 24),
        (5, 'GSE67751', 'whole blood', 0.00036, 6.7,
         'Using 2 independent blood data sets', 'Figure 3D; Supplementary 1', 23),
    ]:
        link = next(r for r in mapping if r[0] == str(dataset))
        q = evidence(path, 'Table 1 dataset ' + str(dataset), ' | '.join(link))
        q2 = paragraph(path, fragment)
        stats = [statistic('p_value', pvalue, 'probability', method='Kruskal-Wallis on age acceleration')]
        if effect is not None:
            stats.append(statistic('age_acceleration', effect, 'years', method='Reported linear-model age acceleration'))
        result.append(record('hivhorvath', gse, tissue, 'Horvath 2013 (353-CpG epigenetic clock)',
                             'HIV-infected subjects', 'HIV-uninfected subjects',
                             'reported_case_control_statistic', 'Table 1; ' + location, [q, q2], stats,
                             'case higher', n_case=cases,
                             n_control=68 if dataset == 4 else None, n_total=int(link[4]),
                             notes=['P value and age difference are distinct source analyses.',
                                    'Dataset 4 tissue labels conflict between Table 1 and Methods.']))
    for dataset, gse, tissue in [(1, 'GSE59457', 'frontal lobe'), (2, 'GSE67749', 'frontal lobe')]:
        q = paragraph(path, 'But we did not observe an accelerated aging effect in the second frontal lobe data set')
        link = next(r for r in mapping if r[0] == str(dataset))
        result.append(record('hivhorvath', gse, tissue, 'Horvath 2013 (353-CpG epigenetic clock)',
                             'HIV-infected subjects', 'HIV-uninfected subjects',
                             'reported_case_control_direction', 'Figures 1H and 2D; Table 1',
                             [evidence(path, 'Table 1 dataset ' + str(dataset), ' | '.join(link)), q],
                             direction='no significant difference'))

    path = cache / 'hivgross-pmc.html'
    for tissue, clock, value, pvalue, relation, fragment, location, ci in [
        ('whole blood', 'Consensus Hannum and Horvath model', 4.9, 1e-8, '<',
         'In contrast, the HIV+ patients had a biological age advancement of 4.9 years', 'Figure 2E-F', [3.4, 7.1]),
        ('neutrophils', 'Hannum model', 2.5, 0.03, '<',
         'In neutrophils, the Hannum model predicted', 'Figure 3E', [0.6, 5.0]),
        ('neutrophils', 'Horvath model', 0.4, 0.05, '>',
         'In neutrophils, the Hannum model predicted', 'Results, validation in purified cells', None),
        ('CD4+ T cells', 'Consensus Hannum and Horvath model', 5.7, 1e-5, '<',
         'In neutrophils, the Hannum model predicted', 'Figure 3F', [3.4, 7.9]),
    ]:
        q = paragraph(path, fragment)
        link = paragraph(path, 'All methylation data for HIV+ and HIV− subjects were deposited')
        stats = [statistic('age_advancement', value, 'years'),
                 statistic('p_value', pvalue, 'probability', relation,
                           method='Student t-test' if tissue == 'whole blood' else None)]
        if ci is not None:
            stats.append(statistic('confidence_interval_95', ci, 'years'))
        result.append(record('hivgross', 'GSE67705', tissue, clock,
                             ('HIV+ subjects on sustained antiretroviral therapy' if tissue == 'whole blood'
                              else 'HIV+ subjects in the sorted-cell validation cohort'),
                             'HIV− matched controls', 'reported_case_control_statistic', location,
                             [q, link], stats,
                             'no significant difference' if relation == '>' else 'case higher',
                             n_case=137 if tissue == 'whole blood' else None,
                             n_control=44 if tissue == 'whole blood' else None,
                             notes=['Whole-blood counts describe recruited subjects, not an inferred post-QC analysis count.',
                                    'The consensus model is not mapped to either component clock.']))

    path = cache / 'obesity-pmc.html'
    for gse, pvalue, total, fragment in [
        ('GSE48325', 0.00068, 62, 'yielding ∼2.2 y of additional DNAm age'),
        ('GSE61258', 0.00012, 79, 'To replicate the BMI-related epigenetic age acceleration in liver')]:
        q = paragraph(path, fragment)
        rows = [r for node in document(path).findall('.//table')
                for r in node.findall('.//tr') if gse in text_of(r)]
        if not rows:
            raise ValueError(f'GEO link absent from obesity Table 1: {gse}')
        link = evidence(path, 'Table 1 row ' + gse, ' | '.join(text_of(c) for c in rows[0]))
        result.append(record('obesity', gse, 'liver', 'Horvath 2013 (353-CpG epigenetic clock)',
                             'Higher body mass index (continuous exposure)',
                             'Lower body mass index (continuous exposure, no binary control arm)',
                             'reported_continuous_association', 'Table 1; Figures 1E and 2E-F', [link, q],
                             [statistic('pearson_r', 0.42, 'dimensionless'),
                              statistic('p_value', pvalue, 'probability', method='Pearson correlation')],
                             n_total=total, association_direction='positive',
                             notes=['Not a reported obese-versus-lean mean difference; no binary arm sizes assigned.']))

    path = cache / 'huntington.xml'
    rows = table(path, 'T2')
    row = next(r for r in rows if r[0] == 'Huntington')
    link = paragraph(path, 'Our novel DNA methylation data have been posted')
    for i, (beta, se, pvalue) in enumerate([(2.06, 0.517, 6.7e-5), (1.704, 0.503, 0.00070), (0.9, 0.486, 0.065)], 1):
        q = evidence(path, f'Table 2 Huntington row, Model {i}',
                     ' | '.join([row[0], row[1 + (i - 1) * 2], row[2 + (i - 1) * 2]]))
        result.append(record('huntington', 'GSE72778', 'multiple brain regions',
                             'Horvath 2013 (353-CpG epigenetic clock)', 'Huntington disease brain samples',
                             'Non-HD brain samples in the multivariate analysis (Table 1 also contains AD samples)',
                             'reported_case_control_statistic', f'Table 2, Model {i}', [link, q],
                             [statistic('beta', beta, 'DNAm years'), statistic('beta_standard_error', se, 'DNAm years'),
                              statistic('p_value', pvalue, 'probability')],
                             'case higher' if i < 3 else 'no significant difference',
                             notes=['The source treats P values as descriptive because brain regions repeat within individuals.',
                                    'Model 1 age-increase text is inconsistent (3.3 and 3.2 years); only the reported coefficient is extracted.',
                                    'Table 1 donor/sample counts are not assumed to equal the Model 2 analysis counts.']))

    path = cache / 'parkinson.xml'
    q = paragraph(path, 'To estimate the actual amount of age acceleration')
    rows = table(path, 'T2')
    covariate = next(r for r in rows if r[0] == 'Measure of Age Acceleration')
    q2 = evidence(path, 'Table 2 AgeAccel model, Measure of Age Acceleration', ' | '.join(covariate[:4]))
    sizes = next(r for r in table(path, 'T1') if r[0] == 'Sample size')
    counts_quote = evidence(path, 'Table 1, Sample size row, Caucasian and Hispanic cohorts', ' | '.join(sizes))
    counts = {'GSE72774': {'n_case': 289, 'n_control': 219},
              'GSE72776': {'n_case': 46, 'n_control': 38}}
    geo_links = []
    for gse in ['GSE72774', 'GSE72776']:
        geo_path = cache / (gse + '-soft')
        text = geo_path.read_text()
        if '!Series_pubmed_id = 26655927' not in text:
            raise ValueError(f'GEO does not cite the PD paper: {gse}')
        geo_links.append(evidence(geo_path, 'GEO SOFT Series_pubmed_id', '!Series_pubmed_id = 26655927'))
    for gse in ['GSE72774', 'GSE72776']:
        result.append(record('parkinson', gse, 'whole blood', 'Horvath 2013 (353-CpG epigenetic clock)',
                             'PD cases from PEG wave 1, Caucasian and Hispanic cohorts pooled',
                             'PEG controls from both cohorts pooled', 'reported_pooled_case_control_statistic',
                             'Results, multivariate linear regression; Table 1', [q, counts_quote] + geo_links,
                             [statistic('age_acceleration', 1.5, 'years', method='Linear model adjusts disease status, age, granulocytes, smoking, ethnicity and sex')],
                             direction='case higher', pooled_gse_ids=['GSE72774', 'GSE72776'],
                             pooling_id='parkinson-pooled-linear-model',
                             per_series_reported_counts=counts,
                             notes=['The pooled effect is not a dataset-specific estimate; do not count the two accession links as independent results.']))
        result.append(record('parkinson', gse, 'whole blood', 'Horvath 2013 (353-CpG epigenetic clock)',
                             'PD cases from PEG wave 1, Caucasian and Hispanic cohorts pooled',
                             'PEG controls from both cohorts pooled', 'reported_pooled_case_control_statistic',
                             'Table 2, AgeAccel logistic model', [q2, counts_quote] + geo_links,
                             [statistic('beta', 0.036, 'log odds of PD per year age acceleration'),
                              statistic('beta_standard_error', 0.017, 'log odds of PD per year age acceleration'),
                              statistic('p_value', 0.037, 'probability', method='Adjusted logistic model with PD as outcome')],
                             direction='case higher', pooled_gse_ids=['GSE72774', 'GSE72776'],
                             pooling_id='parkinson-pooled-logistic-model',
                             per_series_reported_counts=counts,
                             notes=['Outcome is PD, not DNAm age; beta is not an age difference.',
                                    'No per-GEO inference from the pooled regression.']))
    return result


def pdf_excerpt(path, fragment, length=400):
    """Quote an excerpt of the cached pdftotext output and retain its page."""
    for page, content in enumerate(path.read_text().split('\f'), 1):
        text = normalized(content)
        start = text.find(fragment)
        if start >= 0:
            return evidence(path, f'Supplement 1 PDF page {page} (1-based), text extract', text[start:start + length])
    raise ValueError(f'PDF excerpt is absent: {fragment}')


def phenoage_records(cache):
    path = cache / 'phenoage-supplement1.txt'
    result = []
    ds = pdf_excerpt(path, 'We also find that DNAm PhenoAge relates to Down syndrome', 600)
    down_link = table(cache / 'down.xml', 'tbl1')
    for gse, tissue, pvalue, total, cases, controls in [
        ('GSE25395', 'leukocytes', 0.0046, 56, 35, 21),
        ('GSE52588', 'whole blood', 4e-11, 87, 29, None)]:
        link = next(r for r in down_link if gse in r)
        counts = pdf_excerpt(path, 'We previously reported that individuals with Down syndrome', 1000)
        result.append(record('phenoage', gse, tissue, 'DNAm PhenoAge', 'Down syndrome',
                             'Non-DS controls; dataset 2 contains mothers and unaffected siblings',
                             'reported_case_control_statistic', 'Supplement 1, pp. 6 and 8, Down syndrome applications',
                             [ds, counts, evidence(cache / 'down.xml', 'Cited Horvath 2015 Table 1 ' + gse, ' | '.join(link))],
                             [statistic('p_value', pvalue, 'probability')], 'case higher',
                             n_case=cases, n_control=controls, n_total=total,
                             notes=['GSE identity is cross-linked through the explicitly cited Horvath 2015 dataset and total sample count.',
                                    'No unreported PhenoAge mean difference or test method is filled.']))
    link = next(r for r in down_link if 'GSE63347' in r)
    result.append(record('phenoage', 'GSE63347', 'brain', 'DNAm PhenoAge', 'Down syndrome brain samples',
                         'Controls as stated in Supplement 1', 'reported_case_control_direction',
                         'Supplement 1, p. 6, Down syndrome applications',
                         [pdf_excerpt(path, 'PhenoAge did not reveal a significant epigenetic age acceleration effect', 155),
                          evidence(cache / 'down.xml', 'Cited Horvath 2015 Table 1 GSE63347', ' | '.join(link))],
                         direction='no significant difference', n_case=15, n_control=54,
                         notes=['Levine reports 15 DS/54 controls here; the earlier Horvath paper reports 71 total brain samples. Counts are not reconciled by inference.']))
    hiv = pdf_excerpt(path, 'We find that HIV infection is associated', 250)
    hiv_link = table(cache / 'hivhorvath.xml', 'JIV277TB1')
    counts = pdf_excerpt(path, 'The first dat set involved Peripheral blood mononuclear cells', 700)
    for dataset, gse, pvalue, cases, controls in [(4, 'GSE53841', 6e-6, 24, 68), (5, 'GSE67751', 8.6e-6, 23, 69)]:
        link = next(r for r in hiv_link if r[0] == str(dataset))
        result.append(record('phenoage', gse, 'PBMCs (as named in Levine Supplement 1)', 'DNAm PhenoAge',
                             'HIV-infected subjects', 'HIV-uninfected subjects', 'reported_case_control_statistic',
                             'Supplement 1, pp. 6 and 8, HIV applications',
                             [hiv, counts, evidence(cache / 'hivhorvath.xml', 'Cited Horvath 2015 Table 1, dataset ' + str(dataset), ' | '.join(link))],
                             [statistic('p_value', pvalue, 'probability')], 'case higher',
                             n_case=cases, n_control=controls, n_total=92,
                             notes=['GSE and P-value order follow the cited datasets and the explicitly stated case/control counts.',
                                    'No PhenoAge effect size is stated in this passage.']))
    pd = pdf_excerpt(path, "Parkinson's disease status (p=0.028) the n=508 individuals", 220)
    geo_path = cache / 'GSE72774-soft'
    result.append(record('phenoage', 'GSE72774', 'whole blood', 'DNAm PhenoAge',
                         'European-ancestry PD cases in PEG', 'European-ancestry PEG controls',
                         'reported_case_control_statistic', 'Supplement 1, pp. 6 and 8, PD application',
                         [pd, evidence(geo_path, 'Cited paper link, Series_pubmed_id', '!Series_pubmed_id = 26655927')],
                         [statistic('p_value', 0.028, 'probability')], n_case=289, n_control=219, n_total=508,
                         notes=['The passage reports a relationship but does not state its direction; direction remains null.']))
    obesity = pdf_excerpt(path, 'Using the same data, we replicated this finding using the new measure of', 410)
    obesity_path = cache / 'obesity-pmc.html'
    for gse, r, pvalue in [('GSE48325', 0.32, 0.011), ('GSE61258', 0.48, 7.7e-6)]:
        links = [row for node in document(obesity_path).findall('.//table')
                 for row in node.findall('.//tr') if gse in text_of(row)]
        result.append(record('phenoage', gse, 'liver', 'DNAm PhenoAge', 'Higher BMI (continuous exposure)',
                             'Lower BMI (continuous exposure)', 'reported_continuous_association',
                             'Supplement 1, p. 8, Effect of obesity on liver and adipose tissue',
                             [obesity, evidence(obesity_path, 'Cited Horvath 2014 Table 1 ' + gse, ' | '.join(text_of(c) for c in links[0]))],
                             [statistic('pearson_r', r, 'dimensionless'), statistic('p_value', pvalue, 'probability')],
                             association_direction='positive', notes=['Dataset 1/2 mapping follows the explicitly cited Horvath 2014 Table 1.']))
    heading = pdf_excerpt(path, 'Supplementary Table S5F. Body Mass Index versus epigenetic age acceleration in liver', 700)
    for label, clock, beta, se, pvalue, fragment in [
        ('DNAmAge', 'Horvath 2013 (353-CpG epigenetic clock)', 0.15, 0.042, 4.2e-4, 'BMI 0.15 0.042 3.7 4.2e-04'),
        ('DNAmPhenoAge', 'DNAm PhenoAge', 0.13, 0.051, 0.01, 'BMI 0.13 0.051 2.60 1.0e-02')]:
        q = pdf_excerpt(path, fragment, len(fragment))
        result.append(record('phenoage', 'GSE48325', 'liver', clock, 'Higher BMI (continuous exposure)',
                             'Lower BMI (continuous exposure)', 'reported_continuous_association',
                             'Supplement 1, Table S5F liver (second table with this label), pp. 18-19, ' + label,
                             [heading, q], [statistic('beta', beta, 'DNAm years per BMI unit (kg/m²)'),
                                           statistic('beta_standard_error', se, 'DNAm years per BMI unit (kg/m²)'),
                                           statistic('p_value', pvalue, 'probability', method='Age- and sex-adjusted linear model')],
                             n_total=85, association_direction='positive',
                             notes=['Supplement reports n=85 arrays; this differs from n=62 in the 2014 discovery paper and is not reconciled.']))
    return result


def computage_records(cache):
    path = cache / 'hf-meta-pinned.tsv'
    with path.open(newline='') as handle:
        metadata = list(csv.DictReader(handle, delimiter='\t'))
    groups = defaultdict(list)
    for index, row in enumerate(metadata, 2):
        groups[(row['DatasetID'], row['Tissue'], row['Condition'])].append((index, row))
    card = cache / 'hf-readme.md'
    excerpt = ('model ability to predict increased ages in samples from patients with *pre-defined* '
               '**aging-accelerating conditions** (AACs) relative to samples from \nhealthy controls (HC).')
    if excerpt not in card.read_text():
        raise ValueError('ComputAge expected-direction definition changed')
    result = []
    for (gse, tissue, condition), members in sorted(groups.items()):
        if condition == 'HC':
            continue
        index, example = members[0]
        controls = groups.get((gse, tissue, 'HC'), [])
        q = evidence(path, f'TSV line {index}; group DatasetID={gse}, Tissue={tissue}, Condition={condition}',
                     '\t'.join(example.values()))
        q2 = evidence(card, 'Introduction, expected-direction definition', excerpt)
        result.append(record('computage', gse, tissue, None, 'Condition=' + condition,
                             'Condition=HC (healthy controls); may be supplied from other datasets in the benchmark',
                             'literature_expected_direction', 'Pinned Hugging Face computage_bench_meta.tsv', [q, q2],
                             expected_direction='case higher', condition_code=condition,
                             condition_classes=sorted({m['Class'] for _, m in members}),
                             cell_types=sorted({m['CellType'] for _, m in members}),
                             metadata_n_case=len(members), metadata_n_control=len(controls),
                             metadata_n_case_unique_gsm=len({m['SampleID'] for _, m in members}),
                             metadata_n_control_unique_gsm=len({m['SampleID'] for _, m in controls}),
                             control_available_within_series=bool(controls),
                             metadata_count_method='Count metadata lines and distinct SampleID within DatasetID/Tissue/Condition.',
                             notes=['Expected aging direction from benchmark curation, not a published clock statistic.',
                                    'Metadata counts are computed here; reported arm counts remain null.',
                                    'Condition codes are retained without inventing diagnostic thresholds.']))
    return result


def signature_records(cache):
    result = []
    path = cache / 'tysh2023-mmc6.xlsx'
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    sheet = workbook['A (Public datasets)']
    rows = list(sheet.values)
    for index, cells in enumerate(rows[2:], 3):
        if not isinstance(cells[2], str) or not re.fullmatch('GSE[0-9]+', cells[2]):
            continue
        for col in range(3, len(cells), 3):
            if cells[col] is None:
                continue
            signature = rows[0][col]
            coordinates = [openpyxl.utils.get_column_letter(c + 1) + str(index) for c in range(col, col + 3)]
            source_cells = {'Intervention': cells[0], 'Tissue': cells[1], 'GSE ID': cells[2],
                            'Signature': signature, 'NES': cells[col],
                            'P.Value': cells[col + 1], 'P.Adjusted': cells[col + 2]}
            q = evidence(path, sheet.title + ' row ' + str(index) + ', ' + ','.join(coordinates),
                         json.dumps(source_cells, ensure_ascii=False))
            result.append(record('tysh2023', cells[2], cells[1], None, cells[0],
                                 'Control definition not stated in Table S6; consult original series design',
                                 'reported_signature_statistic', 'Table S6A, ' + sheet.title + ', ' + signature, [q],
                                 [statistic('normalized_enrichment_score', cells[col], 'dimensionless signature score'),
                                  statistic('p_value', cells[col + 1], 'probability', method='Permutation test'),
                                  statistic('adjusted_p_value', cells[col + 2], 'probability', method='Benjamini-Hochberg')],
                                 modality='mouse transcriptomics', signature_source=signature,
                                 signature_score_sign='positive' if cells[col] > 0 else 'negative' if cells[col] < 0 else 'zero',
                                 notes=['Signature enrichment is not tAge, predicted chronological age, or a case-control clock effect.',
                                        'No direction on a ClockBase clock is inferred from NES.']))
    workbook.close()
    metadata_path = cache / 'tysh2023-mmc1.xlsx'
    metadata = openpyxl.load_workbook(metadata_path, read_only=True, data_only=True)
    designs = list(metadata['C (RNAseq, CMap compounds)'].values)
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    sheet = workbook['B (Chemical compounds)']
    rows = list(sheet.values)
    for index, cells in enumerate(rows[2:], 3):
        design = next((r for r in designs[1:] if r[0] and r[0].strip().lower() == str(cells[0]).strip().lower()
                       and r[3] == cells[1] and r[12] == 'GSE131868'), None)
        if design is None:
            continue
        control = next(r for r in designs[1:] if r[0] == 'Control' and r[3] == design[3] and r[12] == design[12])
        link = evidence(metadata_path, 'Table S1C, explicit GEO-linked intervention and control rows',
                        json.dumps({'intervention': design, 'control': control}, ensure_ascii=False))
        for col in range(2, len(cells), 3):
            if cells[col] is None:
                continue
            signature = rows[0][col]
            source_cells = {'CMap Compound': cells[0], 'Tissue': cells[1], 'Signature': signature,
                            'NES': cells[col], 'P.Value': cells[col + 1], 'P.Adjusted': cells[col + 2]}
            q = evidence(path, f'{sheet.title} row {index}, columns {col + 1}-{col + 3}',
                         json.dumps(source_cells, ensure_ascii=False))
            result.append(record('tysh2023', design[12], cells[1], None,
                                 f'{design[0]}, {design[1]}, {design[5]} months old, male UM-HET3 mice',
                                 'Male UM-HET3 controls, four months old (Control in Table S1C)',
                                 'reported_signature_statistic', 'Table S6B; Table S1C', [q, link],
                                 [statistic('normalized_enrichment_score', cells[col], 'dimensionless signature score'),
                                  statistic('p_value', cells[col + 1], 'probability', method='Permutation test'),
                                  statistic('adjusted_p_value', cells[col + 2], 'probability', method='Benjamini-Hochberg')],
                                 modality='mouse RNA-seq', n_case=design[4], n_control=control[4],
                                 signature_source=signature,
                                 signature_score_sign='positive' if cells[col] > 0 else 'negative' if cells[col] < 0 else 'zero',
                                 notes=['GEO linkage is explicit in Table S1C. Rows labelled Study are excluded.',
                                        'Signature enrichment is not a tAge clock effect.']))
    workbook.close()
    metadata.close()
    path = cache / 'tysh2019-pmc.html'
    linking = paragraph(path, 'For interventions from publicly available sources (Fig. 7E (upper))')
    specifications = [
        ('GSE21060', 'liver', 'Intraperitoneal IL-6 injection', 'all longevity signatures', 0.025, 'negative',
         'We detected a significant negative association of this intervention with all longevity signatures'),
        ('GSE77082', 'liver', 'Mat1a knockout', 'four of six longevity signatures', 0.02, 'negative',
         'Indeed, our test revealed such a pattern with statistical significance for 4 out of 6 signatures'),
        ('GSE15891', 'liver', '11.5 kPa Po2 hypoxia (11.8% oxygen), 32 days', 'all longevity signatures except rapamycin', 0.034, 'positive',
         'We detected a significant positive association of the response to hypoxia'),
        ('GSE11287', 'liver', 'Hepatocyte-specific conditional Keap1 knockout', 'almost all longevity signatures except rapamycin', 0.0015, 'positive',
         'identified statistically significant positive association with almost all longevity signatures'),
        ('GSE10421', 'liver', 'C57BL/6 males, seven weeks old', 'all longevity signatures except rapamycin', 5.3e-4, 'positive',
         'All longevity signatures except for rapamycin demonstrated a significant positive association with C57BL/6'),
        ('GSE49000', 'liver', 'SIRT1 activator SRT2104 in males', 'CR and common gene patterns', None, 'positive',
         'Accordingly, we detected significant positive associations for these models in males with CR and common gene patterns'),
        ('GSE104234', 'liver', 'Calorie restriction in rhesus monkeys', 'CR signature', None, 'positive',
         'We observed a strong significant association with the CR signature'),
    ]
    for gse, tissue, case, signature, fdr, sign, fragment in specifications:
        q = paragraph(path, fragment)
        stats = [statistic('adjusted_p_value', fdr, 'probability', '<')] if fdr is not None else []
        result.append(record('tysh2019', gse, tissue, None, case,
                             'DBA/2 males, seven weeks old' if gse == 'GSE10421' else 'Original intervention control; not defined in the quoted result',
                             'reported_signature_statistic' if stats else 'reported_signature_direction',
                             'Figure 7E; Results; STAR Methods GEO linkage', [q, linking], stats,
                             modality='rhesus macaque transcriptomics' if gse == 'GSE104234' else 'mouse transcriptomics',
                             signature_source=signature, signature_association_direction=sign,
                             notes=['The adjusted P bound applies to the stated collection of signatures, not to an invented individual clock.',
                                    'No ClockBase clock or case-higher/lower direction is assigned.']))
    return result


def biolearn_records(cache):
    """Retain retrieved clock validations with their non-case-control estimand."""
    path = cache / 'biolearn-43587_2025_987_MOESM3_ESM.xlsx'
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    result = []
    for sheet_name in ['Fig2a', 'Fig2c']:
        sheet = workbook[sheet_name]
        rows = list(sheet.values)
        for index, cells in enumerate(rows[1:], 2):
            values = dict(zip(rows[0], cells))
            q = evidence(path, f'{sheet_name}, row {index}', json.dumps(values, ensure_ascii=False))
            stats = ([statistic('pearson_r', values['Pearson_R'], 'dimensionless')]
                     if sheet_name == 'Fig2a' else
                     [statistic('mae', values['MAE'], 'years'), statistic('rmse', values['RMSE'], 'years')])
            result.append(record('biolearn', values['Dataset'], 'whole blood', values['Clock'], None, None,
                                 'reported_clock_validation', 'Source Data Figure 2, ' + sheet_name, [q], stats,
                                 n_total=values.get('N_Samples'),
                                 notes=['Age prediction accuracy, not a disease-versus-control clock comparison.']))
    workbook.close()
    return result


def bind_atlas(records, atlas_path):
    if not atlas_path.exists():
        raise FileNotFoundError(f'Cannot proceed - required data file not found: {atlas_path}')
    index = defaultdict(list)
    with atlas_path.open(newline='') as handle:
        for row in csv.DictReader(handle):
            index[row['gse_id']].append({key: row[key] for key in
                                       ['atlas_row_ordinal', 'modality', 'clock', 'comparison_name', 'condition_name']})
    for ordinal, row in enumerate(records, 1):
        row['record_id'] = f'curated-{ordinal:05d}'
        matches = index.get(row['gse_id'], [])
        clock_matches = [r for r in matches if r['clock'] == row['clock_clockbase'] and r['modality'] == row['modality']]
        row['atlas'] = {'series_present': bool(matches), 'modalities': sorted({r['modality'] for r in matches}),
                        'row_count': len(matches), 'rows': matches,
                        'series_clock_match': bool(clock_matches),
                        'series_clock_row_ordinals': [r['atlas_row_ordinal'] for r in clock_matches],
                        'comparison_match': None}
        row['direct_case_control_statistic'] = row['ground_truth_type'] == 'reported_case_control_statistic'
        row['notes'].append('Atlas links establish series/clock overlap only; arm, tissue and estimand matching remain unadjudicated.')
    return records


def bind_evidence(records, cache):
    """Resolve every excerpt to cached bytes and the originating GET receipt."""
    for row in records:
        row['source_url'] = ('https://pmc.ncbi.nlm.nih.gov/articles/' + row['pmcid'] + '/'
                             if row['pmcid'] else
                             'https://huggingface.co/datasets/computage/computage_bench/tree/'
                             '3cf52681537b303a5b7e21df03a3ae0a94ed497e')
        for quote in row['quote_evidence']:
            path = cache / quote['file']
            payload = path.read_bytes()
            quote['sha256'] = hashlib.sha256(payload).hexdigest()
            name = path.name
            chain = []
            while not (cache / (name + '.receipt.json')).exists():
                derived = json.loads((cache / (name + '.derived.json')).read_text())
                if hashlib.sha256((cache / name).read_bytes()).hexdigest() != derived['sha256']:
                    raise ValueError(f'Derived evidence byte mismatch: {name}')
                chain.append(derived)
                name = derived['parent_file']
            receipt = json.loads((cache / (name + '.receipt.json')).read_text())
            if receipt.get('error'):
                raise ValueError(f'Evidence points to a failed download: {name}')
            if hashlib.sha256((cache / name).read_bytes()).hexdigest() != receipt['sha256']:
                raise ValueError(f'Download evidence byte mismatch: {name}')
            quote.update(url=receipt['url'], retrieved_at=receipt['retrieved_at'],
                         download_sha256=receipt['sha256'], derivation=chain)
    return records


def manifest(cache, records, atlas_path):
    downloads = [json.loads(path.read_text()) for path in sorted(cache.glob('*.receipt.json'))]
    evidence_files = {q['file'] for row in records for q in row['quote_evidence']}
    evidence_files.update(d['parent_file'] for row in records for q in row['quote_evidence']
                          for d in q['derivation'])
    for entry in downloads:
        path = Path(entry['path'])
        entry['used_as_evidence'] = path.name in evidence_files
        if not path.exists():
            entry['content_validation'] = 'request_failed'
            continue
        payload = path.read_bytes()
        if path.suffix in {'.zip', '.xlsx', '.docx'} and not payload.startswith(b'PK'):
            entry['content_validation'] = 'not_requested_binary_format; cached interstitial or API error'
        elif b'recaptcha/challengepage' in payload or b'POW_CHALLENGE' in payload:
            entry['content_validation'] = 'browser_challenge; not source evidence'
        else:
            entry['content_validation'] = 'cached_response; evidence identity checked during extraction if used'
    statuses = json.loads(Path(__file__).with_name('source_limits_2026_10_05.json').read_text())
    sources = []
    for key, (doi, pmid, pmcid) in SOURCE_DEFINITIONS.items():
        selected = [r for r in records if r['source_id'] == key]
        sources.append({'source_id': key, 'doi': doi, 'pmid': pmid, 'pmcid': pmcid,
                        'records': len(selected), 'atlas_series_records': sum(r['atlas']['series_present'] for r in selected),
                        'atlas_series_clock_records': sum(r['atlas']['series_clock_match'] for r in selected),
                        **statuses[key]})
    derived = [json.loads(path.read_text()) for path in sorted(cache.glob('*.derived.json'))]
    return {'schema_version': 1, 'assembled_at': datetime.now(timezone.utc).isoformat(),
            'artifact_folder': str(cache.parent), 'cache_folder': str(cache), 'atlas': {'path': str(atlas_path),
             'sha256': hashlib.sha256(atlas_path.read_bytes()).hexdigest(),
             'mtime_utc': datetime.fromtimestamp(atlas_path.stat().st_mtime, timezone.utc).isoformat()},
            'sources': sources, 'downloads': downloads, 'derived_files': derived,
            'counts': {'records': len(records), 'by_type': dict(Counter(r['ground_truth_type'] for r in records)),
                       'atlas_series_records': sum(r['atlas']['series_present'] for r in records),
                       'atlas_series_clock_records': sum(r['atlas']['series_clock_match'] for r in records),
                       'direct_case_control_series_clock_records': sum(r['direct_case_control_statistic'] and r['atlas']['series_clock_match'] for r in records)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--cache', type=Path, help='Reuse an existing download cache for a fresh output folder.')
    parser.add_argument('--atlas', type=Path, default=ATLAS)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    try:
        cache = args.cache if args.cache is not None else args.root / 'downloads'
        args.root.mkdir(parents=True, exist_ok=True)
        records = (classic_records(cache) + phenoage_records(cache) + computage_records(cache)
                   + signature_records(cache) + biolearn_records(cache))
        bind_atlas(records, args.atlas)
        bind_evidence(records, cache)
        provenance = manifest(cache, records, args.atlas)
        provenance['artifact_folder'] = str(args.root)
        for name in ['records.jsonl', 'sources.json']:
            if (args.root / name).exists():
                raise FileExistsError(f'Preserve existing output; choose a fresh --root: {args.root / name}')
        (args.root / 'records.jsonl').write_text(''.join(json.dumps(r, ensure_ascii=False, allow_nan=False) + '\n' for r in records))
        (args.root / 'sources.json').write_text(json.dumps(provenance, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
        logging.info('Assembled %s records in %s', len(records), args.root)
        print(json.dumps(provenance['counts'], indent=2))
    except (OSError, ValueError, KeyError, ET.ParseError) as exc:
        logging.error('Curated assembly failed: %s', exc)
        raise


if __name__ == '__main__':
    main()
