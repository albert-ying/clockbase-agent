# Curated published-clock sources, 2026-10-05

The collection contains 360 source records from 14 investigated sources. Of these,
137 reference a series present in the supplied ClockBase atlas, and 28 have an
exact series, modality and clock match. Sixteen of those 28 contain per-series
case-control statistics. These are candidate comparisons: their arms, tissue,
adjustment and estimand still require adjudication before scoring ClockBase.
The 360 records are not 360 independent clock effects.

The requested output folder is
`/Entropy/kying/clockbase-benchmark-2026-10-05/curated/`.
It contains the final `records.jsonl`, `sources.json`, this README, and a
`downloads/` cache. The repository contains copies of the final records and
manifest, extraction code, source limitations, and eight source-bound tests.
Downloaded papers and supplements remain in the requested external folder.
Earlier owned assemblies were preserved in dated subfolders; they are not inputs
to the final collection.

**Source coverage.** Counts describe records, including separate regression
models, separate tests and different signatures. The last column requires exact
clock and modality as well as series membership.

| Source ID and publication | Records | Atlas series | Atlas series + clock |
| --- | ---: | ---: | ---: |
| tysh2023 — Tyshkovskiy, Cell 2023 | 242 | 66 | 0 |
| tysh2019 — Tyshkovskiy, Cell Metabolism 2019 | 7 | 0 | 0 |
| computage — ComputAgeBench metadata | 66 | 35 | 0 |
| biolearn — Ying et al., Nature Aging 2025 version | 8 | 8 | 4 |
| down — Horvath et al., Down syndrome 2015 | 8 | 8 | 8 |
| hivhorvath — Horvath and Levine, HIV 2015 | 6 | 4 | 4 |
| hivgross — Gross et al., HIV 2016 | 4 | 4 | 1 |
| obesity — Horvath et al., liver/BMI 2014 | 2 | 1 | 1 |
| huntington — Horvath et al., Huntington disease 2016 | 3 | 0 | 0 |
| parkinson — Horvath and Ritz, Parkinson disease 2015 | 4 | 4 | 4 |
| werner — Maierhofer et al., Werner syndrome 2017 | 0 | 0 | 0 |
| phenoage — Levine et al., PhenoAge 2018 | 10 | 7 | 6 |
| grimage — Lu et al., GrimAge 2019 | 0 | 0 | 0 |
| pace — Belsky et al., DunedinPACE 2022 | 0 | 0 | 0 |
| **Total** | **360** | **137** | **28** |

There are 89 distinct GEO series in the records, of which 42 occur in the atlas.
The 16 case-control statistic records with matching series/modality/clock cover
12 distinct series/clock pairs. All mapped candidate case-control clocks are
DNA methylation clocks; no retrieved signature effect was equated with tAge.

**Method and evidence.** Primary sources were obtained from PMC/Europe PMC,
publisher attachment servers and Hugging Face. GEO SOFT responses establish the
Parkinson cohorts' link to PMID 26655927. DOI, PMID/PMCID where available, source
location, original clock name, accession, tissue, cases, controls and quotes are
stored per record. `sources.json` records all GET attempts, including failed
requests, redirects, response formats, retrieval timestamps and SHA-256 hashes.
Binary-format checks exclude HTML download interstitials and browser challenges
from evidence. Incorrect early supplement-name assumptions and unrelated search
returns are retained in the cache and explicitly excluded.

Every quote has a file, selector, SHA-256 hash and originating GET receipt.
Paragraph quotations collapse whitespace. Table quotations serialize selected
literal cell values with coordinates or source row labels; they are not invented
prose. PDF excerpts use `pdftotext -layout`; their PDF page and derivation chain
back to the downloaded ZIP are recorded. Derived DOCX visible-text files exclude
Word field instructions. No statistical values are computed from literature
figures or reconstructed from P values. Source values are copied as stated.
Nulls indicate information not supplied or a linkage not established.

ComputAge metadata is pinned to Hugging Face revision
`3cf52681537b303a5b7e21df03a3ae0a94ed497e`
(dataset DOI `10.57967/hf/4347`). Records group metadata by series, tissue and
non-HC condition. Diagnostic codes remain as published. Metadata row counts and
distinct GSM counts are computed from that file and stored in separate
`metadata_n_*` fields, never in reported `n_case`/`n_control`. Absence of healthy
controls in the same series/tissue is flagged. The dataset card's expected
increase in age for aging-accelerating conditions is an expected direction;
it is not a reported effect or significance test.

The atlas input is
`/Entropy/kying/clockbase-verifier-v3-2026-09-27/final-2026-10-05/atlas-v3/atlas_v3_combined.csv`.
Its hash and modification timestamp are in the manifest. Each record includes
all same-series atlas row ordinals, modalities, clocks, comparison names and
condition names. `comparison_match` remains null: neither an accession match
nor a clock match proves that ClockBase analysed the published contrast.

**Record types and mapping rules.**

| `ground_truth_type` | Records | Interpretation |
| --- | ---: | --- |
| reported_case_control_statistic | 24 | Per-series case-control statistic copied from source |
| reported_case_control_direction | 3 | Direction/significance statement without a numerical statistic |
| reported_continuous_association | 6 | BMI correlation or regression, not binary obesity arms |
| reported_pooled_case_control_statistic | 4 | Two pooled PD results linked to both contributing series |
| literature_expected_direction | 66 | ComputAge curated expectation, no reported statistic |
| reported_signature_statistic | 247 | Signature enrichment or association; not an aging clock |
| reported_signature_direction | 2 | Signature association without a numerical statistic |
| reported_clock_validation | 8 | Biolearn chronological-age prediction validation |

Only `reported_case_control_statistic` sets `direct_case_control_statistic=true`.
The strict initial filter is that flag together with `atlas.series_clock_match`.
It yields 16 records, and still requires contrast adjudication. Numerical values
include `name`, `value`, inequality `relation`, `units` and test method when the
method is stated. Requested effect/beta/P/FDR fields are null unless that precise
statistic is reported; an age-acceleration value is not relabelled a raw mean
difference. A non-significant statement is not proof of equivalence. Direction
stays null when the source does not state it.

The explicitly named Horvath 2013 353-CpG model maps to `Horvathv1`; explicitly
named DNAm PhenoAge maps to `PhenoAge`. Exact Biolearn `Horvathv1` and `PhenoAge`
names map directly. Hannum/Horvath consensus, Hannum, Lin and longevity signatures
remain unmapped. Original 2019 GrimAge is not mapped to `GrimAgeV2`. No guessed
mapping to a Hazard/tAge model is made. Species/modality descriptions that do
not exactly match the atlas do not obtain clock matches.

**Examples and limits.**

- Down syndrome, GSE52588, Horvathv1: Table 1 reports P=0.00056; Table 2 reports
  DS beta=3.65564 DNAm years, SE=1.08842, P=0.0012 and age acceleration=4.6 years.
  These are two records with different tests, not two independent cohorts.
- HIV, GSE67705: the source's 4.9-year whole-blood effect uses a consensus model
  and remains unmapped. The separately reported neutrophil Horvath result is
  0.4 years, P>0.05, with no significant difference.
- PhenoAge, GSE72774: the source reports P=0.028 and 289 cases/219 controls but
  does not state direction in the result passage; direction remains null.
- Tyshkovskiy 2023 Table S6 gives NES and permutation/adjusted P values. These
  are signatures, not per-GEO tAge clock ages. GSE131868 is explicitly linked to
  three liver interventions in Table S1C; unlinked rows labelled Study are omitted.

Both Tyshkovskiy papers were verified and relevant tables obtained, but no
per-GEO tAge effect was found in the checked material. Their signature evidence
is retained under separate types. Biolearn's later 2025 paper reports GSE40279
clock validation, retained as validation; its Down Syndrome disease-score
comparison is not an identified aging-clock comparison and is omitted.

The 2017 Werner source reports a result but no verified GEO linkage was obtained;
it contributes no record. GSE131752 is not assigned to that paper by inference.
The original GrimAge main text and all retrieved supplements supplied no
verified GEO case-control result. DunedinPACE's GSE55763 linkage concerns
technical replicates rather than cases and controls. These sources contribute
zero records, with retrieval provenance and reasons preserved. Cohort-level
effects from FHS, WHI, JHS, InCHIANTI or E-Risk are not assigned guessed GEO IDs.

Parkinson pooled effects share `pooling_id` and must be deduplicated before
analysis. PhenoAge GEO links sometimes require the paper's explicit reference,
dataset ordering and sample counts; supporting sources and reasoning are in each
record. Conflicting tissue labels, different cohort sizes between publications,
repeated donor brain samples and adjusted-versus-unadjusted estimands remain
documented rather than resolved by inference. Source-specific limits also appear
in `sources.json` and `source_limits_2026_10_05.json`.

**Reproduction and checks.** Run from the repository root with the existing
Python 3.12 environment. The extractor uses the environment's existing openpyxl
installation; no package was installed. The retriever uses the standard library,
GET requests, a 0.5-second delay, cached error bodies and immutable cache identity
checks. The assembler makes no network requests and refuses existing outputs.

The completed default assembly command produced 360 records, with 137 series
overlaps, 28 series/modality/clock overlaps and 16 case-control statistic overlaps:

```bash
TMPDIR=/Entropy/kying/tmp .venv-pipeline/bin/python -B -m benchmark.published_clock_2026_10_05.curated.assemble_2026_10_05
```

For another run, supply a fresh output folder inside the owned directory and the
existing cache:

```bash
TMPDIR=/Entropy/kying/tmp .venv-pipeline/bin/python -B -m benchmark.published_clock_2026_10_05.curated.assemble_2026_10_05 --root /Entropy/kying/clockbase-benchmark-2026-10-05/curated/offline-rebuild-2026-10-05 --cache /Entropy/kying/clockbase-benchmark-2026-10-05/curated/downloads
```

Tests validate quote byte provenance, copied source values, expected-direction
boundaries, signature/clock separation, consensus mapping, pooled-result handling,
strict atlas joins and refusal to change cached identities. The command passes
all eight tests:

```bash
TMPDIR=/Entropy/kying/tmp .venv-pipeline/bin/python -B -m unittest benchmark.published_clock_2026_10_05.curated.test_curated_2026_10_05
```

`validation_2026_10_05.json` records final file hashes, download and derivation
hash verification, output-copy equality, and byte-identical offline record
reproduction. Its check command is:

```bash
TMPDIR=/Entropy/kying/tmp .venv-pipeline/bin/python -B -m benchmark.published_clock_2026_10_05.curated.verify_2026_10_05 --write-report
```

Omit `--write-report` on subsequent runs to preserve the existing check receipt.
No empirical input, atlas file, virtual environment or file outside
the two owned directories was modified. No git mutation was performed.
