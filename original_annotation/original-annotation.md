# Project Overview
A small pipeline that annotates a user-supplied list of human genes with IDs, function text, psychiatric GWAS hits, Harmonizome psychiatric signals, and PubMed mental-health literature. Outputs both a full CSV for auditing and a formatted Excel file for sharing.

# Repository Layout
- `main.py` - orchestrates the pipeline and writes outputs.
- `modules/hgnc_pull.py` - pulls HGNC/NCBI identifiers, synonyms, summaries, UniProt accessions.
- `modules/function_annotation.py` - gathers function text from NCBI, UniProt, HGNC and builds a combined block.
- `modules/gwas_association.py` - joins precomputed psychiatric GWAS summaries.
- `modules/pubmed_association.py` - queries PubMed for mental-health-related hits and aggregates per gene.
- `modules/harmonizome_association.py` - queries Harmonizome and keeps psychiatric-like disease/phenotype associations.
- `gene_psych_gwas.py` - builds `resources/gwas/gene_psych_gwas.tsv` from the GWAS Catalog TSV.
- `resources/` - supporting data (GWAS Catalog tables, mappings, caches).
- `requirements.txt` - Python dependencies.

# Prerequisites
- Python 3.10+.
- Install deps (recommend a virtual env):
  ```bash
  pip install -r requirements.txt
  ```
- Network access to NCBI Datasets, NCBI E-utilities, UniProt REST, HGNC REST, Harmonizome API.
- GWAS Catalog associations file: `resources/gwas/gwas-catalog-associations_ontology-annotated.tsv` (not included; refresh from GWAS Catalog latest FTP).

# Configuration
`config.json` holds NCBI credentials for PubMed:
```json
{
  "NCBI_EMAIL": "your.email@example.edu",
  "NCBI_API_KEY": "optional-api-key"
}
```
- `NCBI_EMAIL` is required for the PubMed step; without it, PubMed is skipped.
- `NCBI_API_KEY` is optional but increases rate limits.

Input genes: plain text, one symbol or LOC ID per line (lines starting with `#` are ignored), e.g. `genes.txt` or `genes_test.txt`.

# Pipeline (main.py)
1) Read gene list path from CLI or prompt with path completion.
2) HGNC/NCBI pull: resolve symbols/LOC IDs to HGNC/Entrez/Ensembl, gene type, synonyms, UniProt, summaries.
3) Function aggregation: fetch UniProt function comments and HGNC text; build per-source fields and a combined block.
4) Psychiatric GWAS attach: join `resources/gwas/gene_psych_gwas.tsv`; add counts/traits/labels/PMIDs/EFO/study lists.
5) Harmonizome attach: query Harmonizome, keep disease/phenotype datasets, filter by psychiatric keywords, summarize per gene.
6) PubMed attach: use Entrez Gene → PubMed links when Entrez ID exists, else symbol search; filter by mental-health MeSH/text terms; tag genetic-related hits.
7) Outputs: write `annotated_genes_raw.csv` (all columns), build the deliverable view, preview first rows, and write `annotated_genes.xlsx`.

# Module Notes
- `hgnc_pull`: human-only (taxon 9606); adds `annotation_status` and per-row `log`.
- `function_annotation`: prefers UniProt comments; notes duplicate text across sources in the combined block.
- `gwas_association`: merges on `approved_symbol` (or `symbol` fallback); fills empty columns instead of failing when data is missing.
- `harmonizome_association`: uppercases symbols, queries `/download/associations?gene=`, keeps only disease/phenotype datasets and psychiatric-like labels.
- `pubmed_association`: ELink preferred, ESearch fallback; batches EFetch; mental-health detection via MeSH + title keywords; genetic flag via MeSH + regex.
- `gene_psych_gwas`: filters the GWAS Catalog by psychiatric keywords in trait text, explodes multi-gene rows, aggregates per gene into the summary TSV.

# Keyword / Heuristic Lists
<details>
<summary>GWAS psychiatric keywords (gene_psych_gwas.py & Harmonizome alignment)</summary>

- <details>
  <summary>Core / clinical</summary>

  - adhd
  - agoraphobia
  - anxiety
  - asperger
  - autism
  - bipolar
  - conduct disorder
  - depression
  - depressive
  - disruptive behavior
  - dysthymia
  - hyperactivity
  - major depressive
  - mania
  - mental disorder
  - mental or behavioral
  - mental or behavioural
  - mood disorder
  - oppositional defiant
  - panic disorder
  - pervasive developmental
  - personality disorder
  - psychosis
  - psychotic
  - schizophrenia
  - social phobia
  - somatoform
  - stress related
  - stress-related
  - tic disorder
  - tourette
  - unipolar
  </details>

- <details>
  <summary>Substance / behavioral</summary>

  - alcohol dependence
  - alcohol-use disorder
  - binge eating
  - bulimia
  - drug dependence
  - eating disorder
  - nicotine dependence
  - obsessive compulsive
  - obsessive-compulsive
  - ocd
  - post-traumatic stress
  - posttraumatic stress
  - ptsd
  - substance abuse
  - substance use
  - substance-use
  </details>
</details>

<details>
<summary>Harmonizome psychiatric keywords (modules/harmonizome_association.py)</summary>

- <details>
  <summary>Core / clinical</summary>

  - adhd
  - agoraphobia
  - anxiety
  - asperger
  - autism
  - bipolar
  - conduct disorder
  - depression
  - depressive
  - disruptive behavior
  - mania
  - manic
  - mental disorder
  - mental or behavioral
  - mental or behavioural
  - mood disorder
  - neuropsychiatric
  - neuroticism
  - oppositional defiant
  - panic disorder
  - pervasive developmental
  - personality disorder
  - phobia
  - psychiatric
  - psychosis
  - psychotic
  - schizophrenia
  - self harm
  - self-harm
  - social phobia
  - somatoform
  - stress related
  - stress-related
  - suicide
  - suicidal
  - tic disorder
  - tourette
  - unipolar
  </details>

- <details>
  <summary>Substance / behavioral</summary>

  - addiction
  - alcohol dependence
  - alcohol-use disorder
  - anorexia nervosa
  - binge eating
  - bulimia
  - eating disorder
  - nicotine dependence
  - obsessive compulsive
  - obsessive-compulsive
  - ocd
  - post-traumatic stress
  - posttraumatic stress
  - ptsd
  - substance abuse
  - substance use
  - substance-use
  </details>
</details>

<details>
<summary>PubMed mental-health MeSH terms (modules/pubmed_association.py)</summary>

- <details>
  <summary>Core / clinical</summary>

  - Anxiety Disorders
  - Attention Deficit Disorder with Hyperactivity
  - Autism Spectrum Disorder
  - Bipolar Disorder
  - Cyclothymic Disorder
  - Depressive Disorder
  - Depressive Disorder, Major
  - Intellectual Disability
  - Mental Disorders
  - Mental Health
  - Obsessive-Compulsive Disorder
  - Panic Disorder
  - Phobic Disorders
  - Psychotic Disorders
  - Schizoaffective Disorder
  - Schizophrenia
  - Stress Disorders, Post-Traumatic
  </details>

- <details>
  <summary>Substance / behavioral</summary>

  - Alcohol-Related Disorders
  - Cocaine-Related Disorders
  - Opioid-Related Disorders
  - Substance-Related Disorders
  - Suicide
  - Suicidal Ideation
  </details>
</details>

<details>
<summary>PubMed mental-health text keywords (modules/pubmed_association.py)</summary>

- <details>
  <summary>Core / clinical</summary>

  - adhd
  - anxiety
  - asperger
  - autism
  - bipolar
  - depressi
  - mania
  - mental disorder
  - mental health
  - panic disorder
  - post-traumatic
  - psychosis
  - psychotic
  - ptsd
  - schizophrenia
  </details>

- <details>
  <summary>Substance / behavioral</summary>

  - addiction
  - alcohol use
  - alcohol-use
  - attention-deficit
  - substance use
  - substance-use
  - suicide
  - suicidal
  </details>
</details>

<details>
<summary>PubMed genetic flag patterns</summary>

- <details>
  <summary>Genetic / variant terms</summary>

  - cnv
  - copy number variation
  - genetic
  - genome-wide association
  - gwas
  - mutation
  - mutations
  - polymorphism
  - polymorphisms
  - snp
  - variant
  - variants
  </details>

- <details>
  <summary>Genetic MeSH terms</summary>

  - Genetic Predisposition to Disease
  - Genetic Variation
  - Genome-Wide Association Study
  - Genotype
  - Polymorphism, Genetic
  </details>
</details>

# Running the Pipeline
```bash
python main.py path/to/genes.txt
```
Outputs to the project root:
- `annotated_genes_raw.csv` full table with all columns.
- `annotated_genes.xlsx` formatted deliverable.
- Console preview of the first five rows of the deliverable view.

# Refreshing the Psychiatric GWAS Summary
```bash
python gene_psych_gwas.py \
  --assoc-path resources/gwas/gwas-catalog-associations_ontology-annotated.tsv \
  --out resources/gwas/gene_psych_gwas.tsv
```
Optional single-gene check:
```bash
python gene_psych_gwas.py --gene COMT
```

# Data and Methods Notes
- Species fixed to human (taxon 9606) for NCBI calls.
- Symbols are uppercased/trimmed for joins; input order preserved.
- Nullable `Int64` distinguishes zero vs missing; list fields use `; ` separators and render as multi-line in Excel.
- PubMed rate limits follow NCBI guidance (~3 req/sec without API key); adjust `pause` parameters if needed.
- Non-fatal issues land in the per-row `log` instead of aborting the run.
- API-heavy: large gene lists are slow (Harmonizome alone took ~9 minutes for ~75 genes). Prefer small/medium inputs.

# Troubleshooting
- PubMed skipped: add `NCBI_EMAIL` (and optionally `NCBI_API_KEY`) to `config.json`.
- GWAS columns empty: ensure `resources/gwas/gene_psych_gwas.tsv` exists or rebuild it.
- Harmonizome all NA: check network; only psychiatric-like labels are retained by design.
- UniProt function empty: some accessions lack `FUNCTION` comments; combined block will note `n/a`.
- Excel export skipped: install `openpyxl` (preferred) or `XlsxWriter`.

# Deliverables and Audit
- `annotated_genes_raw.csv` keeps all intermediate columns for traceability.
- Excel deliverable prunes source-specific function columns, flattens synonyms, reorders name/symbol, and adds spacing in function text for readability.

# Limitations and Warnings
- PubMed E-utilities gene mapping can be noisy: ELink relies on NCBI’s gene→PubMed links, which can include off-target hits for symbols with aliases/legacy mappings; symbol search (ESearch) is even noisier and can pull unrelated PMIDs.
- Symbol ambiguity: genes with common words, legacy symbols, or LOC placeholders may map incorrectly or not at all; outputs should be manually sanity-checked.
- Stale inputs: `resources/gwas/gene_psych_gwas.tsv` is a snapshot of the GWAS Catalog; rebuild it after catalog updates to stay current.
- Harmonizome coverage is uneven across datasets; some genes or psychiatric labels may be missing entirely, and requests are slow for larger lists.
- No caching across runs; repeated calls will re-hit external APIs.
- This pipeline is not validated for clinical or production use; treat outputs as exploratory and confirm with primary sources.

# License and Attribution
Uses data from NCBI (Datasets, E-utilities), UniProt, HGNC, GWAS Catalog, and Harmonizome. Follow their terms of use.
Personal/non-production tool; use at your own risk, no warranty or support provided.
