# General Workflow
Birds-eye view workflow building `annotated_genes.xlsx`:
1) HGNC/NCBI pull: resolve symbols/LOC IDs to HGNC/Entrez/Ensembl, gene type, synonyms, UniProt, summaries.
1a) Function aggregation: fetch UniProt function comments and HGNC text; build per-source fields and a combined block.
2) gene association via **GWAS Catalog**, **HArmonizome**, **PubMed**, and **DISGENET**
2a) Psychiatric GWAS attach: join `gene_psych_gwas.tsv`; add counts/traits/labels/PMIDs/EFO/study lists.
2b) Harmonizome attach: query Harmonizome, keep disease/phenotype datasets, filter by psychiatric keywords, summarize per gene.
2c) PubMed attach: use Entrez Gene → PubMed links when Entrez ID exists, else symbol search; filter by mental-health MeSH/text terms; tag genetic-related hits.
2d) Received academic API for DISGENET's curated dataset, gathered associations via 
3) Incorporating CpG-disease associations


# Detailed workflows
<h3 style="
        margin: 0px">Gene Annotation</h3>
<details>
<summary>HGNC & NCBI Metadata</summary>

**[NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/)** ([REST API](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/)) & **[HGNC](https://www.genenames.org/)** ([REST API](https://www.genenames.org/help/rest/)):\
unique `ewas_res_groupsig_128.xlsx` gene list → base annotation

1. Input resolution
    - Accepts gene symbols or LOC IDs
    - LOC IDs (e.g., `LOC101926933`) resolved via NCBI Gene ID
    - Symbols resolved via NCBI symbol search (restricted to human)
2. Metadata extraction
    - NCBI: Retrieves HGNC `symbol`, `name`, `hgncID`, `entrezID`, `ensemblID`,  `gene_type`, `synonyms`, UniProt accession as `uniprot`
    - HGNC Fallback: If NCBI lacks an Ensembl ID, queries HGNC using the HGNC ID to fill the gap
    - creates `annotated_genes.xlsx` with above columns
</details>

<details>
<summary>Function Aggregation</summary>

**[UniProt](https://www.uniprot.org/)** ([REST API](https://www.uniprot.org/help/api)); **[HGNC](https://www.genenames.org/)** ([REST API](https://www.genenames.org/help/rest/#!/#tocAnchor-1-1)):\
Gene symbol → enriched annotation

1. Fetching
    - UniProt: Queries API using the `uniprot` accession from NCBI. Extracts "FUNCTION" comments, entry status (Reviewed/Unreviewed), and protein existence evidence (n/a, 1–5)
    - HGNC: Queries API; attempts to find "curator_summary", "gene_group", or "name" to serve as a functional description
2. Aggregation & Deduping
    - Combines text from NCBI (Summary), UniProt, and HGNC
    - Normalizes text (lowercase, whitespace) to detect duplicates automatically
    - Format: Returns a multi-line string labeling the source and collapsing duplicates in `function` column
    - Example:
        ```text
        [SOURCE: NCBI] Protein coding gene...
        [SOURCE: UniProt] same as NCBI
        [SOURCE: HGNC] n/a
        ```
3. Columns Added: `uniprot_review`, `protein_existence`, `function`
</details>

<details>
<summary>GWAS Catalog (Psychiatric)</summary>

**[GWAS Catalog](https://www.ebi.ac.uk/gwas/)** ([Data Downloads](https://www.ebi.ac.uk/gwas/docs/file-downloads)):\
`gwas-catalog-associations_ontology-annotated.tsv` → `misc/gene_psych_gwas.tsv`

1. Pre-processing
    - Filters EBI GWAS Catalog associations for psychiatric keywords (see *Appendix A*)
    - Explodes multi-gene entries (comma-separated MAPPED_GENE)
    - Aggregates traits, studies, and PMIDs per gene
2. Attachment
    - Joins the pre-computed psychiatric GWAS table to the main gene list by gene symbol
    - Logic: Left join on standardized symbol
3. Columns Added: `gwas_assoc_count`, `gwas_traits`, `gwas_labels`, `gwas_pmids`, `gwas_efo_uris`, `gwas_study_accessions`
</details>

<details>
<summary>Harmonizome (Psychiatric)</summary>

**[Harmonizome](https://maayanlab.cloud/Harmonizome/)** ([REST API](https://maayanlab.cloud/Harmonizome/documentation)):\
Gene Symbol → Psychiatric Associations (`annotated_genes.xlsx`)

1. Query all associations for a gene symbol via `download/associations` REST endpoint
2. Filtering
    - Datasets: Restricts to key disease datasets (*CTD, DisGeNET, GAD, GWASdb, DISEASES*)
    - Keywords: Filters associations for psychiatric terms (see *Appendix A*)
3. Columns added: `harmonizome_count`, `harmonizome_terms`, `harmonizome_datasets`
</details>

<details>
<summary>PubMed Literature</summary>

**[NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/)** ([API](https://www.ncbi.nlm.nih.gov/books/NBK25499/)):\
Entrez ID / Symbol → Mental Health Literature Hits (`annotated_genes.xlsx`)

1. Identification
    - Primary: `elink` (Gene ID → PubMed ID) for high-confidence links
    - Secondary: `esearch` (Symbol in Title/Abstract/MeSH) if Entrez ID is missing
2. Filtering
    - Retrieves full XML for identified PMIDs stored in `pubmed_pmid`
    - Mental Health Filter: Keeps papers matching specific MeSH terms or Title keywords (see *Appendix B*), lists term hits in `pubmed_terms`
    - Genetic Tagging: Flags papers as "Genetic" if they contain terms like "Polymorphism", "GWAS", "Variant"
3. Formatting
    - Generates a brief summary string for quick review in `pubmed_brief` column
    - Example:
        ```text
        Schizophrenia genetics 2024 PMID: 12345 Genetic
        Depression study 2023 PMID: 67890 Not Genetic
        ```
4. Columns added: `pubmed_count`, `pubmed_genetic_count`, `pubmed_pmid`, `pubmed_terms`, `pubmed_brief`
</details>

<details>
<summary>DisGeNET</summary>

**`disgenet2r`**:\
gene symbol → Gene-Disease Associations (`data/disgenet_gda.csv`)\
gene symbol → Gene-Evidence Associations (`data/disgenet_gea.csv`)  

1. Filter `disgenet_gea.csv` by column `diseaseClasses_UMLS_ST` == "Mental or Behavioral Dysfunction (T048)"
2. For every gene in `annotated_genes.xlsx`:
    1. `disgenet_psych_diseases` column:
        -   Group `disease_name` values in `disgenet_gea.csv`
        -   Logic:
            - `score` consistency across entries; else "error"
            - Sort unique disease names by score
        -   Format: `[Disease Name], [Score]`
            - Separator: `;` + newline
            - *Example*:
                ```text
                Schizophrenia, 0.7;
                Bipolar Disorder, 0.5
                ```
    2. `disgenet_evidence` column:
        - Separate `disease_name` (rows in `disgenet_gea.csv`)
        - Logic:
            - Polarity: "Positive", "Negative", or "NAPolarity" (for NA)
            - Reference: if `reference_type` == "PMID" → use `reference` ID; else `reference("NA")(source, associationType)`
            - Sort by (1) descending disease score, then (2) polarity ( + > na > - ) and chronological publication
        - Format: `[Disease], [Polarity], [Year], [Ref]`
            - *Example*:
                ```text
                Schizophrenia, Positive, 2006, 2489764;
                Schizophrenia, NAPolarity, 2004, NA(CLINVAR, GeneticVariation);
                Schizophrenia, Negative, 2005, 99348737;
                Bipolar Disorder, Positive, 2015, 9398387
                ```
</details>

<h3 style="
        margin-top: 10px;
        margin-bottom: 0px">Epigene Annotation</h3>
<details>
<summary>EWAS Atlas</summary>

**[NGDC EWAS Atlas](https://ngdc.cncb.ac.cn/ewas/atlas)** ([REST API](https://ngdc.cncb.ac.cn/ewas/api)):\
CpG probe ID → trait associations (`data/ewas_atlas.csv`)

1. Query API per CpG `probeID` in `ewas_res_groupsig_128.xlsx`:
    - Resolve symbols from relatedTranscription as `genes`
    - Capture cpgIsland status as `cpg_island`
    - Flatten `associationList`, containing `trait`, `correlation`, `rank`, `pmid`
        - `correlation`: hypermethylated (`pos`), hypomethylated (`neg`), not reported (`NA`)
        - `rank`: rank in study
2. Query API per `pmid` in resulted associations
    - Count `total_associations`
    - Capture CpG-of-interest p-value as `p`
3. Generate `ewas_atlas.csv`:
    - Logic: One row per trait association; metadata preserved for null-association probes
    - Headers:
        ```text
        cpg, genes, cpg_island, trait, correlation, p, rank, pmid, total_associations
        ```
</details>
<details>
<summary>UCSC Screen</summary>

**[UCSC Genome Browser](https://genome.ucsc.edu/)** ([REST API](https://genome.ucsc.edu/goldenPath/help/api.html)):\
CpG coordinates (`ewas_res_groupsig_128.xlsx`) → Track annotations (`data/ewas_ucsc_annotated.xlsx`)

1. Confirm coordinates resolve to the correct CpG ID via `snpArrayIllumina850k` track
2. **Track Queries**:
    1. Direct Overlap: query `clinvarMain` at exact CpG coordinates to identify clinical significance or mutations (unlikely)
    2. Neighborhood: query `gwasCatalog` within a 5kb window (±5000bp) to identify nearby disease-associated SNPs
3. Generate flattened hits written to Excel file with separate sheets for each track (`clinvarMain`, `gwasCatalog`)
    - Headers determined by resulted json fields
</details>
<br>
<br>

# Interesting findings
- **EWAS Atlas** results for `cg09420691` is only associated with `body mass index (BMI)` trait, despite its related gene, `PRDM15` being associated with `Mental Disorders` at 100% (**DISGENET**) and is the only relevant association from this dataset, suggesting.. what?

Although selection was based solely on genomic proximity, the two closest traits of 765 interestingly ranged from a classical biological phenotype to an educational attainment proxy, namely the highest mathematics course completed (see below).

| cpg_id      | distance | pubMedID | trait                              | region    | genes       |
|-------------|----------|----------|------------------------------------|-----------|-------------|
| cg06941159  | 6        | 30038396 | Highest math class taken (MTAG)    | 16p11.2   | Intergenic  |
| cg20910361  | 8        | 36224396 | Height                             | 4q31.21   |             |

# Considerations/Limitations
- Some keyword captures like "Disease of mental health" were not used during earlier gene annotations, thus Harmonizome extraction is not comprehensive (identified by manual exploration of PRDM15 and the `DISEASES Experimental Gene-Disease Association Evidence Sores 2025` data set)
- Many publicly-available datasets identified during this research (PsyGeNET, MRC-IEU EWAS Catalog, etc.) were difficult to retrieve or not available via the web through any devices/internets i tried accessing them by. Furthermore, datasets I used to build this extraction, like NHGRI-EBI GWAS Catalog's web browser and web endpoints for their downloadable datasets are no longer available, within a few months span of my accessing of them. This highly volatile methodology limits reproducibility and makes me question if I should have attempted a more robust method of navigating API-only datasets like MRC-IEU's OpenGWAS; while I subconsciously filtered databases with UI web-endpoints first to identify extractable variables and validate my programming process at each step, this effectively biased my original datasets I used to capture. Thankfully, I eventually got DISGENET API access, which is a simpler and more methodologically complete extraction and contains sources from previously-unavailable datasets like PsyGeNET.
- I later was able to access MRC-IEU EWAS Catalog. while their browser service wasn't working, I downloaded their results and studies files to explore. I filtered by CpG to identify relevant studies and their effect sizes, then cross referenced the studies to pull in trait data. This process seemed unfulfilling since the resulting data was (1) based traits generalized to entire epigenome studies and (2) the traits were not very telling. See `misc/combined_ewas_catalog.xlsx` if you're interested.
- Some genes in DISGENET are not listed, like `KLRC4-KLRK1` (confirmed via symbol and entrez ID), while some are listed but have no associations (e.g., LINC01500); this is likely due to being forced to use the limited **Curated** dataset through academic plan; the distinction between unlisted and no associations was not captured, rather a single csv, `misc/disgenet_missing_associations.csv`
- there exists a non-zero amount of genes with identifiable psych associations from DISGENET but zero associations from EBI GWAS Catalog or even Harmonizome, indicating imperfect/incomplete capture
- consider playing with the UCSC GWAS Catalog window to narrow screening of neighboring gene-trait associations
- a veriety of resources required further processing to make sense of the data, sources like the [Psychiatric Genomics Consortium](https://pgc.unc.edu) have [downloadable study data](https://pgc.unc.edu/for-researchers/download-results/) from curated genomic papers categorized by psychiatric condition. Further exploration using these data sets could contribute to the accuracy of gene annotations and illuminate potential epigenetic components. 

# *Appendix*

<h3 style="
        margin-top: 10px;
        margin-bottom: 0px">A. Psychiatric Keywords (GWAS & Harmonizome)</h3>
Specific keywords used to filter trait and disease labels across the GWAS Catalog and Harmonizome datasets; matches are case-sensitive and inclusive of substrings
<details>
<summary>Keywords by Category</summary>

```C
# core psych disorders
"schizophrenia",
"psychosis",
"psychotic",
"bipolar",
"mania",
"mood disorder",
"depression",
"depressive",
"major depressive",
"unipolar",
"dysthymia",

# neurodevelopmental / autism / adhd
"autism",
"asperger",
"pervasive developmental",
"adhd",
"attention deficit",
"hyperactivity",
"conduct disorder",
"oppositional defiant",
"disruptive behavior",

# anxiety / ocd / ptsd / stress
"anxiety",
"panic disorder",
"agoraphobia",
"social phobia",
"obsessive-compulsive",
"obsessive compulsive",
"ocd",
"post-traumatic stress",
"posttraumatic stress",
"ptsd",
"stress-related",
"stress related",

# eating / substance
"eating disorder",
"anorexia nervosa",
"bulimia",
"binge eating",
"substance use",
"substance-use",
"substance abuse",
"drug dependence",
"alcohol dependence",
"alcohol-use disorder",
"nicotine dependence",

# other psych/neuropsychiatric baskets
"personality disorder",
"somatoform",
"tic disorder",
"tourette",
"neuropsychiatric",
"mental disorder",
"mental or behavioural",
"mental or behavioral"
```
</details>

<h3 style="
        margin-top: 10px;
        margin-bottom: 0px">B. Mental Health MeSH, Title, & Genetic Terms (PubMed)</h3>
Medical Subject Headings (MeSH) and fallback title keywords used to identify relevant literature from gene queries via NCBI E-utilities and classify as genetic
<details>
<summary>MeSH Terms by Category</summary>

```C
MESH_TERMS: List[str] = [
# core mood disorders
"Depressive Disorder",
"Depressive Disorder, Major",
"Bipolar Disorder",
"Cyclothymic Disorder",
# psychotic disorders
"Schizophrenia",
"Psychotic Disorders",
"Schizoaffective Disorder",
# anxiety, trauma, stress
"Anxiety Disorders",
"Panic Disorder",
"Phobic Disorders",
"Obsessive-Compulsive Disorder",
"Stress Disorders, Post-Traumatic",
# neurodevelopmental
"Autism Spectrum Disorder",
"Attention Deficit Disorder with Hyperactivity",
"Intellectual Disability",
# substance use
"Substance-Related Disorders",
"Alcohol-Related Disorders",
"Opioid-Related Disorders",
"Cocaine-Related Disorders",
# self-harm/suicide
"Suicide",
"Suicidal Ideation",
# broad
"Mental Disorders",
"Mental Health",
]
```
</details>

<details>
<summary>Title Keywords</summary>

```C
TITLE_TERMS: List[str] = [
"schizophrenia",
"bipolar",
"depressi",     # depression, depressive
"mania",
"psychosis",
"psychotic",
"autism",
"asperger",
"adhd",
"attention-deficit",
"anxiety",
"panic disorder",
"ptsd",
"post-traumatic",
"suicide",
"suicidal",
"substance use",
"substance-use",
"addiction",
"alcohol use",
"alcohol-use",
"mental disorder",
"mental health",
]
```
</details>

<details>
<summary>Genetic Keywords</summary>

```C
GENETIC_MESH_TERMS: List[str] = [
"Polymorphism, Genetic",
"Genetic Variation",
"Genome-Wide Association Study",
"Genotype",
"Genetic Predisposition to Disease",
]
```
</details>