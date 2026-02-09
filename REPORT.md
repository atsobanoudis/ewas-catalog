# General Workflow
This project provides a comprehensive bioinformatics pipeline for annotating gene lists derived from epigenetic studies, specifically targeting psychiatric research. The final output is a multi-layered dataset (`annotated_genes.xlsx`), with the `Augmented` sheet serving as the combined genetic-epigenetic annotation.

<h3 style="margin: 0px">1. Core Annotation (Foundation)</h3>

- **Gene Metadata:** Resolving gene identifiers (Symbols, LOC IDs) to standardized HGNC, Entrez, and Ensembl IDs via NCBI Datasets.
- **Functional Synthesis:** Aggregating biological summaries and "FUNCTION" comments from NCBI, UniProt, and HGNC into a unified functional description.

<h3 style="margin: 0px">2. Evidence Integration (Multi-Omic)</h3>

- **Psychiatric Associations:** Identifying gene-disease links across four major evidence streams:
    - **GWAS Catalog:** Significant SNP-trait associations filtered for psychiatric phenotypes.
    - **Harmonizome:** Cross-database (CTD, GAD, etc.) disease associations enriched for mental health keywords.
    - **PubMed Literature:** Automated mining of high-confidence gene-literature links via MeSH terms and title-based genetic tagging.
    - **DisGeNET:** Integration of curated psychiatric disease-gene associations, providing both association scores and supporting evidence.

<h3 style="margin: 0px">3. Epigenetic Context & Augmentation (Final State)</h3>

- **CpG-to-Gene Mapping:** Explicitly linking CpG probe IDs and chromosomal coordinates back to the annotated gene list.
- **EWAS Atlas Integration:** Pulling trait associations directly associated with the CpG probes to provide a primary epigenetic signal.
- **Broad Association Mapping:** Expanding the dataset beyond psychiatric traits to include a broad-spectrum view of all known gene-disease associations (via DisGeNET), allowing for a comprehensive analysis of pleiotropy and non-psychiatric context.



# Detailed workflows
<h3 style="margin: 0px">Gene Annotation</h3>
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
>[Note:]
*The `disgenet_diseases` column represents a broad extraction of all associated traits, bypassing the mental health filter applied to the `disgenet_psych_diseases` column.*
</details>

<h3 style="margin-top: 10px; margin-bottom: 0px">Epigene Annotation</h3>
<details>
<summary>EWAS Atlas</summary>

**[NGDC EWAS Atlas](https://ngdc.cncb.ac.cn/ewas/atlas)** ([REST API](https://ngdc.cncb.ac.cn/ewas/api)):\
CpG probe ID → trait associations (`data/ewas_atlas.csv`)

1. **Query API** per CpG `probeID` in `ewas_res_groupsig_128.xlsx`:
    - Resolve symbols from `relatedTranscription` as `genes`.
    - Capture `cpgIsland` status as `cpg_island`.
    - Flatten `associationList`, containing `trait`, `correlation`, `rank`, `pmid`.
2. **Methodological Refinement** (Post-Capture):
    - **Rank Scoring:** Calculate `rank_score = rank / total_associations` (3 decimal places). If rank is missing, default to `0.000` for output readability.
    - **Correlation Mapping:** Map raw values: `pos` → `hyper`, `neg` → `hypo`, `NA` → `NR`.
    - **Unmapped Gene Gap Analysis:** Cross-reference all genes in `ewas_atlas.csv` against the validated `symbol` and `synonyms` columns of the core annotation.
3. **Data Integration Logic:**
    - **Formatting:** Traits are aggregated into a single cell, separated by `; \n`, in the format: `[trait], [rank_score], [correlation], [pmid]`.
    - **Sorting:** Associations are ordered by `rank_score` descending, then alphabetically by trait.
    - **Column Generation:** Results are partitioned into `ewas_atlas_traits`, `ewas_unmapped_gene` (established symbols), and `ewas_unmapped_regions` (decimal-named clone/LncRNA loci).
4. **Diagnostic Export:** 
    - `unaccounted_genes.csv`: List of uncaptured genes (no decimals) with their original nearest-gene mapping and synonym status.
    - `appendix_c_source.csv`: List of uncaptured genomic regions (decimal-named).
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

<details>
<summary>Augmentation</summary>

**[EWAS Atlas](https://ngdc.cncb.ac.cn/ewas/atlas)** & **[DisGeNET](https://www.disgenet.org/)**:\
`annotated_genes.xlsx` (living file) + `ewas_res_groupsig_128.xlsx` → `annotated_genes.xlsx` (Augmented sheet)

1. **Comprehensive Join**
    - Joins the PI-provided `ewas_res_groupsig_128.xlsx` chromosomal coordinates and probe IDs to the existing annotated gene list.
    - Logic: **Right Join** on the PI's CpG list to ensure every probe provided is retained, even if it lacks a current gene mapping.
    - Key: Primary join on the `input` column to handle gene name discrepancies between raw data and validated symbols.
2. **EWAS Atlas Enrichment**
    - Cross-references CpG probe IDs with the EWAS Atlas dataset to identify additional trait associations and unmapped loci.
    - Logic: Empty associations are left as nulls to maintain clarity in the living document.
3. **Broad-Spectrum DisGeNET**
    - Pulls all gene-disease associations without psychiatric filtering to provide broad phenotypic context.
    - Note: The `disgenet_disease` column contains these broad extractions.
4. **Columns Added**: `cpg`, `cpg_chr`, `cpg_start`, `cpg_end`, `ewas_atlas_traits`, `ewas_unmapped_gene`, `ewas_unmapped_regions`, `disgenet_disease`
</details>


<br>
<br>

# Interesting findings
Because multiple genes are now linked to the same CpG, filtering by CpG can yield interesting, albeit expected results. Genes *UGTA10* and *UGT1A8* are both close enough to *cg00922271* to be listed as a unique gene, and their disease association profiles are expectedly similar from DISGENET data (see below).
| cpg        | symbol  | disgenet_diseases |
|------------|---------|------------------|
| cg00922271 | UGT1A10 | Gilbert Disease, 0.75;<br>Increased bilirubin level (finding), 0.65;<br>Crigler Najjar syndrome, type 1, 0.6;<br>Crigler-Najjar syndrome, 0.55;<br>BILIRUBIN, SERUM LEVEL OF, QUANTITATIVE TRAIT LOCUS 1, 0.4;<br>Crigler Najjar syndrome, type 2, 0.4;<br>GILBERT SYNDROME, SUSCEPTIBILITY TO, 0.4;<br>Lucey-Driscoll syndrome (disorder), 0.4 |
| cg00922271 | UGT1A8  | Diarrhea, 0.5;<br>Gilbert Disease, 0.45;<br>Crigler Najjar syndrome, type 2, 0.4;<br>BILIRUBIN, SERUM LEVEL OF, QUANTITATIVE TRAIT LOCUS 1, 0.4;<br>Crigler-Najjar syndrome, 0.4;<br>Lucey-Driscoll syndrome (disorder), 0.4;<br>Crigler Najjar syndrome, type 1, 0.4;<br>Increased bilirubin level (finding), 0.4;<br>Diarrheal disorder, 0.4;<br>GILBERT SYNDROME, SUSCEPTIBILITY TO, 0.4 |

### Discovery of Distal Gene-CpG Associations (FMN1 Case Study)
A critical finding during the augmentation phase was the identification of gene associations that significantly exceed standard proximity-based windows. A primary example is **cg02255242**, which the EWAS Atlas robustly associates with the gene **FMN1**. 
- **Genomic Context:** FMN1 is located at [[INSERT_DETAIL: FMN1_START]] - [[INSERT_DETAIL: FMN1_END]] on chromosome 15, while cg02255242 resides at 33,128,710 - 33,128,712. This represents a distal association of approximately **70.9kb**, which would be missed by traditional 5kb or 10kb proximity scans.
- **Evidence Profile:** EWAS Atlas reports a [[INSERT_DETAIL: FMN1_PERCENT_BODY]]% body mapping for this CpG-gene pair with a 100% hypermethylated correlation based on current data.
- **Trait Overlap:** FMN1 itself is linked to high-impact traits including **Alzheimer's Disease** (5 associations), **Smoking**, **Mild Cognitive Impairment**, and **Down Syndrome**. This discovery validates the shift from "Nearest Gene" identification to "EWAS-Associated Gene" capture.

### Uncaptured Genes and Genomic Loci
The gap analysis between EWAS Atlas and our core annotation identified several "uncaptured" genes that warrant further exploration.
- **Unaccounted Genes:** We identified [[INSERT_DETAIL: COUNT_UNACCOUNTED]] established symbols present in the EWAS Atlas but missing from our initial mapping (see `unaccounted_genes.csv`).
- **Uncharacterized Loci (Appendix C):** A significant number of probes are associated with uncharacterized long non-coding RNA (lncRNA) loci or genomic regions named using the Human Genome Project clone system (e.g., `AP003039.3`). These represent potential "dark matter" in the psychiatric epigenome.


# Considerations and Methodological Limitations

This research involves the synthesis of multi-omic data from disparate public repositories. While every effort was made to ensure accuracy and completeness, several methodological considerations and limitations must be acknowledged.

### 1. Shift from Proximity to Association
The discovery of distal associations (e.g., FMN1 at 70.9kb) demonstrates that relying solely on "nearest gene" mapping can overlook biologically relevant and validated signals. This highlights a paradigm shift in our approach: capturing genes with identifiable EWAS-associated links to our CpG sites of interest provides a more comprehensive view of the epigenetic landscape than simple genomic proximity.

### 2. Rank-Based Scoring and Contextual Gaps
The `rank_score` metric utilized in the EWAS Atlas integration is limited to associations where an explicit rank was reported in the source study. This necessarily excludes contextual information for traits with large total association counts but unreported ranks (e.g., **cg02255242** and its hypermethylated association with infertility, which lacks a rank among 2,751 associations in PMID 25753583).

### 3. Statistical Context and Trait Broadening
To facilitate a broad-spectrum view, associations were included even when p-values exceeded traditional significance thresholds. In the current augmented set, only [[INSERT_DETAIL: COUNT_SIG_P]] hits were below 1e-04. While many of these "lower significance" traits (e.g., bariatric surgery) may be of less immediate clinical relevance, the decision to broaden the capture could lead to misinterpretation without the relevant statistical context.

### 4. Data Provenance and Snapshot Volatility
The bioinformatics landscape is characterized by high volatility in dataset availability and API accessibility. Results presented here represent a specific point-in-time snapshot. Future efforts to reproduce these results may encounter discrepancies due to database versioning or the deprecation of legacy web services (e.g., the recent shifts in the NHGRI-EBI GWAS Catalog).


<br>

# *Appendix*

<h3 style="margin: 0px">A. Psychiatric Keywords (GWAS & Harmonizome)</h3>
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

<h3 style="margin-top: 10px; margin-bottom: 0px">B. Mental Health MeSH, Title, & Genetic Terms (PubMed)</h3>
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

<h3 style="margin-top: 10px; margin-bottom: 0px">C. Uncharacterized Genomic Regions (EWAS Atlas)</h3>
Detailed list of CpG sites associated with uncharacterized loci, including long non-coding RNA (lncRNA) and genomic regions named using the Human Genome Project clone system (e.g., RP11, AC nomenclature).

| CpG ID | Associated Loci (Atlas) |
|---|---|
| cg07252486 | AP003039.3 |
| cg12790145 | RP11-283I3.2 |
| cg02780130 | RP11-283I3.2 |
| cg25329573 | RP11-283I3.2 |
| cg09911534 | RP11-622C24.2 |
| cg16469117 | AC007092.1 |
| cg17240725 | AC073846.5 |
| cg04070200 | AC007879.5 |
| cg06941159 | RP11-22P6.2 |
| cg13689053 | AP000688.14 |
| cg12658972 | RP11-95P2.3 |
| cg01267120 | AL157871.2, RP11-638I2.6 |
| cg02834909 | RP11-316F12.1 |
| cg05406088 | RP11-321F6.1 |
| cg02282594 | RP11-298I3.1 |
| cg05860956 | RP11-283I3.2 |
| cg01334824 | RP11-100M12.3 |
| cg24468780 | RP11-112J1.2 |
| cg15587955 | RP11-1080G15.1 |
| cg24947255 | XXbac-BPG181B23.6 |
| cg14545602 | AC005498.3 |
| cg23102195 | CTD-2269E23.4 |
| cg11340603 | RP1-78O14.1 |
| cg23184739 | AC110781.3 |
| cg13303179 | RP11-316F12.1 |
| cg10075163 | RP11-283I3.2 |
| cg20188212 | RP11-283I3.2 |
| cg15501526 | RP11-526P5.2 |
| cg05715076 | RP11-543C4.1 |
| cg09535960 | RP11-177H13.2 |
| cg10957166 | AC084018.1, RP11-347I19.7, RP11-347I19.8 |