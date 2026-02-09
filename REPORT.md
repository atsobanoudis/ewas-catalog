---
title: Psychiatric Gene-Epigene Annotation
author: Alex Tsobanoudis
date: 20026-02-08
---
# Psychiatric Gene-Epigene Annotation
## General Workflow
This project provides a comprehensive bioinformatics pipeline for annotating gene lists derived from epigenetic studies, specifically targeting psychiatric research. The final output is a multi-layered dataset (`annotated_genes.xlsx`), with the `augmented` sheet serving as the combined genetic-epigenetic annotation

<h3 style="margin: 0px">1. Core annotation</h3>

\- Gene metadata: resolving gene identifiers (Symbols, LOC IDs) to standardized HGNC, Entrez, and Ensembl IDs via NCBI Datasets.\
\- Functional synthesis: aggregating biological summaries and "FUNCTION" comments from NCBI, UniProt, and HGNC into a unified functional description.

<h3 style="margin: 0px">2. Evidence Integration</h3>

Psychiatric associations: identifying gene-disease links across four major evidence streams\
1\. GWAS Catalog: significant SNP-trait associations filtered for psychiatric phenotypes\
2\. Harmonizome: Cross-database (CTD, GAD, etc.) disease associations enriched for mental health keywords\
3\. PubMed Literature: Automated mining of high-confidence gene-literature links via MeSH terms and title-based genetic tagging\
4\. DisGeNET: Integration of curated psychiatric disease-gene associations, providing both association scores and supporting evidence

<h3 style="margin: 0px">3. Epigenetic Context & Augmentation</h3>

\- CpG-to-Gene mapping: explicitly links CpG probe IDs and chromosomal coordinates back to the annotated gene list\
\- EWAS Atlas integration: pulling trait associations directly associated with the CpG probes to provide a primary epigenetic signal.\
\- Broad association mapping: expanding the dataset beyond psychiatric traits to include a broad-spectrum view of all known gene-disease associations (via DisGeNET), allowing for a comprehensive analysis of pleiotropy and non-psychiatric context.



## Detailed workflows
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
<summary>EWAS Atlas</summary>

**[NGDC EWAS Atlas](https://ngdc.cncb.ac.cn/ewas/atlas)** ([REST API](https://ngdc.cncb.ac.cn/ewas/api)):\
CpG probe ID → trait associations (`data/ewas_atlas.csv`)

1. Query API per CpG `probeID` in `ewas_res_groupsig_128.xlsx`:
    - Resolve symbols from relatedTranscription as `genes`
    - Capture cpgIsland status as `cpg_island`
    - Flatten `associationList`, containing `trait`, `correlation`, `rank`, `pmid`
2. Methodological refining
    - Rank scoring: calculate `rank_score = rank / total_associations` (3 decimal places). If rank is missing, empty value for output readability
    - Correlation (methylation) mapping:
        - `pos` ↔ `hyper`
        - `neg` ↔ `hypo`
        - `NA` ↔ `NR`
    - Unmapped gene gap: cross-reference all genes in `ewas_atlas.csv` against validated `symbol` and `synonyms` columns from `annotated_genes.xlsx` *and* the nearest gene (`unique_gene_name`) from `ewas_res_groupsig_128.xlsx`
3. Data integration logic
    - Formatting: traits are aggregated into a single cell, separated by `; \n`, in the format: `[trait], [rank_score], [correlation], [pmid]`, sorted by descending rank score, then alphabetically by trait
    - Columns: results partitioned into `ewas_atlas_traits`, `ewas_unmapped_gene` (established symbols), and `ewas_unmapped_regions` (decimal-named clone/LncRNA loci; "AP003039.3" vs. "NTM"). Established symbols table can be found in Table 1, while uncharacterized regions are in *Appendix C* Table C.1
    - Example `ewas_atlas_traits` output
        ```
        preterm birth, 0.9, hyper, 1239847;
        allergies, 0.8, hypo, 10923874;
        Alzheimer, 0.8, NR, 019834
        ```
</details>

<details>
<summary>Augmentation</summary>

**[EWAS Atlas](https://ngdc.cncb.ac.cn/ewas/atlas)** & **[DisGeNET](https://www.disgenet.org/)**:\
`annotated_genes.xlsx` → 
```python
annotated_genes.xlsx[
    sheets = "annotated_genes",   # original
             "ewas_atlas",        # copy of data/ewas_atlas.csv
             "augmented",         # NEW: combined sheet
             "ewas_res_groupsig_128"
]
```

1. Comprehensive join
    - Joins the provided `ewas_res_groupsig_128.xlsx` chromosomal coordinates and probe IDs to new augmented sheet
    - Logic: right join provided CpG list to ensure every probe is retained regardless of gene annotation
    - Primary join on the `input` column to handle gene name discrepancies between raw data (nearest genes) and validated symbols (via Entrez → HGNC queries in above section)
2. EWAS Atlas enrichment
    - Cross-references CpG probe IDs with the EWAS Atlas dataset to identify additional trait associations and unmapped genes/loci
    - Logic: empty associations are left as nulls to maintain clarity in the sheet, associations without a rank, and subsequently a rank score, are given a placeholder score of `0.000`
3. Broad-spectrum DisGeNET: pulls all gene-disease associations without psychiatric filtering to provide broad phenotypic context in the `disgenet_disease` column
4. Columns added: `cpg`, `cpg_chr`, `cpg_start`, `cpg_end`, `ewas_atlas_traits`, `ewas_unmapped_gene`, `ewas_unmapped_regions`, `disgenet_disease`
</details>


<br>
<br>

## Findings

<h3 style="margin-top: 10px; margin-bottom: 0px">DISTAL GENE-CPG UNMAPPED ASSOCIATIONS</h3>

A critical finding during the augmentation phase was the identification of CpG-gene associations via EWAS Atlas that were previously unmapped to genes from proximity identification. Established CpG-gene associations (e.g., GSE1, CNTD2) were recorded and displayed in Table 1; uncharacterized long non-coding RNA (lncRNA) loci or genomic regions named using the Human Genome Project clone system (e.g., `AP003039.3`) can be found in *Appendix C*.\
A primary example is *cg02255242*, associated with the gene *FMN1*. The probe *cg02255242* is located at `chr15:33128710–33128712`, approximately 70.9kb from the FMN1 transcription start site (`33057746`). Despite this distance, EWAS Atlas identifies a validated association from a single study characterized by 100% hypermethylation in the gene body. According to EWAS Atlas, FMN1 is broadly implicated in the epigenetic landscape, associated with 42 unique probes and 36 traits; notably, these include high-priority conditions such as Alzheimer’s disease (ranked 3rd; 5 associations) and Mild Cognitive Impairment (tied for 1st; 6 associations).

<p style="margin-top: 20px;"><b>Table 1: Unmapped Gene Against Known Synonyms and Originating Input</b></p>
<table style="border-collapse: collapse; width: 100%; border-top: 1px solid black; border-bottom: 1px solid black;">
  <thead>
    <tr style="border-bottom: 1px solid black;">
      <th style="text-align: left;">CpG</th>
      <th style="text-align: left;">EWAS Atlas Genes</th>
      <th style="text-align: left;">Unmapped Genes</th>
      <th style="text-align: left;">Synonym</th>
      <th style="text-align: left;">Nearest Gene, Original Input</th>
    </tr>
  </thead>
  <tbody>
    <tr><td>cg15035382</td><td>CNTD2</td><td>CNTD2</td><td></td><td>CNTD2</td></tr>
    <tr><td>cg26456563</td><td>CNTD2</td><td>CNTD2</td><td></td><td>CNTD2</td></tr>
    <tr><td>cg10957166</td><td>AC084018.1;RHOF;RP11-347I19.7;RP11-347I19.8</td><td>RHOF</td><td></td><td>LOC338799</td></tr>
    <tr><td>cg24315757</td><td>MAGI1;MAGI1-AS1;MAGI1-IT1</td><td>MAGI1-IT1</td><td></td><td>MAGI1-AS1, MAGI1</td></tr>
    <tr><td>cg10638439</td><td>HCP5;MICA;Y_RNA</td><td>Y_RNA</td><td></td><td>MICA</td></tr>
    <tr><td>cg10260205</td><td>PCBP3;PRED62</td><td>PRED62</td><td></td><td>PCBP3</td></tr>
    <tr><td>cg25857471</td><td>DPCR1;SFTA2</td><td><u>SFTA2</u></td><td></td><td>DPCR1</td></tr>
    <tr><td>cg02255242</td><td>FMN1</td><td><u>FMN1</u></td><td></td><td></td></tr>
    <tr><td>cg16185115</td><td>GSE1</td><td><u>GSE1</u></td><td></td><td></td></tr>
    <tr><td>cg04732357</td><td>ANKRD36C</td><td><u>ANKRD36C</u></td><td></td><td></td></tr>
    <tr><td>cg21088344</td><td>SLFN12L</td><td><u>SLFN12L</u></td><td></td><td></td></tr>
    <tr><td>cg16910670</td><td>NADSYN1</td><td><u>NADSYN1</u></td><td></td><td></td></tr>
    <tr><td>cg13796823</td><td>SAMD3;TMEM200A</td><td><u>SAMD3</u></td><td></td><td></td></tr>
    <tr><td>cg13796823</td><td>SAMD3;TMEM200A</td><td><u>TMEM200A</u></td><td></td><td></td></tr>
  </tbody>
</table>
<p style="font-size: 0.9em;">Genes (<i>Unmapped Genes</i>) extracted from EWAS Atlas by provided CpG ProbeID via REST API were cross-referenced against known gene aliases from previous annotation (<i>Synonym</i>) and initial suspected genes by proximity(<i>Nearest Gene</i>, <i>Original Input</i>). Underlined genes (<i>SFTA2</i>, <i>FMN1</i>, <i>GSE1</i>, <i>ANKRD36C</i>, <i>SLFN12L</i>, <i>NADSYN1</i>, <i>SAMD3</i>, <i>TMEM200A</i>) are not listed in error (<i>CNTD2</i>, <i>RHOF</i>), possible noise (<i>MAGI1-IT1</i>, <i>Y_RNA</i>), or deprecated (<i>PRED62</i>), and thus serve as a potentially significant expansions for exploration. Genes with cross-referenced synonyms accounted for in <code>annotated_genes.xlsx</code> (<i>WBSCR17</i>, <i>DPCR1</i>, <i>PRAMEF23</i>, <i>B3GNTL1</i>, <i>RARS</i>, <i>FAM134B</i>) were removed, while the <i>Synonym</i> variable remained to guide interpretation</p>

<h3 style="margin-top: 10px; margin-bottom: 0px">PSYCHIATRIC AND EPIGENETIC ASSOCIATIONS</h3>

A synthesis of the DisGeNET and EWAS Atlas integrations highlights genes with strong evidence for psychiatric disorders and significant epigenetic trait clustering. Table 2 details genes with the highest confidence scores (Score > 0.5) for psychiatric conditions, led by <i>KCNN3</i> for Schizophrenia. Table 3 illustrates the "pleiotropic load" of these genes, ranking them by the count (traits ≥5) of distinct epigenetic traits associated with their linked CpGs. Notably, <i>BICDL3P</i> and <i>ARAP1</i> exhibit high trait counts, suggesting they may act as broader epigenetic hubs beyond their specific disease associations.

<p style="margin-top: 20px;"><b>Table 2: Top Psychiatric Associations (DisGeNET)</b></p>
<table style="border-collapse: collapse; width: 100%; border-top: 1px solid black; border-bottom: 1px solid black;">
  <thead>
    <tr style="border-bottom: 1px solid black;">
      <th style="text-align: left;">Gene</th>
      <th style="text-align: left;">Disease</th>
      <th style="text-align: left;">Score</th>
    </tr>
  </thead>
  <tbody>
    <tr><td>KCNN3</td><td>Schizophrenia</td><td>0.9</td></tr>
    <tr><td>MAD1L1</td><td>Schizophrenia</td><td>0.8</td></tr>
    <tr><td>VIPR2</td><td>Schizophrenia</td><td>0.75</td></tr>
    <tr><td>JMJD1C</td><td>Intellectual Disability</td><td>0.65</td></tr>
    <tr><td>JMJD1C</td><td>Autism Spectrum Disorders</td><td>0.65</td></tr>
    <tr><td>JMJD1C</td><td>Autistic Disorder</td><td>0.65</td></tr>
    <tr><td>MAGI1</td><td>Bipolar Disorder</td><td>0.65</td></tr>
    <tr><td>THSD7A</td><td>Bipolar Disorder</td><td>0.65</td></tr>
    <tr><td>BCR</td><td>Bipolar Disorder</td><td>0.6</td></tr>
    <tr><td>BCR</td><td>Unipolar Depression</td><td>0.6</td></tr>
    <tr><td>BCR</td><td>Major Depressive Disorder</td><td>0.6</td></tr>
    <tr><td>KCNN3</td><td>Bipolar Disorder</td><td>0.6</td></tr>
    <tr><td>MAGI1</td><td>Schizophrenia</td><td>0.6</td></tr>
    <tr><td>SLC6A12</td><td>Schizophrenia</td><td>0.6</td></tr>
    <tr><td>SYT1</td><td>Neurodevelopmental Disorders</td><td>0.6</td></tr>
  </tbody>
</table>
<p style="font-size: 0.9em;">Genes with the strongest evidence for psychiatric disorders (Score &gt; 0.5) extracted from DisGeNET (<code>data/disgenet_gda.csv</code>). The <i>Score</i> represents the confidence level of the gene-disease association, with higher values indicating stronger evidence from curated sources. <i>KCNN3</i> shows the highest association with <i>Schizophrenia</i>.</p>

<p style="margin-top: 20px;"><b>Table 3: Top Epigenetic Trait Clusters (EWAS Atlas)</b></p>
<table style="border-collapse: collapse; width: 100%; border-top: 1px solid black; border-bottom: 1px solid black;">
  <thead>
    <tr style="border-bottom: 1px solid black;">
      <th style="text-align: left;">CpG</th>
      <th style="text-align: left;">Gene</th>
      <th style="text-align: left;">Trait Count</th>
      <th style="text-align: left;">Top Traits</th>
    </tr>
  </thead>
  <tbody>
    <tr><td>cg17240725</td><td>BICDL3P</td><td>8</td><td>ancestry, 0.874, hypo, 30563547;<br> papillary thyroid carcinoma, 0.378, hypo, 28469731;<br> colorectal laterally spreading tumor, 0.330, hypo, 30183087;<br> papillary thyroid carcinoma, 0.285, hypo, 28938489;<br> leukoaraiosis (LA), 0.241, hypo, 29875652;<br> Alzheimer’s Disease, 0.051, hypo, 37388929;<br> papillary thyroid carcinoma, 0.000, hypo, 31328658;<br> ancestry, 0.000, NR, 31399127</td></tr>
    <tr><td>cg15035382</td><td>ARAP1</td><td>7</td><td>Sex-specific DNA methylation, 0.388, hyper, 37300819;<br> Sex-specific DNA methylation, 0.378, hyper, 37300819;<br> gestational diabetes mellitus, 0.211, hyper, 29596946;<br> sexually dimorphic, 0.161, hyper, 37608375;<br> gender, 0.112, hypo, 26553366;<br> infant sex, 0.088, hypo, 34569420;<br> infant sex, 0.060, hypo, 34569420</td></tr>
    <tr><td>cg09535960</td><td>LOXL2</td><td>6</td><td>alcohol consumption, 0.600, hyper, 27843151;<br> papillary thyroid carcinoma, 0.263, hypo, 28469731;<br> adenoma, 0.224, hypo, 30183087;<br> papillary thyroid carcinoma, 0.203, hypo, 28938489;<br> prostate cancer, 0.188, hyper, 29740534;<br> papillary thyroid carcinoma, 0.000, hypo, 31328658</td></tr>
    <tr><td>cg14906510</td><td></td><td>6</td><td>down syndrome, 0.939, hypo, 29601581;<br> ancestry, 0.927, hyper, 30563547;<br> Alzheimer’s Disease, 0.618, hyper, 37388929;<br> bariatric surgery, 0.442, hyper, 37834223;<br> bariatric surgery, 0.442, hyper, 37834223;<br> Alzheimer's disease (AD), 0.007, hypo, 37388929</td></tr>
    <tr><td>cg05715076</td><td></td><td>5</td><td>Alzheimer’s Disease, 0.608, hyper, 37388929;<br> ancestry, 0.323, hypo, 33108347;<br> mild cognitive impairment, 0.105, hyper, 37388929;<br> ancestry, 0.081, hyper, 30563547;<br> ankylosing spondylitis, 0.066, hyper, 31128893</td></tr>
    <tr><td>cg09506675</td><td>GPR85</td><td>5</td><td>Alzheimer’s Disease, 0.505, hyper, 37388929;<br> perinatally-acquired HIV, 0.487, hyper, 31324826;<br> maternal smoking, 0.361, hyper, 27040690;<br> mild cognitive impairment, 0.174, hyper, 37388929;<br> maternal smoking, 0.072, hyper, 31536415</td></tr>
    <tr><td>cg16744531</td><td>B3GNT3</td><td>5</td><td>perinatally-acquired HIV, 0.993, hyper, 31324826;<br> aging, 0.950, hypo, 29064478;<br> acute myelocytic leukemia (AML), 0.542, hypo, 39052947;<br> respiratory allergies (RA), 0.000, hypo, 26999364;<br> pre- and post-lenalidomide treatm ..., 0.000, NR, 33903952</td></tr>
    <tr><td>cg27031099</td><td></td><td>5</td><td>systemic lupus erythematosus (SLE), 0.862, hypo, 29437559;<br> systemic lupus erythematosus (SLE), 0.804, hypo, 31428085;<br> Parkinson's disease (PD), 0.785, hypo, 28851441;<br> air pollution (NO2), 0.467, hyper, 29410382;<br> smoking, 0.170, hypo, 31552803</td></tr>
  </tbody>
</table>
<p style="font-size: 0.9em;">Genes with the highest number of associated traits in EWAS studies. <i>Trait Count</i> indicates the number of unique traits associated with the gene via CpG probes. <i>Top Traits</i> lists a selection of these traits with their association details (trait, rank score, methylation status, PMID).</p>

<h3 style="margin-top: 10px; margin-bottom: 0px">MISCELLANEOUS FINDINGS</h3>

Although selection was based solely on genomic proximity while exploring the UCSC screening data (`data/ewas_ucsc_annotated.xlsx`), the two closest traits of 765 interestingly ranged from a classical biological phenotype to an educational attainment proxy, namely the highest mathematics course completed (see Table 4).

<p style="margin-top: 20px;"><b>Table 4: Proximal Traits Identified via UCSC Screening</b></p>
<table style="border-collapse: collapse; width: 100%; border-top: 1px solid black; border-bottom: 1px solid black;">
  <thead>
    <tr style="border-bottom: 1px solid black;">
      <th style="text-align: left;">CpG ID</th>
      <th style="text-align: left;">Distance (bp)</th>
      <th style="text-align: left;">PubMed ID</th>
      <th style="text-align: left;">Trait</th>
      <th style="text-align: left;">Region</th>
      <th style="text-align: left;">Genes</th>
    </tr>
  </thead>
  <tbody>
    <tr><td>cg06941159</td><td>6</td><td>30038396</td><td>Highest math class taken (MTAG)</td><td>16p11.2</td><td>Intergenic</td></tr>
    <tr><td>cg20910361</td><td>8</td><td>36224396</td><td>Height</td><td>4q31.21</td><td></td></tr>
  </tbody>
</table>
<p style="font-size: 0.9em;">Top associations identified by genomic proximity within 10bp to the query CpG sites using UCSC screening data (<code>data/ewas_ucsc_annotated.xlsx</code>). Traits range from behavioral proxies (<i>Highest math class taken</i>) to physical phenotypes (<i>Height</i>). <i>Distance</i> represents the offset in base pairs from the CpG site. <i>Genes</i> are noted as <i>Intergenic</i> when the CpG falls outside defined gene bodies within the specified <i>Region</i>.</p>

Because multiple genes are now linked to the same CpG, filtering by CpG in `annotated_genes.xlsx, sheet = 'augmented'` can yield interesting, albeit expected results. Genes *UGTA10* and *UGT1A8* are both close enough to *cg00922271* to be listed as a unique gene, and their disease association profiles are expectedly similar from DISGENET data (see Table 5).
<p style="margin-top: 20px;"><b>Table 5: DisGeNET Disease Profiles for Co-Located Genes</b></p>
<table style="border-collapse: collapse; width: 100%; border-top: 1px solid black; border-bottom: 1px solid black;">
  <thead>
    <tr style="border-bottom: 1px solid black;">
      <th style="text-align: left;">CpG</th>
      <th style="text-align: left;">Symbol</th>
      <th style="text-align: left;">DisGeNET Diseases (Score)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
        <td style="vertical-align: top;">cg00922271</td>
        <td style="vertical-align: top;">UGT1A10</td>
        <td>Gilbert Disease, 0.75;<br>Increased bilirubin level (finding), 0.65;<br>Crigler Najjar syndrome, type 1, 0.6;<br>Crigler-Najjar syndrome, 0.55;<br>BILIRUBIN, SERUM LEVEL OF, QUANTITATIVE TRAIT LOCUS 1, 0.4;<br>Crigler Najjar syndrome, type 2, 0.4;<br>GILBERT SYNDROME, SUSCEPTIBILITY TO, 0.4;<br>Lucey-Driscoll syndrome (disorder), 0.4</td>
    </tr>
    <tr>
        <td style="vertical-align: top;">cg00922271</td>
        <td style="vertical-align: top;">UGT1A8</td>
        <td>Diarrhea, 0.5;<br>Gilbert Disease, 0.45;<br>Crigler Najjar syndrome, type 2, 0.4;<br>BILIRUBIN, SERUM LEVEL OF, QUANTITATIVE TRAIT LOCUS 1, 0.4;<br>Crigler-Najjar syndrome, 0.4;<br>Lucey-Driscoll syndrome (disorder), 0.4;<br>Crigler Najjar syndrome, type 1, 0.4;<br>Increased bilirubin level (finding), 0.4;<br>Diarrheal disorder, 0.4;<br>GILBERT SYNDROME, SUSCEPTIBILITY TO, 0.4</td>
    </tr>
  </tbody>
</table>
<p style="font-size: 0.9em;">Comparison of selected disease association profiles for <i>UGT1A10</i> and <i>UGT1A8</i>, two genes linked to the same CpG probe (<i>cg00922271</i>) in the augmented dataset (<code>annotated_genes.xlsx</code>). Disease associations were retrieved from DisGeNET and ordered by association score. The overlap in phenotypes (e.g., <i>Gilbert Disease</i>, <i>Crigler-Najjar syndrome</i>) highlights the functional similarity of these clustered genes.</p>

## Discussion

This research involves the synthesis of multi-omic data from disparate public repositories. While every effort was made to ensure accuracy and completeness, several methodological considerations and limitations must be acknowledged to guide interpretations.

The bioinformatics landscape is characterized by high volatility in dataset availability and API accessibility. During the course of this project, several key resources (e.g., specific web-endpoints for the NHGRI-EBI GWAS Catalog) became unavailable or shifted to API-only access. This volatility extends to gene identifiers themselves; for instance, the unmapped gene *PRED62* (associated with `cg10260205`) could not be indexed via HGNC or Entrez. Manual tracing through the [ExPheWas Browser](https://exphewas.statgen.org/) identified it as Ensembl ID *ENSG00000268040*, which was revealed to be a deprecated identifier archived in February 2014. Additionally, automated extraction logic occasionally yielded "ghost" associations: `cg10957166` returned a multitude of genes (e.g., *AC084018.1*; *RHOF*; *RP11-347I19.7*; *RP11-347I19.8*) via the EWAS Atlas logic, yet manual validation via the browser confirmed zero associated genes. Consequently, the results presented here represent a specific point-in-time snapshot, and future reproduction efforts may encounter similar discrepancies.

Beyond temporal volatility, the associations identified via the EWAS Atlas, GWAS Catalog, and PubMed mining are primarily descriptive. While a high correlation or frequent literature co-occurrence suggests a potential biological link, these metrics do not imply direct mechanistic causation. Epigenetic signals (CpG methylation) and genomic variants (SNPs) often act in complex regulatory networks where proximity to a gene’s transcription start site does not always equate to functional regulation.

This distinction is underscored by the discovery of unmapped associations, like *FMN1*, which demonstrates that relying on proximity mapping can overlook biologically relevant signals. For example, `cg02255242` is located at `chr15:33128710–33128712`, approximately 70.9kb from the *FMN1* start site (`33057746`). Despite this significant distance, EWAS Atlas contains association of our probe to *FMN1* from 1 study, which showed 100% hypermethylation in the gene body. More importantly, *FMN1* is associated with 42 probes and 36 traits via EWAS Atlas, including high-priority conditions like Alzheimer’s disease (ranked 3rd; 5 associations) and Mild Cognitive Impairment (tied for 1st; 6 associations). This finding effectively transforms our approach from proximity scanning to capturing genes with identifiable, validated EWAS-associated links.

Similar interpretative caution is needed for DisGeNET association scores and the `rank_score` metric utilized in the EWAS Atlas integration. The `rank_score` is restricted to associations where an explicit rank was reported in the source study, which excludes contextual information for traits with large total association counts but unreported ranks. A pertinent example is the aforementioned `cg02255242`: while it has a hypermethylated association with "infertility" in PMID 25753583, it lacks a rank among the study's 2,751 associations, leaving it without a score in our augmented dataset.

These evidentiary constraints are further influenced by statistical thresholds. To facilitate a broad-spectrum view, associations were included even when p-values exceeded traditional significance thresholds; in the current augmented set, only 11 hits were below 1e-04, and only one was below 1e-01. While this decision broadens the capture to include potentially relevant traits (e.g., Alzheimer's, diabetes), it introduces the risk that traits of lower clinical relevance to this study (e.g., bariatric surgery, Gulf War illness) could be misinterpreted without strict statistical context.

Search parameters were also shaped by dataset accessibility and bias. Early gene annotations utilized a specific set of psychiatric keywords (see *Appendix A*). Manual exploration revealed that broader terms (e.g., "Diseases of mental health") might yield additional hits in repositories like the Harmonizome. Furthermore, the decision to rely on HGNC as a fallback mechanism for gene resolution introduced specific extraction errors. For instance, *CNTD2* was incorrectly resolved to *ARAP* (aliases: *CENTD2*, *cnt-d2*) via Entrez, whereas HGNC correctly identifies it as *CCNP* (aliases: *CNTD2*). While HGNC is authoritative for established genes, it fails on loci like *LOC338799* (requiring prefix removal processing) and the aforementioned, deprecated *PRED662*. The decision to prioritize capture via Entrez without a secondary verification step compromised specificity.

Finally, a gap analysis reveals that a non-zero number of genes exhibit significant psychiatric associations in DisGeNET but remain absent from the GWAS Catalog or Harmonizome results. This indicates a potential refinement in earlier extraction logic and that, while the pipeline is effective at capturing well-established links, the "dark matter" of gene-disease associations requires further investigation. Conversely, there exists expected gene-disease associations that are not indexed in DisGeNET. While the curated database serves as a high-fidelity filter, distinguishing established links from low-confidence predictions, this exclusivity introduces a bias against less-characterized or complex genomic entities. Specifically, a subset of genes (`misc/disgenet_missing_associations.csv`) returned no data during the enrichment phase. This void arises from two distinct scenarios: genes that are indexed but lack associations within the specific curated subset, and those that are entirely absent from the DisGeNET index (e.g., *KLRC4-KLRK1*). Crucially, the current extraction logic does not differentiate between these two states, effectively flattening a database gap and a biological null result into a single error class. This data requires further exploration but suggests the pipeline’s reliance on the curated DisGeNET subset may be overly conservative, potentially discarding novel candidates that simply lack the historical literature density required for curation.

**Immediate Data Validation Required**

* Rectify Entrez-HGNC resolution: the reliance on Entrez for gene annotation and resolution (e.g., *CNTD2* resolving to *ARAP*) requires a new verification layer. Future logic must confirm Entrez results against HGNC rather than accepting the first Entrez result as absolute.
* Investigate ghost associations: the anomaly with `cg10957166` (yielding *AC084018.1*, *RHOF*, etc., despite no browser evidence), while seemingly isolated, must be investigated to determine if this is an API parsing error or a database inconsistency

**Future Strategic Improvements**

* Reprocess and annotate newly-identified CpG-gene associations to increase granularity and accuracy of data, and potentially illuminate clinically-relevant psychiatric associations.
* Analyze proximity of other unmapped CpG-gene assocations beyond *FMN1*. Explore reasons they may have been missed and generate more descriptive, locational statistics for mapped and unmapped CpG-gene assocations. 
* Adjusting the 5kb window in the UCSC neighboring GWAS screening to better capture distal regulatory elements.
* Direct processing of raw summary statistics from curated [Psychiatric Genomics Consortium](https://pgc.unc.edu) (PGC) studies to enhance annotation accuracy beyond literature-level metadata. [Downloadable study data](https://pgc.unc.edu/for-researchers/download-results/) are available, categorized by psychiatric condition, but require more complex data processing than I felt was relevant.
<br>

## *Appendix*

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
Detailed list of CpG sites associated with uncharacterized loci, including long non-coding RNA (lncRNA) and genomic regions named using the Human Genome Project clone system (e.g., RP11, AC nomenclature)

<p style="margin-top: 20px;"><b>Table C.1: Uncharacterized Genomic Regions (EWAS Atlas)</b></p>
<table style="border-collapse: collapse; width: 100%; border-top: 1px solid black; border-bottom: 1px solid black;">
  <thead>
    <tr style="border-bottom: 1px solid black;">
      <th style="text-align: left;">CpG ID</th>
      <th style="text-align: left;">Associated Loci (Atlas)</th>
    </tr>
  </thead>
  <tbody>
    <tr><td>cg01267120</td><td>AL157871.2, RP11-638I2.6</td></tr>
    <tr><td>cg01334824</td><td>RP11-100M12.3</td></tr>
    <tr><td>cg02282594</td><td>RP11-298I3.1</td></tr>
    <tr><td>cg02780130</td><td>RP11-283I3.2</td></tr>
    <tr><td>cg02834909</td><td>RP11-316F12.1</td></tr>
    <tr><td>cg04070200</td><td>AC007879.5</td></tr>
    <tr><td>cg05406088</td><td>RP11-321F6.1</td></tr>
    <tr><td>cg05715076</td><td>RP11-543C4.1</td></tr>
    <tr><td>cg05860956</td><td>RP11-283I3.2</td></tr>
    <tr><td>cg06941159</td><td>RP11-22P6.2</td></tr>
    <tr><td>cg07252486</td><td>AP003039.3</td></tr>
    <tr><td>cg09535960</td><td>RP11-177H13.2</td></tr>
    <tr><td>cg09911534</td><td>RP11-622C24.2</td></tr>
    <tr><td>cg10075163</td><td>RP11-283I3.2</td></tr>
    <tr><td>cg10957166</td><td>AC084018.1, RP11-347I19.7, RP11-347I19.8</td></tr>
    <tr><td>cg11340603</td><td>RP1-78O14.1</td></tr>
    <tr><td>cg12658972</td><td>RP11-95P2.3</td></tr>
    <tr><td>cg12790145</td><td>RP11-283I3.2</td></tr>
    <tr><td>cg13303179</td><td>RP11-316F12.1</td></tr>
    <tr><td>cg13689053</td><td>AP000688.14</td></tr>
    <tr><td>cg14545602</td><td>AC005498.3</td></tr>
    <tr><td>cg15501526</td><td>RP11-526P5.2</td></tr>
    <tr><td>cg15587955</td><td>RP11-1080G15.1</td></tr>
    <tr><td>cg16469117</td><td>AC007092.1</td></tr>
    <tr><td>cg17240725</td><td>AC073846.5</td></tr>
    <tr><td>cg20188212</td><td>RP11-283I3.2</td></tr>
    <tr><td>cg23102195</td><td>CTD-2269E23.4</td></tr>
    <tr><td>cg23184739</td><td>AC110781.3</td></tr>
    <tr><td>cg24468780</td><td>RP11-112J1.2</td></tr>
    <tr><td>cg24947255</td><td>XXbac-BPG181B23.6</td></tr>
    <tr><td>cg25329573</td><td>RP11-283I3.2</td></tr>
  </tbody>
</table>
<p style="font-size: 0.9em;">List of CpG sites associated with uncharacterized loci, including long non-coding RNA (lncRNA) and genomic regions named using the Human Genome Project clone system</p>