# AncestryPGx

This repository contains the complete R code and datasets used in the article:  
**_Pharmacogenomics and Genetic Ancestry in Colombia: A Study on All Variant Drug Annotations of PharmGKB_** (under submission).

## Purpose

This study implements a cross-sectional design to investigate the pharmacogenomic landscape of diabetes-related variant annotations in Colombian populations. By integrating variant annotations from PharmGKB (1) with allele frequency data from the Consortium for Genomic Diversity, Ancestry, and Health in Colombia (CÓDIGO) (2), our aim is to assess ancestry-related differences in allele frequencies and to conduct descriptive bibliometric analyses of the associated primary studies.

This repository includes both the full datasets and analysis scripts, enabling full reproducibility of all results.

## Repository Structure

- `main_analysis.R`: Main script with numbered and annotated sections for reproducibility.
- `data/`: Includes:
  - `study_parameters.xlsx`: Study-level metadata, including population sizes and the biogeographical classification of the study cohorts.
  - `var_drug_ann.xlsx`: Full set of variant–drug annotations, which includes associations related to drug dosage, response, metabolism, among others.
  - `SNPs_FQs.csv`: Allele frequencies of each SNP identified.


### Statistical Analysis

- Spearman’s correlation analyses were conducted using allele frequencies from each population. Variants with missing data in any pairwise comparison were excluded.
- For ancestry comparisons, African ancestry was represented by the PLQ population, while the mean frequency across ATQCES, ATQPGC, and CLM defined European ancestry.
- Cases and controls in the original studies were categorized by ancestry (European, Asian, American, African, Other, or Unknown) based on the biogeographical information provided.

### Bibliometric Analysis

- Metadata for the 1,225 studies was retrieved using PubMed IDs through the NCBI Entrez API, Unpaywall API, and Crossref API.
- Journal metrics (H-index and quartile) were obtained via SCImago Journal Rank (SJR) corresponding to each study’s publication year (1999–2024).
- Each study was assigned to a World Health Organization (WHO) region (3) and a World Bank income classification (4) based on the first author’s affiliation.
- Variables included: article title, publication date, number of authors, first author’s country, open/closed access status, citation count, and journal metrics, providing a comprehensive overview for the descriptive bibliometric analysis.

## How to Use

1. Clone or download this repository.
2. Open `main_analysis.R` in RStudio.
3. Ensure the files in the `data/` folder are correctly referenced.
4. Run the script from beginning to end or explore specific numbered sections.

## Reproducibility

> ✅ **This repository is fully reproducible**: All datasets and scripts are included, and no external data download is required.

## Citation

If you use or adapt this code, please cite the associated article (once published) and acknowledge this repository. A DOI will be added upon publication.

## License

This project is licensed under the **MIT License**. You are free to use, modify, and distribute the code with attribution.


## References

1. PharmGKB: https://www.pharmgkb.org  
2. CÓDIGO Database: https://codigo.uniandes.edu.co  
3. World Health Organization (WHO) Regions: https://www.who.int/  
4. World Bank Income Classifications: https://datahelpdesk.worldbank.org/knowledgebase/articles/378833-how-are-the-income-group-thresholds-determined  

