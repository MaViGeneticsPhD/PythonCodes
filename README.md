# Final Project - HORT530
## Maria Victoria Pereira de Souza

# Single Human SNP Analyzer with Polygenic Risk Score Report

---

## Description  
This Python-based project analyzes an individual's SNP (Single Nucleotide Polymorphism) data to identify genetic predispositions 
for various traits and conditions based on published risk allele data from GWAS studies.  

It includes:
- Descriptive statistics on genomic variation  
- Matching of genotypes against known phenotype-associated SNPs  
- Polygenic Risk Score (PRS) calculation  
- Ranked risk contribution summary  
- Automated output report in `.txt` format

> Note: The current PRS implementation uses additive scoring with GWAS-reported effect sizes and 
does **not** perform population-level corrections or statistical significance testing.

---

## What can you expect from this code?
- Data Quality Control: Cleans SNP data (missing values, invalid entries, sex-specific chromosome filtering).  
- Descriptive Stats: SNP counts, distribution per chromosome, and data completeness.  
- Risk Allele Matching: Maps individual genotypes to known phenotype associations.  
- Polygenic Risk Score (PRS) Calculation: Includes dosage-based contribution per trait.  
- Top-Risk Traits Identification: Lists most genetically influenced traits by cumulative PRS score.  
- Graphical Output: in `.png` format.  

---

## Project Structure  

This project contains:
- `README.txt` (this file)  
- `1stats.py` — runs descriptive statistics  
- `PRS.py` — analyzes SNPs and generates PRS report  
- `307383_mavi_seq.csv` — raw SNP data  
- `gwas_catalog.tsv` — Gwas Catalog association data downloaded from https://www.ebi.ac.uk/gwas/

---

## Getting Started  

1. Open a terminal.  
2. Log in to your Purdue Scholar user account.  
3. Navigate to:  
   ```bash
   cd /scratch/scholar/mavi/FinalProjectMavi
   ```

### Environment Setup

Make sure you have Anaconda installed and load it:
```bash
module load anaconda
```

Then, set up the environment:
```bash
conda create --name py3117 python=3.11.7
conda activate py3117
conda install numpy scipy pandas matplotlib 
```

---

## Part 1 - Descriptive Statistics

Run:
```bash
python 1stats.py
```

This script:
- Loads and cleans raw SNP data  
- Outputs cleaned data (`cleaned_snp_data.csv`)  
- Saves a SNP scatter plot as `.png`

Example output:
```
Total SNPs: 60,902  
Chromosomes analyzed: 22  
SNPs per chromosome:  
 - Chromosome 1: 5085  
 - Chromosome 2: 4972  
...
```

---

## Part 2 - SNP-Phenotype Associations & PRS Calculation

Run:
```bash
python PRS.py
```

This script:
- Reads cleaned SNP data  
- Compares individual genotypes to phenotype-associated SNPs  
- Calculates dosage-based PRS  
- Writes detailed report to `polygenic_risk_score_report.txt`  
- Prints Top Risk-Contributing Traits

Example terminal output:
```
Top 20 Risk-Contributing Traits:
----------------------------------------
Type 2 diabetes: 3.2871  
Breast cancer: 2.9504   
...
```

---

## Data Formats

### SNP Input File (`.csv`)
```
RSID,CHROMOSOME,POSITION,GENOTYPE  
rs4916723,1,1857048,CT  
rs1835740,7,2749849,AG  
...
```

### Phenotype-SNP Association File
```
DISEASE/TRAIT,RSID,RISK_ALLELE,EFFECT_SIZE  
ADHD,rs4916723,C,0.67  
Migraine,rs1835740,A,0.51  
...
```

---

## Future Development
- Normalize PRS using population-level allele frequencies  
- Interactive dashboard for trait visualization  
  

---

## References  
1. Aissani, B., Zhang, K. & Wiener, H. Evaluation of GWAS candidate susceptibility loci for uterine leiomyoma in the 
multi-ethnic NIEHS uterine fibroid study. Front Genet 6, 241 (2015).  
2. GWAS Catalog: https://www.ebi.ac.uk/gwas/  
3. Visscher et al., *10 Years of GWAS Discovery*, Nature Reviews Genetics, 2017.

---



