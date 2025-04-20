# Final Project - HORT530 
# Maria Victoria Pereira de Souza

# Single Human SNP analyzer

## Description
This Project is a Python-based code designed to analyze genetic variation in an individual's SNP (Single Nucleotide Polymorphism) data. First, it shows the basic descriptive statistics after the data has been cleaned. Secondly, it matches the individual's genotype against known risk alleles (from publications) for various phenotypes to calculate risk scores. 
Note: these risk scores are biased since they came from an weighted average and not from a statistical test. 
## The code performs:
- Data quality control
- Summary statistics on genetic variation
- Genotype frequency analysis
- Graphs (.png) for visualization
- SNP matching against known phenotype associations (some are self-reported by the individual, others are expectations) 
- Risk score calculations for various traits (phenotypes)

This project enables individuals to understand their genetic predispositions to various traits and conditions based on their raw SNP data (such as that provided by consumer genetic testing services).

## Features
- Data Cleaning: Pre-processes raw SNP data files for analysis
- Phenotype Association: Maps SNPs to associated phenotypes based on scientific literature
- Risk Assessment: Calculates weighted risk scores based on matching risk alleles
- Reporting: Generates detailed reports of findings with risk levels
- Supports CSV-based input
- Performance Optimization: Efficiently processes large SNP datasets (60k+ SNPs)

#######################################################################################################################################################

FIRST PART - DESCRIPTIVE STATISTICS

```bash
python 1stats.py
```
This will analyze the raw SNP data file, generate a cleaned SNP data file in csv format and some graphics as well.

### Sample output

Total number of SNPs: 60902
Number of unique chromosomes: 22
SNPs per chromosome:
  Chromosome 1: 5085 SNPs
  Chromosome 2: 4972 SNPs
  Chromosome 3: 4293 SNPs
  Chromosome 4: 3812 SNPs

##########################################################################################################################################################

SECOND PART - SNPS X PHENOTYPES ASSOCIATIONS

```bash
python SNPs_analyzer.py
```
This will analyze the default SNP data file (`cleaned_snp_data.csv`) against the phenotype associations defined in `combined_phenotypes_rsid_effect.csv`.

Sample Output

Detailed Phenotype-SNP Analysis:
Phenotype SNP Genotype Risk Allele Effect Size Match
--------------------------------------------------------------------------------
ADHD:
  rs6537401 Not found GA 0.945 N/A
  rs4916723 AA AC 0.918 Yes
  rs77960 Not found GA 0.929 N/A

...

============================================================
PHENOTYPE RISK SUMMARY:
------------------------------------------------------------
Phenotype Risk Score
------------------------------------------------------------
brown_eye 100.0% (High)
HbF 22.4% (Low)
colorectal_cancer 17.0% (Low)
ADHD 14.5% (Low)

Run time duration: 0.23 seconds


To analyze your own data:

1. Prepare your SNP data in CSV format with at least columns for RSID and GENOTYPE
2. Create a phenotype associations file with columns for phenotype, rsID, Genotype (risk allele), and effect_size
3. Run the analysis


### Data Format

SNP Data CSV Format

RSID,CHROMOSOME,POSITION,GENOTYPE
rs4916723,1,1857048,CT
rs1835740,7,2749849,AG
rs12913832,21,987689,GG
...


Phenotype Associations CSV Format
```
phenotype,rsID,Genotype,effect_size
ADHD,rs4916723,C,0.67
ADHD,rs11420276,G,0.22
Migraine,rs1835740,A,0.51
...
```

### Future development plans (based on current limitations) include:

- Expanded Phenotype Database: Include more phenotypes and associations from recent studies
- API Integration: Connect with scientific databases for up-to-date phenotype associations
- Population Comparison: Compare individual results with population averages


### References

1. Aissani, B., Zhang, K. & Wiener, H. Evaluation of GWAS candidate susceptibility loci for uterine leiomyoma in the multi-ethnic NIEHS uterine fibroid study. Front Genet 6, 241 (2015).
2. Day, F. et al. Large-scale genome-wide meta-analysis of polycystic ovary syndrome suggests shared genetic architecture for different diagnosis criteria. PLoS Genet 14, e1007813 (2018).
3. Demontis, D. et al. Genome-wide analyses of attention deficit hyperactivity disorder identify 27 risk loci, refine the genetic architecture, and implicate several cognitive domains. Nat Genet 55, 198–208 (2023).
4. Hautakangas, H. et al. Genome-wide analysis of 102,084 migraine cases identifies 123 risk loci and subtype-specific risk alleles. Nat Genet 54, 152–160 (2022).
5. Liu, F. et al. Meta-analysis of genome-wide association studies identifies 8 novel loci involved in shape variation of human head hair. Human Molecular Genetics 27, 559–575 (2018).
6. Zhang, H. et al. Genome-wide association study identifies 32 novel breast cancer susceptibility loci from overall and subtype-specific analyses. Nat Genet 52, 572–581 (2020).
7. Genome-wide analysis in over 1 million individuals of European ancestry yields improved polygenic risk scores for blood pressure traits | Nature Genetics. https://www.nature.com/articles/s41588-024-01714-w.
8. Genome-wide association analyses identify 143 risk variants and putative regulatory mechanisms for type 2 diabetes | Nature Communications. https://www.nature.com/articles/s41467-018-04951-w.
9. Genome-wide association study of major anxiety disorders in 122,341 European-ancestry cases identifies 58 loci and highlights GABAergic signaling. PubMed Central (PMC) https://pmc.ncbi.nlm.nih.gov/articles/PMC11245051/.
10. Table 1 Summary results for the new colorectal cancer risk loci in Europeans.
11. Table 2 Association results of detected risk variants.
12. Liu, F. et al. Meta-analysis of genome-wide association studies identifies 8 novel loci involved in shape variation of human head hair. Hum Mol Genet (2018) doi:10.1093/hmg/ddx416.

