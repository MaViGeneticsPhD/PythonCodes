#!/usr/bin/python3

import pandas as pd
from difflib import get_close_matches

# Load GWAS catalog
gwas_df = pd.read_csv("gwas_catalog.tsv", sep='\t', low_memory=False)

# Split 'STRONGEST SNP-RISK ALLELE' into RSID and RISK_ALLELE
gwas_df[['RSID', 'RISK_ALLELE']] = gwas_df['STRONGEST SNP-RISK ALLELE'].str.extract(r'(rs\d+)-([ACGT])')

# Choose effect size column
effect_column = 'OR or BETA' if 'OR or BETA' in gwas_df.columns else 'P-VALUE (TEXT)'

# Lowercase disease trait for easier matching
gwas_df['DISEASE/TRAIT_clean'] = gwas_df['DISEASE/TRAIT'].str.lower().str.strip()

# Disease list provided
disease_list = """Type 2 diabetes
Parkinson's disease
Migraine
Melanomas of skin
Colon cancer
Lymphoid leukemia
Asthma
Hypothyroidism
Thyroid cancer
Cancer of bronchus 
Malignant neoplasm of female breast
Herpes simplex
Insomnia
Tuberculosis
Herpes zoster
Sarcoidosis
Endometriosis
Infectious mononucleosis
Viral hepatitis B
Anorexia nervosa
Cholangitis
Pneumonia
Barrett's esophagus
Hodgkin's disease
Systemic lupus erythematosus
Cancer of nasopharynx
Ulcerative colitis
Inflammatory bowel disease and other gastroenteritis and colitis
Vitiligo
Rheumatoid arthritis
Hyperglyceridemia
Cancer of prostate
Psoriasis
Squamous cell carcinoma
Myocardial infarction
Pancreatic cancer
Cancer of stomach
Colorectal cancer
Cataplexy and narcolepsy
Malignant neoplasm of uterus
Multiple sclerosis
Abdominal aortic aneurysm
Polycystic ovaries
Coronary atherosclerosis
Schizophrenia
Basal cell carcinoma
Alzheimer's disease
Posttraumatic stress disorder
Psychosis
Uterine leiomyoma
Excessive or frequent menstruation
Primary pulmonary hypertension
Congestive heart failure
Intestinal infection
Calculus of kidney
Cervical cancer
Allergic rhinitis
Polymyositis
Keloid scar
Type 1 diabetes
Obesity
Leprosy
Essential tremor
Osteoporosis
Celiac disease
Malignant neoplasm of ovary
Vascular dementia
Giant cell arteritis
Malignant neoplasm of bladder
Dermatomyositis
Multiple myeloma
ADHD
Cerebral aneurysm
Graves' disease
Wegener's granulomatosis
Psoriatic arthropathy
Chronic sinusitis
Nasal polyps
Major depressive disorder
Attention deficit hyperactivity disorder
Cardiac dysrhythmias
Cataract
Proteinuria
Hematuria
Cancer of esophagus
Malignant neoplasm of liver, primary
Urinary calculus
Cirrhosis of liver without mention of alcohol
Glaucoma
Viral hepatitis C
Diverticulitis
Chronic lymphocytic thyroiditis
Rheumatic disease of the heart valves
Hypertension
Alcoholism
Osteitis deformans [Paget's disease of bone]
Osteoarthrosis
Myeloproliferative disease
Decreased white blood cell count
Inguinal hernia
Anxiety disorders
Weight
Hair type
Hair color
Eye color
Skin color
Blood type
Sleep duration
Acne
hirsutism
allergies
lactose intolerance
vitamin d levels""".strip().split('\n')

# Lowercase for matching
disease_list_clean = [d.lower().strip() for d in disease_list]

# Optional: trait mapping (only add mappings that matter in your context)
trait_mapping = {
    "melanomas of skin": "melanomas of skin (phecode 172.11)",
    "lymphoid leukemia": "lymphocytic leukemia",
    "cancer of bronchus": "cancer of bronchus; lung (phecode 165.1)",
    "malignant neoplasm of female breast": "malignant neoplasm of female breast (phecode 174.11)",
    "herpes zoster": "herpes zoster (phecode 53)",
    "infectious mononucleosis": "mononucleosis",
    "viral hepatitis b": "hepatitis b",
    "cholangitis": "cholecystitis",
    "cancer of nasopharynx": "cancer of oropharynx (phecode 149.1)",
    "hyperglyceridemia": "hypertriglyceridemia",
    "polycystic ovaries": "polycystic ovary syndrome",
    "posttraumatic stress disorder": "post-traumatic stress disorder",
    "cirrhosis of liver without mention of alcohol": "cirrhosis of liver (phecode 571.51)",
    "psoriatic arthropathy": "psoriatic arthritis",
    "hair type": "hair shape",
    "allergies": "allergic disease",
    "lactose intolerance": "lactose tolerance test",
}

# Apply mapping
gwas_df['DISEASE/TRAIT_clean'] = gwas_df['DISEASE/TRAIT_clean'].replace(trait_mapping)

# Filter GWAS
filtered_gwas = gwas_df[gwas_df['DISEASE/TRAIT_clean'].isin(disease_list_clean)]

# Save filtered GWAS
filtered_gwas[['RSID', 'RISK_ALLELE', effect_column, 'DISEASE/TRAIT']].to_csv("gwas_catalog_filtered.csv", index=False)

# Load cleaned SNP data (must include columns RSID and GENOTYPE)
snp_df = pd.read_csv("cleaned_snp_data.csv")

# Load GWAS filtered file
gwas_filtered = pd.read_csv("gwas_catalog_filtered.csv")

# Merge by RSID
merged = pd.merge(snp_df, gwas_filtered, on='RSID')

# Define dosage calculator
def get_dosage(genotype, risk_allele):
    return genotype.upper().count(risk_allele)

# Calculate dosage and PRS component
merged['Dosage'] = merged.apply(lambda row: get_dosage(row['GENOTYPE'], row['RISK_ALLELE']), axis=1)
merged['PRS_Component'] = merged['Dosage'] * merged[effect_column]

# Compute PRS
total_prs_score = merged['PRS_Component'].sum()

# Save detailed breakdown
breakdown_cols = ['RSID', 'DISEASE/TRAIT', 'GENOTYPE', 'RISK_ALLELE', 'Dosage', effect_column, 'PRS_Component']
merged[breakdown_cols].to_csv("prs_trait_breakdown.csv", index=False)


# Generate report text
with open("polygenic_risk_score_report.txt", "w") as report:
    report.write("Polygenic Risk Score Report\n")
    report.write("===========================\n\n")
    report.write(f"Total Polygenic Risk Score (PRS): {total_prs_score:.4f}\n\n")
    report.write("Trait-Level Contributions:\n")
    report.write("----------------------------------------\n")
    for _, row in merged.iterrows():
        report.write(f"RSID: {row['RSID']}\n")
        report.write(f"Trait: {row['DISEASE/TRAIT']}\n")
        report.write(f"Genotype: {row['GENOTYPE']} | Risk Allele: {row['RISK_ALLELE']}\n")
        report.write(f"Dosage: {row['Dosage']} | Effect Size: {row[effect_column]} | Contribution to PRS: {row['PRS_Component']:.4f}\n")
        report.write("----------------------------------------\n")

print("PRS computation complete. Outputs:")
print("- Polygenic Risk Score:", round(total_prs_score, 4))
print("- Report file: polygenic_risk_score_report.txt")
print("- Breakdown (table) file: prs_trait_breakdown.csv")

# Summarize PRS contribution by trait
trait_prs_summary = merged.groupby('DISEASE/TRAIT')['PRS_Component'].sum().sort_values(ascending=False)

# what sre the top risk-contributing traits
top_n = 20
print(f"\nTop {top_n} Risk-Contributing Traits:")
print("----------------------------------------")
for trait, score in trait_prs_summary.head(top_n).items():
    print(f"{trait}: {score:.4f}")

# Append to report
with open("polygenic_risk_score_report.txt", "a") as report:
    report.write("\nTop Risk-Contributing Traits:\n")
    report.write("========================================\n")
    for trait, score in trait_prs_summary.head(top_n).items():
        report.write(f"{trait}: {score:.4f}\n")

