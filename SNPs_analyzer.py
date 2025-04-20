#!/usr/bin/python

import csv                 # for reading and writing csv files
import time                # for tracking execution time
from collections import defaultdict        #dictionary that provides default values for missing keys

# Track run time
start_time = time.time()

# Read SNP data cleaned
with open('cleaned_snp_data.csv', 'r') as f:
    reader = csv.reader(f)                       # creates a reader object
    header = next(reader)                        # reads the first row as the header
    rsid_idx = header.index('RSID')              # find the indext position of the rsID column
    cleaned_data = list(reader)                  # reads all remaining rows into a list

# Create a dictionary of SNP data
genotype_idx = header.index('GENOTYPE')
snp_dict = {}                                    # initialize a dictionary to store SNP data
for row in cleaned_data:                         # loop to iterate through each row in the cleaned_data
    rsid = row[rsid_idx]
    genotype = row[genotype_idx]
    snp_dict[rsid] = genotype                    # adds the rsID (key) and genotype (value) to the dictionary

# Read phenotype data from CSV file
phenotype_snps = defaultdict(lambda: {"snps": []})         # defaultdict that will return {"snps": []} for any missing key
with open('combined_phenotypes_rsid_effect.csv', 'r') as f:
    reader = csv.reader(f)
    header = next(reader)

    # Determine column indices from header
    phenotype_idx = header.index('Phenotype')
    rsid_idx = header.index('rsID')
    risk_allele_idx = header.index('Genotype')
    effect_size_idx = header.index('effect_size')

    # Read each row and organize by phenotype
    for row in reader:
        phenotype = row[phenotype_idx]
        rsid = row[rsid_idx]
        risk_allele = row[risk_allele_idx]

        # Convert effect size to float
        try:
            effect_size = float(row[effect_size_idx])
        except (ValueError, IndexError):
            effect_size = 0.0  # Set to zero if not available or not numeric

        # Add this SNP to the phenotype's list
    # each SNP is stored as a dictionary with rsID, risk_allele, and effect_size keys
        phenotype_snps[phenotype]["snps"].append({
            "rsID": rsid,
            "risk_allele": risk_allele,
            "effect_size": effect_size
        })

# Print the phenotype-genotype table with detailed SNP information
print("\nDetailed Phenotype-SNP Analysis:")
print(f"Phenotype SNP Genotype Risk_Allele Effect_Size Match")
print("-" * 80)

# Track total risk score by phenotype
phenotype_risk_scores = {}              # open an empty dictionary to store risk scores for each phenotype

for phenotype, data in phenotype_snps.items():
    total_found_effect = 0
    total_possible_effect = 0                         # initializes counters for tracking effects and matches
    matches_found = 0

    print(f"\n{phenotype}:")

# Loop through each SNP associated with the current phenotype
# extract the information and add the effect size to the total possible effects
    for snp_info in data["snps"]:
        rsid = snp_info["rsID"]
        risk_allele = snp_info["risk_allele"]
        effect_size = snp_info["effect_size"]

        total_possible_effect += effect_size

        # Get the genotype from my data
        if rsid in snp_dict:
            genotype = snp_dict[rsid]

            # Check for any matching alleles between risk_allele and genotype
            has_risk_allele = False

            # condition for single allele
            if len(risk_allele) == 1:
                has_risk_allele = risk_allele in genotype
            # condition for 2 alleles
            else:
                for allele in risk_allele:
                    if allele in genotype:
                        has_risk_allele = True
                        break

            match_status = "Yes" if has_risk_allele else "No"

           # If there's a match, add the effect to the risk score
            if has_risk_allele:
                matches_found += 1
                total_found_effect += effect_size

            print(f"  {rsid} {genotype} {risk_allele} {effect_size} {match_status}")
        else:
            print(f"  {rsid} Not found {risk_allele} {effect_size} N/A")

    # Calculate overall risk percentage for this phenotype
    if total_possible_effect > 0:
        risk_percentage = (total_found_effect / total_possible_effect) * 100
        phenotype_risk_scores[phenotype] = risk_percentage

        print(f"  Summary: {matches_found}/{len(data['snps'])} SNPs matched, Risk Score: {risk_percentage:.1f}%")
    else:
        print(f"  Summary: No valid SNPs found")

# Display overall summary
print("\n" + "=" * 60)
print("PHENOTYPE RISK SUMMARY:")
print("-" * 60)
print(f"Phenotype Risk Score")
print("-" * 60)

# Sorts phenotypes by risk score in descending order
for phenotype, score in sorted(phenotype_risk_scores.items(), key=lambda x: x[1], reverse=True):
    risk_level = "High" if score > 60 else "Medium" if score > 30 else "Low"
    print(f"{phenotype} {score:.1f}% ({risk_level})")

# Track run time
end_time = time.time()
print('\nRun time duration: {:.2f} seconds'.format(end_time - start_time))

