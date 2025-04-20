#!/urs/bin/python3

# Descriptive statistics for SNP data 

filename = '307383_mavi_seq.csv'

# Read the file
with open(filename, 'r') as file:
    lines = file.readlines()

# Process header and data
header = lines[0].strip().split(',')
data = [line.strip().split(',') for line in lines[1:]]

# Identify column indices
chrom_idx = header.index('CHROMOSOME')
genotype_idx = header.index('GENOTYPE')

# Clean the data
cleaned_data = []
valid_bases = {'A', 'C', 'G', 'T'}

for row in data:
    chrom = row[chrom_idx]
    genotype = row[genotype_idx]

    # Skip unwanted chromosomes and invalid genotypes
    if chrom in ['X', 'Y', 'MT']:
        continue
    if genotype == '--':
        continue
    if not all(base in valid_bases for base in genotype):
        continue

    cleaned_data.append(row)

# Descriptive statistics
total_snps = len(cleaned_data)
chromosomes = set()
snps_per_chromosome = {}

for row in cleaned_data:
    chrom = row[chrom_idx]
    chromosomes.add(chrom)
    if chrom not in snps_per_chromosome:
        snps_per_chromosome[chrom] = 0
    snps_per_chromosome[chrom] += 1

# Print results
print(f"Total number of SNPs: {total_snps}")
print(f"Number of unique chromosomes: {len(chromosomes)}")
print("SNPs per chromosome:")

# Convert chromosome keys to integers for sorting
for chrom in sorted(snps_per_chromosome, key=lambda x: int(x)):
    print(f"  Chromosome {chrom}: {snps_per_chromosome[chrom]} SNPs")

import matplotlib.pyplot as plt

# count Homozygous and Heterozygous
homozygous_count = 0
heterozygous_count = 0

for row in cleaned_data:
    genotype = row[genotype_idx]
    if len(genotype) == 2:
        if genotype[0] == genotype[1]:
            homozygous_count += 1
        else:
            heterozygous_count += 1

print(f"\nGenotype Counts:")
print(f"  Homozygous: {homozygous_count}")
print(f"  Heterozygous: {heterozygous_count}")

# Scatter Plot of SNP positions
# Prepare data: x = chromosome (numeric), y = position
chrom_pos = []
chrom_list = []
position_idx = header.index('POSITION')

for row in cleaned_data:
    chrom = row[chrom_idx]
    if chrom.isdigit():  # Only include numeric chromosomes
        chrom_list.append(int(chrom))
        chrom_pos.append(int(row[position_idx]))

# Plotting
plt.figure(figsize=(12, 6))
plt.scatter(chrom_list, chrom_pos, s=1, alpha=0.5)
plt.xlabel('Chromosome')
plt.ylabel('SNP Position')
plt.title('Scatter Plot of SNP Positions per Chromosome')
plt.xticks(range(1, 23))
plt.grid(True)
plt.tight_layout()
plt.savefig('snp_scatter_plot.png', dpi=300)

import csv

# Save cleaned data to a new CSV file
with open('cleaned_snp_data.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)        # Write the header
    writer.writerows(cleaned_data) # Write the cleaned SNP data rows

print("Cleaned data saved to 'cleaned_snp_data.csv'")
