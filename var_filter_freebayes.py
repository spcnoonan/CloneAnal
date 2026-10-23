import sys

def filter_wey_wis_variants(input_vcf, output_vcf, wey_wis_samples):
    # Open input VCF file for reading
    with open(input_vcf, 'r') as infile:
        # Open output VCF file for writing
        with open(output_vcf, 'w') as outfile:
            for line in infile:
                if line.startswith('#'):
                    # Write header lines to the output VCF file
                    outfile.write(line)
                else:
                    # Split the line into fields
                    fields = line.strip().split('\t')

                    # Extract genotypes for Wey and Wis samples
                    wey_genotype = fields[wey_wis_samples[0]].split(':')[0]
                    wis_genotype = fields[wey_wis_samples[1]].split(':')[0]

                    # Check if either the reference allele ('0') or the alternative allele ('1') is present in both Wey and Wis
                    if ('0' in wey_genotype and '0' in wis_genotype) or ('1' in wey_genotype and '1' in wis_genotype) or ('2' in wey_genotype or '2' in wis_genotype) or ('3' in wey_genotype or '3' in wis_genotype) or ('4' in wey_genotype or '4' in wis_genotype) or ('5' in wey_genotype or '5' in wis_genotype) or ('.' in wey_genotype or '.' in wis_genotype):
                        continue  # Skip this variant if either allele is present in both Wey and Wis

                    # Write the variant to the output VCF file
                    outfile.write(line)

# Replace 'input.vcf' and 'output.vcf' with your input and output file names
# Replace [9, 10] with the column indices of 'Wey' and 'Wis' samples (0-based index)
filter_wey_wis_variants(sys.argv[1], sys.argv[2], [9, 10])