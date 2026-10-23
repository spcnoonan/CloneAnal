import sys

def split_vcf_by_genotype(vcf_file, wey_output, wis_output):
    # Open the input and output files
    with open(vcf_file, 'r') as vcf, open(wey_output, 'w') as wey_out, open(wis_output, 'w') as wis_out:
        for line in vcf:
            # Write the header lines (those starting with '#') to both files
            if line.startswith('#'):
                wey_out.write(line)
                wis_out.write(line)
            else:
                # Split the line by tabs to get the fields
                fields = line.strip().split('\t')

                # Extract the genotype information for columns 27 and 28
                sample_9_data = fields[9]
                sample_10_data = fields[10]

                # Extract the GT field from both samples
                gt_9 = sample_9_data.split(':')[0]  # GT field of column 10
                gt_10 = sample_10_data.split(':')[0]  # GT field of column 11

                # If GT field of column 27 contains a '1', write to wey_vars.vcf
                if '1' in gt_9:
                    wis_out.write(line)

                # If GT field of column 28 contains a '1', write to wis_vars.vcf
                if '1' in gt_10:
                    wey_out.write(line)

# Example usage
vcf_file = sys.argv[1]
wey_output = 'WeyVars.vcf'
wis_output = 'WisVars.vcf'

split_vcf_by_genotype(vcf_file, wey_output, wis_output)