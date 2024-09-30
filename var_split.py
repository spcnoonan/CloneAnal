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
                sample_27_data = fields[27]
                sample_28_data = fields[28]

                # Extract the GT field from both samples
                gt_27 = sample_27_data.split(':')[0]  # GT field of column 27
                gt_28 = sample_28_data.split(':')[0]  # GT field of column 28

                # If GT field of column 27 contains a '1', write to wey_vars.vcf
                if '1' in gt_27:
                    wey_out.write(line)

                # If GT field of column 28 contains a '1', write to wis_vars.vcf
                if '1' in gt_28:
                    wis_out.write(line)

# Example usage
vcf_file = sys.argv[1]
wey_output = 'WeyVars.vcf'
wis_output = 'WisVars.vcf'

split_vcf_by_genotype(vcf_file, wey_output, wis_output)