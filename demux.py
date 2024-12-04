import subprocess

# demuxlet can be installed with:
# conda install bioconda::demuxlet

def run_demuxlet(input_bam, reference_vcf, output_prefix):
    demuxlet_cmd = [
        'demuxlet',
        '--sam', input_bam,
        '--vcf', reference_vcf,
        '--field', 'GT',
        '--out', output_prefix,
    ]
    subprocess.run(demuxlet_cmd, check=True)
    print("Demuxlet completed.")

path_vcf = '/home/ehb/Desktop/scRNA_analysis/CellRanger_and_VCF_Data/vcf_data/valid_8patient.vcf'

# path_bam = '/home/ehb/Desktop/scRNA_analysis/CellRanger_and_VCF_Data/GEX/gex.bam'
# path_output_with_prefix = '/home/ehb/Desktop/scRNA_analysis/CellRanger_and_VCF_Data/GEX/gex'

path_bam = '/home/ehb/Desktop/scRNA_analysis/CellRanger_and_VCF_Data/FB/fb.bam'
path_output_with_prefix = '/home/ehb/Desktop/scRNA_analysis/CellRanger_and_VCF_Data/FB/fb'


run_demuxlet(path_bam, path_vcf, path_output_with_prefix)
