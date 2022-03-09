#stat_summary.py
#author: Noah A. Legall
#date created: June 13th 2019
#date finished:
#dates edited:
#purpose: creation of alignment stats after alignment pipeline is complete.

import sys # use to access arguments
import os # use in order to call commands from terminal script is called in
import re # regular expressions.
import subprocess as sp

#this program will create the stats files that data will be extracted from.
#if there are 15 samples studied then there will be 15 rows in the outputted csv.

#inputs
r1 = sys.argv[1]
r2 = sys.argv[2]
#trim_r1 = sys.argv[3]
#trim_r2 = sys.argv[4]
bam = sys.argv[3]

# Determine statistics revolving around SNP quality.
#vcf = sys.argv[6]

# FASTQ statistics
print("generating fastq statistics...")
sample_name = re.sub(".trimmed_R1.fastq","",r1)
r1_output = re.sub(".trimmed_R1.fastq",".R1_stats.txt",r1)
r2_output = re.sub(".trimmed_R2.fastq",".R2_stats.txt",r2)
generate_fastq_stats = "fastx_quality_stats -i {} -o {}"

os.system(generate_fastq_stats.format(r1,r1_output))
os.system(generate_fastq_stats.format(r2,r2_output))
R1_size = sp.check_output("ls -Ll {} | cut -f 8 -d ' '".format(r1_output),shell=True).decode('ascii')
R2_size = sp.check_output("ls -Ll {} | cut -f 8 -d ' '".format(r2_output),shell=True).decode('ascii')
Q_ave_R1 = sp.check_output("cat {} | sed \"1d\" | awk '{{sum+=$6}} END{{print sum/NR}}'".format(r1_output),shell=True).decode('ascii')
Q_ave_R2 = sp.check_output("cat {} | sed \"1d\" | awk '{{sum+=$6}} END{{print sum/NR}}'".format(r2_output),shell=True).decode('ascii')
R1_ave_read_length = sp.check_output("awk '{{if(NR%4==2) {{count++; bases += length}} }} END{{print bases/count}}' {}".format(r1),shell=True).decode('ascii')
R2_ave_read_length = sp.check_output("awk '{{if(NR%4==2) {{count++; bases += length}} }} END{{print bases/count}}' {}".format(r2),shell=True).decode('ascii')
fastq_stats = "{},{},{},{},{},{},{}".format(sample_name,R1_size.strip(),R2_size.strip(),Q_ave_R1.strip(),Q_ave_R2.strip(),R1_ave_read_length.strip(),R2_ave_read_length.strip())

# Mapping statistics
print("generating mapping statistics...")
os.system("samtools index {}".format(bam))
ave_coverage = sp.check_output("samtools depth {} | awk '{{sum+=$3}} END {{ print sum/NR}}'".format(bam),shell=True).decode('ascii').strip()
total_reads = sp.check_output("samtools flagstat {} | awk -F '[ ]' 'NR==1 {{print $1}}'".format(bam),shell=True).decode('ascii').strip()
mapped_reads = sp.check_output("samtools flagstat {} | awk -F '[ ]' 'NR==7 {{print $1}}'".format(bam),shell=True).decode('ascii').strip()
unmapped_reads = int(total_reads) - int(mapped_reads)
mapping_stats = "{},{},{}".format(mapped_reads,ave_coverage,unmapped_reads)

# VCF statistics 
# rtg-tools


print("{},{}".format(fastq_stats,mapping_stats))
