from Bio import SeqIO
import sys

chrom_ranges = open("chrom_ranges.txt", 'w')
for record in SeqIO.parse(sys.argv[1], "fasta"):
    chrom = record.id
    total_len = len(record.seq)
    min_number = 0
    step = 100000
    if step < total_len:
        for chunk in range(min_number, total_len, step)[1:]:
            chrom_ranges.write("{}:{}-{}\n".format(chrom, min_number, chunk))
            min_number = chunk
    chrom_ranges.write("{}:{}-{}\n".format(chrom, min_number, total_len))
chrom_ranges.close()

    