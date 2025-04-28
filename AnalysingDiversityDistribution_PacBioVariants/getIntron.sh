
REF_DIR=/Volumes/Alter/LHISI/References

for GENE in $(cat protsOfInterest.txt)
do

# Get whole gene bed
grep ${GENE} protsOfInterest.gff | awk '$3 == "gene" {print $1"\t"$4"\t"$5"\t"substr($9,4,8)}' > gene_temp.bed

# Get exons and UTRs
grep ${GENE} protsOfInterest.gff | awk '$3 ~ "UTR" || $3 == "exon" {print $1"\t"$4"\t"$5"\t"substr($9,4,8)}' > exonutr_temp.bed

# Subtract
bedtools subtract -a gene_temp.bed -b exonutr_temp.bed >> intron.bed

done

rm *temp.bed
