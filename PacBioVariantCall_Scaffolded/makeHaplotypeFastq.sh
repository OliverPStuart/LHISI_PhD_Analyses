
file=PacBio_MaxCov_Table.txt.candidateHaplotypes.txt

# For each unique contig in the file
for contig in $(cut -f1 ${file} | uniq)
do

  # Make empty header
  samtools view -H ${contig}:*bam > temp1

  # For each unique phase set on that contig in the file
  for phase in $(awk "/${contig}\t/" ${file} | cut -f2)
  do

    #And for both haplotypes of that phase set
    for hap in 1 2
    do

      # View bamfile with this contig in name
      samtools view -@ 2 ${contig}:*bam | \
      # Select reads with matching haplotype tag and phase tag
      awk "/HP:i:${hap}\tPS:i:${phase}/" > temp2
      # Convert to fastq for assembly
      cat temp1 temp2 | samtools fastq - > ${contig}.phase${phase}.hap${hap}.fastq
      gzip ${contig}.phase${phase}.hap${hap}.fastq

      rm temp2

    done
  done

    rm temp1

    # Move all fastqs into contig specific directory
    mkdir ${contig}Fastqs
    mv ${contig}*fastq.gz ${contig}Fastqs

done
