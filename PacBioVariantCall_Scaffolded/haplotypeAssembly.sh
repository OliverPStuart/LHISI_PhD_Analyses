
# For all lines, one per phase set
for LINE in $(awk "{print NR}" PacBio_MaxCov_Table.txt.candidateHaplotypes.txt)
do
  # Get phase set info to input to wtdbg2
  CTG=$(awk "NR==${LINE}" PacBio_MaxCov_Table.txt.candidateHaplotypes.txt | cut -f1)
  PHASE=$(awk "NR==${LINE}" PacBio_MaxCov_Table.txt.candidateHaplotypes.txt | cut -f2)
  LEN=$(awk "NR==${LINE}" PacBio_MaxCov_Table.txt.candidateHaplotypes.txt | cut -f4)

  # Make output dir
  mkdir ${CTG}phase${PHASE}assemblies

  # Now for both haplotypes
  for HAP in 1 2
  do

    # Assemble with wtdbg2 pipeline, setting estimated genome size to estimated phase block length
    # This is the difference between ID of the phase block (the first base) and the position of the last variant detected in it
    /home/oliver/wtdbg2/wtdbg2.pl \
    -t 8 -x sq -g ${LEN} \
    -o ${CTG}phase${PHASE}hap${HAP} \
    ${CTG}Fastqs/${CTG}.phase${PHASE}.hap${HAP}.fastq.gz

    # Get that shit out of here
    mv ${CTG}phase${PHASE}hap${HAP}* ${CTG}phase${PHASE}assemblies

  done
done



# Now a separate loop to get the info from them

echo "CTG PHASE HAP ASSEMBLY LEN PREDLEN" > AssemblyStats.txt

for LINE in $(awk "{print NR}" PacBio_MaxCov_Table.txt.candidateHaplotypes.txt)
do

  # Get phase set info
  CTG=$(awk "NR==${LINE}" PacBio_MaxCov_Table.txt.candidateHaplotypes.txt | cut -f1)
  PHASE=$(awk "NR==${LINE}" PacBio_MaxCov_Table.txt.candidateHaplotypes.txt | cut -f2)
  PREDLEN=$(awk "NR==${LINE}" PacBio_MaxCov_Table.txt.candidateHaplotypes.txt | cut -f4)

  # For both haplotypes
  for HAP in 1 2
  do

  N=$(awk "/>/" ${CTG}phase${PHASE}assemblies/${CTG}phase${PHASE}hap${HAP}.cns.fa | wc -l)
  grep ">" ${CTG}phase${PHASE}assemblies/${CTG}phase${PHASE}hap${HAP}.cns.fa > ctgs

  # For all contigs in the cns.fa for that haplotype
    for ASSEMBLY in $(seq ${N})
    do

      LEN=$(awk "NR == ${ASSEMBLY}" ctgs | cut -d "=" -f2)

      # Get contig length and name
      echo "${CTG} ${PHASE} ${HAP} ${ASSEMBLY} ${LEN} ${PREDLEN}"  >> AssemblyStats.txt

    done
  done
done
