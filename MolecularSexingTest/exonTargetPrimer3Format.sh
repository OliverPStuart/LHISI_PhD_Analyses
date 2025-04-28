
rm Primer3InputFile

REF_DIR=/Volumes/Alter/LHISI/References

lines=$(cat ExonTargetCoords.txt | wc -l | sed -e 's/^[ \t]*//')

for LINE in $(seq 1 ${lines})
do

coord=$(awk -v line=${LINE} 'NR==line {print $1}' ExonTargetCoords.txt)
name=$(awk -v line=${LINE} 'NR==line {print $2}' ExonTargetCoords.txt)
samtools faidx ${REF_DIR}/LHISI_Scaffold_Assembly.fasta ${coord} > tmp.fa
seq=$(bioawk -c fastx '{print $seq}' tmp.fa)

echo "SEQUENCE_ID=${name}" >> Primer3InputFile
echo "SEQUENCE_TEMPLATE=${seq}" >> Primer3InputFile
echo "PRIMER_TASK=generic" >> Primer3InputFile
echo "PRIMER_OPT_SIZE=30" >> Primer3InputFile
echo "PRIMER_MIN_SIZE=28" >> Primer3InputFile
echo "PRIMER_MAX_SIZE=33" >> Primer3InputFile
echo "PRIMER_MAX_NS_ACCEPTED=0" >> Primer3InputFile
echo "PRIMER_MIN_TM=66" >> Primer3InputFile
echo "PRIMER_MAX_TM=74" >> Primer3InputFile
echo "PRIMER_OPT_TM=70" >> Primer3InputFile
echo "=" >> Primer3InputFile

done
