
REF_DIR=/Volumes/Alter/LHISI/References

lines=$(cat ExonTargetCoords.txt | wc -l | sed -e 's/^[ \t]*//')

for LINE in $(seq 1 ${lines})
do

coord=$(awk -v line=${LINE} 'NR==line {print $1}' ExonTargetCoords.txt)
name=$(awk -v line=${LINE} 'NR==line {print $2}' ExonTargetCoords.txt)

echo ">"${name} >> tmp
samtools faidx ${REF_DIR}/LHISI_Scaffold_Assembly.fasta ${coord} >> tmp

done

grep -v "Scaffold" tmp > ExonTargets.fa
rm tmp
