# This script runs the half_sib_het model with
# replication and combines outputs into a nice table.
# The script only needs only one input: the number
# of replicates desired, and outputs a timestamped
# output table.

DATE=$(date '+%d%m%Y_%H%M')
REPS=$1

echo "rep,ind1,ind2" > ${DATE}_half_sib_het_results.txt

for rep in $( seq 1 $REPS )
do

slim -d rep=$rep half_sib_het.slim | \
tail -1 >> ${DATE}_half_sib_het_results.txt

done
