# This script runs the nonWF_demography.slim model with
# replication and combines outputs into a nice containing
# folder. The script only needs only one input: the
# number of replicates desired, and outputs a timestamped
# output folder.

DATE=$(date '+%d%m%Y_%H%M')
REPS=$1

for rep in $( seq 1 $REPS )
do

slim -d rep=$rep nonWF_demography.slim &> ${DATE}.log

done

mkdir $DATE
mv *txt $DATE
