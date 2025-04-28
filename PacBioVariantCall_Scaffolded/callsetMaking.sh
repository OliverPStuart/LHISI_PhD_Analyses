
# Make lists of allele balances
sh tableMaking.sh Scaffold_PacBio.vcf
mawk '{print $1"\t"$2"\t"($7/($7+$8))}' Scaffold_PacBio.vcf.table > siteBals
mawk '$3 < 0.2 || $3 > 0.8 {print $1"\t"$2}' siteBals > offBal_liberal
mawk '$3 < 0.3 || $3 > 0.7 {print $1"\t"$2}' siteBals > offBal_strict

# Make list of sites that do not pass and are below default site quality threshold
# Also check for homozygotes here
mawk '!/#/' Scaffold_PacBio.vcf | grep "\t1|1:\|\t1/1:\|\t0|0:\|\t0/0:" | cut -f1,2 > homs

# Combine into lists for filtering
cat offBal_liberal homs | sort | uniq > badSites_liberal
cat offBal_strict homs | sort | uniq > badSites_strict
rm offBal* homs

# PASS filter, QUALITY > 136, GQ > 20, AB < 0.8, > 0.2, minDP = 100, maxDP = 180
vcftools --vcf Scaffold_PacBio.vcf \
--not-chr Scaffold_4 \
--minGQ 20 \
--minQ 136 \
--remove-filtered-all \
--max-missing 1 \
--exclude-positions badSites_liberal \
--minDP 100 \
--maxDP 180 \
--recode \
--recode-INFO-all \
--out Scaffold_PacBio_Filt1

# PASS filter, QUALITY > 136, GQ > 30, AB < 0.8, > 0.2, minDP = 100, maxDP = 180
vcftools --vcf Scaffold_PacBio.vcf \
--not-chr Scaffold_4 \
--minGQ 30 \
--minQ 136 \
--remove-filtered-all \
--max-missing 1 \
--exclude-positions badSites_liberal \
--minDP 100 \
--maxDP 180 \
--recode \
--recode-INFO-all \
--out Scaffold_PacBio_Filt2

# PASS filter, QUALITY > 136, GQ > 30, AB < 0.7, > 0.2, minDP = 100, maxDP = 180
vcftools --vcf Scaffold_PacBio.vcf \
--not-chr Scaffold_4 \
--minGQ 30 \
--minQ 136 \
--remove-filtered-all \
--max-missing 1 \
--exclude-positions badSites_strict \
--minDP 100 \
--maxDP 180 \
--recode \
--recode-INFO-all \
--out Scaffold_PacBio_Filt3

# Also make a stringently filtered callset for Scaffold_4, it may be useful later
vcftools --vcf Scaffold_PacBio.vcf \
--chr Scaffold_4 \
--minGQ 30 \
--minQ 78 \
--remove-filtered-all \
--max-missing 1 \
--exclude-positions badSites_liberal \
--minDP 40 \
--maxDP 120 \
--recode \
--recode-INFO-all \
--out Scaffold_4_PacBio_Filt

rename 's/\.recode//g' *  ; rm *log bals*
