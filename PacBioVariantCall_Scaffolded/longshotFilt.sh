
VCF=$1
depth=$2

# Select only sites with quality greater than the median read depth
# As per the longshot paper
# Remove genotype with genotype quality < 30 just as a hard filter
mawk '!/#/' ${VCF} | cut -f1,2,6,7 | awk -v depth=${depth} '$3 >= depth && $4 == "PASS"' > goodFilteredSites.txt

# Select only these sites
vcftools --vcf ${VCF} --positions goodFilteredSites.txt --minGQ 30 --max-missing 1 --recode --recode-INFO-all --out ${VCF}.QUAL_FILT

rm goodFilteredSites.txt *log
rename 's/\.recode//g' *
