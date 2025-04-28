## Windowed heterozygosity, custom script
## Required bedtools and bioawk
## Calls an Rscript with tidyverse

SCRIPT_DIR=/Volumes/Alter/LHISI/Scripts

VCF=${1}
REF=${2}
WIN=${3}

## Make the windows file

bioawk -c fastx 'length($seq) > 1000000 {print $name"\t"length($seq)}' ${REF} > lengths
bedtools makewindows -g lengths -w ${WIN} > windows

## Get list of variants as bedfile

mawk -v win=${WIN} '!/#/ {print $1"\t"$2-1"\t"$2}' ${VCF} > sites

## Now intersect them, counting overlaps in the windows with the sites
## And then calculating H(o) in the window size

bedtools intersect -a windows -b sites -c > counts

## Now launch Rscript to plot

Rscript ${SCRIPT_DIR}/hetPlotting.R ${VCF}

rm sites lengths windows counts
