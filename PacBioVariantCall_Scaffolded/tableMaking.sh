
file=$1

mawk '!/#/' ${file} | cut -f1,2,4,5,6 > temp1
mawk '!/#/' ${file} | cut -f8 | cut -d ";" -f1 | sed "s/DP=//g" > temp2
mawk '!/#/' ${file} | cut -f8 | cut -d ";" -f2 | sed "s/AC=//g" | tr "," "\t" > temp3
mawk '!/#/' ${file} | cut -f8 | cut -d ";" -f3 | sed "s/AM=//g" | tr "," "\t" > temp4
mawk '!/#/' ${file} | cut -f10 | cut -d ":" -f1 > temp5
mawk '!/#/' ${file} | cut -f10 | cut -d ":" -f3 > temp6
mawk '!/#/' ${file} | cut -f10 | cut -d ":" -f4 > temp7

paste temp1 temp2 temp3 temp4 temp5 temp6 temp7 > ${file}.table

rm temp*
