root="/Data/zoc/result/10X-count/PBMC/10X-VDJ-human/5RNA"
ls $root | grep -v -E "DME|GVHD|PDR" > sample.txt
cat sample.txt | cut -d - -f 1 > disease.txt
paste sample.txt disease.txt > config.txt

