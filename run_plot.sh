Rscript=/home/ye/Software/Python/location/envs/10X/bin/Rscript
dataset=$1

$Rscript plot.R --dataset $dataset  --markers "CD38" "CXCR5" "CCR7" "CD28" "HLA-DRB1" "LILRB1" "CD86" "CCR6" "ICOS" "CXCR3" "CD3D" "CD8A" "CD4" "IL7R"

#$Rscript classifyCells_v2.R --dataset $dataset

#$Rscript StreamAnalysis.R --dataset $dataset
