#!/bin/bash

# This is the bash script to run lefse and generate bar plots using two input table (both tab separated)
# This script has to be run in a conda environment that has lefse and all the dependencies installed
# 1. pheno table that already transposed to horizontal position
# 2. cts table that have sample names correctly sorted according to the pheno table

pheno=$1
cts=$2
outdir=$3

mkdir -p $outdir

cat $pheno $cts > $outdir/lefse_input_table.tsv

format_input.py $outdir/lefse_input_table.tsv  $outdir/lefse_input_table.in -c 1 -u 2 -o 1000000

run_lefse.py $outdir/lefse_input_table.in $outdir/lefse_input_table.res

plot_res.py $outdir/lefse_input_table.res  $outdir/lefse_res_bar.png  --feature_font_size 4 --width 10 --dpi 300

plot_cladogram.py $outdir/lefse_input_table.res $outdir/lefse_res_cladogram.png --dpi 300 --format png --label_font_size 4
