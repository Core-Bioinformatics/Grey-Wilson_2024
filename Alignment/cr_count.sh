#/bin/bash

out_dir=''
raw_dir=''

for p in $(cat rawSamplenames.txt);
do
echo "${p}"
echo "${p}" >> cr.log
mkdir -p ${out_dir}/${p}
~/software/cellranger/cellranger-7.2.0/cellranger \
count --id=${p} \
--transcriptome=~/software/cellranger/refdata-gex-GRCh38-2020-A/ \
--fastqs=${raw_dir}/${p}/ \
--localcores=100 \
--include-introns true \
--output-dir=${out_dir}/${p} \
--localmem=200
done
