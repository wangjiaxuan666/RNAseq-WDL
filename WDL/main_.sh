#!/bin/bash  
annot_gtf=/jdfssz1/ST_HEALTH/P21Z10200N0092/wangjiaxuan/refer/index/hisat2/mu/Mus_musculus.GRCm38.84.gtf
python=/jdfssz1/ST_HEALTH/P21Z10200N0092/wangjiaxuan/biosoft/miniconda3/bin/python
collect_columns=/jdfssz1/ST_HEALTH/P21Z10200N0092/wangjiaxuan/biosoft/miniconda3/bin/collect-columns
samtools=/jdfssz1/ST_HEALTH/P21Z10200N0092/wangjiaxuan/biosoft/miniconda3/envs/samtools/bin/samtools
stringtie=/jdfssz1/ST_HEALTH/P21Z10200N0092/wangjiaxuan/biosoft/miniconda3/envs/RNAseq/bin/stringtie
preDE=/jdfssz1/ST_HEALTH/P21Z10200N0092/wangjiaxuan/workflow/RNAseq/util/preDE.py
read_length=100

cd 5.genome_mapping/
ls *.sorted.bam | while read bam
do
  sample=${bam%%.*}
  echo "${stringtie} -e -p 8 -G ${annot_gtf} -o ${sample}.result.gtf -e -A ${sample}.gene_abundances.tsv ${bam}" > ${sample}_main.sh
  python /jdfssz1/ST_HEALTH/P21Z10200N0092/wangjiaxuan/script/bin/qsub.cpython-38.pyc -s 2 -c 8 -l 8 -g 5g -r ${sample}_main.sh
  echo -e "${sample}\t${sample}.result.gtf" >> exp_gft.list
done

${python} ${preDE} -i exp_gft.list -l ${read_length}
fpkm=$(ls *.gene_abundances.tsv)
${collect_columns} gene_fpkm_matrix.xls ${fpkm} \
-c 7 \
-H \
-a ref_gene_id gene_name \
-g ${annot_gtf} \
-S

${collect_columns} gene_tpm_matrix.xls ${fpkm} \
-c 8 \
-H \
-a ref_gene_id gene_name \
-g ${annot_gtf} \
-S