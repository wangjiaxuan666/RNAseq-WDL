#!/bin/bash  
annot_gtf=/data/wangjiaxuan/refer/index/hisat2/mu/Mus_musculus.GRCm38.84.gtf
python=/data/wangjiaxuan/biosoft/miniconda3/bin/python
collect_columns=/data/wangjiaxuan/biosoft/miniconda3/bin/collect-columns
samtools=/data/wangjiaxuan/biosoft/miniconda3/envs/samtools/bin/samtools
stringtie=/data/wangjiaxuan/biosoft/miniconda3/envs/RNAseq/bin/stringtie
preDE=/data/wangjiaxuan/workflow/RNAseq/util/preDE.py
read_length=100
matedata_tsv=/hwfssz5/ST_HEALTH/P21Z10200N0092/tongyihan/SARS-CoV-2/3.RNAseq_again/info/sample_input_path.tsv
rrna_index=/data/wangjiaxuan/biosoft/rRNA_Data/BOWTIE2_index/rRNA
genome_index=/data/wangjiaxuan/refer/index/hisat2/mu/genome_snp_tran  
annot_gtf=/data/wangjiaxuan/refer/index/hisat2/mu/Mus_musculus.GRCm38.84.gtf
mapping_software=/data/wangjiaxuan/biosoft/hisat2-2.2.0/hisat2

# cd 4.rtRNA_remove/
# cat ${matedata_tsv} | while read group sample raw_fq1 raw_fq2
# do
#   rr_fq1=${sample}_rRNAremoved.fq.1.gz
#   rr_fq2=${sample}_rRNAremoved.fq.2.gz
  
#   echo "
#   ${mapping_software} \
#   -p 4 \
#   --dta \
#   -x ${genome_index} \
#   -1 ${rr_fq1} \
#   -2 ${rr_fq2} \
#   -S ${sample}_tmp.sam

#   ${samtools} sort ${sample}_tmp.sam -O bam -@ 4 -o ../5.genome_mapping/${sample}.sorted.bam && rm ${sample}_tmp.sam " >  ${sample}_genome_mapping.sh
#   python /data/wangjiaxuan/script/bin/qsub.cpython-38.pyc -s 2 -c 8 -l 8 -g 5g -r ${sample}_genome_mapping.sh
# done

cd 5.genome_mapping/
# ls *.sorted.bam | while read bam
# do
#   sample=${bam%%.*}
#   echo "${stringtie} -e -p 8 -G ${annot_gtf} -o ${sample}.result.gtf -e -A ${sample}.gene_abundances.tsv ${bam}" > ${sample}_main.sh
#   python /data/wangjiaxuan/script/bin/qsub.cpython-38.pyc -s 2 -c 8 -l 8 -g 5g -r ${sample}_main.sh
#   echo -e "${sample}\t${sample}.result.gtf" >> exp_gft.list
# done

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