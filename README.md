# Bulk RNAseq 分析流程

本流程是基于WDL构建的转录组分析流程，采用的分析软件为：
1. urmap比对去除核糖体
2. fastp质控
3. Hisat2比对基因组
4. stringtie组装新转录本并定量


后续待完善的是：
1. 用DEseq2做差异分析
2. 用clusterprofile做腹肌分析
3. 用WGCNA做基因权重共表达网络分析

# 流程图

![image-20220309111821786](https://pic-1259340288.cos.ap-guangzhou.myqcloud.com/img/202203091118221646795902hMsqYAimage-20220309111821786.png)

# 运行

## 输入文件

输入文件是json格式，需要提供对应文件和软件的路径。 这个文件是非常重要的输入文件，其中RNAseq.matedata_tsv对应的表格必须是制表符分割的tsv文件，包括四列，
第一列是分组，第二列是样本名称，第三列是read1的fq文件，第四列是read2的fq文件。

> 注意要严格按照规范，整理输入文件，否则必报错或者输出结果不对。

```json
{
    "RNAseq.matedata_tsv" : "./1.pre_info/sample_input_path.tsv",# 输入文件的表格，
    "RNAseq.urmap_ufi" : "~/biosoft/rRNA_Data/URMAP_index/rRNA.dna.8.ufi", # 去核糖体的序列index，不用修改
    "RNAseq.genome_index":"./chrX_data/indexes/chrX_tran",  # hisat2构建的全基因组索引
    "RNAseq.annot_gtf" : "./chrX_data/genes/chrX.gtf",# 基因组的注释文件gtf
    "RNAseq.read_length" : "150", # read 长度，一般而言，华大是100，illumina是150
    #------------下面全是软件的文件途径---------
    "RNAseq.python" : "~/biosoft/miniconda3/bin/python",
    "RNAseq.collect_columns" : "~/biosoft/miniconda3/bin/collect-columns",
    "RNAseq.samtools" : "~/biosoft/miniconda3/envs/samtools/bin/samtools",
    "RNAseq.stringtie" : "~/biosoft/miniconda3/envs/RNAseq/bin/stringtie",
    "RNAseq.gffcompare" : "~/biosoft/miniconda3/envs/SSR/bin/gffcompare",
    "RNAseq.preDE" : "WDL/preDE.py",
    "RNAseq.urmap" : "~/biosoft/urmap/bin/urmap",
    "RNAseq.fastp" : "~/biosoft/miniconda3/envs/RNAseq/bin/fastp",
    "RNAseq.mapping_software" : "~/biosoft/hisat2-2.2.0/hisat2",
    "RNAseq.pigz" : "~/biosoft/miniconda3/envs/RNAseq/bin/pigz"
}
```

## 运行命令

前期准备：
1）下载test数据集合：

```bash
mkdir test
cd test
bash download_test_data.sh
```
根据要求的matedata_tsv要求的格式整理表格，修改`WDL/input.RNAseq.json`

2) 运行需要有cromwell.jar文件，所以需要先下载。

```bash
bash download_cromwell.sh
```

3) 然后运行：

```
./run_RNAseq -i WDL/input.RNAseq.json
```

4) 等待结果出来

## 输出结果

1. 质控过滤的fq
2. 去除核糖体reads的fq
3. hisat2 比对文件 bam
4. 重构转录组 gtf 
5. fpkm_matrix; tpm_matrix; count_matrix