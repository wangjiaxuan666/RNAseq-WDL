version development

workflow RNAseq {
    input {
        File matedata_tsv
        String pigz
        String fastp
        String bowtie2
        String rrna_index
        String mapping_software
        String genome_index
        String samtools
        String annot_gtf
        String stringtie
        String gffcompare
        Int read_length
        String python
        File preDE
        String collect_columns
    }
    
    Array[Array[String]] sample_tsv = read_tsv(matedata_tsv)

    scatter ( sample_sheet in sample_tsv ){
        call getinfo {
            input : sample_sheet = sample_sheet
        }

        call fastp_qc {
            input :
                fastp = fastp,
                read1 = getinfo.read1,
                read2 = getinfo.read2,
                sample = getinfo.sample
        }

        call remove_rtRNA {
            input :
                bowtie2 = bowtie2,
                read1 = fastp_qc.filter_fq1,
                read2 = fastp_qc.filter_fq2,
                rrna_index = rrna_index,
                sample = getinfo.sample,
                samtools  = samtools
        }

        call genome_map {
            input :
            mapping_software = mapping_software,
            genome_index = genome_index,
            sample = getinfo.sample,
            read1 = remove_rtRNA.non_rrna_read1,
            read2 = remove_rtRNA.non_rrna_read2,
            samtools  = samtools
        }

        call stringtie_reGTF {
            input :
            sample = getinfo.sample,
            stringtie = stringtie,
            annot_gtf = annot_gtf,
            bam = genome_map.bam
        }
    }

    call stringtie_mergeGTF {
        input :
        stringtie = stringtie,
        gtf = stringtie_reGTF.gtf,
        annot_gtf = annot_gtf,
        gffcompare = gffcompare
    }

    scatter ( i in range(length(getinfo.sample))){
        call stringtie_quantification1 {
            input :
            bam = "${genome_map.bam[i]}",
            gtf = stringtie_mergeGTF.gtf,
            sample = "${getinfo.sample[i]}",
            stringtie = stringtie
        }
    }

    call get_expression {
        input :
            python = python,
            preDE = preDE,
            sample = getinfo.sample,
            count_s = stringtie_quantification1.count_s,
            read_length = read_length,
            collect_columns = collect_columns,
            fpkm = stringtie_quantification1.fpkm,
            gtf = stringtie_mergeGTF.gtf
    }

    output {
      Array[String] sample = getinfo.sample
      Array[File] bam = genome_map.bam
      Array[File] nrrna_read1 = remove_rtRNA.non_rrna_read1
      Array[File] nrrna_read2 = remove_rtRNA.non_rrna_read2
      Array[File] filter_read1 = fastp_qc.filter_fq1
      Array[File] filter_read2 = fastp_qc.filter_fq2
      File gtf = stringtie_mergeGTF.gtf
      File count_matrix = get_expression.count
      File fpkm_matrix = get_expression.fpkm
      File tpm_matrix = get_expression.tpm
    }
}

task getinfo {
    input {
        Array[String] sample_sheet
        String sample_group = sample_sheet[0]
        String sample_name = sample_sheet[1]
        String sample_read1 = sample_sheet[2]
        String sample_read2 = sample_sheet[3]
    }
    command {
        echo "yes,I get the info"
    }

    output {
        String sample  = sample_name
        String group = sample_group
        String read1 = sample_read1
        String read2 = sample_read2
    }

}

task fastp_qc {
    input {
        String fastp
        String read1
        String read2
        String sample
    }

    command {
        ${fastp} \
        -i ${read1} \
        -I ${read2} \
        -o ${sample}.filter.R1.fq.gz \
        -O ${sample}.filter.R2.fq.gz
    }

    output {
        File filter_fq1 = "${sample}.filter.R1.fq.gz"
	    File filter_fq2 = "${sample}.filter.R2.fq.gz"
        File filter_report_json = "fastp.json"
        File filter_report_html = "fastp.html"
    }
}

task remove_rtRNA {
    input {
        String bowtie2
        File read1
        File read2
        String rrna_index
        String sample
        String samtools
    }
    
    command {
        # ${bowtie2} \
        # -ufi ${rrna_index} \
        # -map2  ${read1} \
        # -reverse ${read2} \
        # -samout /dev/stdout | \
        # ${samtools} fastq \
        # -f 13 \
        # -1 ${sample}.filter.non_rrna.R1.fq.gz \
        # -2 ${sample}.filter.non_rrna.R2.fq.gz \
        # -@ 4

        ${bowtie2} \
        --very-sensitive-local \
        --no-unal \
        -I 1 \
        -X 1000 \
        -p 6 \
        -x ${rrna_index} \
        -1 ${read1} \
        -2 ${read2} \
        -S ${sample}.rRNA_mapping.sam \
        --un-conc-gz ${sample}_rRNAremoved.fq.gz 2>${sample}_Map2rRNAStat.xls 
        
        ${samtools} view -S -b -o ${sample}_Map2rRNA.bam ${sample}.rRNA_mapping.sam 
        rm ${sample}.rRNA_mapping.sam 
    }

    output {
        File non_rrna_read1 = "${sample}_rRNAremoved.fq.1.gz"
        File non_rrna_read2 = "${sample}_rRNAremoved.fq.2.gz"
    }
}


task genome_map {
    input {
       String mapping_software
       String genome_index
       File read1
       File read2
       String sample
       String samtools
    }

    command {
        ${mapping_software} \
        -p 4 \
        --dta \
        -x ${genome_index} \
        -1 ${read1} \
        -2 ${read2} \
        -S tmp.sam 
        ${samtools} sort tmp.sam -O bam -@ 4 -o ${sample}.sorted.bam && rm tmp.sam
    }

    output {
        File bam = "${sample}.sorted.bam"
    }
}

task stringtie_reGTF {
    input {
        String stringtie
        File bam
        String annot_gtf
        String sample
    }

    command {
        ${stringtie} \
        ${bam} \
        -l ${sample} \
        -p 4 \
        -G ${annot_gtf} \
        -o ${sample}_stringtie.gtf
    }

    output {
        File gtf = "${sample}_stringtie.gtf"
    }
}

task stringtie_mergeGTF {
    input {
        Array[String] gtf
        String annot_gtf
        String stringtie
        String gffcompare
    }

    command {
        ${stringtie} \
        --merge \
        ${write_lines(gtf)} \
        -p 8 \
        -G ${annot_gtf} \
        -o stringtie_merge.gtf

        ${gffcompare} -r ${annot_gtf} -G -o stringtie_merge stringtie_merge.gtf
    }

    output {
        File gtf = "stringtie_merge.gtf"
        File stats = "stringtie_merge.stats"
    }
}

task stringtie_quantification1 {
    input {
        String bam
        String gtf
        String stringtie
        String sample
    }

    command {
        ${stringtie} -e \
        -p 8 \
        -G ${gtf} \
        -o ${sample}.result.gtf \
        -A ${sample}.gene_abundances.tsv \
        ${bam}
    }

    output {
        File fpkm = "${sample}.gene_abundances.tsv"
        File count_s = "${sample}.result.gtf"
    }
    
}

task get_expression {
    input {
        String python
        File preDE
        Int read_length
        Array[File] count_s
        Array[String] sample
        Array[Array[String]] count_s_list = transpose([sample, count_s])
        String collect_columns
        Array[String]+ fpkm
        File gtf
    }

    command {

        ${python} ${preDE} -i ${write_tsv(count_s_list)} -l ${read_length}
        
        ${collect_columns} gene_fpkm_matrix.xls ${sep = ' ' fpkm} \
        -c 7 \
        -H \
        -a ref_gene_id gene_name \
        -g ${gtf} \
        -n ${sep = ' ' sample} \
        -S

        ${collect_columns} gene_tpm_matrix.xls ${sep = ' ' fpkm} \
        -c 8 \
        -H \
        -a ref_gene_id gene_name \
        -g ${gtf} \
        -n ${sep = ' ' sample} \
        -S
    }

    output {
        File count = "gene_count_matrix.csv"
        File fpkm = "gene_fpkm_matrix.xls"
        File tpm = "gene_tpm_matrix.xls"
    }
}