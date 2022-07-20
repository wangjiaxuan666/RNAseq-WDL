rm(list = ls())
#---------
require(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)
library(DESeq2)
require(xbox)
#---------
files = list.files("7.exp/")
setwd("7.exp/")
#fpkm_data = read_tsv(grep("fpkm",files,value = T))
count_data = read_csv(grep("count",files,value = T))
tpm_data = read_tsv(grep("tpm",files,value = T))
#---------
count_data = count_data %>% separate(gene_id,into = c("feature","gene_name"),sep = "\\|")
count_data = count_data %>% filter(!grepl("MSTRG\\.\\d",feature)) %>% select(-c("feature"))
tpm_data = tpm_data %>% filter(!grepl("MSTRG\\.\\d",feature)) %>% select(-c("feature","ref_gene_id"))
#fpkm_data = fpkm_data %>% filter(!grepl("MSTRG\\.\\d",feature)) %>% select(-c("feature","ref_gene_id"))
setwd("../")
#---------
counts = count_data
fpkms = tpm_data
colnames(fpkms) = colnames(counts)
#---------
group = read_csv("group3.csv")
fpkms = fpkms[!duplicated(fpkms$gene_name),] %>% column_to_rownames("gene_name")
counts = counts[!duplicated(counts$gene_name),] %>% column_to_rownames("gene_name")
if(all(rownames(counts) == rownames(fpkms))){
  print("all gene rowname is same")
}else{
  print("not the same")
}
#------------------
if(!dir.exists("./result")) dir.create("./result")
outdir = "./result/1.All_exp"
xbox::chdir(outdir)
write.csv(counts,"All_gene_exp_count.csv")
write.csv(fpkms,"All_gene_exp_tpm_or_fpkm.csv")

# 设置随机种子
set.seed(20211108)
# 设置差异阈值
threshold_adj_p = 0.05
threshold_fc = 1

# 输入数据
# 表达矩阵
# 分组数据,第一列是样本名称,第二列是分组名称
class_data = group
colnames(class_data) = c("id","group")
class_data %>% distinct() %>% group_by(id) %>% dplyr::count() %>%filter(n >1)
class_data = class_data %>% distinct() %>% column_to_rownames("id")

vs = unique(class_data$group)
class_data$group = factor(class_data$group,levels = vs)

#----------find the intersection----------
setdiff(colnames(fpkms),colnames(counts))
setdiff(colnames(counts),colnames(fpkms))
setdiff(colnames(counts),rownames(class_data))
setdiff(rownames(class_data),colnames(fpkms))
#----check the sample------
class_data = arrange(class_data,group)
exp_data = counts[,rownames(class_data)]
fpkm_data = fpkms[,rownames(class_data)]

rm(group,counts,fpkms)
#++++++++++++++++++++++++++++++
# 差异分析部分的目录
#++++++++++++++++++++++++++++++
# if(!dir.exists("./result")) dir.create("./result")
if(!dir.exists("../2.DEG")) dir.create("../2.DEG")
xbox::chdir("../2.DEG")

require(autopca)
pca_data = pca_data_tidy(as.data.frame(fpkm_data))
pca(pca_data,rename = "diy",sample_group = class_data,add_ploy = "ellipse",point_size = 0.5,display_name = T)
ggsave(paste("The_fpkm_all_gene_sample_group_PCA.pdf",sep = ""),
       width = 8,
       height = 6,
       dpi = 300)
#------------
resname = "group_COB2_vs_CLN2"
vs = unlist(str_split(sub("group_","",resname),"_vs_"))
group_A = grep(vs[1],class_data$group)
group_B = grep(vs[2],class_data$group)
pca_data = pca_data_tidy(as.data.frame(fpkm_data[,c(group_A,group_B)]))
test = data.frame(row.names = rownames(class_data)[c(group_A,group_B)],group = class_data[c(group_A,group_B),])
pca(pca_data,rename = "diy",sample_group = test,add_ploy = "polygon",point_size = 3,display_name = T)
ggsave(paste(vs[1],"_vs_",vs[2],"PCA.pdf",sep = ""),
       width = 8,
       height = 6,
       dpi = 300)
#-----------
group_list=(class_data %>% group_by(group) %>% dplyr::count())[[1]]
group_list = as.vector(group_list)
num=(class_data %>% group_by(group) %>% dplyr::count())[[2]]

Batch_Deseq_differnece<-function(exprSet,group,num,save_dir){
  # exprSet = exp_data
  # group=group_list
  # num=num
  # save_dir = "test"
  ##create a folder 
  save_dir<-paste0(save_dir,"/")
  dir.create(save_dir)
  ## creat a group
  group_list= factor(rep(group,num))
  group_list
  colData=data.frame(row.names = colnames(exprSet),
                     group=group_list)
  
  for (i in 1:length(group)){
    name=unique(group)[i]
    print(name)
    colData$group<-relevel(colData$group,ref=name)
    dds=DESeq2::DESeqDataSetFromMatrix(countData = exprSet,
                                       colData = colData,
                                       design = ~group) 
    dds <- dds[ rowSums(DESeq2::counts(dds)) > 3, ]
    dds <- DESeq2::DESeq(dds)
    for (j in 2:length(DESeq2::resultsNames(dds))){
      
      resname=DESeq2::resultsNames(dds)[j]
      
      res = DESeq2::results(dds, name=resname)
      res_lfc <- lfcShrink(dds, coef=j, res=res, type="apeglm")
      
      #summary(res_lfc)
      summary(res)
      #res_lfc = res
      save_dir_MA=paste0(save_dir,"/",resname)
      dir.create(save_dir_MA)
      write.csv(res,paste0(save_dir_MA,"/",resname,"_res.csv"))
      write.csv(res_lfc,paste0(save_dir_MA,"/",resname,"_reslfc.csv"))
      png(paste0(save_dir_MA,"/",resname,"_MA.png"),width=600*3,height=3*600,res=72*3) 
      plotMA(res, ylim=c(-3,3),main=paste0(resname," MA"))
      dev.off()
      png(paste0(save_dir_MA,"/",resname,"_MAlfc.png"),width=600*3,height=3*600,res=72*3) 
      xlim <- c(1,1e5); ylim<-c(-3,3)
      plotMA( res_lfc, xlim=xlim, ylim=ylim, main=paste0(resname," apeglm"))
      dev.off()
      #-------
      # diff_gene_ds2 <- subset(res_lfc, padj < threshold_adj_p  & abs(log2FoldChange) > threshold_fc)
      # as.data.frame(diff_gene_ds2) -> diff_result
      #-------
      # res_lfc = read.csv("result/DEG_old/group_VLN2_vs_VOB2/group_VLN2_vs_VOB2_reslfc.csv",row.names = 1)
      # resname = "group_VLN2_vs_VOB2"
      # save_dir_MA = "result/DEG_old/group_VLN2_vs_VOB2/"
      vs = unlist(str_split(sub("group_","",resname),"_vs_"))
      group_A = grep(vs[1],class_data$group)
      group_B = grep(vs[2],class_data$group)
      pca_data = pca_data_tidy(as.data.frame(fpkm_data[,c(group_A,group_B)]))
      test = data.frame(row.names = rownames(class_data)[c(group_A,group_B)],group = class_data[c(group_A,group_B),])
      pca(pca_data,rename = "diy",sample_group = test,add_ploy = "ellipse",point_size = 3,display_name = T)
      ggsave(paste(save_dir_MA,"/",vs[1],"_vs_",vs[2],"PCA.pdf",sep = ""),
             width = 8,
             height = 6,
             dpi = 300)
      diff_gene_fc = col_ttest(fpkm_data+0.0001,g1 = group_A ,g2 = group_B,adjust = F,p = F)
      # get the gene pass pvalue filter
      diff_gene_ds2_p <- subset(res_lfc, padj < threshold_adj_p)
      as.data.frame(diff_gene_ds2_p) -> diff_gene_ds2_p
      cat("Have",nrow(diff_gene_ds2_p),"Gene pass the P <" ,threshold_adj_p)
      # get the gene pass log2fc filter
      diff_gene_ds2_fc <- subset(diff_gene_fc, abs(log2fc) > threshold_fc)
      cat("Have",nrow(diff_gene_ds2_fc),"Gene pass the |log2fc| >" ,threshold_fc)
      #--------------------------
      diff_id = intersect(rownames(diff_gene_ds2_fc),rownames(diff_gene_ds2_p))
      diff_result = cbind(diff_gene_ds2_fc[diff_id,c(group_A,group_B,ncol(diff_gene_ds2_fc))],diff_gene_ds2_p[diff_id,c("pvalue","padj")])
      all_result = cbind(diff_gene_fc[rownames(res_lfc),c(group_A,group_B,ncol(diff_gene_fc))],res_lfc[,c("pvalue","padj")])
      write.csv(diff_result, paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc.csv"))
      write.csv(exprSet[,c(group_A,group_B)],paste0(save_dir_MA,"/","all_exp_",resname,"_count.csv"))
      write.csv(fpkm_data[,c(group_A,group_B)],paste0(save_dir_MA,"/","all_exp_",resname,"_fpkm.csv"))
      if(nrow(diff_result) > 3){
        select_gene_number = 20
        heatmap_color = colorRampPalette(c("#44397d","white","#da2a46"),bias = 1)(100)
        #heatmap_color =  colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(20)
        result_id = c("log2fc","pvalue","padj")
        h_fc_id = diff_result[,result_id] %>% filter(log2fc > 0) %>% top_n(log2fc,n = select_gene_number) %>%rownames() %>% rev()
        l_fc_id = diff_result[,result_id] %>% filter(log2fc < 0) %>% top_n(log2fc,n = -select_gene_number) %>%rownames() %>% rev()
        pheatmap::pheatmap(log10(exp_data[c(h_fc_id,l_fc_id),c(group_A,group_B)]+1),
                           color = heatmap_color,
                           scale = "row",
                           annotation_col = class_data,
                           annotation_row = data.frame(row.names = c(h_fc_id,l_fc_id),Class= c(rep("UP",length(h_fc_id)),rep("DOWN",length(l_fc_id)))),
                           cluster_rows = F,
                           cluster_cols = F,
                           fontsize_row = 1,
                           angle_col = 0,
                           show_colnames = T,
                           filename = paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_pheatmap.pdf"),
                           width = 8,
                           height = 6)
        write.csv(exp_data[c(h_fc_id,l_fc_id),], paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_pheatmap_data.csv"))
        diff_result[h_fc_id,result_id] %>%
          rownames_to_column("id") %>%
          mutate(id = fct_reorder(id, log2fc)) %>%
          ggplot(aes(x = id,y = 2^log2fc,fill = -log10(padj)))+
          geom_col()+
          geom_text(aes(label = id),size = 1)+
          scale_fill_gradient(low = "#f1ad64",high = "#da2a46")+
          theme_void()+
          coord_polar()+
          ggtitle("The TOP20 Up Gene")->p1
        
        diff_result[l_fc_id,result_id] %>%
          rownames_to_column("id") %>%
          mutate(id = fct_reorder(id, log2fc)) %>%
          ggplot(aes(x = id,y = 2^-log2fc,fill = -log10(padj)))+
          geom_col()+
          geom_text(aes(label = id),size = 1)+
          scale_fill_gradient(low = "#95bbce",high = "#44397d")+
          theme_void()+
          coord_polar()+
          ggtitle("The TOP20 Down Gene")->p2
        require(patchwork)
        p1+p2
        ggsave(paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_luoxuan.pdf"),width = 16,height = 6,dpi = 300)
        #-----------------------------------------
        # 火山图
        vl_data = all_result[,result_id]
        vl_data$label <- as.factor(
          ifelse(vl_data$padj < threshold_adj_p & abs(vl_data$log2fc) > threshold_fc,
                 ifelse(vl_data$log2fc > threshold_fc ,'Up','Down'   ),'Stable'
          )
        )
        
        for_label <- vl_data[c(h_fc_id,l_fc_id),]
        for_label = for_label[-3,]
        
        plot2 = ggplot(data=vl_data, aes(x=log2fc, y=-log10(padj), colour=label)) +
          scale_color_manual(values=c("#44397d", "grey","#da2a46"))+
          geom_point(alpha=0.4, size=1.2) +
          theme_bw(base_size = 12, base_family = "Times") +
          geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
          geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
          theme(legend.position="right",
                panel.grid=element_blank(),
                legend.title = element_blank(),
                legend.text= element_text(face="bold", color="black",family = "Times", size=8),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(face="bold", color="black", size=12),
                axis.text.y = element_text(face="bold",  color="black", size=12),
                axis.title.x = element_text(face="bold", color="black", size=12),
                axis.title.y = element_text(face="bold",color="black", size=12))+
          labs(x="log2 (fold change)",y="-log10 (p-value)",title="Volcano picture of DEG")
        plot2+geom_point(size = 0, data = for_label, aes(colour=label, fill=label)) +
          ggrepel::geom_label_repel(
            aes(label = rownames(for_label)),
            max.overlaps = 100,
            data = for_label,
            show.legend = F,
            label.size = NA,
            fill = alpha(c("white"),0))
        ggsave(paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_huoshan.pdf"),width = 8,height = 6,dpi = 300)
        #------
        
        upid = diff_result %>% filter(log2fc > 0) %>% rownames()
        downid = diff_result %>% filter(log2fc < 0) %>% rownames()

        gene1 = bitr(downid,fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
        ego <- enrichGO(gene          = gene1$ENTREZID,
                        OrgDb         = org.Mm.eg.db,
                        keyType = "ENTREZID",
                        ont = "all",
                        pvalueCutoff = 0.5)
        if(!is.null(ego)){
          dotplot(ego)
          ggplot2::ggsave(paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_down_enrichment.pdf"),
                          width = 10,height = 8,dpi = 300)
          write.csv(ego@result,paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_down_enrichment_data.csv"),quote = FALSE)
        }

        gene1 = bitr(upid,fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
        ego <- enrichGO(gene          = gene1$ENTREZID,
                        OrgDb         = org.Mm.eg.db,
                        keyType = "ENTREZID",
                        ont = "all",
                        pvalueCutoff = 0.5)
        if(!is.null(ego)){
          dotplot(ego)
          ggplot2::ggsave(paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_up_enrichment.pdf"),
                          width = 10,height = 8,dpi = 300)
          write.csv(ego@result,paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_up_enrichment_data.csv"),quote = FALSE)
        }
      }
    }
  }
}
Batch_Deseq_differnece(exp_data,group=group_list,num,save_dir = "./")

# #---------NO STANTARD WORKFOLW-----------
# Batch_Deseq_differnece<-function(exprSet,group,num,save_dir){
#   # exprSet = exp_data
#   # group=group_list
#   # num=num
#   # save_dir = "test"
#   ##create a folder 
#   # save_dir<-paste0(save_dir,"/")
#   # dir.create(save_dir)
#   ## creat a group
#   group_list= factor(rep(group,num))
#   group_list
#   colData=data.frame(row.names = colnames(exprSet),
#                      group=group_list)
  
#   for (i in 1:length(group)){
#     name=unique(group)[i]
#     print(name)
#     colData$group<-relevel(colData$group,ref=name)
#     dds=DESeq2::DESeqDataSetFromMatrix(countData = exprSet,
#                                        colData = colData,
#                                        design = ~group) 
#     dds <- dds[ rowSums(DESeq2::counts(dds)) > 3, ]
#     dds <- DESeq2::DESeq(dds)
#     for (j in 2:length(DESeq2::resultsNames(dds))){
      
#       resname=DESeq2::resultsNames(dds)[j]
      
#       res = DESeq2::results(dds, name=resname)
#       res_lfc <- lfcShrink(dds, coef=j, res=res, type="apeglm")
      
#       #summary(res_lfc)
#       summary(res)
#       #res_lfc = res
#       save_dir_MA=paste0(save_dir,"/",resname)
#       dir.create(save_dir_MA)
#       write.csv(res,paste0(save_dir_MA,"/",resname,"_res.csv"))
#       write.csv(res_lfc,paste0(save_dir_MA,"/",resname,"_reslfc.csv"))
#       png(paste0(save_dir_MA,"/",resname,"_MA.png"),width=600*3,height=3*600,res=72*3) 
#       plotMA(res, ylim=c(-3,3),main=paste0(resname," MA"))
#       dev.off()
#       png(paste0(save_dir_MA,"/",resname,"_MAlfc.png"),width=600*3,height=3*600,res=72*3) 
#       xlim <- c(1,1e5); ylim<-c(-3,3)
#       plotMA( res_lfc, xlim=xlim, ylim=ylim, main=paste0(resname," apeglm"))
#       dev.off()
#       #-------
#       # diff_gene_ds2 <- subset(res_lfc, padj < threshold_adj_p  & abs(log2FoldChange) > threshold_fc)
#       # as.data.frame(diff_gene_ds2) -> diff_result
#       #-------
#       # res_lfc = read.csv("result/DEG_old/group_VLN2_vs_VOB2/group_VLN2_vs_VOB2_reslfc.csv",row.names = 1)
#       # resname = "group_VLN2_vs_VOB2"
#       # save_dir_MA = "result/DEG_old/group_VLN2_vs_VOB2/"
#       vs = unlist(str_split(sub("group_","",resname),"_vs_"))
#       group_A = grep(vs[1],class_data$group)
#       group_B = grep(vs[2],class_data$group)
#       pca_data = pca_data_tidy(as.data.frame(fpkm_data[,c(group_A,group_B)]))
#       test = data.frame(row.names = rownames(class_data)[c(group_A,group_B)],group = class_data[c(group_A,group_B),])
#       pca(pca_data,rename = "diy",sample_group = test,add_ploy = "ellipse",point_size = 3,display_name = T)
#       ggsave(paste(save_dir_MA,"/",vs[1],"_vs_",vs[2],"PCA.pdf",sep = ""),
#              width = 8,
#              height = 6,
#              dpi = 300)
#       diff_gene_fc = col_ttest(fpkm_data+0.0001,g1 = group_A ,g2 = group_B,adjust = F,p = F)
#       # get the gene pass pvalue filter
#       # ---change---
#       res_lfc$padj = qvalue::qvalue(res_lfc$pvalue)[["qvalues"]]
#       write.csv(res_lfc,paste0(save_dir_MA,"/",resname,"_reslfc.csv"))
#       #=============
#       diff_gene_ds2_p <- subset(res_lfc, padj < threshold_adj_p)
#       as.data.frame(diff_gene_ds2_p) -> diff_gene_ds2_p
#       cat("Have",nrow(diff_gene_ds2_p),"Gene pass the P <" ,threshold_adj_p)
#       # get the gene pass log2fc filter
#       diff_gene_ds2_fc <- subset(diff_gene_fc, abs(log2fc) > threshold_fc)
#       cat("Have",nrow(diff_gene_ds2_fc),"Gene pass the |log2fc| >" ,threshold_fc)
#       #--------------------------
#       diff_id = intersect(rownames(diff_gene_ds2_fc),rownames(diff_gene_ds2_p))
#       diff_result = cbind(diff_gene_ds2_fc[diff_id,c(group_A,group_B,ncol(diff_gene_ds2_fc))],diff_gene_ds2_p[diff_id,c("pvalue","padj")])
#       all_result = cbind(diff_gene_fc[rownames(res_lfc),c(group_A,group_B,ncol(diff_gene_fc))],res_lfc[,c("pvalue","padj")])
#       write.csv(diff_result, paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc.csv"))
#       write.csv(exprSet[,c(group_A,group_B)],paste0(save_dir_MA,"/","all_exp_",resname,"_count.csv"))
#       write.csv(fpkm_data[,c(group_A,group_B)],paste0(save_dir_MA,"/","all_exp_",resname,"_fpkm.csv"))
#       if(nrow(diff_result) > 3){
#         select_gene_number = 20
#         heatmap_color = colorRampPalette(c("#44397d","white","#da2a46"),bias = 1)(100)
#         #heatmap_color =  colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(20)
#         result_id = c("log2fc","pvalue","padj")
#         h_fc_id = diff_result[,result_id] %>% filter(log2fc > 0) %>% top_n(log2fc,n = select_gene_number) %>%rownames() %>% rev()
#         l_fc_id = diff_result[,result_id] %>% filter(log2fc < 0) %>% top_n(log2fc,n = -select_gene_number) %>%rownames() %>% rev()
#         pheatmap::pheatmap(log10(exp_data[c(h_fc_id,l_fc_id),c(group_A,group_B)]+1),
#                            color = heatmap_color,
#                            scale = "row",
#                            annotation_col = class_data,
#                            annotation_row = data.frame(row.names = c(h_fc_id,l_fc_id),Class= c(rep("UP",length(h_fc_id)),rep("DOWN",length(l_fc_id)))),
#                            cluster_rows = F,
#                            cluster_cols = F,
#                            fontsize_row = 1,
#                            angle_col = 0,
#                            show_colnames = T,
#                            filename = paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_pheatmap.pdf"),
#                            width = 8,
#                            height = 6)
#         write.csv(exp_data[c(h_fc_id,l_fc_id),], paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_pheatmap_data.csv"))
#         diff_result[h_fc_id,result_id] %>%
#           rownames_to_column("id") %>%
#           mutate(id = fct_reorder(id, log2fc)) %>%
#           ggplot(aes(x = id,y = 2^log2fc,fill = -log10(padj)))+
#           geom_col()+
#           geom_text(aes(label = id),size = 1)+
#           scale_fill_gradient(low = "#f1ad64",high = "#da2a46")+
#           theme_void()+
#           coord_polar()+
#           ggtitle("The TOP20 Up Gene")->p1
        
#         diff_result[l_fc_id,result_id] %>%
#           rownames_to_column("id") %>%
#           mutate(id = fct_reorder(id, log2fc)) %>%
#           ggplot(aes(x = id,y = 2^-log2fc,fill = -log10(padj)))+
#           geom_col()+
#           geom_text(aes(label = id),size = 1)+
#           scale_fill_gradient(low = "#95bbce",high = "#44397d")+
#           theme_void()+
#           coord_polar()+
#           ggtitle("The TOP20 Down Gene")->p2
#         require(patchwork)
#         p1+p2
#         ggsave(paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_luoxuan.pdf"),width = 16,height = 6,dpi = 300)
#         #-----------------------------------------
#         # 火山图
#         vl_data = all_result[,result_id]
#         vl_data$label <- as.factor(
#           ifelse(vl_data$padj < threshold_adj_p & abs(vl_data$log2fc) > threshold_fc,
#                  ifelse(vl_data$log2fc > threshold_fc ,'Up','Down'   ),'Stable'
#           )
#         )
        
#         for_label <- vl_data[c(h_fc_id,l_fc_id),]
#         for_label = for_label[-3,]
        
#         plot2 = ggplot(data=vl_data, aes(x=log2fc, y=-log10(padj), colour=label)) +
#           scale_color_manual(values=c("#44397d", "grey","#da2a46"))+
#           geom_point(alpha=0.4, size=1.2) +
#           theme_bw(base_size = 12, base_family = "Times") +
#           geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
#           geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
#           theme(legend.position="right",
#                 panel.grid=element_blank(),
#                 legend.title = element_blank(),
#                 legend.text= element_text(face="bold", color="black",family = "Times", size=8),
#                 plot.title = element_text(hjust = 0.5),
#                 axis.text.x = element_text(face="bold", color="black", size=12),
#                 axis.text.y = element_text(face="bold",  color="black", size=12),
#                 axis.title.x = element_text(face="bold", color="black", size=12),
#                 axis.title.y = element_text(face="bold",color="black", size=12))+
#           labs(x="log2 (fold change)",y="-log10 (p-value)",title="Volcano picture of DEG")
#         plot2+geom_point(size = 0, data = for_label, aes(colour=label, fill=label)) +
#           ggrepel::geom_label_repel(
#             aes(label = rownames(for_label)),
#             max.overlaps = 100,
#             data = for_label,
#             show.legend = F,
#             label.size = NA,
#             fill = alpha(c("white"),0))
#         ggsave(paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_huoshan.pdf"),width = 8,height = 6,dpi = 300)
#         #------
        
#         upid = diff_result %>% filter(log2fc > 0) %>% rownames()
#         downid = diff_result %>% filter(log2fc < 0) %>% rownames()
        
#         gene1 = bitr(downid,fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
#         ego <- enrichGO(gene          = gene1$ENTREZID,
#                         OrgDb         = org.Mm.eg.db,
#                         keyType = "ENTREZID",
#                         ont = "all",
#                         pvalueCutoff = 0.5)
#         if(!is.null(ego)){
#           dotplot(ego)
#           ggplot2::ggsave(paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_down_enrichment.pdf"),
#                           width = 10,height = 8,dpi = 300)
#           write.csv(ego@result,paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_down_enrichment_data.csv"),quote = FALSE)
#         }
        
#         gene1 = bitr(upid,fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
#         ego <- enrichGO(gene          = gene1$ENTREZID,
#                         OrgDb         = org.Mm.eg.db,
#                         keyType = "ENTREZID",
#                         ont = "all",
#                         pvalueCutoff = 0.5)
#         if(!is.null(ego)){
#           dotplot(ego)
#           ggplot2::ggsave(paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_up_enrichment.pdf"),
#                           width = 10,height = 8,dpi = 300)
#           write.csv(ego@result,paste0(save_dir_MA,"/","diff_exp_",resname,"_reslfc_up_enrichment_data.csv"),quote = FALSE)
#         }
#       }
#     }
#   }
# }
# Batch_Deseq_differnece(exp_data,group=group_list,num,save_dir = "./")


# #-------------------
# data <- read_csv("group_CLN2_vs_COB2/all_exp_group_CLN2_vs_COB2_fpkm.csv")
