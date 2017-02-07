#Dr Wargo  second project

library(Packzhu)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)

allunique<-function(v){
  if (any(duplicated(v))){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

endsWith<-function(str,patten){
  ends<-paste(patten,"$",sep='')
  return(grepl(ends,str))
}
get_data<-function(...){
  
  file=file.choose()
  if (endsWith(file,"txt")){
    data<-read.table(file,...)
  }
  else {
    data<-read.csv(file,...)
  }
  
  return(data)
}

selectbycolumnin<-function(df,col,vec){
  col<-substitute(col)
  col<-deparse(col)
  col<-df[,col]
  selected<-sapply(col,function(x) x %in% vec)
  return(df[selected,,drop=F])
}

n_top_genes_with_most_mad<-function(df,n=1000){
  mads<-apply(df,1,mad)
  mads<-sort(mads,decreasing=T)
  data<-df[names(mads),]
  data<-data[1:n,]
  return(data)
}

triplevenn<-function(df1,df2,df3,venn_names,higher=FALSE,lower=FALSE){
  if(higher){
    df1<-df1[df1$log2FoldChange<0,]
    df2<-df2[df1$log2FoldChange<0,]
    df3<-df1[df3$log2FoldChange<0,]
  }
  if(lower){
    df1<-df1[df1$log2FoldChange>0,]
    df2<-df2[df1$log2FoldChange>0,]
    df3<-df1[df3$log2FoldChange>0,]
  }
  lists<-list(rownames(df1),rownames(df2),rownames(df3))
  names(lists)=venn_names
  area1=length(lists[[1]])
  area2=length(lists[[2]])
  area3=length(lists[[3]])
  n12=length(intersect(lists[[1]],lists[[2]]))
  n23=length(intersect(lists[[2]],lists[[3]]))
  n13=length(intersect(lists[[1]],lists[[3]]))
  n123=length(intersect(intersect(lists[[1]],lists[[2]]),lists[[3]]))
  category=names(lists)
  col=c("green","blue","red")
  grid.newpage()
  draw.triple.venn(area1=area1,area2=area2,area3=area3,n12=n12,n23=n23,n13=n13,n123=n123,
                   category=category,col=col,cat.col=col)
}






# directory is C:\Users\Woodmanlab\Desktop\Wargo II
setwd("C:/Users/Woodmanlab/Desktop/Wargo II")

# read raw count data from 0409_rawcounts.txt
expression<-get_data(header=T,row.names=1,stringsAsFactors=F,sep="")

#clean expression data
ensembl_id<-sapply(expression$gene_id,function(x) unname(unlist(strsplit(x,"\\."))[1]))
expression$gene_id<-unname(ensembl_id)
expression<-aggregate(.~gene_id+gene_name,data=expression,sum)
rownames(expression)<-expression$gene_id
expression<-expression[,-1]
expression$gene_name<-sapply(expression$gene_name,function(x) trim(x))


#read protein coding RNA list from protein coding RNA.csv
pcRNA_genes<-get_data(header=T,stringsAsFactors=F)
#genes_with_ensembl<-as.character(pcRNA_genes$Ensembl_ID[pcRNA_genes$Ensembl_ID!=""])
#genes_wo_ensembl<-as.character(pcRNA_genes$Symbol[pcRNA_genes$Ensembl_ID==""])

#extract protein coding RNA data,clean and save to protein coding rna expression.csv
selected_pcRNA<-sapply(expression$gene_name,function(x) x%in%pcRNA_genes$Symbol)
pcRNA_expression<-expression[selected_pcRNA,]
pcRNA_expression<-aggregate(.~gene_name,data=pcRNA_expression,sum)
rownames(pcRNA_expression)<-pcRNA_expression$gene_name
pcRNA_expression<-pcRNA_expression[,-1]
colnames(pcRNA_expression)<-sapply(colnames(pcRNA_expression),function(x) gsub("X","",x))
write.csv(pcRNA_expression,"protein coding rna expression.csv")

#read profile data from 0409_sampleinfo.txt
profile<-get_data(header=T,row.names=1,stringsAsFactors=F,sep="")
profile<-t(profile)
profile<-as.data.frame(profile)
rownames(profile)<-sapply(rownames(profile),function(x) gsub("X","",x))

profile$treatment<-factor(paste(profile$patient,profile$timepoint,sep="_"))
profile$path="pCR"
profile$path[profile$patient==96]<-"NR"
profile$path[profile$patient==80]<-"unknown"
profile$path[profile$patient==97]<-"non_pCR"
profile$path[profile$patient==109]<-"non_pCR"
profile$path[profile$patient==129]<-"non_pCR"
profile$path[profile$patient==163]<-"non_pCR"

pCR_profile_on_pre<-profile[(profile$timepoint!="SURG")&(profile$path=="pCR"),]
pCR_profile_on_pre$timepoint<-factor(as.character(pCR_profile_on_pre$timepoint))
pCR_data<-pcRNA_expression[,rownames(pCR_profile_on_pre)]
non_pCR_profile_on_pre<-profile[(profile$timepoint!="SURG")&(profile$path=="non_pCR"),]
non_pCR_profile_on_pre$timepoint<-factor(as.character(non_pCR_profile_on_pre$timepoint))
non_pCR_data<-pcRNA_expression[,rownames(non_pCR_profile_on_pre)]

patients_included<-c("97","109","163","85","145","153","172")

profile_included<-selectbycolumnin(profile,patient,patients_included)

pcRNA_expression_included<-pcRNA_expression[,rownames(profile_included)]

included_dds<-DESeqDataSetFromMatrix(countData = pcRNA_expression_included,
                                     colData=profile_included,
                                     design=~timepoint)
included_dds$sample<-profile_included$treatment
included_dds$run<-paste0("run",1:nrow(profile_included))
included_ddsColl<-collapseReplicates(included_dds,included_dds$sample,included_dds$run)
included_dds<-DESeq(included_ddsColl)
included_rlog<-rlog(included_dds)
included_log_data<-assay(included_rlog)
data<-n_top_genes_with_most_mad(included_log_data,n=1000)
profile_coll<-colData(included_dds)
ha<-HeatmapAnnotation(df=profile_coll[,c(3,6),drop=F])
Heatmap(t(scale(t(data))),cluster_columns = T,top_annotation = ha,name="Expression",row_names_gp = gpar(fontsize=1))
for(an in colnames(profile_coll[,c(3,6),drop=F])) {
  decorate_annotation(an, {
    grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
  })
}

#comparision with pre and surg of #109,97 and 163
nonpCR_profile_included<-profile_included[profile_included$patient %in% c("97","109","163"),]
nonpCR_profile_included<-nonpCR_profile_included[-c(1,2),]
nonpCR_profile_included$patient<-as.character(nonpCR_profile_included$patient)
nonpCR_profile_included$timepoint<-as.character(nonpCR_profile_included$timepoint)
nonpCR_pcRNA_expression_included<-pcRNA_expression[,rownames(nonpCR_profile_included)]
nonpCR_included_dds<-DESeqDataSetFromMatrix(countData = nonpCR_pcRNA_expression_included,
                                     colData=nonpCR_profile_included,
                                     design=~patient+timepoint)
nonpCR_included_dds$sample<-nonpCR_profile_included$treatment
nonpCR_included_dds$run<-paste0("run",1:nrow(nonpCR_profile_included))
nonpCR_included_ddsColl<-collapseReplicates(nonpCR_included_dds,nonpCR_included_dds$sample,nonpCR_included_dds$run)
nonpCR_included_dds<-DESeq(nonpCR_included_ddsColl)
nonpCR_results<-results(nonpCR_included_dds,contrast=c("timepoint","PRE","SURG"))
nonpCR_results<-nonpCR_results[order(nonpCR_results$padj),]
nonpCR_results<-nonpCR_results[complete.cases(nonpCR_results),]
nonpCR_results<-nonpCR_results[(nonpCR_results$padj<=0.05)&(abs(nonpCR_results$log2FoldChange)>=1),]
write.csv(nonpCR_results,"comparision between pre and surg in 97_109_163 nonpCR.csv")
nonpCR_included_rlog<-rlog(nonpCR_included_dds)
nonpCR_included_log_data<-assay(nonpCR_included_rlog)
nonpCR_profile_coll<-colData(nonpCR_included_dds)
heatmap_data<-nonpCR_included_log_data[rownames(nonpCR_results),]
write.csv(heatmap_data,"data for comparision between pre and surg in 97_109_163 nonpCR.csv")
ha<-HeatmapAnnotation(df=nonpCR_profile_coll[,c(3,6),drop=F],col=list(timepoint=c("PRE"="red","SURG"="blue"),path=c("non_pCR"="orange","pCR"="green")))
Heatmap(t(scale(t(heatmap_data))),col=colorRamp2(c(-2,0,2),c("green","black","red")),cluster_columns = T,top_annotation = ha,name="Expression",row_names_gp = gpar(fontsize=3))
for(an in colnames(nonpCR_profile_coll[,c(3,6),drop=F])) {
  decorate_annotation(an, {
    grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
  })
}

data_t<-t(scale(t(heatmap_data)))

#comparision with pre and surg of #109,97 and 163
# also include 109on
nonpCR_profile_included<-profile_included[profile_included$patient %in% c("97","109","163"),]
nonpCR_profile_included$timepoint[c(1,2)]<-"SURG"
nonpCR_profile_included$patient<-as.character(nonpCR_profile_included$patient)
nonpCR_profile_included$timepoint<-as.character(nonpCR_profile_included$timepoint)
nonpCR_pcRNA_expression_included<-pcRNA_expression[,rownames(nonpCR_profile_included)]
nonpCR_included_dds<-DESeqDataSetFromMatrix(countData = nonpCR_pcRNA_expression_included,
                                            colData=nonpCR_profile_included,
                                            design=~patient+timepoint)
nonpCR_included_dds$sample<-nonpCR_profile_included$treatment
nonpCR_included_dds$run<-paste0("run",1:nrow(nonpCR_profile_included))
nonpCR_included_ddsColl<-collapseReplicates(nonpCR_included_dds,nonpCR_included_dds$sample,nonpCR_included_dds$run)
nonpCR_included_dds<-DESeq(nonpCR_included_ddsColl)
nonpCR_results<-results(nonpCR_included_dds,contrast=c("timepoint","PRE","SURG"))
nonpCR_results<-nonpCR_results[order(nonpCR_results$padj),]
nonpCR_results<-nonpCR_results[complete.cases(nonpCR_results),]
nonpCR_results_2<-nonpCR_results[(nonpCR_results$padj<=0.05)&(abs(nonpCR_results$log2FoldChange)>=1),]
write.csv(nonpCR_results,"comparision between pre and surg in 97_109_163 and 109_on nonpCR.csv")
nonpCR_included_rlog<-rlog(nonpCR_included_dds)
nonpCR_included_log_data<-assay(nonpCR_included_rlog)
nonpCR_profile_coll<-colData(nonpCR_included_dds)
nonpCR_profile_coll$timepoint<-as.character(nonpCR_profile_coll$timepoint)
nonpCR_profile_coll$timepoint[1]<-"ON"
heatmap_data<-nonpCR_included_log_data[rownames(nonpCR_results),]
write.csv(heatmap_data,"data for comparision between pre and surg in 97_109_163 and 109_on nonpCR.csv")
ha<-HeatmapAnnotation(df=nonpCR_profile_coll[,c(3,6),drop=F],col=list(timepoint=c("PRE"="red","SURG"="blue","ON"="purple"),path=c("non_pCR"="orange","pCR"="green")))
Heatmap(t(scale(t(heatmap_data))),col=colorRamp2(c(-2,0,2),c("green","black","red")),cluster_columns = T,top_annotation = ha,name="Expression",row_names_gp = gpar(fontsize=3))
for(an in colnames(nonpCR_profile_coll[,c(3,6),drop=F])) {
  decorate_annotation(an, {
    grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
  })
}

#comparision with pre and on of #85,145,153,172
pCR_profile_included<-profile_included[profile_included$patient %in% c("85","145","153","172"),]
pCR_profile_included<-pCR_profile_included[pCR_profile_included$timepoint!="SURG",]
pCR_profile_included$patient<-as.character(pCR_profile_included$patient)
pCR_profile_included$timepoint<-as.character(pCR_profile_included$timepoint)
pCR_pcRNA_expression_included<-pcRNA_expression[,rownames(pCR_profile_included)]
pCR_included_dds<-DESeqDataSetFromMatrix(countData = pCR_pcRNA_expression_included,
                                            colData=pCR_profile_included,
                                            design=~patient+timepoint)
pCR_included_dds$sample<-pCR_profile_included$treatment
pCR_included_dds$run<-paste0("run",1:nrow(pCR_profile_included))
pCR_included_ddsColl<-collapseReplicates(pCR_included_dds,pCR_included_dds$sample,pCR_included_dds$run)
pCR_included_dds<-DESeq(pCR_included_ddsColl)
pCR_results<-results(pCR_included_dds,contrast=c("timepoint","PRE","ON"))
pCR_results<-pCR_results[order(pCR_results$padj),]
pCR_results<-pCR_results[complete.cases(pCR_results),]
pCR_results_1<-pCR_results[(pCR_results$padj<=0.05)&(abs(pCR_results$log2FoldChange)>=1),]
write.csv(pCR_results,"comparision between pre and on in 85_145_153_172 pCR.csv")
pCR_included_rlog<-rlog(pCR_included_dds)
pCR_included_log_data<-assay(pCR_included_rlog)
pCR_profile_coll<-colData(pCR_included_dds)
heatmap_data<-pCR_included_log_data[rownames(pCR_results),]
write.csv(heatmap_data,"data for comparision between pre and on in 85_145_153_172 pCR.csv")
ha<-HeatmapAnnotation(df=pCR_profile_coll[,c(3,6),drop=F],col=list(timepoint=c("PRE"="red","SURG"="blue","ON"="purple"),path=c("non_pCR"="orange","pCR"="green")))
Heatmap(t(scale(t(heatmap_data))),col=colorRamp2(c(-2,0,2),c("green","black","red")),cluster_columns = T,top_annotation = ha,name="Expression",row_names_gp = gpar(fontsize=1))
for(an in colnames(pCR_profile_coll[,c(3,6),drop=F])) {
  decorate_annotation(an, {
    grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
  })
}

#comparision with pre and surg of #85,145,153,172
pCR_profile_included<-profile_included[profile_included$patient %in% c("85","145","153","172"),]
pCR_profile_included<-pCR_profile_included[pCR_profile_included$timepoint!="ON",]
pCR_profile_included$patient<-as.character(pCR_profile_included$patient)
pCR_profile_included$timepoint<-as.character(pCR_profile_included$timepoint)
pCR_pcRNA_expression_included<-pcRNA_expression[,rownames(pCR_profile_included)]
pCR_included_dds<-DESeqDataSetFromMatrix(countData = pCR_pcRNA_expression_included,
                                         colData=pCR_profile_included,
                                         design=~patient+timepoint)
pCR_included_dds$sample<-pCR_profile_included$treatment
pCR_included_dds$run<-paste0("run",1:nrow(pCR_profile_included))
pCR_included_ddsColl<-collapseReplicates(pCR_included_dds,pCR_included_dds$sample,pCR_included_dds$run)
pCR_included_dds<-DESeq(pCR_included_ddsColl)
pCR_results<-results(pCR_included_dds,contrast=c("timepoint","PRE","SURG"))
pCR_results<-pCR_results[order(pCR_results$padj),]
pCR_results<-pCR_results[complete.cases(pCR_results),]
pCR_results_2<-pCR_results[(pCR_results$padj<=0.05)&(abs(pCR_results$log2FoldChange)>=1),]
write.csv(pCR_results,"comparision between pre and surg in 85_145_153_172 pCR.csv")
pCR_included_rlog<-rlog(pCR_included_dds)
pCR_included_log_data<-assay(pCR_included_rlog)
pCR_profile_coll<-colData(pCR_included_dds)
heatmap_data<-pCR_included_log_data[rownames(pCR_results),]
write.csv(heatmap_data,"data for comparision between pre and surg in 85_145_153_172 pCR.csv")
ha<-HeatmapAnnotation(df=pCR_profile_coll[,c(3,6),drop=F],col=list(timepoint=c("PRE"="red","SURG"="blue","ON"="purple"),path=c("non_pCR"="orange","pCR"="green")))
Heatmap(t(scale(t(heatmap_data))),col=colorRamp2(c(-2,0,2),c("green","black","red")),cluster_columns = T,top_annotation = ha,name="Expression",row_names_gp = gpar(fontsize=1))
for(an in colnames(pCR_profile_coll[,c(3,6),drop=F])) {
  decorate_annotation(an, {
    grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
  })
}

#draw venndiagram 
#non_pCR,pre and surg of #109,97 and 163
#pCR,pre and on of #85,145,153,172
#pCR,pre and surg of #85,145,153,172

triplevenn(nonpCR_results_1,pCR_results_1,pCR_results_2,c("pre_surg_in_non_pCR","pre_on_in_pCR"," pre_surg_in_pCR"))
triplevenn(nonpCR_results_2,pCR_results_1,pCR_results_2,c("pre_surg_in_non_pCR","pre_on_in_pCR"," pre_surg_in_pCR"))
triplevenn(nonpCR_results_2,pCR_results_1,pCR_results_2,
           c("pre_surg_in_non_pCR","pre_on_in_pCR"," pre_surg_in_pCR"),
           higher=T)




pCR_dds<-DESeqDataSetFromMatrix(countData = pCR_data,
                                colData=pCR_profile_on_pre,
                                design=~timepoint)
pCR_dds$sample<-pCR_profile_on_pre$treatment
pCR_dds$run<-paste0("run",1:nrow(pCR_profile_on_pre))
pCR_ddsColl<-collapseReplicates(pCR_dds,pCR_dds$sample,pCR_dds$run)
pCR_dds<-DESeq(ddsColl)
pCR_result<-results(pCR_dds)
pCR_result<-pCR_result[order(pCR_result$padj),]
pCR_selected_results<-pCR_result[complete.cases(pCR_result),]
sig_genes_pCR_on_pre<-rownames(pCR_selected_results[(abs(pCR_selected_results$log2FoldChange)>=1)&(pCR_selected_results$padj<=0.05),])

pCR_log<-rlog(pCR_ddsColl)
pCR_log_data<-assay(pCR_log)
pCR_log_data_selected<-pCR_log_data[sig_genes_pCR_on_pre,]
coll_profile<-colData(pCR_ddsColl)
coll_profile<-coll_profile[order(coll_profile$timepoint),]
pCR_log_data_selected<-pCR_log_data_selected[,rownames(coll_profile)]
ha<-HeatmapAnnotation(df=coll_profile[,3,drop=F])
Heatmap(t(scale(t(pCR_log_data_selected))),cluster_columns = F,top_annotation = ha,name="Expression",row_names_gp = gpar(fontsize=1))
write.csv(pCR_result,"comparision between pre vs on in pCR groups.csv")
write.csv(pCR_log_data_selected,"data for heatmap between pre vs on in pCR groups.csv")





non_pCR_dds<-DESeqDataSetFromMatrix(countData = non_pCR_data,
                                colData=non_pCR_profile_on_pre,
                                design=~timepoint)
non_pCR_dds$sample<-non_pCR_profile_on_pre$treatment
non_pCR_dds$run<-paste0("run",1:nrow(non_pCR_profile_on_pre))
non_pCR_ddsColl<-collapseReplicates(non_pCR_dds,non_pCR_dds$sample,non_pCR_dds$run)
non_pCR_dds<-DESeq(non_pCR_ddsColl)

non_pCR_result<-results(non_pCR_dds)
non_pCR_result<-non_pCR_result[order(non_pCR_result$padj),]
non_pCR_selected_results<-non_pCR_result[complete.cases(non_pCR_result),]
sig_genes_non_pCR_on_pre<-rownames(non_pCR_selected_results[(abs(non_pCR_selected_results$log2FoldChange)>=1)&(non_pCR_selected_results$padj<=0.05),])

non_pCR_log<-rlog(non_pCR_ddsColl)
non_pCR_log_data<-assay(non_pCR_log)
non_pCR_log_data_selected<-non_pCR_log_data[sig_genes_non_pCR_on_pre,]
coll_profile<-colData(non_pCR_ddsColl)
coll_profile<-coll_profile[order(coll_profile$timepoint),]
non_pCR_log_data_selected<-non_pCR_log_data_selected[,rownames(coll_profile)]
ha<-HeatmapAnnotation(df=coll_profile[,3,drop=F])
Heatmap(t(scale(t(non_pCR_log_data_selected))),cluster_columns = F,top_annotation = ha,name="Expression",row_names_gp = gpar(fontsize=1))
write.csv(non_pCR_result,"comparision between pre vs on in non_pCR groups.csv")
write.csv(non_pCR_log_data_selected,"data for heatmap between pre vs on in non_pCR groups.csv")


#make dds object using pcRNA_expression and profile

dds<-DESeqDataSetFromMatrix(countData = pcRNA_expression,
                            colData=profile,
                            design=~timepoint)
dds$sample<-profile$treatment
dds$run<-paste0("run",1:nrow(profile))

ddsColl<-collapseReplicates(dds,dds$sample,dds$run)

coll_profile<-colData(ddsColl)

#prepare heatmap data
heat_data<-rlog(ddsColl)
heatmap_data<-assay(heat_data)
write.csv(heatmap_data,"normalized expression using rlog.csv")
heatmap_data<-read.csv("normalized expression using rlog.csv",header=T,row.names=1)
colnames(heatmap_data)<-sapply(colnames(heatmap_data),function(x) gsub("X","",x))
mads<-sort(apply(heatmap_data,1,mad),decreasing=T)
top_genes_with_most_variation<-names(mads)[1:3000]




#do comparisions
dds_r<-DESeq(ddsColl)
res_on_pre<-results(dds_r,contrast=c("timepoint","ON","PRE"))
res_on_pre<-res_on_pre[order(res_on_pre$padj),]
write.csv(res_on_pre,"comparision results on vs pre.csv")
res_on_pre<-read.csv("comparision results on vs pre.csv",header=T,row.names=1)
res_on_pre<-res_on_pre[complete.cases(res_on_pre),]
sig_genes_on_pre<-rownames(res_on_pre[(abs(res_on_pre$log2FoldChange)>=1)&(res_on_pre$padj<=0.05),])
data<-heatmap_data[sig_genes_on_pre,]
profile_on_pre<-coll_profile[coll_profile$timepoint!="SURG",]
profile_on_pre<-profile_on_pre[order(profile_on_pre$timepoint),]
ha<-HeatmapAnnotation(df=profile_on_pre[,3,drop=F])
Heatmap(t(scale(t(data[,rownames(profile_on_pre)]))),cluster_columns = F,top_annotation = ha,name="Expression",row_names_gp = gpar(fontsize=1))
data_for_heatmap_on_vs_pre<-data[,rownames(profile_on_pre)]
write.csv(data_for_heatmap_on_vs_pre,"data for heatmap of different genes between on and pre.csv")


pCR_profile<-profile_on_pre[which(profile_on_pre$path=="pCR",arr.ind=T),]
non_pCR_profile<-profile_on_pre[which(profile_on_pre$path=="non_pCR",arr.ind=T),]
path_profile<-rbind(pCR_profile,non_pCR_profile)

#comparision of on and pre in pCR group
data<-heatmap_data[,rownames(pCR_profile)]











res_on_surg<-results(dds_r,contrast=c("timepoint","ON","SURG"))
res_on_surg<-res_on_surg[order(res_on_surg$padj),]
write.csv(res_on_surg,"comparision results on vs surg.csv")
res_on_surg<-read.csv("comparision results on vs surg.csv",header=T,row.names=1)
res_on_surg<-res_on_surg[complete.cases(res_on_surg),]
sig_genes_on_surg<-rownames(res_on_surg[(abs(res_on_surg$log2FoldChange)>=1)&(res_on_surg$padj<=0.05),])
data<-heatmap_data[sig_genes_on_surg,]
profile_on_surg<-coll_profile[coll_profile$timepoint!="PRE",]
profile_on_surg<-profile_on_surg[order(profile_on_surg$timepoint),]
data_for_heatmap_on_vs_surg<-data[,rownames(profile_on_surg)]
ha<-HeatmapAnnotation(df=profile_on_surg[,2:3,drop=F])
Heatmap(t(scale(t(data[,rownames(profile_on_surg)]))),cluster_columns = F,top_annotation = ha,name="Expression",row_names_gp = gpar(fontsize=1))
write.csv(data_for_heatmap_on_vs_surg,"data for heatmap of different genes between on and surg.csv")

res_pre_surg<-results(dds_r,contrast=c("timepoint","PRE","SURG"))
res_pre_surg<-res_pre_surg[order(res_pre_surg$padj),]
write.csv(res_pre_surg,"comparision results pre vs surg.csv")
res_pre_surg<-read.csv("comparision results pre vs surg.csv",header=T,row.names=1)
res_pre_surg<-res_pre_surg[complete.cases(res_pre_surg),]
sig_genes_pre_surg<-rownames(res_pre_surg[(abs(res_pre_surg$log2FoldChange)>=1)&(res_pre_surg$padj<=0.05),])
data<-heatmap_data[sig_genes_pre_surg,]
profile_pre_surg<-coll_profile[coll_profile$timepoint!="ON",]
profile_pre_surg<-profile_pre_surg[order(profile_pre_surg$timepoint),]
data_for_heatmap_pre_vs_surg<-data[,rownames(profile_pre_surg)]
write.csv(data_for_heatmap_pre_vs_surg,"data for heatmap of different genes between pre and surg.csv")
ha<-HeatmapAnnotation(df=profile_pre_surg[,2:3,drop=F])
Heatmap(t(scale(t(data[,rownames(profile_pre_surg)]))),cluster_columns = F,top_annotation = ha,name="Expression",row_names_gp = gpar(fontsize=1))
Heatmap(t(scale(t(data[,rownames(profile_pre_surg)]))),cluster_columns = T,top_annotation = ha,name="Expression",row_names_gp = gpar(fontsize=1))
top_genes_with_most_variation<-top_n_vars_genes(heatmap_data,n=1500)

#unsupervised cluster
data<-heatmap_data[top_genes_with_most_variation,]
profile_unsuper<-coll_profile
ha<-HeatmapAnnotation(df=profile_on_pre[,2:3,drop=F])
Heatmap(t(scale(t(data))),cluster_columns = T,top_annotation = ha,name="Expression",row_names_gp = gpar(fontsize=1))


