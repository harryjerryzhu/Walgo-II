##normalize long non coding gene name
lnc_list<-get_data(header=T,row.names=1,stringsAsFactors=F,sep="\t")
colnames(lnc_list)<-sapply(colnames(lnc_list),function(x) gsub("\\.","_",x))
lnc_expression<-selectbycolumnin(expression,gene_name,lnc_list$Approved_Symbol)
lnc_expression<-aggregate(.~gene_name,data=lnc_expression,sum)
rownames(lnc_expression)<-lnc_expression$gene_name
lnc_expression<-lnc_expression[,-1]
write.csv(lnc_expression,"lnc_expression.csv")

lnc_dds<-DESeqDataSetFromMatrix(countData = lnc_expression,
                                     colData=profile,
                                     design=~timepoint)
sizeFactors(lnc_dds)<-sizeFactors(SFs)
lnc_dds <- estimateDispersions(lnc_dds)
lnc_dds<-nbinomWaldTest(lnc_dds)

lnc_dds$sample<-profile$treatment
lnc_dds$run<-paste0("run",1:nrow(profile))
lnc_ddsColl<-collapseReplicates(lnc_dds,lnc_dds$sample,lnc_dds$run)
lnc_ddsColl <- lnc_ddsColl[rowSums(counts(lnc_ddsColl)) > 1, ]
lnc_dds<-DESeq(lnc_ddsColl)
lnc_rlog<-varianceStabilizingTransformation(lnc_dds,blind=F)
lnc_log_data<-assay(lnc_rlog)
write.csv(lnc_log_data,"lnc_rlog_data.csv")


data<-n_top_genes_with_most_mad(lnc_log_data,n=100)
profile_coll<-colData(lnc_dds)

ha<-HeatmapAnnotation(df=profile_coll[,c(3,6),drop=F])
Heatmap(t(scale(t(data))),cluster_columns = T,top_annotation = ha,name="Expression",row_names_gp = gpar(fontsize=1))
