library(Packzhu)

# read total expression data, and clean the expression
expression<-get_data(header=T,row.names=1,stringsAsFactors=F,sep="")
ensembl_id<-sapply(expression$gene_id,function(x) unname(unlist(strsplit(x,"\\."))[1]))
expression$gene_id<-unname(ensembl_id)
expression<-aggregate(.~gene_id+gene_name,data=expression,sum)
rownames(expression)<-expression$gene_id
expression<-expression[,-1]
expression$gene_name<-sapply(expression$gene_name,function(x) trim(x))
colnames(expression)<-sapply(colnames(expression),function(x) gsub("X","",x))
total_expression<-expression[,-1]


extract<-function(df,col,pat){
  col<-deparse(substitute(col))
  selected<-sapply(df[,col],function(x) grepl(pat,x,ignore.case = T))
  return(df[selected,])
}

longestcommonsubstring<-function(str1,str2){
  str1_chars<-unlist(strsplit(str1,""))
  pre_l_substring=0
  substring=""
  for (i in 1:length(str1_chars)){
    l_substring=1
    pat<-str1_chars[i]
    
    if (grepl(pat,str2)){
      start=i+1
      for(j in start:length(str1_chars)){
       
        pat<-paste(pat,str1_chars[j],sep="")
        
        if (grepl(pat,str2)){
          l_substring=l_substring+1
        }
        else{
          break()
        }
        
      }
    }
    if(l_substring>=pre_l_substring){
      pre_l_substring<-l_substring
      substring=pat
    }
  }
  return(substring)
}

norm_gene_name<-function(expression,col,microRNA_list){
  col<-deparse(substitute(col))
  for(i in seq_len(nrow(expression))){
    if(!(expression[,col][i] %in% microRNA_list$Approved_Symbol)){
      if(any(grepl(expression[,col][i],microRNA_list$Synonyms))){
        expression[,col][i]<-microRNA_list$Approved_Symbol[grepl(expression[,col][i],microRNA_list$Synonyms)]
      }
    else{
      
      com_string_l<-sapply(microRNA_list$Approved_Name,function(x) nchar(longestcommonsubstring(expression[,col][i],x)))
      ind<-which.max(com_string_l)
      expression[,col][i]<-microRNA_list$Approved_Symbol[ind]
      
    }
    }
  }
  return(expression)
}
##normalize microRNA gene name
microRNA_list<-get_data(header=T,row.names=1,stringsAsFactors=F,sep="\t")
colnames(microRNA_list)<-sapply(colnames(microRNA_list),function(x) gsub("\\.","_",x))
microRNA_expression<-extract(expression,gene_name,"MIR")
microRNA_expression<-norm_gene_name(microRNA_expression,gene_name,microRNA_list)
microRNA_expression<-aggregate(.~gene_name,data=microRNA_expression,sum)
rownames(microRNA_expression)<-microRNA_expression$gene_name
microRNA_expression<-microRNA_expression[,-1]
write.csv(microRNA_expression,"microRNA_expression.csv")


##prepare profile data
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

total_reading_dds<-DESeqDataSetFromMatrix(countData = total_expression,
                                                        colData=profile,
                                                        design=~timepoint)
SFs<-estimateSizeFactors(total_reading_dds)

microRNA_dds<-DESeqDataSetFromMatrix(countData = microRNA_expression,
                                     colData=profile,
                                     design=~timepoint)
sizeFactors(microRNA_dds)<-sizeFactors(SFs)
microRNA_dds <- estimateDispersions(microRNA_dds)
microRNA_dds<-nbinomWaldTest(microRNA_dds)

microRNA_dds$sample<-profile$treatment
microRNA_dds$run<-paste0("run",1:nrow(profile))
microRNA_ddsColl<-collapseReplicates(microRNA_dds,microRNA_dds$sample,microRNA_dds$run)
microRNA_ddsColl <- microRNA_ddsColl[rowSums(counts(microRNA_ddsColl)) > 1, ]
microRNA_dds<-DESeq(microRNA_ddsColl)
microRNA_rlog<-varianceStabilizingTransformation(microRNA_dds,blind=F)
microRNA_log_data<-assay(microRNA_rlog)
write.csv(microRNA_log_data,"microRNA_rlog_data.csv")
write.csv(profile_coll,"collase replicated profile.csv")

data<-n_top_genes_with_most_mad(microRNA_log_data,n=300)
profile_coll<-colData(microRNA_dds)

copyhalist<-function(sourcehalist){
  targetlist<-list()
  for (listname in slotNames(sourcehalist)){
    targetlist[[listname]]<-slot(sourcehalist,listname)
  }
  return(targetlist)
}
hAnnotation<-function(df){
  ha<-copyhalist(HeatmapAnnotation(df))
  
  ha$profile<-df
  return(ha)
}

ha<-HeatmapAnnotation(df=profile_coll[,c(3,6),drop=F])
ha<-hAnnotation(df=profile_coll[,c(3,6),drop=F])
Heatmap(t(scale(t(data))),cluster_columns = T,top_annotation = ha,name="Expression",row_names_gp = gpar(fontsize=1))
for(an in colnames(profile_coll[,c(3,6),drop=F])) {
  decorate_annotation(an, {
    grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
  })
}
