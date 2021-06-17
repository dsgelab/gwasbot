#### Read phenotype file ####
library(data.table)
library(dplyr)
library(jsonlite)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(BiocManager)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) #BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(bumphunter) #BiocManager::install("bumphunter")
library(org.Hs.eg.db) #BiocManager::install("org.Hs.eg.db")


pheno <- fread("data/finngen_R5_pheno_n.tsv")

wrapper <- function(x, ...) 
{
		paste(strwrap(x, width=80), collapse = "\n")
}


## This is later for closest gene annotation
genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg38.knownGene)




set.seed(123)
RES <- NULL
for (i in pheno$phenocode[pheno$num_gw_significant>0])
{
    system(paste0("gsutil cp gs://finngen-public-data-r5/summary_stats/finngen_R5_" ,i,".gz ."))
  
	df1 <- fread(cmd=paste0("gunzip -c finngen_R5_",i,".gz"), header=T, sep="\t", select=c('#chrom', 'pos', 'pval')) 
	
	if(min(df1$pval,na.rm=T)<0.00000005)
	{
	  df <- df1 %>% rename("#chrom" = "chrom") %>%
	    filter(is.finite(pval), pval > 0) %>%
	    mutate(chrom=gsub('X', '23', chrom), pos=as.integer(pos), pval_t=-log10(pval),indic=if_else(pval < 5e-8,1,2), odd=as.numeric(chrom) %% 2) %>% mutate(chromnum=as.numeric(chrom))
	  
	  df$variant <- as.character(1:nrow(df))
	  
	  posmin <- tapply(df$pos,df$chromnum, min)
	  posmax <- tapply(df$pos,df$chromnum, max)
	  posshift <- head(c(0,cumsum(as.numeric(posmax))),-1)
	  names(posshift) <- names(posmin)
	  
	  for (k in unique(df$chrom))
	  {
	    df$pos_new[df$chrom==k] <-  df$pos[df$chrom==k] + posshift[names(posshift) == k]
	  }
	  
	  dfmsplit <- split(df, df$chrom)
	  xbreaks <- sapply(dfmsplit,function(x) x$pos_new[length(x$pos_new)/2])
	  
	  df_manhattan <- df %>%
	    filter(pval <= 0.01)
	  
	  ymax <- ifelse(max(df_manhattan$pval_t) < -log10(5e-8), -log10(5e-8), max(df_manhattan$pval_t)) + 0.2
	  
	  s1 <- wrapper(pheno$name[pheno$phenocode==i])
	  case_or_not <- paste0(s1,"\nN. cases=",pheno$num_cases[pheno$phenocode==i],"; N. controls=",pheno$num_controls[pheno$phenocode==i])
	  
	  
	  # Closest gene
	  labeldf <- NULL
	  k <- 0
	  dftemp <- df_manhattan
	  while (k < 5 & (min(dftemp$pval) < 0.00000005))
	  {
	    indmin <- which(df_manhattan$variant==dftemp$variant[which.min(dftemp$pval)])
	    posmin <- dftemp$pos[which.min(dftemp$pval)]
	    chromin <- dftemp$chrom[which.min(dftemp$pval)]
	    closestgene <- matchGenes(makeGRangesFromDataFrame(data.frame(chr=paste0("chr",chromin),start=posmin,end=posmin+1)),genes, type="any")$name
	    dftemp <- dftemp[!(dftemp$chrom==chromin & dftemp$pos>(posmin-2000000) & dftemp$pos<(posmin+2000000)),]
	    k <- k + 1
	    
	    labeldf <- rbind(labeldf,c(indmin,closestgene))
	  }
	  
	  
	  df_manhattan$label <- ""
	  df_manhattan$label[as.numeric(labeldf[,1])] <- labeldf[,2]
	  
	  
	  outp <- ggplot(df_manhattan, aes(x = pos_new,y = pval_t)) +
	    geom_point(aes(size=as.factor(indic),
	                   colour=as.factor(odd))) +
	    scale_x_continuous(breaks = xbreaks, labels = names(xbreaks), expand = c(0,0)) +
	    scale_y_continuous(expand = c(0,0), limits = c(-log10(0.01),ymax), breaks=c(seq(2,round(ymax),length.out=7)),labels=c(round(seq(2,round(ymax),length.out=7),0)))+
	    expand_limits(x = 23.3) +
	    guides(colour = FALSE,alpha=FALSE, size=FALSE, fill=FALSE) +
	    labs(x = "chromosome", y = expression(-log[10](italic(P)))) + 
	    theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "red"), plot.title = element_text(hjust = 0.5,face ="bold")) +
	    geom_hline(aes(yintercept= -log10(0.00000005)),colour = "red", lwd=0.6, linetype = 5)  + 
	    scale_colour_manual(values = c("darkorange2","darkmagenta")) + 
	    scale_size_manual(values=c(1,0.4)) +
	    ggtitle(case_or_not) + geom_text_repel(aes(label=label))
	  
	   ggsave(paste0("data/manhattan_FINNGEN/",i,"_MF.png"), width = 12, height = 6, dpi = 200)
	  
	   RES <- rbind(RES,cbind(pheno[pheno$phenocode==i,c("phenocode","name","num_cases","num_controls")],pheno$path_bucket[pheno$phenocode==i], paste0(i,"_MF.png"), paste0("http://r5.finngen.fi/pheno/",i),paste0("https://r5.risteys.finngen.fi/phenocode/",i)))
	}

    system(paste0("rm finngen_R5_",i,".gz"))

    print(which(i==pheno$phenocode))
}

df <- data.frame(RES)
colnames(df) <- c("phenocode","description","num_cases","num_controls","file","pheweb_link","risteys_link")

write(toJSON(df, pretty=TRUE),file="data/images_finngen.js")
