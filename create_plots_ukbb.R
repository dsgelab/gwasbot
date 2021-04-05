#### Read phenotype file ####
library(data.table)
library(dplyr)
library(jsonlite)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(gdalUtils)
library(BiocManager)
library(TxDb.Hsapiens.UCSC.hg19.knownGene) #BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(bumphunter) #BiocManager::install("bumphunter")
library(org.Hs.eg.db) #BiocManager::install("org.Hs.eg.db")


wrapper <- function(x, ...) 
{
		paste(strwrap(x, width=80), collapse = "\n")
}

h2 <- fread("data/ukb31063_h2_topline.02Oct2019.tsv", sep="\t", header=T)
# This should be the right set of phenotypes to include
h2 <- h2[h2$confidence %in% c("high","medium") & h2$h2_sig %in% c("z4","z7"),]

# file containing genetic correlations
load("data/geno_correlation_sig.Rdata")
# This is instead the traits that we are including
h2 <- h2[gsub("_irnt","",h2$phenotype) %in% geno_corr_df$p1,]

manifest <- fread("data/Pan_UKBB_manifest_2020-09-07.csv", sep=",", header=T)
#colnames(manifest)[1] <- "phenotype"
#colnames(manifest)[6] <- "dropbox"

## Add link to ukbiobank ##
#h2$ukbl <- ifelse(grepl('^([0-9])',h2$phenotype),sapply(strsplit(h2$phenotype,"_"),"[[",1),41202)

## This is later for closest gene annotation
genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

manifest$phenocode2 <- ifelse(manifest$coding!="" & manifest$phenocode != manifest$coding,paste0(manifest$phenocode,"_",manifest$coding),ifelse(manifest$modifier!="",paste0(manifest$phenocode,"_",manifest$modifier),manifest$phenocode))

h2 <- h2[!h2$phenotype %in% c("129_irnt","20015_irnt") & h2$phenotype %in% manifest$phenocode2,]


set.seed(123)
RES <- NULL
for (i in h2$phenotype)
{

    maniget <- manifest %>% filter(phenocode2==i, pheno_sex=="both_sexes")

	system(maniget$wget)
	original_filename <- strsplit(maniget$wget,"/")[[1]][5]
	new_filename <- paste0(remove_file_extension(original_filename),".gz")
	system(paste0("mv ",original_filename," ",new_filename))

    rowsnames <- strsplit(system(paste0("zcat < ",new_filename, " | head -1"),intern = TRUE),"\t")[[1]]
    loconfname <- rowsnames[grepl("low_confidence",rowsnames)]
    loconfname <- ifelse("low_confidence_EUR" %in% loconfname,"low_confidence_EUR",loconfname[1])
    pvalname <- rowsnames[grepl("pval_",rowsnames)][1]

	df <- fread(cmd=paste0("zcat < ",new_filename), header=T, sep="\t", select=c('chr','pos','ref','alt',loconfname, pvalname)) 

    colnames(df)[grepl(loconfname,colnames(df))] <- "low_confidence_EUR"
    colnames(df)[grepl(pvalname,colnames(df))] <- "pval_meta"
        
    df <- df %>% filter(is.finite(pval_meta), pval_meta > 0, low_confidence_EUR==FALSE)  %>%
        mutate(chr=gsub('X', '23', chr), pos=as.integer(pos), pval_t=-log10(pval_meta),indic=if_else(pval_meta < 5e-8,1,2), odd=as.numeric(chr) %% 2) %>% mutate(chromnum=as.numeric(chr))


    posmin <- tapply(df$pos,df$chromnum, min)
    posmax <- tapply(df$pos,df$chromnum, max)
    posshift <- head(c(0,cumsum(as.numeric(posmax))),-1)
    names(posshift) <- names(posmin)

    for (k in unique(df$chr))
    {
        df$pos_new[df$chr==k] <-  df$pos[df$chr==k] + posshift[names(posshift) == k]
    }

    dfmsplit <- split(df, df$chr)
    xbreaks <- sapply(dfmsplit,function(x) x$pos_new[length(x$pos_new)/2])

    df_manhattan <- df %>%
       filter(pval_meta <= 0.01)

    ymax <- ifelse(max(df_manhattan$pval_t) < -log10(5e-8), -log10(5e-8), max(df_manhattan$pval_t)) + 0.2

    s1 <- wrapper(maniget$description)
    case_or_not <- ifelse(is.na(maniget$n_controls_EUR),paste0(s1,"\nN=",maniget$n_cases_full_cohort_both_sexes,"\n",maniget$pops),paste0(s1,"\nN. cases=",maniget$n_cases_full_cohort_both_sexes,"; N. controls=",sum(maniget$n_controls_AFR,maniget$n_controls_AMR,maniget$n_controls_CSA,maniget$n_controls_EAS,maniget$n_controls_EUR,maniget$n_controls_MID,na.rm=T),"\n",maniget$pops))
  
    # Closest gene
    labeldf <- NULL
    k <- 0
    df_manhattan$variant <- paste0(df_manhattan$chr,df_manhattan$pos,df_manhattan$ref,df_manhattan$alt)
    dftemp <- df_manhattan
    while (k < 5 & (min(dftemp$pval_meta) < 0.00000005))
    {
      indmin <- which(df_manhattan$variant==dftemp$variant[which.min(dftemp$pval_meta)])
      posmin <- dftemp$pos[which.min(dftemp$pval_meta)]
      chromin <- dftemp$chr[which.min(dftemp$pval_meta)]
      closestgene <- matchGenes(makeGRangesFromDataFrame(data.frame(chr=paste0("chr",chromin),start=posmin,end=posmin+1)),genes, type="any")$name
      dftemp <- dftemp[!(dftemp$chr==chromin & dftemp$pos>(posmin-1000000) & dftemp$pos<(posmin+1000000)),]
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

    ggsave(paste0("data/manhattan_UKBB_trans/",i,"_trans_MF.png"), width = 12, height = 6, dpi = 200)

    system(paste0("rm ",new_filename))

    #RES <- rbind(RES,cbind(h2[h2$pheno==i,c("description","n","n_controls","n_cases","ukbl")],maniget$dropbox,paste0("plot_ukbb/",i,"_MF.png")))

    print(which(i==h2$phenotype))
}

#df <- data.frame(RES)
#colnames(df) <- c("description","n","n_controls","n_cases","ukbl","dropbox","file")

#write(toJSON(df, pretty=TRUE),file="images_ukbb.js")


#gcloud compute scp manifest_201807.tsv anachiner:/home/aganna/

#gcloud compute scp anachiner:/home/aganna/*.png /Users/andreaganna/Documents/Work/Post_doc/twitter_bot/manhattan/