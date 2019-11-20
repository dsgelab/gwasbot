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

h2 <- fread("ukb31063_h2_topline.02Oct2019.tsv", sep="\t", header=T)
# This should be the right set of phenotypes to include
#h2 <- h2[h2$confidence %in% c("high","medium") & h2$h2_sig %in% c("z4","z7"),]

# file containing genetic correlations
load("geno_correlation_sig.Rdata")
# This is instead the traits that we are including
h2 <- h2[gsub("_irnt","",h2$phenotype) %in% geno_corr_df$p1,]

manifest <- fread("Manifest_201807.csv", sep=",", header=T)
colnames(manifest)[1] <- "phenotype"
colnames(manifest)[6] <- "dropbox"

## Add link to ukbiobank ##
h2$ukbl <- ifelse(grepl('^([0-9])',h2$phenotype),sapply(strsplit(h2$phenotype,"_"),"[[",1),41202)

## This is later for closest gene annotation
genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

set.seed(123)
RES <- NULL
for (i in h2$phenotype)
{

	maniget <- manifest %>% filter(phenotype==i, Sex=="both_sexes")

	system(maniget$dropbox)
	original_filename <- strsplit(maniget$dropbox," ")[[1]][4]
	new_filename <- paste0(remove_file_extension(original_filename),".gz")
	system(paste0("mv ",original_filename," ",new_filename))

	
	df <- fread(cmd=paste0("gzcat ",new_filename), header=T, sep="\t", select=c('variant', 'low_confidence_variant', 'pval')) %>%
        filter(is.finite(pval), pval > 0, low_confidence_variant==FALSE) %>%
        separate(variant, c('chrom', 'pos', 'ref', 'alt'), sep=':', remove=FALSE) %>%
        mutate(chrom=gsub('X', '23', chrom), pos=as.integer(pos), pval_t=-log10(pval),indic=if_else(pval < 5e-8,1,2), odd=as.numeric(chrom) %% 2) %>%
        mutate(chromnum=as.numeric(chrom))

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

    s1 <- wrapper(h2$description[h2$pheno==i])
    case_or_not <- ifelse(is.na(h2$n_cases[h2$pheno==i]),paste0(s1,"\nN=",h2$n[h2$pheno==i]),paste0(s1,"\nN. cases=",h2$n_cases[h2$pheno==i],"; N. controls=",h2$n_controls[h2$pheno==i]))
  
    # Closest gene
    labeldf <- NULL
    k <- 0
    dftemp <- df_manhattan
    while (k < 3 & (min(dftemp$pval) < 0.00000005))
    {
      indmin <- which(df_manhattan$variant==dftemp$variant[which.min(dftemp$pval)])
      posmin <- dftemp$pos[which.min(dftemp$pval)]
      chromin <- dftemp$chrom[which.min(dftemp$pval)]
      closestgene <- matchGenes(makeGRangesFromDataFrame(data.frame(chr=paste0("chr",chromin),start=posmin,end=posmin+1,strand="+")),genes, type="any")$name
      dftemp <- dftemp[!(dftemp$chrom==chromin & dftemp$pos>(posmin-1000000) & dftemp$pos<(posmin+1000000)),]
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

    ggsave(paste0(i,"_MF.png"), width = 12, height = 6, dpi = 200)

    system(paste0("rm ",new_filename))

    RES <- rbind(RES,cbind(h2[h2$pheno==i,c("description","n","n_controls","n_cases","ukbl")],maniget$dropbox,paste0(i,"_MF.png")))

    print(which(i==h2$phenotype))
}

df <- data.frame(RES)
colnames(df) <- c("description","n","n_controls","n_cases","ukbl","dropbox","file")

write(toJSON(df, pretty=TRUE),file="images_ukbb.js")


#gcloud compute scp manifest_201807.tsv anachiner:/home/aganna/

#gcloud compute scp anachiner:/home/aganna/*.png /Users/andreaganna/Documents/Work/Post_doc/twitter_bot/manhattan/