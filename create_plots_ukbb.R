#### Read phenotype file ####
library(data.table)
library(dplyr)
library(jsonlite)
library(tidyr)
library(ggplot2)


system("gsutil cp gs://ukb-mega-gwas-results/round2/additive-tsvs/phenotypes.both_sexes.tsv.bgz .")
system("zcat phenotypes.both_sexes.tsv.bgz > phenotypes.both_sexes.tsv")


wrapper <- function(x, ...) 
{
		paste(strwrap(x, width=80), collapse = "\n")
}


pheno <- fread("phenotypes.both_sexes.tsv", sep="\t", header=T)
pheno <- pheno[pheno$n_cases>500 | is.na(pheno$n_cases),]
manifest <- fread("manifest_201807.tsv", sep="\t", header=T)
colnames(manifest)[1] <- "phenotype"
colnames(manifest)[6] <- "dropbox"

## Add link to ukbiobank ##
pheno$ukbl <- ifelse(grepl('^([0-9])',pheno$phenotype),sapply(strsplit(pheno$phenotype,"_"),"[[",1),41202)


set.seed(123)
RES <- NULL
for (i in pheno$phenotype)
{

	maniget <- manifest %>% filter(phenotype==i, Sex=="both_sexes")

	system(paste0("gsutil cp gs://ukb-mega-gwas-results/round2/additive-tsvs/",i,".gwas.imputed_v3.both_sexes.tsv.bgz ."))

	df <- fread(input=paste0("zcat ",i,".gwas.imputed_v3.both_sexes.tsv.bgz"), header=T, sep="\t", select=c('variant', 'low_confidence_variant', 'pval')) %>%
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

    s1 <- wrapper(pheno$description[pheno$pheno==i])
    case_or_not <- ifelse(is.na(pheno$n_cases[pheno$pheno==i]),paste0(s1,"\nN=",pheno$n_non_missing[pheno$pheno==i]),paste0(s1,"\nN. cases=",pheno$n_cases[pheno$pheno==i],"; N. controls=",pheno$n_controls[pheno$pheno==i]))

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
    	ggtitle(case_or_not)


    ggsave(paste0(i,"_MF.png"), width = 12, height = 6, dpi = 200)

    system(paste0("rm ",i,".gwas.imputed_v3.both_sexes.tsv.bgz"))

    RES <- rbind(RES,cbind(pheno[pheno$pheno==i,c("description","n_non_missing","n_controls","n_cases","ukbl")],maniget$dropbox,paste0(i,"_MF.png")))

    print(which(i==pheno$phenotype))
}

df <- data.frame(RES)
colnames(df) <- c("description","n_non_missing","n_controls","n_cases","ukbl","dropbox","file")

write(toJSON(df, pretty=TRUE),file="images_ukbb.js")



#gcloud compute scp manifest_201807.tsv anachiner:/home/aganna/

#gcloud compute scp anachiner:/home/aganna/*.png /Users/andreaganna/Documents/Work/Post_doc/twitter_bot/manhattan/