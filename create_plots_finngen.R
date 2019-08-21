#### Read phenotype file ####
library(data.table)
library(dplyr)
library(jsonlite)
library(tidyr)
library(ggplot2)


info <- read_json("R2_pheno.json")

pheno <- data.frame(phenotype=unlist(lapply(info, `[`, "phenocode")),
                 phenostring=unlist(lapply(info, `[`, "phenostring")),
                 ncases=unlist(lapply(info, `[`, "num_cases")),
                 ncontrol=unlist(lapply(info, `[`, "num_controls")),
                 gw_sign=unlist(lapply(info, `[`, "num_gw_significant")))

# Keep only endpoints with GW-significants hits
pheno_sign <- pheno[pheno$gw_sign > 0,]

wrapper <- function(x, ...) 
{
		paste(strwrap(x, width=80), collapse = "\n")
}



set.seed(123)
RES <- NULL
for (i in pheno_sign$phenotype)
{
  system(paste0("gsutil cp gs://finngen-production-library-green/R2/summary_stats/release/",i,".gz ."))
	      
	df <- fread(input=paste0("gunzip -c ",i,".gz"), header=T, sep="\t", select=c('#chrom', 'pos', 'pval')) %>%
	      rename("chrom" = "#chrom") %>%
        filter(is.finite(pval), pval > 0) %>%
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

    s1 <- wrapper(pheno_sign$phenostring[pheno_sign$phenotype==i])
    case_or_not <- paste0(s1,"\nN. cases=",pheno_sign$ncases[pheno_sign$phenotype==i],"; N. controls=",pheno_sign$ncontrol[pheno_sign$phenotype==i])

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

    system(paste0("rm ",i,".gz"))

    RES <- rbind(RES,cbind(pheno_sign[pheno_sign$phenotype==i,c("phenostring","ncontrol","ncases")],paste0("gs://finngen-production-library-green/R2/summary_stats/release/",i,".gz"), paste0(i,"_MF.png")))

    print(which(i==pheno_sign$phenotype))
}

df <- data.frame(RES)
colnames(df) <- c("description","n_cases","n_controls","link","file")

write(toJSON(df, pretty=TRUE),file="images_finngen.js")



#gcloud compute scp manifest_201807.tsv anachiner:/home/aganna/

#gcloud compute scp anachiner:/home/aganna/*.png /Users/andreaganna/Documents/Work/Post_doc/twitter_bot/manhattan/