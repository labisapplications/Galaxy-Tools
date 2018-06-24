## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("pathview")
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("gage")

#Import gage library
library(gage)

#Define the file directory using "/"
Directory <- "F:/BIOINFORMATICA_Doutorado/Proteômica/Pathview/"

#Set the system directory as file directory
setwd(Directory)

#Read file inside the new system directory and save inside the Gene Set "gs" object_
#_inside demo object
demo.gs=readList("kegg_pathway.gmt")

#Show DataFrame "gs"
demo.gs

#Set the system directory ti the file folder, read the file into a data frame (df),
#reaname the firt colum to ID.
setwd("F:/BIOINFORMATICA_Doutorado/Proteômica/Experimentos/EXP11e12_ok/Estatistica11e12/Novo_ScriptR/R_EXP11_15jun")
df_exp11 <-read.delim("Proteins For Anova.txt",  
                       header= TRUE, 
                       na.strings = 'NaN')
names(df_exp11)[1] <- "ID"

setwd("F:/BIOINFORMATICA_Doutorado/Proteômica/Experimentos/EXP11e12_ok/Estatistica11e12/Novo_ScriptR/R_EXP12_15jun")
df_exp12_N <-read.delim("Control_Fixed_EXP12_sem5e7.txt",  
                      header= TRUE, 
                      na.strings = 'NaN')
names(df_exp12_N)[1] <- "ID"

#Merge two data frames (datasets) by ID - EXP-N and EXP+N
df_total_exp11and12 <- merge(df_exp11, df_exp12_N, by="ID", all = TRUE)

#Set the row names of df with values of first column
rownames(df_total_exp11and12) <- df_total_exp11and12[,1]

#Exclude the first column
df_total_exp11and12 <- df_total_exp11and12[,-1]

#Get the column names as a vector
cn=colnames(df_total_exp11and12)

#Find the column index of group name
control=grep('LFQ.intensity.12_TP',cn, ignore.case =T)
treat = grep('LFQ.intensity.11_TP',cn, ignore.case =T)
#Test with one Time Point (TP)
#control=grep('LFQ.intensity.12_TP6',cn, ignore.case =T)
#treat = grep('LFQ.intensity.11_TP6',cn, ignore.case =T)

#Detecting metabolic pathways that change - usign GAGE
#Unsuccessfully
teste_gage <- gage(df_total_exp11and12, 
                   use.stouffer=TRUE, 
                   rank.test = T, 
                   compare="unpaired", 
                   same.dir = T, 
                   gsets = demo.gs, 
                   ref= control, 
                   samp=treat,
                   FDR.adj = T)

