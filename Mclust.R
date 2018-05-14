#Select the folder where your input file is located
directory = "G:/BIOINFORMATICA_Doutorado/Proteômica/Experimentos/EXP11_MAS_CentralAnalitica_20042018/Estatística_exp11"
setwd(directory)

#Load the file that contains the dataset in a variable
filename = "Matrix_ANOVA1way_pvalue_631prot_mediana.txt"
proteins = read.csv(filename, sep = "\t", header=TRUE) #now.names = 2?
selected_proteins = proteins[3:633, 1:10]
id_proteins = proteins[3:633, 29]
ids = matrix(, ncol=1, nrow=631)
j = 1
for(i in id_proteins){
  id_proteins_aux = strsplit(i, ";")
  unlist(id_proteins_aux)
  ids[j] = id_proteins_aux[[1]][1]
  j = j + 1
}

library("mclust")
mod1 = Mclust(selected_proteins, G = 2:30)
summary(mod1)
result = matrix(ids, ncol=2, nrow=631)
for(k in 1:631){
  result[k,2] = mod1$classification[k]
}
write.csv(result, "testest3.csv")
