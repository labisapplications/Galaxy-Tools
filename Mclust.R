###
# title: ""
# author: "Diego M. Riano-Pachon, Rodrigo R. Dorado Goitia"
# date: ""
# output: pdf
###

###
# Description
###

rm(list=ls())

#Load required packages
library(qvalue)
library(ggplot2)
library(mclust)
library(reshape2)
###

# Set Functions

# Set The directory to start the execution and read a csv file
# directoryName String Name of the directory.
# filename String Name of the csv file.
# Rodrigo Dorado
readCsvDirectory<-function(directoryName, filename) {
  setwd(directoryName)
  fileData <- read.csv(filename,
                       sep = "\t", 
                       header=TRUE, 
                       na.strings = 'NaN',
                       dec='.',
                       stringsAsFactors=FALSE)
  return(fileData[-1:-2,])
}

# Print the dim, columns name and the first 5 rows of a variable
# dataToInform  Object Name of the varaible.
# Rodrigo Dorado
printInfo <- function(dataToInform) {
  print('Length')
  print(length(dataToInform))
  print('---------------')
  print('Dim')
  print(dim(dataToInform))
  print('***************')
  print('Column Names')
  print(colnames(dataToInform))
  print('_______________')
  print('Firsts rows')
  print(head(dataToInform))
  print('_______________')
}

## Put the id of protein as row name, get only the data of the experiments
## proteins matrix All the proteins data
## column_id int The number of the column where the is the id
## init_column int The number of the colmun of the initial experiment
## number_exp int The number of experiments
## Rodrigo Dorado 
getDataofInfo <- function(proteins, column_id, init_column, number_exp) {
  selected_columns <- proteins[1:nrow(proteins), c(init_column:(number_exp * 3))]
  id_proteins_column <- proteins[1:nrow(proteins), column_id]
  ids <- matrix(, ncol = 1, nrow = nrow(proteins))
  j <- 1
  for(i in id_proteins_column){
    id_proteins_aux <- strsplit(i, ";")
    unlist(id_proteins_aux)
    ids[j] <- id_proteins_aux[[1]][1]
    j <- j + 1
  }
  row.names(selected_columns) <- ids
  return(selected_columns)
}

## Verify if the protein contain at least one group with all the values different to NaN
## proteins matrix All the proteins data
## number_exp int The number of experiments
## num_div int The number of divisions per experiment
## Rodrigo Dorado
VarifyAComplteGroup <- function(proteins, number_exp, num_div) {
  caso1 <- 0
  caso2 <- 0
  caso3 <- 0
  caso4 <- 0
  protein_erase <- c()
  col_names <- colnames(proteins)
  ## Get every protein
  count <- 0
  for (protein in 1:nrow(proteins)) {
    ok_protein <- FALSE
    count_p <- 0
    ## Move by a experiment
    for (result_exp in seq(1, (number_exp * num_div), by = num_div)) {
      ok_group <- TRUE
      pos_NaN <- c()
      pos_valid <- c()
      ## Get the name of the experiment
      column_name <- getNameColumn(col_names, result_exp)
      ## Move into the divisions
      for (div_pos in 1:num_div) {
        pos <- result_exp + (div_pos - 1)
        ## check if is NaN
        if(is.na(proteins[protein, pos])){
          ok_group <- FALSE
          pos_NaN <- c(pos_NaN, pos)
        } else {
          pos_valid <- c(pos_valid, pos)
        }
      }
      ##Check if in the group, if the three are valid
      if (ok_group) {
        ok_protein <- TRUE
      }
      switch_var <- length(pos_NaN)
      allPos <- c(pos_valid, pos_NaN)
      ## Execute conditions
      if (switch_var == 0) {
        resultValue <- getTheMedian(proteins, protein, allPos)
        caso1 <- caso1 + 1
      }
      if (switch_var == 1) {
        proteins[protein, pos_NaN] <- mean(as.numeric(proteins[protein, pos_valid]))
        resultValue <- getTheMedian(proteins, protein, allPos)
        caso2 <- caso2 + 1
      }
      if (switch_var == 2) {
        proteins[protein, pos_NaN] <- 0
        ## resultValue <- 0
        resultValue <- proteins[protein, pos_valid]
        caso3 <- caso3 + 1
        count_p <- count_p + 1 
      }
      if (switch_var == 3) {
        proteins[protein, pos_NaN] <- 0
        resultValue <- 0
        caso4 <- caso4 + 1
      }
      ## Get the unique value
      proteins[protein, column_name] <- resultValue 
    }
    if(count_p > 0) {
      count <- count + 1
    }
    ## Check if the protein conatin at least one group valid
    if (!ok_protein) {
      protein_erase <- c(protein_erase, protein)
    }
  }
  ## return the list
  print(caso1)
  print(caso2)
  print(caso3)
  print(caso4)
  print(count)
  if(length(protein_erase) > 0) {
    return(proteins[-protein_erase,])
  }else{
    return(proteins)
  }
}

## Get the column Name of the group
## col_names array All the column names
## pos int The number of column to change
## Rodrigo Dorado
getNameColumn <- function(col_names, pos) {
  return(getNameOfEspColumn(col_names[pos]))
}

## Split the name of the group
## col_names String The Id of the protein
## Rodrigo Dorado
getNameOfEspColumn <- function(col_names) {
  return(strsplit(strsplit(col_names, "[.]")[[1]][3], "[_]")[[1]][1])
}

## Get the median of a group
## proteins matrix All the proteins
## protein int The actual protein
## positions int The number of experiment
## Rodrigo Dorado
getTheMedian <- function(proteins, protein, positions) {
  return(median(as.numeric(proteins[protein, positions])))
}

## Merge two data frames by the id
## proteins Data Frame The first data frame
## p_values Data Frame The second data frame
mergeById <- function(proteins, p_values) {
  proteins$auxiliar <- row.names(proteins)
  p_values$auxiliar <- row.names(p_values)
  proteins <- merge(proteins, p_values, by.x = 'auxiliar', by.y = 'auxiliar')
  row.names(proteins) <- proteins[, 'auxiliar']
  return(proteins[, -1])
}

## Get the p values of the ANOVA
## proteins Data Frame All the data of the experiments
## Rodrigo Dorado
calculateANOVA <- function(proteins) {
  all_p_values <- c()
  for (i in 1:nrow(proteins)) { #nrow(proteins)
    groupNames <- c()
    for (j in colnames(proteins)) {
      groupNames <- c(groupNames, getNameOfEspColumn(j)) 
    }
    MToANOVA <- data.frame(groupNames, as.numeric(t(proteins[i, ])))
    colnames(MToANOVA) <- c('GROUP', 'DATA')
    aov1 <- aov(DATA ~ GROUP, data = MToANOVA)
    all_p_values <- c(all_p_values, summary(aov1)[[1]][["Pr(>F)"]][1])
  }
  p_values <- data.frame(all_p_values)
  row.names(p_values) <- row.names(proteins)
  return(p_values)
}

## Select the data that has less than 0.05 the p_value; Replace NA by 0
## selected_columns_proteins Matrix the selected columns of the data
## Rodrigo Dorado
selectSignificanProteinsExpression <- function(selected_columns_proteins) {
  protein_Qvalues <- qvalue(selected_columns_proteins$all_p_values, lambda=0)
  selected_columns_proteins$QValue <- protein_Qvalues$qvalues
  selectedSignificantExpr <- selected_columns_proteins[which(selected_columns_proteins$QValue < 0.05), c(1:10)]
  all_row_names <- row.names(selectedSignificantExpr)
  selectedSignificantExpr <- as.data.frame(sapply(selectedSignificantExpr, as.numeric))
  row.names(selectedSignificantExpr) <- all_row_names
  ##selectedSignificantExpr[is.na(selectedSignificantExpr)] <- 0
  return(selectedSignificantExpr)
}

## Get Z Values of the proteins Expression
## selectedSignificantExpr Matrix The values of the selected proteins
## Rodrigo Dorado
getZValueSgnificantProteinExpression <- function(selectedSignificantExpr) {
  return(as.data.frame(t(scale(t(selectedSignificantExpr)))))
}

## Execute cluster function, save the result in a csv file. Note: add the condition of max clusters
## ZValues Matrix The z values to cluster
## numOfClusters int The max number of clusters
## filename String The name of the csv fil
## Rodrigo Dorado
executeMClust <- function(ZValues, filename, numOfClusters = 60) {
  print("Clustering")
  print(numOfClusters)
  clusterResult <- Mclust(ZValues, G = 2:numOfClusters)
  print(summary(clusterResult))
  Num_clusters <- max(clusterResult$classification, na.rm=TRUE)
  if (Num_clusters >= numOfClusters) {
    newNumOfClusters <- numOfClusters + 30
    executeMClust(ZValues, filename, newNumOfClusters)
  } else {
    clusters <- data.frame(GeneID = names(clusterResult$classification), Cluster = clusterResult$classification, row.names = NULL)
    table(clusters$Cluster)
    write.csv(clusters, filename, row.names = FALSE, quote = FALSE)
    return(clusters) 
  }
}

## Prepare data structure
## ZValues Matrix the Z Values of the proteins selected
## clusters Matrix The clusters of the proteins
## Rodrigo Dorado
getZvaluesMelt <- function(ZValues, clusters) {
  ZValues_melt <- ZValues
  ZValues_melt$GeneID <- rownames(ZValues_melt)
  ZValues_melt <- merge(ZValues_melt, clusters, by.x = 'GeneID', by.y = 'GeneID')
  ZValues_melt <- melt(ZValues_melt, id.vars = c('GeneID','Cluster'))
  colnames(ZValues_melt) <- c('GeneID', 'Cluster', 'Time', 'RelativeExpression')
  ZValues_melt$Cluster <- as.factor(ZValues_melt$Cluster)
  ##
  ZValues_melt$Time2[ZValues_melt$Time == 'TP5'] <- -48
  ZValues_melt$Time2[ZValues_melt$Time == 'TP6'] <- -24
  ZValues_melt$Time2[ZValues_melt$Time == 'TP7'] <- 0
  ZValues_melt$Time2[ZValues_melt$Time == 'TP9'] <- 1
  ZValues_melt$Time2[ZValues_melt$Time == 'TP10'] <- 2
  ZValues_melt$Time2[ZValues_melt$Time == 'TP11'] <- 4
  ZValues_melt$Time2[ZValues_melt$Time == 'TP12'] <- 6
  ZValues_melt$Time2[ZValues_melt$Time == 'TP14'] <- 24
  ZValues_melt$Time2[ZValues_melt$Time == 'TP15'] <- 48
  ZValues_melt$Time2[ZValues_melt$Time == 'TP16'] <- 72
  ##
  return(ZValues_melt)
}

## Test the cluster function with diferent sizes
## zvalues data Frame The Z values of the data to get clustered
## number_times array The sizes of the test
## Rodrigo Dorado
testCluster <- function (zvalues, number_times) {
  for (i in number_times) {
    print("***********")
    print(i)
    cluster_result <- Mclust(zvalues, G = 2:i)
    print(summary(cluster_result))
    print("-----------")
  }
}


## Set the color of the background of the image
## Rodrigo dorado
setBackGroundColors <- function() {
  rects <- as.data.frame(matrix(ncol = 3,nrow=2))
  colnames(rects) <- c('xstart','xend','col')
  rects[1,] <- c(-49,0,'gray0')
  rects[2,] <- c(0,73,'white')
  rects$xend <- as.numeric(rects$xend)
  rects$xstart <- as.numeric(rects$xstart)
  return(rects)
}

createFolder <- function(mainDir, folder) {
  if (file.exists(folder)){
    setwd(file.path(mainDir, folder))
  } else {
    dir.create(file.path(mainDir, folder))
    setwd(file.path(mainDir, folder))
  }
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## Save all the plots in a pdf file
## name String The name of the file
## data Array All the plots
## Rodrigo Dorado
savePlotPDF <- function(name, data) {
  pdf(name)
  invisible(lapply(data, print))
  dev.off()
}

directory <- "C:/Users/Rodrigo Dorado/Desktop/Git/Galaxy-Tools"
filename <- "log2_categorizada_1344prot.txt"
project_name <- "last_prpject_june_6"
id_column <- 37
first_experiment <- 0
Num_Experiments <- 10
Num_Divisions <- 3

Column_division <- Num_Experiments * Num_Divisions
All_proteins <- readCsvDirectory(directory, filename)
printInfo(All_proteins)
createFolder(directory, project_name)
proteins_complete <- getDataofInfo(All_proteins, id_column, first_experiment, Num_Experiments)
printInfo(proteins_complete)
proteins_New_column <- VarifyAComplteGroup(proteins_complete, Num_Experiments, Num_Divisions)
printInfo(proteins_New_column)
proteins_for_ANOVA <- proteins_New_column[,1:Column_division]
proteins <- proteins_New_column[,(Column_division + 1):(Column_division + Num_Experiments)]
printInfo(proteins_for_ANOVA)
printInfo(proteins)
write.csv(proteins_for_ANOVA, "Proteins For Anova.csv",row.names = TRUE, quote = FALSE)
write.csv(proteins, "Proteins With Conditions.csv",row.names = TRUE, quote = FALSE)
p_values <- calculateANOVA(proteins_for_ANOVA)
printInfo(p_values)
proteins <- mergeById(proteins, p_values)
printInfo(proteins)
selectedSignificantExpr <- selectSignificanProteinsExpression(proteins)
printInfo(selectedSignificantExpr)
Sig_Pro_Exp_ZValues <- getZValueSgnificantProteinExpression(selectedSignificantExpr)
printInfo(Sig_Pro_Exp_ZValues)
clusters <- executeMClust(Sig_Pro_Exp_ZValues, 'clusters.csv')
###
# testCluster(Sig_Pro_Exp_ZValues, c(30, 45, 60, 35, 39))
###
printInfo(clusters)
Sig_Pro_Exp_ZValues_melt <- getZvaluesMelt(Sig_Pro_Exp_ZValues, clusters)
rects <- setBackGroundColors()


##plot the result. note: make function(s)
plots_withBackground <- list()
plots_withOutBackground <- list()

for (clusterNumber in 1:max(clusters$Cluster)){
  averageProfile <- as.data.frame(matrix(ncol = 10,nrow=1))
  colnames(averageProfile) <- c('TP5','TP6','TP7','TP9','TP10','TP11','TP12','TP14','TP15','TP16')
  averageProfile$TP5 <- mean(selectedSignificantExpr[row.names(selectedSignificantExpr) %in% clusters[which(clusters$Cluster == clusterNumber),"GeneID"],"TP5"])
  averageProfile$TP6 <- mean(selectedSignificantExpr[row.names(selectedSignificantExpr) %in% clusters[which(clusters$Cluster == clusterNumber),"GeneID"],"TP6"])
  averageProfile$TP7 <- mean(selectedSignificantExpr[row.names(selectedSignificantExpr) %in% clusters[which(clusters$Cluster == clusterNumber),"GeneID"],"TP7"])
  averageProfile$TP9 <- mean(selectedSignificantExpr[row.names(selectedSignificantExpr) %in% clusters[which(clusters$Cluster == clusterNumber),"GeneID"],"TP9"])
  averageProfile$TP10 <- mean(selectedSignificantExpr[row.names(selectedSignificantExpr) %in% clusters[which(clusters$Cluster == clusterNumber),"GeneID"],"TP10"])
  averageProfile$TP11 <- mean(selectedSignificantExpr[row.names(selectedSignificantExpr) %in% clusters[which(clusters$Cluster == clusterNumber),"GeneID"],"TP11"])
  averageProfile$TP12 <- mean(selectedSignificantExpr[row.names(selectedSignificantExpr) %in% clusters[which(clusters$Cluster == clusterNumber),"GeneID"],"TP12"])
  averageProfile$TP14 <- mean(selectedSignificantExpr[row.names(selectedSignificantExpr) %in% clusters[which(clusters$Cluster == clusterNumber),"GeneID"],"TP14"])
  averageProfile$TP15 <- mean(selectedSignificantExpr[row.names(selectedSignificantExpr) %in% clusters[which(clusters$Cluster == clusterNumber),"GeneID"],"TP15"])
  averageProfile$TP16 <- mean(selectedSignificantExpr[row.names(selectedSignificantExpr) %in% clusters[which(clusters$Cluster == clusterNumber),"GeneID"],"TP16"])
  
  zval_DEGenes_cluster <- as.data.frame(t(scale(t(averageProfile))))
  zval_DEGenes_cluster_melt <- melt(zval_DEGenes_cluster)
  colnames(zval_DEGenes_cluster_melt) <- c("Time","RelativeExpression")
  zval_DEGenes_cluster_melt$GeneID <- "Average"
  
  zval_DEGenes_cluster_melt$Time2[zval_DEGenes_cluster_melt$Time == 'TP5'] <- -48
  zval_DEGenes_cluster_melt$Time2[zval_DEGenes_cluster_melt$Time == 'TP6'] <- -24
  zval_DEGenes_cluster_melt$Time2[zval_DEGenes_cluster_melt$Time == 'TP7'] <- 0
  zval_DEGenes_cluster_melt$Time2[zval_DEGenes_cluster_melt$Time == 'TP9'] <- 1
  zval_DEGenes_cluster_melt$Time2[zval_DEGenes_cluster_melt$Time == 'TP10'] <- 2
  zval_DEGenes_cluster_melt$Time2[zval_DEGenes_cluster_melt$Time == 'TP11'] <- 4
  zval_DEGenes_cluster_melt$Time2[zval_DEGenes_cluster_melt$Time == 'TP12'] <- 6
  zval_DEGenes_cluster_melt$Time2[zval_DEGenes_cluster_melt$Time == 'TP14'] <- 24
  zval_DEGenes_cluster_melt$Time2[zval_DEGenes_cluster_melt$Time == 'TP15'] <- 48
  zval_DEGenes_cluster_melt$Time2[zval_DEGenes_cluster_melt$Time == 'TP16'] <- 72
  
  numberGenes = length(unique(clusters[which(clusters$Cluster == clusterNumber),"GeneID"]))
  
  outFilePlot1 = paste("cluster1",clusterNumber,"graphWithoutBackground","pdf",sep=".")
  outFilePlot2 = paste("cluster1",clusterNumber,"graphWithBackground","pdf",sep=".")
  
  p <- ggplot(Sig_Pro_Exp_ZValues_melt[ which(Sig_Pro_Exp_ZValues_melt$Cluster == clusterNumber),], aes(x = Time2, y = RelativeExpression, group = GeneID)) + 
    geom_line(col = "lightgray") + 
    geom_line(data = zval_DEGenes_cluster_melt, aes(x = Time2, y = RelativeExpression)) +
    theme_bw() +
    xlab("Time point") +
    ylab("Relative expression") +
    ggtitle(paste("Cluster 1", clusterNumber, sep = ".")) +
    annotate("text", x = 3, y = 1, label = paste(numberGenes,"proteins", sep = " "))
  ##ggsave(filename = outFilePlot1, plot = p)
  plots_withOutBackground[[clusterNumber]] <- p
  p2 <- ggplot() + 
    geom_line(data = Sig_Pro_Exp_ZValues_melt[which(Sig_Pro_Exp_ZValues_melt$Cluster == clusterNumber),], aes(x = Time2, y = RelativeExpression, group = GeneID), col = "lightgray") +
    geom_line(data = zval_DEGenes_cluster_melt, aes(x = Time2, y = RelativeExpression, group = 1), col = "black") +
    geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col), alpha = 0.4, show.legend = F) +
    scale_fill_manual(values = c('#cccccc',"#FFFFFF")) +
    theme_bw() +
    xlab("Time point") +
    ylab("Relative expression") +
    ggtitle(paste("Cluster 1", clusterNumber, sep = ".")) +
    annotate("text", x = 3, y = 1, label = paste(numberGenes, "proteins", sep = " "))
  ##ggsave(filename = outFilePlot2, plot = p2)
  plots_withBackground[[clusterNumber]] <- p2
}
savePlotPDF("withOutBackground.pdf",plots_withOutBackground)
savePlotPDF("withBackground.pdf",plots_withBackground)
