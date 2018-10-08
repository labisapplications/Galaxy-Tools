#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("Wrong number of parameters!\n\n")
  cat("Usage:\n")
  cat("    Rscript main.R <input data csv> <output file>\n\n")
  quit(save = "no", status = 1)
}

source('iota/IOTA.R')
library('igraph', warn.conflicts = FALSE)

filename <- args[1]
output <- args[2]
cat(paste("Processing file", filename, "\n"))

raw_data <- read.csv(filename, sep = ';')

# Transpose data
data <- raw_data[, -1]
rownames(data) <- raw_data[, 1]

# Normalize data
tmp <- data - apply(data, 1, min, na.rm = TRUE)
data <- tmp / apply(tmp, 1, max, na.rm = TRUE)

# Run IOTA
a <- IOTAsigned(data, method = 'sqrt')
diag(a) <- 0  # Ignore selfcorrelation

threshold <- 0.95

up.regulated <- graph.adjacency(a > threshold)
down.regulated <- graph.adjacency(a < -threshold)

result <- rbind(
                cbind(get.edgelist(up.regulated), "+"),
                cbind(get.edgelist(down.regulated), "-")
                )

write.table(result, file = output,
            quote = FALSE, sep = '\t',
            col.names = FALSE, row.names = FALSE)
