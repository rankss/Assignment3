# Load libraries
library(tidyverse)
library(rentrez)
library(Biostrings)
library(stringi)
library(dplyr)

# Set working directory, for myself only
setwd("C:/Users/benya/OneDrive - University of Guelph/Courses/BINF6210/Assignment2")

# Installing tensorflow/keras

# Uncomment the below 2 lines if needing to download tensorflow from scratch
# install.packages("remotes")
# remotes::install_github("rstudio/tensorflow")

# Uncomment the below 2 lines if needing to install python (git must be in PATH)
# library(reticulate)
# reticulate:install_python()

library(tensorflow)
library(keras)
install_keras() # Need to run initially before the rest of the code

# NCBI nuccore search results (Data was too large so I downloaded instead)
entrez_search_fetch <- function(db, term, retmax, rettype) {
    search <- entrez_search(db = db, term = term, retmax = retmax, use_history = TRUE)
    fetch <- entrez_fetch(db = db, web_history = search$web_history, rettype = rettype)
    return(fetch)
}

# Commented out so they were not ran
entrez_search_fetch(db = "nuccore", term = "Homo sapiens[ORGN] AND biomol_mrna[PROP] AND refseq[filter] AND 700:750[SLEN]", retmax = 1000, rettype = "fasta") %>%
     write("mRNA_sequences.fasta", sep = "\n")

entrez_search_fetch(db = "nuccore", term = "Homo sapiens[ORGN] AND biomol_rrna[PROP] AND refseq[filter]", retmax = 1000, rettype = "fasta") %>%
    write("rRNA_sequences.fasta", sep = "\n")

entrez_search_fetch(db = "nuccore", term = "animals[filter] AND biomol_trna[PROP]", retmax = 1000, rettype = "fasta") %>%
    write("tRNA_sequences.fasta", sep = "\n")

# Preprocessing of data into working set for training neural network
# Concept: Split sequence into profiles of 3-mer frequency probability as input into network

# Mapping of integer to types of RNA
RNAMap <- c("mRNA", "rRNA", "tRNA")
RNAMapFunc <- function(x) {
    if (x[1]) {
        return(RNAMap[1])
    } else if (x[2]) {
        return(RNAMap[2])
    } else if (x[3]) {
        return(RNAMap[3])
    } 
    return("Other")
}

# Create function to process all 3 data frames similarly
readFastaToDataFrame <- function(filename, type) {
    # Read fasta file, obtain trinucleotide frequency and add necessary columns to data frame
    DNA <- readDNAStringSet(filename)
    df <- trinucleotideFrequency(DNA, as.prob = TRUE) %>% as.data.frame()
    df$type <- rep(type, each = nrow(df))
    df$names <- names(DNA)
    df$accession <- word(df$names, 1L)
    return(df)
}

# Function that convert to one-hot encoding matrix and then into data frame
typeToOneHot <- function(df) {
    df$mRNA <- as.integer(df$type == "mRNA")
    df$rRNA <- as.integer(df$type == "rRNA")
    df$tRNA <- as.integer(df$type == "tRNA")
    return(df)
}

# Read fasta sequence into data frame
mRNA.df <- readFastaToDataFrame("mRNA_sequences.fasta", "mRNA")
rRNA.df <- readFastaToDataFrame("rRNA_sequences.fasta", "rRNA")
tRNA.df <- readFastaToDataFrame("tRNA_sequences.fasta", "tRNA")

# Summarize the average of 3-mers
mRNA.summary <- summarize_all(subset(mRNA.df, select = c(1:4^3)), mean)
rRNA.summary <- summarize_all(subset(rRNA.df, select = c(1:4^3)), mean)
tRNA.summary <- summarize_all(subset(tRNA.df, select = c(1:4^3)), mean)

# Concatenate the summary and produce barplot
df.summary <- rbind(mRNA.summary, rRNA.summary, tRNA.summary)
rownames(df.summary) <- c("mRNA", "rRNA", "tRNA")
par(las=2)
barplot(as.matrix(df.summary), ylab = "Normalized Frequency", legend = TRUE, args.legend = list(cex = 0.5))

# Combine and scramble data frame, reset indices
RNA.df <- rbind(mRNA.df, rRNA.df, tRNA.df)
RNA.df <- RNA.df[sample(nrow(RNA.df)), ]
RNA.df <- typeToOneHot(RNA.df)
row.names(RNA.df) <- NULL

# Set seed for reproducibility for random training/testing set selection
# set.seed(42)

# Determine training/testing set size (default 90% training)
train.size = round(nrow(RNA.df) * 0.9)

# Split dataframe into training/testing set
train.df <- RNA.df[1:train.size, ]
test.df <- RNA.df[-(1:train.size), ]

# Split training/testing dataframe into x (3-mer profile) and y (one hot) matrices
train.x <- subset(train.df, select = c(1:4^3)) %>% as.matrix()
train.y <- subset(train.df, select = RNAMap) %>% as.matrix()
test.x <- subset(test.df, select = c(1:4^3)) %>% as.matrix()
test.y <- subset(test.df, select = RNAMap) %>% as.matrix()

# Define the simple neural network sequential model
model <- keras_model_sequential() %>%
    layer_dense(units = ncol(train.x), activation = "relu", input_shape = ncol(train.x)) %>%
    layer_dense(units = 512, activation = "relu") %>%
    layer_dense(units = 1024, activation = "relu") %>%
    layer_dense(units = ncol(train.y), activation = "softmax") %>%
    compile(loss = "categorical_crossentropy", optimizer = "adam", metrics = "accuracy")

# Obtain and plot the history of training
history <- fit(model, train.x, train.y, epochs = 10, batch_size = 128, validation_split = 0.1)
plot(history)

# Evaluate and report the results of the model
report <- function(model, x) {
    result <- data.frame(round(predict(model, x)))
    result <- apply(result, 1, RNAMapFunc)
    return(result)
}

evaluate(model, test.x, test.y, verbose = 0)
result <- cbind(subset(test.df, select = c("accession", "type")), report(model, test.x))
colnames(result) <- c("accession", "type", "predict")

# Check which sequences produced different than true results
result.diff <- result[result$type != result$predict, ]
result.diff

# Function to easily test model accuracy on other fasta files
modelAccuracy <- function(filename) {
    # Read and process file, retain only mRNA, rRNA, and tRNA sequences
    DNA <- readDNAStringSet(filename)
    df <- trinucleotideFrequency(DNA, as.prob = TRUE) %>% as.data.frame()
    df$names <- names(DNA)
    # Extract RNA type via regex, if the type is not in the title, it is likely an mRNA
    df$type <- stri_extract_first_regex(df$names, "[:lower:]{1,3}RNA") %>% replace_na("mRNA")
    df <- df[df$type == "mRNA" | df$type == "rRNA" | df$type == "tRNA", ]
    df <- typeToOneHot(df)
    
    # Split data frame into x (input) and y (output)
    df.x <- subset(df, select = c(1:4^3)) %>% as.matrix()
    df.y <- subset(df, select = RNAMap) %>% as.matrix()
    
    # Evaluate and report results of model
    print(evaluate(model, df.x, df.y, verbose = 0))
    df.result <- cbind(subset(df, select = c("names", "type")), report(model, df.x))
    colnames(df.result) <- c("names", "type", "predict")
    df.diff <- df.result[df.result$type != df.result$predict, ]
    
    #return(df.diff)
}

# Testing on other data sets
modelAccuracy("genome.fasta") # Downloaded from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz

modelAccuracy("tRNA_testing.fasta") # Downloaded from http://gtrnadb.ucsc.edu/Hsapi19/hg19-tRNAs.fa
