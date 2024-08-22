install.packages("ape")    
install.packages("pegas")   
install.packages("adegenet") 
install.packages("poppr")
install.packages("picante") 
library(ape)
library(pegas)
library(adegenet)
library(poppr)
library(picante)

# PD
nj_tree <- read.tree("/Users/juanjuan/new/cleaned_7ktree.tre")
plot(tree)
sequence_ids <- tree$tip.label
comm_matrix <- matrix(1, nrow = 1, ncol = length(sequence_ids))
colnames(comm_matrix) <- sequence_ids
rownames(comm_matrix) <- "Sample1"
pd_value <- pd(comm_matrix, tree)
print(pd_value)

#Bayesian fit
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("Biostrings")
}
library(Biostrings)
library(rstan)

fasta_file <- "/Users/juanjuan/new/4k.fasta"
fasta_data <- readDNAStringSet(fasta_file)
print(fasta_data)
sequence_counts <- table(as.character(fasta_data))
print(sequence_counts)
species_count <- as.numeric(sequence_counts) 
data_list <- list(
  N = length(species_count),
  y = species_count
)
fit <- stan(model_code = model_code, data = data_list, iter = 2000, chains = 4)
print(fit)
stan_plot(fit, pars = "lambda")
stan_trace(fit, pars = "lambda")
summary(fit)
print(fit)

#shannon and simpson

feature_table <- read.csv("/Users/juanjuan/new/otu_counts.tsv", sep="\t", header=TRUE,)
otu_counts <- aggregate(. ~ OTU_ID, data=feature_table, sum)
total_counts <- sum(feature_table[, 2])
head(otu_counts)
shannon_index <- diversity(otu_counts[, -1], index="shannon")
simpson_index <- diversity(otu_counts[, -1], index="simpson")

#Chao1

chao1_estimate <- estimateR(feature_table[, 2])
chao1 <- chao1_estimate["S.chao1"]
print(paste("estimated species richness:", chao1))

