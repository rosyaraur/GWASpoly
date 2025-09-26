# GWASpoly
# R package for performing GWAS in polyploids 
# Example how we can use this for working with copy number variation 
# To simulate copy number variant (CNV) data and run a GWAS using GWASpoly 
# in R, you’ll need to adapt the input format to reflect CNV states and treat them as numeric genotypes. 
# Step 1: Simulate CNV Data
# CNVs can be represented as integer copy numbers (e.g., 0 to 4). 
# Here's how to simulate a dataset:
# Parameters
num_individuals <- 100
num_cnvs <- 500

# Simulate CNV states (0 to 4 copies)
set.seed(123)
cnv_matrix <- matrix(sample(0:4, num_individuals * num_cnvs, replace = TRUE), 
                     nrow = num_cnvs, ncol = num_individuals)

# Create marker info
marker_info <- data.frame(
  Marker = paste0("CNV", 1:num_cnvs),
  Chrom = sample(1:5, num_cnvs, replace = TRUE),
  Position = sample(1:1e6, num_cnvs)
)

# Combine into GWASpoly format
cnv_df <- cbind(marker_info, as.data.frame(cnv_matrix))
colnames(cnv_df)[-(1:3)] <- paste0("Ind", 1:num_individuals)

# Save to CSV
write.csv(cnv_df, "cnv_genotype.csv", row.names = FALSE)

# Step 2: Simulate Phenotype Data
# Simulate a quantitative trait
phenotype <- data.frame(
  Sample = paste0("Ind", 1:num_individuals),
  Trait1 = rnorm(num_individuals, mean = 10, sd = 2)
)

write.csv(phenotype, "cnv_phenotype.csv", row.names = FALSE)

# Step 3: Run GWASpoly
library(GWASpoly)

# Load genotype and phenotype
geno <- read.GWASpoly("cnv_genotype.csv", format = "numeric", ploidy = 4)
pheno <- read.table("cnv_phenotype.csv", header = TRUE, sep = ",")

# Attach phenotype
geno <- set.pheno(geno, pheno)

# Set model and run GWAS
geno <- set.model(geno, models = c("additive"))
results <- run.GWASpoly(geno)

# View results
head(results$Trait1$additive)

# Plot Results
plot(results, trait = "Trait1", model = "additive", type = "manhattan")
plot(results, trait = "Trait1", model = "additive", type = "qq")

# GWASpoly treats numeric genotypes flexibly, so CNV states (0–4) can be modeled as additive effects.
# Simulate CNV Genotype Data with 4 Populations
set.seed(123)
num_individuals <- 200
num_cnvs <- 500
num_pops <- 4

# Assign individuals to populations
pop_labels <- rep(1:num_pops, each = num_individuals / num_pops)

# Simulate population-specific CNV frequencies
cnv_freqs <- matrix(runif(num_cnvs * num_pops, min = 0.1, max = 0.9), nrow = num_cnvs)

# Generate CNV states (0–4 copies) per individual based on population
cnv_matrix <- matrix(0, nrow = num_cnvs, ncol = num_individuals)
for (i in 1:num_individuals) {
  pop <- pop_labels[i]
  cnv_matrix[, i] <- rbinom(num_cnvs, size = 4, prob = cnv_freqs[, pop])
}

# Format for GWASpoly
marker_info <- data.frame(
  Marker = paste0("CNV", 1:num_cnvs),
  Chrom = sample(1:5, num_cnvs, replace = TRUE),
  Position = sample(1:1e6, num_cnvs)
)
cnv_df <- cbind(marker_info, as.data.frame(cnv_matrix))
colnames(cnv_df)[-(1:3)] <- paste0("Ind", 1:num_individuals)
write.csv(cnv_df, "cnv_genotype.csv", row.names = FALSE)

# Step 2: Simulate Traits with CNV Effects
# Simulate 3 traits with population-specific CNV effects
num_traits <- 3
trait_matrix <- matrix(0, nrow = num_individuals, ncol = num_traits)

for (t in 1:num_traits) {
  qtl_indices <- sample(1:num_cnvs, 10)
  effects <- matrix(rnorm(10 * num_pops), nrow = 10)
  
  for (i in 1:num_individuals) {
    pop <- pop_labels[i]
    cnv_values <- cnv_matrix[qtl_indices, i]
    trait_matrix[i, t] <- sum(cnv_values * effects[, pop]) + rnorm(1)
  }
}

phenotype <- data.frame(
  Sample = paste0("Ind", 1:num_individuals),
  Trait1 = trait_matrix[, 1],
  Trait2 = trait_matrix[, 2],
  Trait3 = trait_matrix[, 3]
)
write.csv(phenotype, "cnv_phenotype.csv", row.names = FALSE)

# Step 3: Run GWAS with GWASpoly
library(GWASpoly)

# Load data
geno <- read.GWASpoly("cnv_genotype.csv", format = "numeric", ploidy = 4)
pheno <- read.table("cnv_phenotype.csv", header = TRUE, sep = ",")

# Attach phenotype
geno <- set.pheno(geno, pheno)

# Optional: add population structure (e.g., PCA)
geno <- set.structure(geno, n.PC = 3)

# Set model and run GWAS
geno <- set.model(geno, models = c("additive"))
results <- run.GWASpoly(geno)

# View results
head(results$Trait1$additive)

# Step 4: Visualize Results
plot(results, trait = "Trait1", model = "additive", type = "manhattan")
plot(results, trait = "Trait2", model = "additive", type = "qq")

# GWASpoly supports several models including:

# "general": full genotype model

# "additive": linear dosage effect

# "dominant": dominant allele effect

# "diplo-general" and "diplo-additive": for diploid-only comparisons
library(GWASpoly)

geno <- read.GWASpoly("cnv_genotype.csv", format = "numeric", ploidy = 4)
pheno <- read.table("cnv_phenotype.csv", header = TRUE, sep = ",")
geno <- set.pheno(geno, pheno)

models <- c("general", "additive", "dominant", "diplo-general", "diplo-additive")
geno <- set.model(geno, models = models)

# Example: Trait1
head(results$Trait1$general)
head(results$Trait1$additive)
head(results$Trait1$dominant)

# 5. Visualize Manhattan and QQ Plots
for (model in models) {
  plot(results, trait = "Trait1", model = model, type = "manhattan")
  plot(results, trait = "Trait1", model = model, type = "qq")
}



