
install.packages(c("adegenet", "poppr", "dplyr", "hierfstat", "reshape2", "ggplot2", "RColorBrewer", "scales", "kableExtra", "data.table", "magic", "flextable", "pegas"))
install.packages("haplo.stats")
install.packages("trio")
install.packages("arsenal")
install.packages("readxl")
# Load the haplo.stats package
library(haplo.stats)
library(trio)
library(readxl)

ls(name="package:haplo.stats")
help(haplo.em)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("trio")



geno2<- read.csv("D:/grad/Haploviewd/GH_Dokki-Mandara.xlsx")
lobster = read.csv("D:/grad/popgene/Lobster_SNP_Genotypes.csv")
# Read in the genetic data as a haplo.stats "geno" object
geno <- read.pedfile("D:/grad/Haploviewd/SNP_Dokki.ped")
# Read in the genetic data as a haplo.stats "geno" object
geno <- read.pedfile(file = "D:/grad/Haploviewd/SNP_Dokki.txt")


# Read in the genetic data as a haplo.stats "geno" object
geno <- read_excel("D:/grad/Haploview/GH_Dokki-Mandara.xlsx")
geno2 <- read_xls("D:/grad/Haploview/Haploview-input_N_Ansari-Pour.xls")
typeof(geno2)

haplo_data <- haplo.em(geno2)

# Perform haplotype analysis using the EM algorithm
hap_em <- haplo.em(geno)
# Print the haplotype frequencies
print(hap_em$frequencies)
# Perform haplotype association testing using the score test
hap_assoc <- haplo.score(hap_em, trait)
# Print the haplotype association results
print(hap_assoc$results)
# Visualize the LD patterns using the plot.ld function
plot.ld(hap_em)
# Define haplotype blocks using the block.generation function
blocks <- block.generation(hap_em)
# Print the haplotype block definitions
print(blocks$blocks)
geno