# Load required packages
library(gdata)      # for reading .xls files
library(haplo.stats)
library(adegenet)

library(genetics)
library(hierfstat)
data(hla.demo)
typeof(hla.demo)
attach(hla.demo)


# Read in the Haploview-input_N_Ansari-Pour.xls file
geno <- hla.demo[,c(17,18,21:24)]
haploview_file <- read_xls("D:/grad/Haploview/Haploview-input_N_Ansari-Pour.xls")
genotype_data <- t(haploview_file[,3:ncol(haploview_file)])



genotype_data <- genotype(genotype_data , sep = "")
genotype_data_DF <- makeGenotypes(genotype_data)
# Compute LD table for all 3 genotypes
ldt <- LD(genotype_data_DF)

# display the results
print(ldt)                               # textual display
LDtable(ldt)  




geno_id <- df2genind(geno , sep = "")
geno_light <- makeGenotypes(geno)
geno_mx <- as.matrix(geno_light)
LD_plot <- glPlot(geno_mx , col=NULL, legend=TRUE, posi="bottomleft")



data <- makeGenotypes(genotype_data)
data <- as.matrix(data)
ldt <- LD(data)
plot(ldt, digits=2, marker=19)



label <-c("DQB","DRB","B")
keep <- !apply(is.na(geno) | geno==0, 1, any)

save.em.keep <- haplo.em(geno)
save.em.keep2 <- haplo.em(genotype_data)

score.gaus <- haplo.score(resp, geno, locus.label=label, 
                          trait.type = "gaussian")

plot.haplo.score(score.gaus)


anotherType<- glGeno(genotype_data, sep="")


geneId_hap <- df2genind(genotype_data , sep = "")
#geneId_hap <- tab(geneId_hap , NA.method = "mean")
genLight_haplo <- as.genlight(genotype_data , ncode=2)
# Replace missing values in the genlight object with 0
genLight_haplo <- tab(geneId_hap)
LD_plot <- glPlot(genLight_haplo , col=NULL, legend=TRUE, posi="bottomleft")

# warning: output will not exactly match

print.haplo.em(save.em.keep2)
# Calculate pairwise LD between SNPs in the haplotype object
ld_data <- gl.ld.haplotype(save.em.keep2)

# Generate an LD plot from the pairwise LD data
pdf("ld_plot.pdf") # Open a PDF file to save the plot
plot.ld(ld_data) # Generate the LD plot
dev.off() # Close the PDF file