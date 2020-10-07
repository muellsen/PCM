### Prepare data
library(phyloseq)
library(Matrix)
library(microbiome)

#### Load soil data ####
pHdata <- import_biom("./pH_soil/238_otu_table.biom") # load OTU table and taxonomy table
sample <- read.delim("./pH_soil/88soils_metadata.txt", sep="\t", header=TRUE, row.names=1) # load sample data
sample_data(pHdata) <- sample # assign metadata to phyloseq object
# 
# > pHdata
# otu_table(pHdata)   OTU Table:         [ 7396 taxa and 89 samples ]
# sample_data(pHdata) Sample Data:       [ 89 samples by 69 sample variables ]
# tax_table(pHdata)   Taxonomy Table:    [ 7396 taxa by 7 taxonomic ranks ]
#
# Filter OTUs with fewer than 100 reads
# Ref: https://biocore.github.io/gneiss/docs/v0.4.0/tutorials/qiime2/88soils-qiime2-tutorial.html
pHdata_filt <- prune_taxa(taxa_sums(pHdata) > 100, pHdata)
# Sanity check: are the 116 taxa the same ones that were used in the Balance Tree paper?
library(readxl)
table_paper <- read_excel("./inline-supplementary-material-2.xlsx")
all.equal(sort(as.numeric(rownames(otu_table(pHdata_filt)))), sort(table_paper$`#OTU_ID`)) # Passed 
# Drop Sample 103.BB1 which has counts for only one OTU
nOTUperSample <- colSums(otu_table(pHdata_filt) != 0)
pHdata_filt <- prune_samples(nOTUperSample > 1, pHdata_filt)
y <- sample_data(pHdata_filt)$ph
X <- t(otu_table(pHdata_filt)@.Data)

####
# Replicate mean PH for OTUs as were used in the Balance Tree paper. DID NOT SUCCEED.
# OTU mean PH formula: Eq (3) of https://msystems.asm.org/content/2/1/e00162-16
# Ref: https://biocore.github.io/gneiss/docs/v0.4.0/tutorials/qiime2/88soils-qiime2-tutorial.html
# Ref: https://biocore.github.io/gneiss/docs/v0.4.0/tutorials/python/88soils-python-tutorial.html
#otu_tb <- otu_table(pHdata_filt)
#otu_tb[otu_tb == 0] <- 1 # pseudo count for zero abundance
#X <- otu_tb %*% diag(1/colSums(otu_tb)) # proportions of OTU in every sample (in column)
#X <- diag(1/rowSums(X)) %*% X # proportions of samples containing every OTU (in row)
#table_sample <- sample_data(ag_filt)
#mean_ph <- X %*% table_sample$ph
#range(mean_ph[match(table_paper$`#OTU_ID`, rownames(otu_table(ag_filt)))] - table_paper$mean_ph)
####

# Write data as csv files for use in MATLAB
write.csv(y,file = "pHData.csv")
write.csv(X,file = "soilOTUData.csv")
write.csv(tax_table(pHdata_filt),file = "pHTaxTable.csv")

## Write taxonomy table (gives weird format in MATLAB)
#write_phyloseq(pHdata_filt, 'TAXONOMY')




