###################################
# Practicum 2: Microarrays Part 2 #
###################################


## Set up your working directory
#####################

working_dir <- "C:\\Users\\ruben\\OneDrive - UGent\\AHTA"
setwd(working_dir)


## Install packages
#####################

# BiocManager::install('minfi', update = F)
# BiocManager::install('ENmix', update = F)
# BiocManager::install('IlluminaHumanMethylationEPICanno.ilm10b4.hg19', update = F)
# BiocManager::install('IlluminaHumanMethylationEPICmanifest', update = F)
# BiocManager::install('dplyr', update = F)
# BiocManager::install('tidyr', update = F)

library(minfi) #standaard voor illumina methylatie-arrays
library(ENmix) #Zie dat rlang boven 1.3.0 blijft om te werken, preprocessing
library(limma) #lin-modellen voor dif analyse
library(dplyr) #data-frames manipuleren
library(tidyr) #data-frames manipuleren

## Import Data
####################

## Load the phenotypic annotation data
pheno <- read.table("InfiniumEPIC_Data/Phenotypic_Data_Infinium.txt", header = T, sep = "\t")

## Load the Infinium data
fibromethyl_RG <- read.metharray.exp("InfiniumEPIC_Data", extended = TRUE)

## Match the sample names to the phenotypic sample IDs (first 10 characters)
sampleNames(fibromethyl_RG) <- substr(sampleNames(fibromethyl_RG), 1, 10)

## Make sure that they are in the same order!
all(sampleNames(fibromethyl_RG) == pheno$SampleID) 

## Add the phenotypic annotation to the fibromethyl_RG object
rownames(pheno) <- pheno$SampleID
pData(fibromethyl_RG) <- DataFrame(pheno)

## Take a look at your data
fibromethyl_RG

class(fibromethyl_RG) # object type (RGChannelSetExtended)
dim(fibromethyl_RG) # rows: probes, cols: samples
sampleNames(fibromethyl_RG) # or colnames(fibromethyl_RG),  sample IDs
head(featureNames(fibromethyl_RG)) # or rownames(fibromethyl_RG), probe IDs
pData(fibromethyl_RG) # phenotypic data
head(getAnnotation(fibromethyl_RG)) # probe/gene annotation data
annotation(fibromethyl_RG) # type of microarray
assays(fibromethyl_RG) # available numeric data

head(getRed(fibromethyl_RG)) # red intensities
head(getGreen(fibromethyl_RG)) # green intensities
head(getNBeads(fibromethyl_RG)) # number of beads on the array

## Subset data
fibromethyl_RG <- fibromethyl_RG[, pData(fibromethyl_RG)$Origin %in% c("NAF", "CAF")]

## Change sample names
sampleNames(fibromethyl_RG) <- paste(sampleNames(fibromethyl_RG), pData(fibromethyl_RG)$Origin, sep = "_")

## Take a look at the subset of the data
sampleNames(fibromethyl_RG)
dim(fibromethyl_RG)
pData(fibromethyl_RG)

## Convert RGChannelSetExtended to fibromethyl_N (no preprocessing)
fibromethyl_raw <- preprocessRaw(fibromethyl_RG)

## Take a look at the fibromethyl_N
fibromethyl_raw
class(fibromethyl_raw) # object type (fibromethyl_N)
dim(fibromethyl_raw) # rows: CpGs, cols: samples
sampleNames(fibromethyl_raw) # or colnames(fibromethyl_raw),  sample IDs
head(featureNames(fibromethyl_raw)) # or rownames(fibromethyl_raw), CpG IDs
pData(fibromethyl_raw) # phenotypic data
head(getAnnotation(fibromethyl_raw)) # probe/gene annotation data
annotation(fibromethyl_raw) # type of microarray
assays(fibromethyl_raw) # available numeric data

head(getMeth(fibromethyl_raw)) # methylated intensities
head(getUnmeth(fibromethyl_raw)) # unmethylated intensities

## Function to calculate M values
getM <- function(fibromethyl_N, offset = 1) {
  M <- log2((getMeth(fibromethyl_N) + offset) / (getUnmeth(fibromethyl_N) + offset))
  return(M)
}

## Take a look at the raw M- and beta-values
head(getB(fibromethyl_raw, offset = 100)) # offset is 100 by default
head(getM(fibromethyl_raw, offset = 1))  # offset is 1 by default


## Quality Control on raw data
####################


# Raw data QC plotting function

rawQCplot <- function(rgset, fibromethyl_N, filename = "RawQC.jpg")
{
  jpeg(filename = filename, height = 1000, width = 1200)
  par(mfrow = c(3, 2)) # subplots: 3 rows & 2 columns
  
  # Plot red and green intensities together (use RGChannelSet)
  multifreqpoly(
    mat = log2(getGreen(rgset) + 1),
    col = "green",
    main = "Raw Red & Green Intensities",
    xlab = "Log2 Intensity",
    legend = NULL # remove legend
  )
  multifreqpoly(
    mat = log2(getRed(rgset) + 1),
    col = "red",
    append = TRUE, # use append to combine 2 plots
    legend = list(
      x = "left", # change legend position
      legend = c("Green", "Red"), # change legend values
      fill = c("green", "red") # change legend color
    )
  )
  
  # Plot methylated and unmethylated intensities together (use fibromethyl_N)
  multifreqpoly(
    mat = log2(getMeth(fibromethyl_N) + 1),
    col = "orange",
    main = "Raw Meth & Unmeth Intensities",
    xlab = "Log2 Intensity",
    legend = NULL
  )
  multifreqpoly(
    mat = log2(getUnmeth(fibromethyl_N) + 1),
    col = "blue",
    append = TRUE,
    legend = list(
      x = "left",
      legend = c("Methylated", "Unmethylated"),
      fill = c("orange", "blue")
    )
  )
  
  # Plot summed meth and unmeth intensities
  multifreqpoly(
    mat = log2(getMeth(fibromethyl_N) + getUnmeth(fibromethyl_N) + 1),
    main = "Raw Sum of Meth & Unmeth Intensities",
    xlab = "Log2 Intensity",
    legend = list(
      x = "left",
      legend = colnames(fibromethyl_N),
      fill = 1:ncol(fibromethyl_N)
    )
  )
  
  # Plot combined meth and unmeth intensities
  multifreqpoly(
    mat = log2(rbind(getMeth(fibromethyl_N) + 1, getUnmeth(fibromethyl_N) + 1)),
    main = "Raw Combined Meth & Unmeth Intensities",
    xlab = "Log2 Intensity",
    legend = list(
      x = "left",
      legend = colnames(fibromethyl_N),
      fill = 1:ncol(fibromethyl_N)
    )
  )
  
  # Plot beta-values
  multifreqpoly(
    mat = getB(fibromethyl_N),
    main = "Raw Beta-values",
    xlab = "Beta-value",
    legend = list(
      x = "top",
      legend = colnames(fibromethyl_N),
      fill = 1:ncol(fibromethyl_N)
    )
  )
  
  # Plot M-values
  multifreqpoly(
    mat = getM(fibromethyl_N),
    main = "Raw M-values",
    xlab = "M-value",
    legend = list(
      x = "left",
      legend = colnames(fibromethyl_N),
      fill = 1:ncol(fibromethyl_N)
    )
  )
  
  dev.off() # close the plot
}

# Run the function
rawQCplot(rgset = fibromethyl_RG, fibromethyl_N = fibromethyl_raw, filename = "RawQC.jpg")

## Filtering

# Quality control function
fibromethyl_QC <- getQC(fibromethyl_raw)

# NA probes
NA_CpGs <- is.na(fibromethyl_QC)

## Preprocessing
####################

## Background correction & dye-bias correction & filtering

?preprocessENmix

fibromethyl_BC <- preprocessENmix(fibromethyl_RG, bgParaEst = "oob", dyeCorr="RELIC", QCinfo=NULL, exQCsample=TRUE,
                                  exQCcpg=TRUE, exSample=NULL, exCpG=NULL, nCores = 2)

## Normalisation

?norm.quantile

fibromethyl_N <- norm.quantile(fibromethyl_BC, method = "quantile1")


## Quality Control on preprocessed data
####################

# Preprocessed data QC function
preproQCplot <- function(fibromethyl_N, filename = "PreProQC.jpg") 
  {
  # Create a combined quality control plot
  jpeg(filename = filename, height = 1000, width = 1200)
  par(mfrow = c(3, 2)) # subplots: 3 rows & 2 columns
  
  # Which probes are green and red and TII
  color <- getAnnotation(fibromethyl_N)$Color
  green <- color == "Grn"
  red <- color == "Red"
  t2 <- color == ""
  
  # Plot red and green intensities together
  multifreqpoly(
    mat = log2(rbind(
      getMeth(fibromethyl_N)[green, ] + 1,
      getUnmeth(fibromethyl_N)[(green | t2), ] + 1
    )),
    col = "green",
    main = "Preprocessed Red & Green Intensities",
    xlab = "Log2 Intensity",
    ylim = c(0, 35000),
    legend = NULL # remove legend
  )
  multifreqpoly(
    mat = log2(rbind(
      getMeth(fibromethyl_N)[(red | t2), ] + 1,
      getUnmeth(fibromethyl_N)[red, ] + 1
    )),
    col = "red",
    append = TRUE, # use append to combine 2 plots
    legend = list(
      x = "left", # change legend position
      legend = c("Green", "Red"), # change legend values
      fill = c("green", "red") # change legend color
    )
  )
  
  
  # Plot methylated and unmethylated intensities together (use fibromethyl_N)
  multifreqpoly(
    mat = log2(getMeth(fibromethyl_N) + 1),
    col = "orange",
    main = "Preprocessed Meth & Unmeth Intensities",
    xlab = "Log2 Intensity",
    ylim = c(0, 35000),
    legend = NULL
  )
  multifreqpoly(
    mat = log2(getUnmeth(fibromethyl_N) + 1),
    col = "blue",
    append = TRUE,
    legend = list(
      x = "left",
      legend = c("Methylated", "Unmethylated"),
      fill = c("orange", "blue")
    )
  )
  
  # Plot summed meth and unmeth intensities
  multifreqpoly(
    mat = log2(getMeth(fibromethyl_N) + getUnmeth(fibromethyl_N) + 1),
    main = "Preprocessed Sum of Meth & Unmeth Intensities",
    xlab = "Log2 Intensity",
    legend = list(
      x = "left",
      legend = colnames(fibromethyl_N),
      fill = 1:ncol(fibromethyl_N)
    )
  )
  
  # Plot combined meth and unmeth intensities
  multifreqpoly(
    mat = log2(rbind(getMeth(fibromethyl_N) + 1, getUnmeth(fibromethyl_N)) + 1),
    main = "Preprocessed Combined Meth & Unmeth Intensities",
    xlab = "Log2 Intensity",
    legend = list(
      x = "left",
      legend = colnames(fibromethyl_N),
      fill = 1:ncol(fibromethyl_N)
    )
  )
  
  # Plot beta-values
  multifreqpoly(
    mat = getB(fibromethyl_N),
    main = "Preprocessed Beta-values",
    xlab = "Beta-value",
    legend = list(
      x = "top",
      legend = colnames(fibromethyl_N),
      fill = 1:ncol(fibromethyl_N)
    )
  )
  
  # Plot M-values
  multifreqpoly(
    mat = getM(fibromethyl_N),
    main = "Preprocessed M-values",
    xlab = "M-value",
    legend = list(
      x = "left",
      legend = colnames(fibromethyl_N),
      fill = 1:ncol(fibromethyl_N)
    )
  )
  
  dev.off() # close the plot
  
}

# Run the function
preproQCplot(fibromethyl_N, filename = "PreProQC.jpg")


## Differential methylation analysis
####################

## Turn categorical variables into factors
celltype <- factor(pData(fibromethyl_RG)$Origin, levels = c("NAF","CAF"))
# define your levels with the control group first, then the case group!!

## Fit a linear model on the M-values
design <- model.matrix(~ celltype)
colnames(design) <- c("Intercept","CAFvsNAF")
fit <- lmFit(getM(fibromethyl_N), design)

## Limma: apply eBayes moderation
fit2 <- eBayes(fit)

# Extract DM results and check the number of signigicant CpGs
LIMMAout <- topTable(fit2, coef = "CAFvsNAF", number = Inf, adjust.method = "BH")
sum(LIMMAout$adj.P.Val < 0.05)
summary(decideTests(fit2))

## QC of results

# Volcano plot
plot(LIMMAout$logFC, -log10(LIMMAout$P.Value),
     pch = ".", col = ifelse(LIMMAout$adj.P.Val < 0.05, "red", "black"),
     xlab = "logFC (M-value difference CAF - NAF)", ylab = "-log10(p-value)",
     main = "Volcano plot")

# MA plot
plot(LIMMAout$AveExpr, LIMMAout$logFC,
     pch = ".", col = ifelse(LIMMAout$adj.P.Val < 0.05, "red", "black"),
     xlab = "Average M-value", ylab = "logFC", main = "MA plot")
abline(h = 0, col = "blue")

# p-value histogram
hist(LIMMAout$P.Value, breaks = 50, main = "p-value histogram", xlab = "p-value")

## Check methylation for top results

# Check M-values for top results
top_cpgs <- rownames(LIMMAout)[1:10]
round(getM(fibromethyl_N)[top_cpgs, ], 2)

# Check beta-values for top results
round(getB(fibromethyl_N)[top_cpgs, ], 2)

# Add the mean betas of the CAF and NAF group, and their absolute difference
betas <- getB(fibromethyl_N)[rownames(LIMMAout), ] # same row order as LIMMAout
LIMMAout$MeanBeta_CAF <- rowMeans(betas[, celltype == "CAF"])
LIMMAout$MeanBeta_NAF <- rowMeans(betas[, celltype == "NAF"])
LIMMAout$AbsDeltaBeta <- abs(LIMMAout$MeanBeta_CAF - LIMMAout$MeanBeta_NAF)

## Annotation
#################### 

## Load annotation
anno <- getAnnotation(fibromethyl_N)
head(anno)

# ... (r. 399-404 stay as they are)



## Check whether all CpGs are present in both sets and that there are no multiple annotations per CpGs
dim(anno)
dim(LIMMAout)
all(rownames(LIMMAout) %in% anno$Name)
all(anno$Name %in% rownames(LIMMAout))
sum(duplicated(anno$Name)) # 0 => no multiples

## Get matching indices based on CpG names
ix <- match(rownames(LIMMAout), anno$Name)
# Note: match only returns the first match, not all matches (if multiple annotations are present)

## Add interesting annotation columns to limma output
colnames(anno)
LIMMAout$Chromosome <- anno$chr[ix]
LIMMAout$Position <- anno$pos[ix]
LIMMAout$Gene_ID <- anno$UCSC_RefGene_Name[ix]
LIMMAout$CpG_Feature <- anno$UCSC_RefGene_Group[ix]
LIMMAout$CpG_Island <- anno$Relation_to_Island[ix]
head(LIMMAout)


## Gene-level analysis
####################

## Summarize data on the gene level

# remove CpGs with no gene annotation
LIMMAout_gene <- LIMMAout[LIMMAout$Gene_ID != "" & !is.na(LIMMAout$Gene_ID), ]

# Save CpG rownames (otherwise these are removed in the next step)
LIMMAout_gene$CpG <- rownames(LIMMAout_gene)

# Separate CpGs with multiple genes/features (this takes a while to run)
LIMMAout_gene <- LIMMAout_gene %>% 
  separate_rows(Gene_ID, CpG_Feature, CpG_Island, sep = ";") %>%
  distinct() %>%
  group_by(CpG, Gene_ID) %>%
  mutate(CpG_Feature = paste(unique(CpG_Feature), collapse = ";"),
         CpG_Island = paste(unique(CpG_Island), collapse = ";")) %>%
  ungroup() %>%
  distinct() %>%
  as.data.frame()

# Select the most interesting columns
LIMMAout_gene <- LIMMAout_gene[, c(15, 12, 13, 14, 1, 5, 7, 8, 9)]

head(LIMMAout_gene)

## Number of DM genes
sig_gene <- LIMMAout_gene[LIMMAout_gene$adj.P.Val < 0.05, ]
length(unique(sig_gene$Gene_ID))

# Stricter: also require a biologically relevant beta difference
sig_gene_strict <- sig_gene[sig_gene$AbsDeltaBeta > 0.2, ]
length(unique(sig_gene_strict$Gene_ID))

## Look for multiple CpGs in the same gene undergoing similar methylation differences
gene_summary <- sig_gene %>%
  mutate(Direction = ifelse(logFC > 0, "Hyper_in_CAF", "Hypo_in_CAF")) %>%
  group_by(Gene_ID, Direction) %>%
  summarise(n_CpGs = n(),
            mean_logFC = mean(logFC),
            mean_AbsDeltaBeta = mean(AbsDeltaBeta),
            .groups = "drop") %>%
  filter(n_CpGs >= 3) %>%
  arrange(desc(n_CpGs)) %>%
  as.data.frame()
head(gene_summary, 20)

## Look for multiple CpGs in promoter regions undergoing similar methylation differences
promoter_summary <- sig_gene %>%
  filter(grepl("TSS200|TSS1500", CpG_Feature)) %>%
  mutate(Direction = ifelse(logFC > 0, "Hyper_in_CAF", "Hypo_in_CAF")) %>%
  group_by(Gene_ID, Direction) %>%
  summarise(n_CpGs = n(),
            mean_logFC = mean(logFC),
            mean_AbsDeltaBeta = mean(AbsDeltaBeta),
            .groups = "drop") %>%
  filter(n_CpGs >= 3) %>%
  arrange(desc(n_CpGs)) %>%
  as.data.frame()
head(promoter_summary, 20)
