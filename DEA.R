# Load the library

library("DESeq2")
library("tximeta")
library("ggplot2")
library("SummarizedExperiment") 

# Import data on R
# coldata <- data.frame (tissue = tissuedata, runaccesion = runaccesion, individual = individualname, age = age, replica = replicanumber, libraryprepDate = libraryprepDate, lab = laboratoryName, seqdate = sequencingDate, laboratory = laboratoryName, batch = batchnumber)

dir <- '/home/agalab/Scrivania/RNA_seq_trial'
coldata <- read.csv(file.path(dir, "coldata.csv"))
coldata
coldata$files <- file.path(dir, "quants", coldata$tissue, "quant.sf")
all(file.exists(coldata$files)
coldata    

# Tximeta doesn not contain the horse so we need to manually import it
indexDir <- file.path(dir, "EqIndex") # my salmon index
fastaFTP <- c('ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/863/925/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_rna.fna.gz')
gtfPath <- file.path(dir,"GCF_002863925.1_EquCab3.0_genomic.gff") 


tmp <- tempdir()
jsonFile <- file.path(tmp, paste0(basename(indexDir), ".json"))
makeLinkedTxome(indexDir=indexDir,
                source="NCBI", organism="Equus caballus",
                release="EquCab3.0", genome="equCab3",
                fasta=fastaFTP, gtf=gtfPath,
                jsonFile=jsonFile)

se <- tximeta(coldata = coldata, type = "salmon", dropInfReps = TRUE)
seg <- summarizeToGene(se)
dds <- DESeqDataSet(seg, design = ~ tissue + replica)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
library("BiocParallel")
register(MulticoreParam(6))
dds <- DESeq(dds, parallel = TRUE)
plotMA(results(dds))
resultsNames (dds)
resLFC <- lfcShrink(dds, contrast="", type="apeglm")
resLFC
