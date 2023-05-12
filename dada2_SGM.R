# Author: Benítez Domínguez, Belén
# Date: 11/05/2023
# Project: TFM: Wine yeast community assembly: effects of regionality, agricultural 
# management, and fermentative conditions on their composition and transcriptional regulation


# Section: 5.5. DNA EXTRACTION SAMPLES AND AMPLICON-SEQ ANALYSIS (SGM)

library(dada2)
library(ShortRead)
library(Biostrings)

#load("E:/TFM/ITS/SGM/dada2/dada2_SGM.RData")


# Set path to directory containing R1 and R2 fastq files.
path <- 'E:/TFM/ITS/SGM/dada2/input' # This particular dataset is already demultiplexed

# Define forward and reverse read files.
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE)) # forward reads
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE)) # reverse reads


# Define forward and reverse primer sequences.
FWD <- "TCCTCCGCTTATTGATATGC"  ## forward primer sequence
REV <- "GTGARTCATCGAATCTTTG"  ## reverse primer sequence

# Function to generate all orientations of a given sequence.
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString)) 
}

# Get all possible orientations of the forward and reverse primers.
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


# Pre-filter the sequences to remove those with Ns
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) 
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filtN <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, verbose = TRUE) # multithreading is set to default (= FALSE), as it is not supported by Windows

# Function to count the number of reads containing a given primer sequence.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Count the number of reads containing each orientation of the forward and reverse primers in the filtered reads.
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[2]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[2]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[2]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[2]]))

primerHits(FWD, fnFs.filtN[[2]]) 

# Run cutadapt to remove primers from reads.
cutadapt <- "E:/R/cutadapt.exe" 
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# As primers are not right oriented, we search for the forward REV primers in R1, and forward FWD primers in R2.
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-g", REV, "-G", FWD,  #No hemos puesto los .RC pq luego no sale el plot, pero se entiende q falta un -a FWD.RC (R1) y un REV.RC (R2)
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) #input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names from filenames
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name)) 


# Quality profile of the forward (cutFs) and reverse (cutRs) reads
plotQualityProfile(cutFs[1:10]) 
plotQualityProfile(cutRs[1:3])

# Filter and trim
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, truncLen = c(220,125), rm.phix = TRUE, compress = TRUE, multithread = FALSE) 

# Visualise the estimated error rates as a sanity check
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
plotErrors(errR, nominalQ = TRUE)


# Perform DADA2 denoising
dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = FALSE)


# Merge paired-end reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE, minOverlap = 5) # Overlap=5 por??


# Construct sequence table
seqtab <- makeSequenceTable(mergers)
hist(nchar(getSequences(seqtab)))
row.names(seqtab) <- sapply(strsplit(row.names(seqtab), "_"), `[`, 1)
row.names(seqtab) <- gsub("NGS025-21-ITS2-", "ITS-", row.names(seqtab))
seqtab


### REMOVE CHIMERAS ###

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) 


# Track the number of reads throughout the process
getN <- function(x) sum(getUniques(x))
track <- cbind.data.frame(filtN[,1], out[,2], sapply(dadaFs, getN), 
                          sapply(dadaRs, getN), sapply(mergers, getN), 
                          rowSums(seqtab.nochim))


colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names


# Assign taxonomy to the sequence table
taxa.cut <- assignTaxonomy(seqtab.nochim, 
                           "C:/Users/Windows/Downloads/1E25CA4CC30A31C2E2B8CB2C89824C83D080A7F5A62E6263A0E95B37C6628067/sh_general_release_dynamic_29.11.2022.fasta")

# Add a new column with the percentage of non-chimeric sequences
track <- cbind.data.frame(track, perc = track[,6]*100/track[,1])

# Save the objects to disk
# save.image("E:/TFM/ITS/SGM/dada2/dada2_SGM.RData")
saveRDS(taxa.cut, "E:/TFM/ITS/SGM/dada2/output/tax_SGM.rds")
saveRDS(seqtab.nochim, "E:/TFM/ITS/SGM/dada2/output/ASV_SGM.rds")
write.table(track, "E:/TFM/ITS/SGM/dada2/output/track_SGM.txt", sep = "\t", dec = ",")
