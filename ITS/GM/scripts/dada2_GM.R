# Author: Benítez Domínguez, Belén
# Date: 11/05/2023
# Project: TFM: Wine yeast community assembly: effects of regionality, agricultural 
                # management, and fermentative conditions on their composition and transcriptional regulation


# Section: 5.5. DNA EXTRACTION SAMPLES AND AMPLICON-SEQ ANALYSIS (GM)


library(dada2)
library(ShortRead)
library(Biostrings)

# Set path to directory containing R1 and R2 fastq files.
path <- "E:/ITS_reads"  # This particular dataset is already demultiplexed

# Define forward and reverse read files.
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE)) # forward reads
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE)) # reverse reads

# Define forward and reverse primer sequences.
FWD <- "TCCTCCGCTTATTGATATGC"
REV <- "GTGARTCATCGAATCTTTG"

# Function to generate all orientations of a given sequence.
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# Get all possible orientations of the forward and reverse primers.
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# Filter reads for quality and Ns.
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) 
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filtN <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, verbose = TRUE) # multithreading is set to default (= FALSE), as it is not supported by Windows

# Function to count the number of reads containing a given primer sequence.
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
# Count the number of reads containing each orientation of the forward and reverse primers in the filtered reads.
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# Run cutadapt to remove primers from reads.
cutadapt <- "E:/R/cutadapt.exe" 
system2(cutadapt, args = "--version") 
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-g", REV, "-G", FWD, 
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i],
                             fnFs.filtN[i], fnRs.filtN[i]))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Get the forward and reverse fastq filenames
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names from filenames
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))

# Plot quality profile of first 2 forward and random 3 reverse fastq files
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[sample(1:length(cutRs), 3)])

# Filter and trim reads
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, truncLen = c(220,125), rm.phix = TRUE, compress = TRUE, multithread = FALSE)  # on windows, set multithread = FALSE

# Learn error rates
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)

# Plot error rates
plotErrors(errR, nominalQ = TRUE)

# Perform DADA2 denoising
dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = FALSE)

# Merge paired-end reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE, minOverlap = 5)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
hist(nchar(getSequences(seqtab)))

# Rename sequence table rows
row.names(seqtab) <- sapply(strsplit(row.names(seqtab), "_"), `[`, 1)
row.names(seqtab) <- gsub("NGS025-21-ITS2-", "ITS-", row.names(seqtab))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track the number of reads throughout the process
getN <- function(x) sum(getUniques(x))
track <- cbind.data.frame(filtN[,1], out[,2], sapply(dadaFs, getN), 
                          sapply(dadaRs, getN), sapply(mergers, getN), 
                          rowSums(seqtab.nochim))


# Assign column names and row names to track data frame
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names

# Assign taxonomy to the sequence table
taxa.cut <- assignTaxonomy(seqtab.nochim, 
                           "E:/Database/sh_general_release_10.05.2021")

# Add a new column with the percentage of non-chimeric sequences
track <- cbind.data.frame(track, perc = track[,6]*100/track[,1])

# Save the objects to disk
# save.image("dada2_GM.RData")
saveRDS(taxa.cut, "E:/Resultados_provisionales/tax_GM.rds")
saveRDS(seqtab.nochim, "E:/Resultados_provisionales/ASV_GM.rds")
write.table(track, "E:/Resultados_provisionales/track_GM.txt", sep = "\t", dec = ",")
