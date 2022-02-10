# script for extracting and counting integration sites from sam files. launched from 04_extractsites.sh 

# get file name from command line (or info only)
commandline <- commandArgs(trailingOnly = TRUE)

## setup ####

# load required libraries and functions
suppressMessages(suppressWarnings(library(Biostrings)))
suppressMessages(suppressWarnings(library(tidyverse)))

whichbit <- function(flag, bit) {
  # test imput
  if (!is.vector(flag)) {stop("input must be a vector")}
  # if (length(flag) > 1) {stop("input must be a single flag")}
  
  if (length(flag) == 1){
    # test flag to find read
    flaginbits <- intToBits(flag)
    res <- ifelse(flaginbits[bit] == 01, TRUE, FALSE)
  }
  
  # if multiplt flags are given, loop on same function
  if(length(flag) > 1){
    res <- logical(length(flag))
    for(i in 1:length(flag)){
      res[i] <- whichbit(flag[i], bit)
    }
  }
  
  # report res back
  return(res)
}

parse_cigar_neg <- function(cigar, seqlen){
  if (length(cigar) == 1) {
    # adjust seqlen based on deletions, insertions
    # inspired by script in http://www.bioinformatics.babraham.ac.uk/projects/bismark/deduplicate_bismark
    ops <- unlist(str_split(cigar, "\\d+")) # store the operations
    lens <- unlist(str_split(cigar, "\\D")) # store the length per operation
    ops <- ops[ops != ""]
    lens <- as.integer(lens[lens != ""]) # remove blanks from ends
    
    longcigar <- tibble(ops, lens)
    
    correctionfactor = seqlen + 
      sum(longcigar$lens[longcigar$ops == "D"]) - 
      sum(longcigar$lens[longcigar$ops == "I"]) - 
      parse_cigar_pos(cigar)
  } else{
    correctionfactor <- numeric(length(cigar))
    for (i in 1:length(cigar)) {
      correctionfactor[i] <- parse_cigar_neg(cigar[i], seqlen[i])
    }
  }
  return(correctionfactor)
}

parse_cigar_pos <- function(cigar){
  # function for calculating the correction factor required for negative stranded reads. 
  # NOTE: This is different from the fusion read version of this function
  # output: the correction factor (numeric) to use.
  if (length(cigar) == 1) {
    if (grepl("^\\d+M", cigar)) {
      correctionfactor <- 0
    } else if (grepl("^\\d+H\\d+M", cigar)) {
      correctionfactor <- 0
    } else if (grepl("^\\d+S\\d+M", cigar)) {
      correctionfactor <- as.numeric(gsub("^(\\d+)S\\d+M.*", "\\1", cigar))
    } else if (grepl("^\\d+I\\d+M", cigar)) {
      correctionfactor <- 0
    } else {
      stop("cigar not planned for ", cigar)
    }
  }
  
  if (length(cigar) > 1) {
    correctionfactor <- numeric(length(cigar))
    for (i in 1:length(cigar)) {
      correctionfactor[i] <- parse_cigar_pos(cigar[i])
    }
  }
  
  return(correctionfactor)
}

extractread <- function(df){
  
  if (!require(tidyverse)) {
    stop("tidyverse is required to run this function")
  } else {
    
    df <- df %>% 
      as_tibble() %>%
      filter(grepl("^AS", X14), 
             grepl("^XS", X15)) %>%
      mutate(
        seqlen = nchar(X10),
        chr = factor(X3, levels = c(1:22, "X", "Y")), 
        ort = ifelse(whichbit(X2, 5), "-", "+"),
        cigarcorfact = ifelse(ort == "-", 
                              parse_cigar_neg(X6, seqlen), 
                              parse_cigar_pos(X6)),
        seq = ifelse(ort == "-", 
                     as.character(reverseComplement(DNAStringSet(X10))), 
                     X10),
        pos = ifelse(ort == "-", 
                     X4 + cigarcorfact - 1, 
                     X4 - cigarcorfact), 
        AS = as.numeric(gsub(".+:", "", X14)), 
        XS = as.numeric(gsub(".+:", "", X15)) 
      ) %>%
      filter(AS >= 36) %>%
      select(coor = X1, chr, pos, ort, seq, mapq = X5, seqlen, 
             cigar = X6, AS, XS) 
    
  }
  
  return(df)
  
}

bamname <- commandline[1]

fbname <- gsub(".bam", "", basename(bamname))
dbname <- dirname(bamname)

filename1 <- paste0(fbname, "_R1.sam")
filename2 <- paste0(fbname, "_R2.sam")

outputsisters <- paste0(fbname, ".sisters")
outputsites <- paste0(fbname, ".sites")
outdir <- gsub("04_align", "05_extract_sites", dbname)

cat("\n\nfiles are:", filename1, "and", filename2, "\n")
cat("outputs are:", outputsisters, "and", outputsites, "\n")

## run starts here ####
cat("\n----\nreading in files...\n")
timestamp()

r1 <- read.table(file.path(dbname, filename1),
                 col.names = paste0("X", 1:16),
                 fill = TRUE,
                 quote = "",
                 stringsAsFactors = FALSE,
                 comment.char = "")

r2 <- read.table(file.path(dbname, filename2),
                 col.names = paste0("X", 1:16),
                 fill = TRUE,
                 quote = "",
                 stringsAsFactors = FALSE,
                 comment.char = "")

# Test whether there was a problem with the reading of the quality string
if (length(which(nchar(r1$X11) != nchar(r1$X10))) > 0) {stop("reads not read in correctly!\n"); quit(status = 123)}
if (length(which(nchar(r2$X11) != nchar(r2$X10))) > 0) {stop("reads not read in correctly!\n"); quit(status = 123)}

orig.r1 <- nrow(r1)
orig.r2 <- nrow(r2)

cat("\n----\nextracting files...\n")
timestamp()
r1 <- extractread(r1)
r2 <- extractread(r2)

ext.r1 <- nrow(r1)
ext.r2 <- nrow(r2)


# exit if no reads remain
if (nrow(r1) == 0 | nrow(r2) == 0) {
  tibble(bamname, orig.r1, orig.r2, ext.r1, ext.r2, ext.joined = NA) %>% 
    write_tsv(path = file.path(outdir, "extractionquant.txt"), append = TRUE)
  cat("\n", paste("counts:", bamname, orig.r1, orig.r2, 
                  ext.r1, ext.r2, "NA", sep = "\t"), "\n")
  cat("no reads!\n")
  quit()
  } else {
  cat("\n----\nfiles extracted...\n")
  timestamp()
}

joined <- inner_join(r1, r2, by = "coor", suffix = c("1", "2")) %>%
  mutate(length = abs(pos2 - pos1) + 1)

ext.joined <- nrow(joined)

cat("\n----\nsummerising files...\n")
timestamp()

joined %>%
  group_by(chr1, pos1, ort1, chr2, pos2, ort2, length) %>%
  arrange(desc(seqlen1)) %>%
  summarise(totaldupl = n(), minseqlen1 = min(seqlen1), minmapq1 = min(mapq1), minmapq2 = min(mapq2),
            minAS1 = min(AS1), minAS2 = min(AS2), seqexample = first(seq1)) %>%
  ungroup() %T>%
  write_tsv(file.path(outdir, outputsisters)) %>%
  group_by(chr1, pos1, ort1) %>%
  arrange(desc(nchar(seqexample))) %>%
  summarise(totalshear = n(), totaldupl = sum(totaldupl), minseqlen1 = min(minseqlen1), minmapq1 = min(minmapq1),
            minmapq2 = min(minmapq2), minAS1 = min(minAS1), minAS2 = min(minAS2), seqexample = first(seqexample)) %>%
  ungroup() %>%
  arrange(desc(totalshear), desc(totaldupl)) %>%
  write_tsv(file.path(outdir, outputsites))

cat("\n----\ndone!\n")
timestamp()

tibble(bamname, orig.r1, orig.r2, ext.r1, ext.r2, ext.joined) %>% 
  write_tsv(path = file.path(outdir, "extractionquant.txt"), append = TRUE)
cat("\n", paste("counts:", bamname, orig.r1, orig.r2, 
                ext.r1, ext.r2, ext.joined, sep = "\t"), "\n")
