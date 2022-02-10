#@ script for cleaning integration site lists from extracted data
#@ using agrepl function to find approximate matches in sequence - the same site mapped to multiple locations
#@ using SonicLength (Berry et al,  DOI: 10.1093/bioinformatics/bts004) to estimate abundance from shear site information
#@ script can be quite memory intensive, it is highly suggested to use hpc parallel system to process. here shown using IBM LSF 

#@ metadata: required to define patient id. used to separate situation where there is/isnt expectation of overlap in integration site
#@ (It is advisable to avoid this situation by sequencing samples from the same patient (eg timepoints, replicates) on different lanes)

## setup ####
# report softwares used to user
print(sessionInfo())


# define job parameters
winsize = 5000
jobid <- Sys.getenv('LSB_JOBID')
jobindex <- as.numeric(Sys.getenv('LSB_JOBINDEX'))


# define paths
datapath <- <path/to/data>
scriptpath <- <path/to/scripts> # requires additional util script, twincount.R
tmpdirpath <- <path/to/tmp> # if using parallelization, path for temporary directories.
sitepath <- file.path(datapath, "05_extract_sites") # input dir
cleansitepath <- file.path(datapath, "06_clean_sites") # output dir
metapath <- file.path(datapath, "metadata") # metadata table is required to define which samples are from the same patient


# load required libraries and functions
library(sonicLength)
library(tidyverse)
findtwins <- function(sitestoclean, i){
  twins <- tibble(V1 = i, V2 = sitestoclean$rn[sitestoclean$rn > i] ) %>%
    inner_join(sitestoclean, by = c("V1" = "rn")) %>%
    inner_join(sitestoclean, by = c("V2" = "rn")) %>%
    filter(keep.x, keep.y) %>%
    mutate(subseq.x = str_sub(seqexample.x, 1, 20), 
           subseq.y = str_sub(seqexample.y, 1, 20)) %>%
    mutate(twin = agrepl(unique(subseq.x), subseq.y) | 
             (chr1.x == chr1.y & 
                pos1.x <= pos1.y + 2 & pos1.x >= pos1.y - 2 & 
                ort1.x == ort1.y)) %>%
    filter(twin) %>% 
    select(V1, V2, subseq.x, subseq.y, totaldupl.x, totaldupl.y, 
           totalshear.x, totalshear.y, tag.x, tag.y, ptcode.x, ptcode.y) 
  
  return(twins)
}
cleansites <- function(sitestoclean, alltwins = NULL){
  cat("\n")
  
  for (i in 1:max(sitestoclean$rn)) {
    if (i %% 100 == 0) {cat(".")}
    if (i %% 10000 == 0) {cat(" ", as.character(Sys.time()), "\n")}
    if (!i %in% sitestoclean$rn[sitestoclean$keep]) {next}
    if (i >= max(sitestoclean$rn[sitestoclean$keep])) {break}
    
    if (is.null(alltwins)) {
      twins <- findtwins(sitestoclean, i) %>%
        mutate(decl = case_when(
          tag.x == tag.y ~ TRUE,                 # classical twins
          grepl("JKT|jurkat|control", ptcode.x, ignore.case = TRUE) ~ TRUE,   # always declut against controls
          grepl("JKT|jurkat|control", ptcode.y, ignore.case = TRUE) ~ TRUE, 
          tag.x != tag.y & ptcode.x != ptcode.y ~ TRUE,                       # declut when dif patients
          tag.x != tag.y & ptcode.x == ptcode.y ~ FALSE,                      # do not declut when same patients
          TRUE ~ NA
        )) %>%
        filter(decl)
    }else{
      twins <- alltwins[[i]]
      twins <- twins[twins$V2 %in% sitestoclean$rn[sitestoclean$keep], ]
    }
    
    removetwins <- 
      unique(c(twins$V1[twins$totaldupl.x == twins$totaldupl.y & 
                          twins$totalshear.x == twins$totalshear.y], 
               twins$V2))
    
    if (length(removetwins) != 0) {
      sitestoclean$keep[sitestoclean$rn %in% removetwins] <- FALSE
    }
  }
  
  return(sitestoclean)
}
wrapestabun <- function(x) {
  safeestabun <- safely(estAbund)
  y <- safeestabun(x$cloneid, x$length, min.length = min(x$length)) 
  if (is.null(y$error)) {
    aft <- enframe(y$result$theta, name = "cloneid", value = "sisters")
  } else if (is.null(y$result)) {
    cat("\nsonicLength fail, shearsites kept as sisters, total clones", 
        length(unique(x$cloneid)), ", total shear sites", nrow(x), 
        ", total reads", sum(x$totaldupl), "\n\n") 
    aft <- tibble(cloneid = as.character(x$cloneid)) %>%
      count(cloneid) %>%
      rename(sisters = n)
  } else {
    stop("ERROR check wrapestabun")
  }
  return(aft)
}


# select flowcell lane to process, It is recommended to use parallelization and select from a lane file.   
laneline <- read_tsv(file = file.path(metapath, "lanelist.txt"), 
                     col_types = "cccc") %>%
  slice(ifelse(jobindex != 0, jobindex, 99))

datapath <- file.path(sitepath, laneline$assembly, laneline$virus)
outpath <- file.path(cleansitepath, laneline$assembly, laneline$virus)
flowcell <- laneline$fc
samplelane <- laneline$lane


# load metadata
metabyfc <- read_tsv(file.path(metapath, "allmeta.txt"), 
                     col_types = "cccciciclicccccccdd") %>%
  filter(seqexpid == flowcell, lane == as.numeric(gsub("L00", "", samplelane))) 
barcodes <- read_tsv(file.path(metapath, "barcodes.txt"), 
                     col_types = c("cci"))


# choose files to read
samplestodo <- list.files(path = datapath, 
                          pattern = paste0(flowcell, ".+", samplelane, "\\.sites"))


cat(" datapath is ", datapath, "\n")
cat(" flowcell is ", flowcell, "\n")
cat(" lane is ", samplelane, "\n")
cat(" there are ", length(samplestodo), "samples, including", samplestodo[1], "\n")


## Use C.Berry's sonicLength package to estimate abundance from lengths distribution ####
cat("\n----\nestimating site abundance (sonicLength)...\n")
timestamp()

sisterlist <- list()

for (i in 1:length(samplestodo)) {
  filename <- gsub(".sites", ".sisters", samplestodo)[i]
  sisterlist[[filename]] <- read_tsv(file.path(datapath, filename),
                                     col_types = c("cdccdcidiiiiic")) %>%
    unite("cloneid", chr1, pos1, ort1) %>% 
    mutate(cloneid = as.factor(cloneid)) 
}

allabun <- lapply(sisterlist, wrapestabun) %>%
  bind_rows( .id = "filename") %>%
  separate(filename, into = c(NA, "tag", "lane", NA)) %>%
  separate(cloneid, c("chr1", "pos1", "ort1"), sep = "_") %>%
  mutate(pos1 = as.numeric(pos1))


## set sites for cleaning ####
cat("\n----\nloading site data...\n")
timestamp()


# load all sites from lane to one table
sitelist <- list()

for (i in 1:length(samplestodo)) {
  filename <- samplestodo[i]
  sitelist[[filename]] <- read_tsv(file.path(datapath, filename),
                                   col_types = c("cdciiiiiiic"))
}

sitelist <- bind_rows(sitelist, .id = "filename") %>%
  separate(filename, into = c(NA, "tag", "lane", NA)) %>%
  left_join(select(metabyfc, ptcode, index), by = c("tag" = "index")) %>%
  left_join(allabun, by = c("lane", "tag", "chr1", "pos1", "ort1")) %>%
  arrange(desc(sisters), desc(totalshear), desc(totaldupl)) %>%
  mutate(rn = row_number(), keep = TRUE)

if(nrow(allabun == nrow(sitelist))){
  cat("\nThere are ", nrow(sitelist), " UIS on this lane\n")
}else{
  stop("PROBLEM WITH SITE COUNTS - DOESNT MATCH")
}



# check if metadata available for all samples
napt <- which(is.na(sitelist$ptcode))
if (length(napt) > 0) {
  cat("not all barcodes found in metadata, VUs", 
      barcodes$Description[barcodes$Index %in% unique(sitelist$tag[napt])], 
      "set as controls\n")
  sitelist$ptcode[napt] <- paste("control", sitelist$tag[napt], sep = "-")
}




## clean sites ####

cat("\n----\ncleaning site data...\n")
timestamp()

# declutter twins - if > 50000 UIS, to be split and parallelised
if (nrow(sitelist) <= 50000) {
  
  siteclean <- cleansites(sitelist)
  
}else{

  # setup temp dir
  tmpdirname <- file.path(tmpdirpath, jobid, jobindex)
  dir.create(file.path(tmpdirpath, jobid))
  dir.create(tmpdirname)
  save(sitelist, file = file.path(tmpdirname, "all.Rdata"))
  
  
  # launch sub-jobs (to be modified based on which system is used. here using IBM LSF)
  parts <- ceiling(nrow(sitelist) / winsize)
  command <- paste0('bsub -J "detwin', jobid, '-', jobindex, '[1-', parts, ']%', 
                    parts, '" -n 4 -M 6000 -R "rusage[mem=6000]" -o ', scriptpath, 
                    '/logs/detwinpart', jobid, '-', jobindex, '.%J.%I.out -e ', 
                    scriptpath, '/logs/detwinpart', jobid, '-', jobindex, 
                    '.%J.%I.err  "/usr/bin/Rscript --vanilla  ', scriptpath, 
                    'util/twincount.R ', winsize, ' ', tmpdirname, '"')
  callout <- system(command, intern = TRUE)
  cat("\n", callout, "\n")
  
  
  # hold job until all subjobs are completed
  partsleft <- parts - length(list.files(path = tmpdirname, pattern = "part"))
  
  while (partsleft > 0) {
    newpartsleft <- parts - length(list.files(path = tmpdirname, pattern = "part"))
    if (newpartsleft != partsleft) {
      cat("\n", newpartsleft, "parts left", as.character(Sys.time()), "\n")
      partsleft <- newpartsleft
    }
    cat("|")
    Sys.sleep(30 * partsleft)
  }
  
  
  # load all jobs back
  cat("\n----\ncollating twin site data...\n")
  timestamp()
  
  reslist <- list()
  
  for (part in 1:parts) {
    filename <- file.path(tmpdirname, paste0("part", part, ".Rdata"))
    cat(filename, "\n")
    if (file.exists(filename)) {
      load(filename)
      reslist <- c(reslist, res)
    }else{
      cat("Check files, not all present!\n")
      break()
    }
  }
  
  if (length(reslist) != nrow(sitelist)) {
    cat("Something is wrong - length of list does not match! reslist is ",
        length(reslist), "sitelist is", nrow(sitelist), "\n")
  }
  
  
  # clean sites
  cat("\n----\ncleaning site data...\n")
  timestamp()
  
  siteclean <- cleansites(sitelist, reslist)

}

cat("\n----\nprepare output...\n")
timestamp()

#  prepare output
siteclean <- siteclean[siteclean$keep, ]

cat("\n----\ndeclutter result -", nrow(siteclean), "/", nrow(sitelist), "\n")

# write output to individual files
for (tag.i in unique(siteclean$tag)) {
  filename.i <- paste(flowcell, tag.i, samplelane, "clean.sites", sep = "_")
  datatable <- siteclean %>%
    filter(tag == tag.i) %>%
    select(-keep, -rn )
  write_tsv(datatable, path = file.path(outpath, filename.i))
}
