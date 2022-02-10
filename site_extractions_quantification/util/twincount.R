#@ utility script launched from 05_decluttertwins.R.
#@ script identifies similarity in sequence using agrepl, identifing cases where the same integration sites were mapped to >1 positions due to small changes in sequence.

# get variables from system
commandline <- commandArgs(trailingOnly = TRUE)


# load required libraries and functions
library(tidyverse)
library(parallel)

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

findtwinssection <- function(x){
  i <- unique(x$rn)
  
  if (i >= max(sitelist$rn[sitelist$keep])) {return(NULL)}
  
  twins <- findtwins(sitelist, i) %>%
    mutate(decl = case_when(
      tag.x == tag.y ~ TRUE,                 # classical twins
      grepl("JKT|jurkat|control", ptcode.x, ignore.case = TRUE) ~ TRUE,   # always declut against controls
      grepl("JKT|jurkat|control", ptcode.y, ignore.case = TRUE) ~ TRUE,
      tag.x != tag.y & ptcode.x != ptcode.y ~ TRUE,                       # declut when dif patients
      tag.x != tag.y & ptcode.x == ptcode.y ~ FALSE,                      # do not declut when same patients
      TRUE ~ NA
    )) %>%
    filter(decl)
  
  return(twins)
} 


# set variables
winsize <- as.numeric(commandline[1])
tmpdirused <- commandline[2]
part <- as.numeric(Sys.getenv('LSB_JOBINDEX'))

filename <- file.path(tmpdirused, paste0("part", part, ".Rdata"))


# load required data
load(file.path(tmpdirused, "all.Rdata"))
ls()


# set parameters 
fromnum <- (part - 1) * winsize + 1
tonum <- ifelse(part * winsize > nrow(sitelist), nrow(sitelist), part * winsize)


# Report to user
cat("\n----\ncleaning site data...\n")
timestamp()

cat("\nfile name is", filename, "\n")
cat("from: ", fromnum, "\tto: ", tonum, "\n")


# Process filepart
cl <- makeCluster(4)

clusterEvalQ(cl, {library(tidyverse)})
clusterExport(cl, varlist = c("sitelist", "findtwins"))

res <- parLapply(cl, split(sitelist, sitelist$rn)[fromnum:tonum], findtwinssection)

stopCluster(cl)


save(res, file = filename)





# projectpath <- "~/lmpcr/"
# datapath <- file.path(projectpath, "data/05_extract_sites/original/HIV1")
# metapath <- file.path(projectpath, "metadata")
# projectpath <- "~/Dropbox/Projects/integrationsite/"
# datapath <- file.path(projectpath, "imresult/05_extract_sites/original/HIV1")
# metapath <- file.path(projectpath, "data/metadata")
# library(profvis)
# library(parallel)
# library(microbenchmark)

# flowcell <- "2014R5"
# samplelane <- "L002"
# virus <- gsub(".+\\/", "", datapath)

# load metadata
# metabyfc <- read_tsv(file.path(metapath, "allmeta.txt"), 
#                      col_types = "cccciciclicccccccdd") %>%
#   filter(seqexpid == flowcell, lane == as.numeric(gsub("L00", "", samplelane))) 
# barcodes <- read_tsv(file.path(metapath, "barcodes.txt"), col_types = c("cci"))
# 

# # choose files to read
# allsitefiles <- list.files(path = datapath, pattern = paste0(flowcell, ".+", samplelane, "\\.sites"))
# samplestodo <- grep(samplelane, allsitefiles, value = TRUE)

# cat("datapath is ", datapath, "\n")
# cat(" flowcell is ", flowcell, "\n")
# cat(" lane is ", samplelane, "\n")
# cat(" there are ", length(samplestodo), 
#     "samples, including", samplestodo[1], "\n")


# load all sites from lane to one table
# sitelist <- list()
# for (i in 1:length(samplestodo)) {
#   cat(".")
#   filename <- samplestodo[i]
#   tmpobj <- read_tsv(file.path(datapath, filename),
#                      col_types = c("cdciiiiiiic"))
#   sitelist[[filename]] <- tmpobj
# }
# cat("\n")
# 
# sitelist <- bind_rows(sitelist, .id = "filename") %>%
#   separate(filename, into = c(NA, "tag", "lane", NA)) %>%
#   arrange(desc(totalshear), desc(totaldupl)) %>%
#   left_join(select(metabyfc, ptcode, index), by = c("tag" = "index")) %>%
#   mutate(rn0 = row_number(), keep = TRUE) 
# 

# sitelist <- sitelist[1:50000, ]
# 
# dfsmall <- sitelist %>% mutate(rn = row_number())




# 
# cleansites <- function(sitestoclean, alltwins = NULL, part = NULL){
#   cat("\n")
#   
#   for (i in 1:max(sitestoclean$rn)) {
#     if (i %% 100 == 0) {cat(".")}
#     if (i %% 10000 == 0) {cat("\n")}
#     # if(i == 151){cat(i, "\n"); return(sitestoclean); break}
#     if (!i %in% sitestoclean$rn[sitestoclean$keep]) {next}
#     if (i >= max(sitestoclean$rn[sitestoclean$keep])) {break}
#     
#     if (is.null(alltwins)) {
#       twins <- findtwins(sitestoclean, i) %>%
#         mutate(decl = case_when(
#           tag.x == tag.y ~ TRUE,                 # classical twins
#           grepl("JKT|jurkat|control", ptcode.x, ignore.case = TRUE) ~ TRUE,   # always declut against controls
#           grepl("JKT|jurkat|control", ptcode.y, ignore.case = TRUE) ~ TRUE, 
#           tag.x != tag.y & ptcode.x != ptcode.y ~ TRUE,                       # declut when dif patients
#           tag.x != tag.y & ptcode.x == ptcode.y ~ FALSE,                      # do not declut when same patients
#           TRUE ~ NA
#         )) %>%
#         filter(decl)
#     }else{
#       twins <- alltwins[alltwins$V1 == i &
#                           alltwins$V2 %in% sitestoclean$rn[sitestoclean$keep], ]
#     }
#     
#     removetwins <- 
#       unique(c(twins$V1[twins$totaldupl.x == twins$totaldupl.y & 
#                           twins$totalshear.x == twins$totalshear.y], 
#                twins$V2))
#     
#     if (length(removetwins) != 0) {
#       sitestoclean$keep[sitestoclean$rn %in% removetwins] <- FALSE
#     }
#     
#     # if (!sitestoclean$keep[sitestoclean$rn == 2736]){cat(i, " 2736 excluded")}
#     
#   }
#   
#   return(sitestoclean)
# }
# 
# loopfindtwins <- function(sitestoclean, part = NULL, wins = 5000){
#   
#   fromnum <- (part - 1) * wins + 1
#   tonum <- part * wins
#   cat("from: ", fromnum, "\tto: ", tonum, "\n")
#   
#   twins <- list()
#   
#   for (i in fromnum:tonum) {
#     j <- i - (part - 1) * wins
#     if (i %% 100 == 0) {cat(".")}
#     if (i %% 10000 == 0) {cat("\n")}
#     # if(i == 151){cat(i, "\n"); return(sitestoclean); break}
#     if (!i %in% sitestoclean$rn[sitestoclean$keep]) {next}
#     if (i >= max(sitestoclean$rn[sitestoclean$keep])) {break}
#     
#     twins[[j]] <- findtwins(sitestoclean, i) %>%
#       mutate(decl = case_when(
#         tag.x == tag.y ~ TRUE,                 # classical twins
#         grepl("JKT|jurkat|control", ptcode.x, ignore.case = TRUE) ~ TRUE,   # always declut against controls
#         grepl("JKT|jurkat|control", ptcode.y, ignore.case = TRUE) ~ TRUE, 
#         tag.x != tag.y & ptcode.x != ptcode.y ~ TRUE,                       # declut when dif patients
#         tag.x != tag.y & ptcode.x == ptcode.y ~ FALSE,                      # do not declut when same patients
#         TRUE ~ NA
#       )) %>%
#       filter(decl)
#   }
#   return(twins)
# }
# 
