library(AlphaSimR)
library(dplyr)
library(tidyr)
library(readr)
library(furrr)

args = commandArgs(trailingOnly=TRUE)

set.seed(345575433)
plan(multiprocess, workers = 5)
as.character(seq(100,111)) %>%
  furrr::future_map(~runMacs(nInd = 5223,
						nChr = 10,
						segSites = 20000,
						species = "CATTLE") %>%
						saveRDS(paste0("founderhaps/founderpop", .x, ".RDS")))



#chrom = as.numeric(args[1])
# RAN_snpcounts =
#   read_table("RAN.report.frq",
#            col_types = cols(SNP = col_character())) %>%
#   filter(MAF != 1, MAF != 0) %>%
#   group_by(CHR) %>%
#   count() %>%
#   mutate(snps = round(n/50))

# founderHaps =
#   runMacs(nInd = 7601,
#           nChr = 29,
#           segSites = c(43436, 37095, 35509, 33168, 33023, 33416, 32360, 28007, 28733, 30171, 31392, 24578, 20758, 22051, 26231, 23527, 21740, 21715, 20480, 21201, 19968, 18296, 16453, 17116, 13922, 15603, 13579, 12913, 14951),
#           species = "CATTLE")

# saveRDS(founderHaps, "founderhaps/RAN_founderHaps.truesnps.7601.rds")
