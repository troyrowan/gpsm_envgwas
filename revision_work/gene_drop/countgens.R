library(tidyr)
library(dplyr)
library(readr)
library(furrr)
library(purrr)
library(pedigree)

# ran_ped =
#   read_csv("revision_work/gene_drop/RAAAPedigreeFile.csv") %>%
#   select(contains("reg")) %>%
#   distinct()
# 
# #Rename stupidly named columns lol
# colnames(ran_ped) <- c("id", "sire", "dam", "pgs", "pgd", "mgs", "mgd", "sss", "ssd", "sds", "sdd", "dss", "dsd", "dds", "ddd")
# 
# #Make pedigree into a long version
# ran_longped =
#   rbind(select(ran_ped, id, sire, dam),
#         select(ran_ped, id = sire, sire = pgs, dam = pgd),
#         select(ran_ped, id = dam, sire = mgs, dam = mgd),
#         select(ran_ped, id = pgs, sire = sss, dam = ssd),
#         select(ran_ped, id = pgd, sire = sds, dam = sdd),
#         select(ran_ped, id = mgs, sire = dss, dam = dsd),
#         select(ran_ped, id = mgd, sire = dds, dam = ddd)) %>%
#   distinct() %>%
#   mutate(sire = case_when(sire %in% .$id ~ sire,
#                           TRUE ~ 0),
#          dam = case_when(dam %in% .$id ~ dam,
#                          TRUE ~ 0)) %>%
#   mutate(sire = replace_na(sire, 0),
#          dam = replace_na(dam, 0))
ped = 
  read_tsv("ran_longped.txt") %>% pedigree::orderPed()

bind_cols(ped, as.data.frame(countGen(ped))) %>% 
  write_tsv("ran_gennumbers.txt")


id <- 1:5
dam <- c(0,0,1,1,4)
sire <- c(0,0,2,2,3)
ped <- data.frame(id,sire,dam)
countGen(ped)
bind_cols(ped, as.data.frame(countGen(ped)))
