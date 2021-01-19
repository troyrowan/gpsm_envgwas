library(dplyr)
library(tidyr)
library(readr)
library(pedigree)

redangus =
  read_csv("/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus_phenotypes.csv")

RAN = read_csv("All_Animals_SireDam.csv", 
               col_names = c("id", "sire", "dam"),
               skip = 1)

library(pedigree)

tf = 
  RAN %>% 
  mutate(datavec = case_when(id %in% redangus$regisno ~ TRUE,
                             TRUE ~ FALSE)) %>% 
  .$datavec

RAN = pedigree::orderPed(RAN)

trimmed_RAN = trimPed(RAN, tf)

write_csv(trimmed_RAN, "Trimmed_RAN_FullPed.csv")
