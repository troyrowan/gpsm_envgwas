library(AlphaSimR)
library(tidyr)
library(dplyr)
library(readr)
library(furrr)
library(purrr)

args <- commandArgs(TRUE)


founderHaps = args[1]
run = args [2]
#chrom = args[2]

#Reads 3 gen pedigree given to us by Red Angus Association
# ran_ped =
#   read_csv("RAAAPedigreeFile.csv") %>%
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
redangus =
  read_csv("/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus_phenotypes.csv")

gelbvieh =
  read_csv("/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/gelbvieh_phenotypes.csv")

#Pedigrees need to be in long format
ran_longped = read_csv("Genotyped_FullPed.csv")
gel_longped = read_csv("Gelbvieh_Pedigree.csv", col_names = c("iid", "reg", "sire", "dam", "dob"), skip = 1)
print("Read pedigrees")
#This reads in a previously created founder haplotype file
#Same founder haplotypes are used for each iteration
founders = readRDS(founderHaps)
SP = SimParam$new(founders)
SP$setTrackPed(TRUE)
#SP$setGender("no")
SP$addSnpChip(nSnpPerChr = 20000) #Each segregating site should be genotyped.
print("Params Set")

#I believe AlphaSimR automatically sorts pedigree and assigns founders based on individuals w/ sires/dams = 0
#Genotypes for founder population created with this command:
# set.seed(345575433)
# plan(multiprocess, workers = 5) #Multithreading for speedup
# as.character(seq(100,111)) %>%	#Ten different founder populations, named "founderpop100.RDS", etc.
#   furrr::future_map(~runMacs(nInd = 5223, #Number of individuals in pedigree generation zero
#                     nChr = 10, #Ten chromosomes
#                     segSites = 20000, #20K segregating sites per chromosome
#                     species = "CATTLE") %>% #Demographic history of cattle dictates starting Ne and
#                     saveRDS(paste0("founderhaps/founderpop", .x, ".RDS")))

ran_pop =
  newPop(founders)
gel_pop =
  newPop(founders)
print("New population made")

#Pedigree cross is big AlphaSimR function for genedrops
#Allows haplotypes from founder individuals to "drop" down pedigree in same way they would in real world
ran_output =
  pedigreeCross(ran_pop, #Population created above
                id = ran_longped$id, #
                father = ran_longped$sire,
                mother = ran_longped$dam)

gel_output =
  pedigreeCross(gel_pop, #Population created above
                id = gel_longped$reg, #
                father = gel_longped$sire,
                mother = gel_longped$dam)
print("Pedigree cross ran")
#Names SNPs so that we can iterate through chromosomes in subsequent steps
#Columns of "genotypes" dataframe renamed to reflect the chr:pos naming
#There's got to be an easier way to do this
snps = data.frame(
  rs = paste0("SNP_", 1:sum(200000)),
  pos = SP$snpChips[[1]]@lociLoc,
  chr = rep(1:10, each = 20000)) %>%
  mutate(snp = paste(chr, pos, sep = ":"))
print("Created SNP dataframe")

#Creates annotation file for GEMMA
snps %>%
  select(rs, pos, chr) %>%
  write_csv(paste0("simulations/run", run, "/RAN.genedrop.run", run, ".snps.txt"), col_names = FALSE)
snps %>%
  select(rs, pos, chr) %>%
  write_csv(paste0("simulations/run", run, "/GEL.genedrop.run", run, ".snps.txt"), col_names = FALSE)
print("wrote SNP File")

#Pulls out the genotypes for individuals who were genotyped in initial Red Angus dataset
#This gets rid of extra genotypes for simulated pedigree indivuals (makes GRM creation easier)
ran_genotypes =
  data.frame(ran_output@id,
             stringsAsFactors = FALSE) %>%
  cbind(ran_longped) %>%
  left_join(read_csv("red_angus_phenotypes.csv") %>%
              select(id = regisno, age)) %>%
  distinct() %>%
  cbind(pullSnpGeno(ran_output)) %>%
  filter(id %in% redangus$regisno)
print("Created genotype ped dataframe")

colnames(ran_genotypes) <- c("ran_output.id", "id", "sire", "dam", "age", snps$snp)

#Writes out phenotype file for GEMMA
ran_genotypes %>%
  select(age) %>%
  write_csv(paste0("simulations/run", run, "/RAN.genedrop.run", run, ".phenotypes.txt"), col_names = FALSE)
gel_genotypes %>%
  select(age) %>%
  write_csv(paste0("simulations/run", run, "/GEL.genedrop.run", run, ".phenotypes.txt"), col_names = FALSE)

print("Wrote phenotypes")
#Genotypes should be the only thing that changes between runs
#Phenotypes, SNP annotations remain the same
#This may be the lone step that we need to break into chromosomes
#Runs multi-threaded (helps w/ inverting big matrix)
# plan(multiprocess)
# options(future.globals.maxSize= 1250000000)
#For 1 in user specificied number of chromosomes

#Writing the genotype file in MGF (mean genotype file) format for GEMMA input
seq(1,10) %>%
  as.character() %>%
  #furrr::future_map(
  purrr::map(
                ~select(ran_genotypes, starts_with(paste0(.x, ":"))) %>% #Will select appropriately namedcolumns based on
      t() %>%
      as.data.frame() %>%
      mutate(rs = rownames(.),
             a1 = "a1",
             a2 = "a2") %>%
      select(rs, a1, a2, everything())) %>%
  purrr::reduce(bind_rows) %>%
  write_csv(gzfile(paste0("simulations/run", run, "/RAN.genedrop.run", run, ".mgf.gz")),
            col_names = FALSE)

seq(1,10) %>%
  as.character() %>%
  #furrr::future_map(
  purrr::map(
    ~select(gel_genotypes, starts_with(paste0(.x, ":"))) %>% #Will select appropriately namedcolumns based on
      t() %>%
      as.data.frame() %>%
      mutate(rs = rownames(.),
             a1 = "a1",
             a2 = "a2") %>%
      select(rs, a1, a2, everything())) %>%
  purrr::reduce(bind_rows) %>%
  write_csv(gzfile(paste0("simulations/run", run, "/GEL.genedrop.run", run, ".mgf.gz")),
            col_names = FALSE)
