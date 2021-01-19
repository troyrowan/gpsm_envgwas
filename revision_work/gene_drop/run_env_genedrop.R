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
print("Read Phenotypes")

ran_longped = read_csv("/data/tnr343/gpsm_sims/gene_drop/Genotyped_FullPed.csv")
print("Read pedigree")
#This reads in a previously created founder haplotype file
#Same founder haplotypes are used for each iteration
founders = readRDS(founderHaps)
SP = SimParam$new(founders)
SP$setTrackPed(TRUE)
SP$setGender("no")
SP$addSnpChip(nSnpPerChr = 20000)
print("Set new population parameters")
#Performs gene drop based on long pedigree
#I believe AlphaSimR automatically sorts pedigree and assigns founders based on individuals w/ sires/dams = 0
ran_pop = founders
#  newPop(founders)
print("Created new population")

ran_output =
  pedigreeCross(ran_pop,
                id = ran_longped$id,
                father = ran_longped$sire,
                mother = ran_longped$dam)
print("Performed pedigree cross")
#Names SNPs so that we can iterate through chromosomes in subsequent steps
#Columns of "genotypes" dataframe renamed to reflect the chr:pos naming
#There's got to be an easier way to do this
snps = data.frame(
  rs = paste0("SNP_", 1:sum(200000)),
  pos = SP$snpChips[[1]]@lociLoc,
  chr = rep(1:10, each = 20000)) %>%
  mutate(snp = paste(chr, pos, sep = ":"))

# snps = data.frame(
#   rs = paste0("SNP_", 1:sum(c(43436, 37095, 35509, 33168, 33023, 33416, 32360, 28007, 28733, 30171, 31392, 24578, 20758, 22051, 26231, 23527, 21740, 21715, 20480, 21201, 19968, 18296, 16453, 17116, 13922, 15603, 13579, 12913, 14951))),
#   pos = SP$snpChips[[1]]@lociLoc,
#   chr = c(rep(1, SP$snpChips[[1]]@lociPerChr[1]),
#                    rep(2, SP$snpChips[[1]]@lociPerChr[2]),
#                    rep(3, SP$snpChips[[1]]@lociPerChr[3]),
#                    rep(4, SP$snpChips[[1]]@lociPerChr[4]),
#                    rep(5, SP$snpChips[[1]]@lociPerChr[5]),
#                    rep(6, SP$snpChips[[1]]@lociPerChr[6]),
#                    rep(7, SP$snpChips[[1]]@lociPerChr[7]),
#                    rep(8, SP$snpChips[[1]]@lociPerChr[8]),
#                    rep(9, SP$snpChips[[1]]@lociPerChr[9]),
#                    rep(10, SP$snpChips[[1]]@lociPerChr[10]),
#                    rep(11, SP$snpChips[[1]]@lociPerChr[11]),
#                    rep(12, SP$snpChips[[1]]@lociPerChr[12]),
#                    rep(13, SP$snpChips[[1]]@lociPerChr[13]),
#                    rep(14, SP$snpChips[[1]]@lociPerChr[14]),
#                    rep(15, SP$snpChips[[1]]@lociPerChr[15]),
#                    rep(16, SP$snpChips[[1]]@lociPerChr[16]),
#                    rep(17, SP$snpChips[[1]]@lociPerChr[17]),
#                    rep(18, SP$snpChips[[1]]@lociPerChr[18]),
#                    rep(19, SP$snpChips[[1]]@lociPerChr[19]),
#                    rep(20, SP$snpChips[[1]]@lociPerChr[20]),
#                    rep(21, SP$snpChips[[1]]@lociPerChr[21]),
#                    rep(22, SP$snpChips[[1]]@lociPerChr[22]),
#                    rep(23, SP$snpChips[[1]]@lociPerChr[23]),
#                    rep(24, SP$snpChips[[1]]@lociPerChr[24]),
#                    rep(25, SP$snpChips[[1]]@lociPerChr[25]),
#                    rep(26, SP$snpChips[[1]]@lociPerChr[26]),
#                    rep(27, SP$snpChips[[1]]@lociPerChr[27]),
#                    rep(28, SP$snpChips[[1]]@lociPerChr[28]),
#                    rep(29, SP$snpChips[[1]]@lociPerChr[29]))) %>%
#   mutate(snp = paste(chr, pos, sep = ":"))

#Creates annotation file for GEMMA
snps %>%
  select(rs, pos, chr) %>%
  write_csv(paste0("simulations/envgwas/run", run, "/RAN.genedrop.run", run, ".snps.txt"), col_names = FALSE)


#Pulls out the genotypes for individuals who were genotyped in initial Red Angus dataset
#This gets rid of extra genotypes for simulated pedigree indivuals (makes GRM creation easier)
genotypes =
  data.frame(ran_output@id,
             stringsAsFactors = FALSE) %>%
  cbind(ran_longped) %>%
  left_join(read_csv("red_angus_phenotypes.csv") %>%
              select(id = regisno, meantemp, precip, elev)) %>%
  distinct() %>%
  cbind(pullSnpGeno(ran_output)) %>%
  filter(id %in% redangus$regisno)


colnames(genotypes) <- c("ran_output.id", "id", "sire", "dam", "meantemp", "precip", "elev", snps$snp)


genotypes %>%
  select(meantemp, precip, elev) %>%
  write_csv(paste0("simulations/envgwas/run", run, "/RAN.genedrop.run", run, ".env_phenotypes.txt"), col_names = FALSE)

#Genotypes should be the only thing that changes between runs
#Phenotypes, SNP annotations remain the same
#This may be the lone step that we need to break into chromosomes
#Runs multi-threaded (helps w/ inverting big matrix)
# plan(multiprocess)
# options(future.globals.maxSize= 1250000000)
#For 1 in user specificied number of chromosomes
seq(1,10) %>%
  as.character() %>%
  #furrr::future_map(
  purrr::map(
		~select(genotypes, starts_with(paste0(.x, ":"))) %>% #Will select appropriately namedcolumns based on
      t() %>%
      as.data.frame() %>%
      mutate(rs = rownames(.),
             a1 = "a1",
             a2 = "a2") %>%
      select(rs, a1, a2, everything())) %>%
  purrr::reduce(bind_rows) %>%
  write_csv(paste0("simulations/envgwas/run", run, "/RAN.genedrop.run", run, ".mgf"),
            col_names = FALSE)
