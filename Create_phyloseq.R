#####Finalizing a bunch of shit#####

#install.packages('phytools')

library(phyloseq)
library(tidyverse)
library(phangorn)
library(DECIPHER)
library(dada2)
library(microViz)
library(ranacapa)
library(phytools)
library(dada2)

setwd("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023")

big = readRDS("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023/BIG_seqtab.nochim.rds")
otro = readRDS("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023/OTRO_seqtab.nochim.rds")

mergetab <- mergeSequenceTables(big, otro) # unnamed arguments assumed to be sequence tables
mertab_collapse = collapseNoMismatch(mergetab)

taxa <- assignTaxonomy(mertab_collapse, minBoot = 80, "~/Dropbox (UFL)/Laughinghouse_Lab/CyanoSeq/v_1.2.1/CyanoSeq1.2.1_BLCC_SILVA138.1_dada2.fastq.gz", multithread = T, verbose = F)


#Creating PhyloSeq object
ps <- phyloseq(otu_table(mergetab, taxa_are_rows=FALSE), 
               tax_table(taxa))

#Load in metadata
METADATA <- read.csv("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023/data/LakeO_Bac_Meta_w_DBHYDRO_limitation.csv", header = T, row.names = 1)

#Bind metadata to phyloseq object
ps@sam_data <- sample_data(METADATA)

ps_2 <- microbiomeutilities::add_refseq(ps,tag="ASV")

saveRDS(ps_2, "~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/physeq_objects/SEED_total_Apr24.rds")


#####Subsetting Lake O data#####
#So now we have a phyloseq object with only the taxa were interested in
#Samples data is added, the ASVs are stored in the object

#subsetting LAKE-ONLY phyloseq object for the time being
physeq_Lake = subset_samples(ps_2, `WaterBody` == "Lake")
physeq_Lake

saveRDS(physeq_Lake, "~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/physeq_objects/LakeO_SEED_total_Apr24.rds")
physeq_Lake = readRDS("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/physeq_objects/LakeO_SEED_total_Apr24.rds")

#Filter out low abundance taxa
physeq_Lake2 = tax_filter(physeq_Lake,
                          #prev_detection_threshold = .2 , #ASV's must occur 2% of total realtavie abundance
                          min_prevalence = 0.1, #must occur in at least TWO samples we have 49 samples, 2/49 = 0.04
                          min_total_abundance = 100 #THE ASV MUST OCCUR AT LEAST 100 TIMES IN TOTAL ACROSS ALL SAMPLES
)
physeq_Lake2 #This is bac + cyano

#Filtering out all the shit
#This is bac + cyano
physeq = physeq_Lake2 %>%
  subset_taxa((Class!="Chloroplast") | is.na(Order)) %>%
  subset_taxa((Class !="Cyanobacteriia") | is.na(Class)) %>%
  subset_taxa((Kingdom != "Eukaryota") | is.na(Kingdom)) %>%
  subset_taxa((Family != "Mitochondria") | is.na(Family)) %>%
  subset_taxa((Phylum != "NA") | is.na(Phylum)) %>%
  subset_taxa((Kingdom != "Archaea") | is.na(Kingdom))
#Reduction from 60541 to 4048 ASVs
physeq

#There are some weird names i dont like and needed to fix based on the tree, this is doing that lol

# Change taxa name of ASVX to "name"
taxa <- tax_table(physeq)
taxa["ASV1176", "Genus"] <- "Coelosphaerium"
taxa["ASV339", "Genus"] <- "Microcystaceae Cluster 1"
taxa["ASV2762", "Genus"] <- "Microcystaceae Cluster 1"
taxa["ASV751", "Genus"] <- "Microcystaceae Cluster 1"
taxa["ASV2779", "Genus"] <- "Merismopedia"
taxa["ASV2173", "Genus"] <- "Microcystaceae Cluster 2"
taxa["ASV1987", "Genus"] <- "Amphiheterocytum"
taxa["ASV176", "Genus"] <- "Aphanizomenonaceae Cluster 1"
taxa["ASV3311", "Order"] <- "Pseudanabaenales"
taxa["ASV3311", "Family"] <- "Pseudanabaenaceae"
taxa["ASV3311", "Genus"] <- "Pseudanabaena"
taxa["ASV5289", "Family"] <- "Scytonemataceae"
taxa["ASV5289", "Genus"] <- "Neocylindrospermum"
taxa["ASV125", "Order"] <- "Leptolyngbyales"
taxa["ASV2225", "Family"] <- "Leptolyngbyaceae"
taxa["ASV125", "Family"] <- "Leptolyngbyaceae"
taxa["ASV2225", "Genus"] <- "Leptolyngbyaceae Cluster 1"
taxa["ASV125", "Genus"] <- "Leptolyngbyaceae Cluster 1"
taxa["ASV170", "Genus"] <- "Nodosilinea"
taxa["ASV5103", "Genus"] <- "Nodosilinea"
taxa["ASV669", "Genus"] <- "Nodosilinea"
taxa["ASV793", "Genus"] <- "Nodosilinea"
taxa["ASV437", "Genus"] <- "Nodosilinea"
taxa["ASV1838", "Order"] <- "Synechococcales"
taxa["ASV378", "Order"] <- "Synechococcales"
taxa["ASV1838", "Family"] <- "Synechococcaceae"
taxa["ASV378", "Family"] <- "Synechococcaceae"
taxa["ASV1838", "Genus"] <- "Synechococcaceae Cluster 1"
taxa["ASV378", "Genus"] <- "Synechococcaceae Cluster 1"
taxa["ASV5121", "Genus"] <- "Cyanobium"
taxa["ASV5788", "Genus"] <- "Cyanobium"
taxa["ASV1300", "Genus"] <- "Cyanobium"
taxa["ASV2009", "Genus"] <- "Cyanobium"
taxa["ASV400", "Genus"] <- "Cyanobium"
taxa["ASV1423", "Genus"] <- "Cyanobium"
taxa["ASV3299", "Genus"] <- "Cyanobium"
taxa["ASV2713", "Genus"] <- "Cyanobium"
taxa["ASV2519", "Genus"] <- "Cyanobium"
taxa["ASV1716", "Genus"] <- "Cyanobium"
taxa["ASV1173", "Genus"] <- "Cyanobium"
taxa["ASV1604", "Genus"] <- "Cyanobium"
taxa["ASV1973", "Genus"] <- "Cyanobium"
taxa["ASV3527", "Genus"] <- "Cyanobium"
taxa["ASV152", "Genus"] <- "Cyanobium"
taxa["ASV90", "Genus"] <- "Regnicoccus"
taxa["ASV3125", "Genus"] <- "Regnicoccus"
taxa["ASV65", "Genus"] <- "Regnicoccus"
taxa["ASV1425", "Genus"] <- "Vulcanococcus"
taxa["ASV947", "Genus"] <- "Vulcanococcus"
taxa["ASV379", "Genus"] <- "Prochlorococcaceae_XX"
taxa["ASV2807", "Genus"] <- "Prochlorococcaceae_XX"
taxa["ASV2073", "Genus"] <- "Cyanobium"
taxa["ASV51", "Genus"] <- "Cyanobium"
taxa["ASV292", "Genus"] <- "Lacustricoccus"
taxa["ASV683", "Genus"] <- "Lacustricoccus"
taxa["ASV1974", "Genus"] <- "Lacustricoccus"
taxa["ASV41", "Genus"] <- "Lacustricoccus"
taxa["ASV87", "Genus"] <- "Vulcanococcus"
taxa["ASV3259", "Genus"] <- "Vulcanococcus"
taxa["ASV3433", "Genus"] <- "Vulcanococcus"
taxa["ASV2428", "Genus"] <- "Vulcanococcus"
taxa["ASV4820", "Genus"] <- "Vulcanococcus"
taxa["ASV11", "Genus"] <- "Vulcanococcus"
taxa["ASV322", "Genus"] <- "Vulcanococcus"
taxa["ASV1616", "Genus"] <- "Vulcanococcus"
taxa["ASV4592", "Genus"] <- "Vulcanococcus"
taxa["ASV2957", "Genus"] <- "Vulcanococcus"
taxa["ASV20", "Genus"] <- "Vulcanococcus"
taxa["ASV2399", "Genus"] <- "Vulcanococcus"
taxa["ASV67", "Genus"] <- "Vulcanococcus"
taxa["ASV904", "Genus"] <- "Vulcanococcus"
taxa["ASV35", "Genus"] <- "Vulcanococcus"

tax_table(physeq) <- taxa




saveRDS(physeq, "~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/physeq_objects/Physeq_fullApr24.rds")



#Filtering so just Cyanos
physeq_Cyano = physeq %>% subset_taxa(Class == "Cyanophyceae")
x = sample_sums(physeq_Cyano) 
y = taxa_sums(physeq_Cyano)
y %>% mean(.)

fasta = physeq_Lake3@refseq
Biostrings::writeXStringSet(fasta, "~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/ALL_CyanoASV.fasta")


#Read in newick file
x = read.tree(file="~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/cyanoasvtree.txt")
#Bind tree to physeq object
physeq_Cyano@phy_tree <- phy_tree(x)

#Looking at the tree
plot_tree(physeq_Cyano,  color = "Order", label.tips = "taxa_names")


saveRDS(physeq_Cyano, "~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/physeq_objects/Physeq_CyanoApr24.rds")

physeq = readRDS("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/physeq_objects/Physeq_CyanoApr24.rds")


metadata <- read.csv("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023/data/LakeO_Bac_Meta_w_DBHYDRO_limitation.csv", header = T, row.names = 1)  %>% filter(WaterBody == "Lake")


metadata_2 = metadata %>% 
  dplyr::select(-c('SO4_mg.L_X',"CollectionYear","DATE","CollectionDate","DATE","DO_mg.L_X", 'B_mgL', "TN_TP",#"tot_cya",
                   "NA_mg.L_X","K_mg.L_X","CA_mg.L_X","MG_mg.L_X","CHLA.LC_ug.L_X","PC_ugL","Chl_ugL", 
                   "MC_dimethylLR_ugL","MC_980LR_deriv_ugL", "tot_cya", "Salinity_state", "NOX_mg.L_X", "TOTAL_N_mg.L_X", "TN_TP_X",
                   "Anatoxin_a_ugL","Euglenaphycin_ugL",'Anabaenopeptins',
                   'Time', 'Station_Id','TSI_Chla',
                   'TSI_TN','TSI_TP','TSI_TN2','TSI_TP2','SPCOND_uS.cm_x','TURB_NTU_X','T.SUS.SD_mg.L_X','NH4_mg.L_X',
                   'OPO4_mg.L_X','TPO4_mg.L_X','ALKALNTY_mg.L_X','TSI','Cyanopeptide','WaterBody','Secchi_Depth'
  )) %>%
  mutate("DIN" = (NOx_mgL+NH3_mgL)) %>%
  mutate("DIN_DIP" = (NOx_mgL+NH3_mgL)/TRP_mgL)%>%
  mutate_all(funs(ifelse(is.infinite(.), median(., na.rm = TRUE), .)))

physeq@sam_data <- sample_data(metadata_2)

saveRDS(physeq, "~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/physeq_objects/Physeq_CyanoMay1.rds")







