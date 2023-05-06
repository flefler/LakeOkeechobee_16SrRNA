#####LAKE OKEECHOBEE BACTERIAl primers#####
#https://schuyler-smith.github.io/phylosmith/index.html
#https://david-barnett.github.io/microViz/reference/ps_calc_dominant.html
#https://rdrr.io/github/microbiome/microbiome/man/core_members.html

######## PACKAGES & SUCH ###########

#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
#install.packages('ggthemes', dependencies = TRUE)
#install.packages('gg_ordiplot')
#install.packages("remotes")
#remotes::install_github("microbiome/microbiome")
#remotes::install_github("gauravsk/ranacapa")
#devtools::install_github("david-barnett/microViz")
#install.packages("remotes")
#remotes::install_github("jfq3/ggordiplots")
#install.packages("GUniFrac")
#devtools::install_github('microsud/microbiomeutilities')
#install.packages("eulerr")
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(patchwork)
library(pairwiseAdonis)
library(microbiomeutilities)
library(GUniFrac)
library(eulerr)
library(microViz)
library(gclus)
library(ggrepel)
library(ggforce)
library(tibble)
library(BiodiversityR)# also loads vegan
library(ranacapa)
library(microbiome)
library(tidyverse)
library(phyloseq)
library(ggsci)
library(ggthemes)
library(ade4)
library(adegraphics)
library(adespatial)
library(vegan3d)
library(MASS)
library(ellipse)
library(FactoMineR)
library(rrcov)
library(missMDA)
library(Amelia)
library(BiodiversityR)
library(ggordiplots)
library(indicspecies)







#####Determining periods of nutrient limitation#####
library(tidyverse)
library(lubridate)

#read in report
REPORT <- read.csv("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023/LakeO_Sites_NP.csv")

#Remove those weird samples
report = REPORT %>% filter(Sample != "FCEB")

#Filter down to appropriate times
x1 = with(report, report[(Date >= "20190620" & Date <= "20201102"),])

#Remove the NAs
report1 = x1 %>% drop_na(TOTAL.N)

report2 = report1 %>% select(-Sample)

#Convert to numeric
report2.5 = report2 |> filter(!is.na(as.numeric(TOTAL.N))) 
report2.75 = report2.5 |> filter(!is.na(as.numeric(TPO4))) 
report2.75$TOTAL.N = as.numeric(report2.75$TOTAL.N)
report2.75$TPO4 = as.numeric(report2.75$TPO4)

#Make new variable
report3 = report2.75 %>% mutate(TNTP = TOTAL.N / TPO4) 

#Convert to date
report3$Date <- ymd(report3$Date)

#Add limitation data
report4 = report3 %>%
  mutate(limitation = case_when(TNTP <= 9 ~ 'N',
                                TNTP > 9 & TNTP < 23 ~ 'Co',
                                TNTP >= 23 ~ 'P'))


#Plot
plot <- ggplot(data=report4, aes(x=Date, y=TNTP, group = Station, colour = Station)) + 
  geom_line() +
  labs(y= "TN:TP", x = "Date") +
  geom_point() +
  facet_wrap( ~ Station, ncol=3) 
plot 


#These are sample dates!!!
dates_vline <- as.Date(c("2019-08-20", "2019-09-17","2019-10-23","2019-12-03","2020-01-14","2020-02-19","2020-04-07","2020-06-09","2020-07-15","2020-09-02"))


#add to plot
plot + 
  geom_vline(xintercept = as.numeric(dates_vline), col = "black", lwd = 0.2) + 
  geom_hline(yintercept = c(9,23), col = "red", lwd = 0.5)


#Restart R and clear environment
rm(list=ls())
.rs.restartR()



#####Load phyloseq object and metadata#####

setwd("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023")

#Reads in total physeq object
physeq_tot = readRDS("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/physeq_objects/Physeq_fullApr24.rds") 
#x = sample_sums(physeq_Lake2) 
#y = taxa_sums(physeq)
#x %>% mean(.)

#Reads in Cyano physeq object
physeq_bar = readRDS("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/physeq_objects/Physeq_CyanoMay1.rds") 
#x = sample_sums(physeq_Lake2) 
#y = taxa_sums(physeq)
#x %>% mean(.)


#These samples suck
physeq <- prune_samples(sample_names(physeq_bar) != "LO2", physeq_bar)
physeq <- prune_samples(sample_names(physeq) != "LO6", physeq)


physeq_alpha = readRDS("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023/physeq_objects/LakeO_Cyano_forFD_Jan12.rds")


#View tree, very importatnt to ensure that the tree isnt fucky and polyphletic. 
#This WILL alter your results, assuming you are leveraging these data

#p_p = physeq %>% subset_taxa(Family =="Prochlorococcaceae")
#plot_tree(p_p,  color = "Genus", label.tips = "taxa_names")
#plot_tree(physeq_bac,  color = "Phylum")

metadata_magic <- read.csv("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023/data/LakeO_Bac_Meta_w_DBHYDRO_limitation.csv", header = T, row.names = 1) %>% tibble::rownames_to_column(var = "SampleID") %>% 
  filter(WaterBody == "Lake") %>% dplyr::filter(SampleID != "LO2") %>% dplyr::filter(SampleID != "LO6")

metadata_magic_L <- read.csv("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023/data/LakeO_Bac_Meta_w_DBHYDRO_limitation.csv", header = T, row.names = 1) %>% tibble::rownames_to_column(var = "SampleID") %>% 
  filter(WaterBody == "Lake") %>% dplyr::filter(SampleID != "LO2") %>% dplyr::filter(SampleID != "LO6") %>% filter(Limitation != "P")



metadata <- read.csv("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023/data/LakeO_Bac_Meta_w_DBHYDRO_limitation.csv", header = T, row.names = 1) %>% tibble::rownames_to_column(var = "SampleID") %>% filter(WaterBody == "Lake")
#metadata = as.data.frame(physeq@sam_data)

#These keep reading in as characters, idk why
metadata$TSI = as.numeric(metadata$TSI)
metadata$TN_TP = as.numeric(metadata$TN_TP)
metadata$tot_cya = as.numeric(metadata$tot_cya)
metadata$TN_TP_X = as.numeric(metadata$TN_TP_X)
metadata$Outing = as.factor(metadata$Outing)
metadata$Site = as.factor(metadata$Site)
metadata$Location = as.factor(metadata$Location)
metadata$Limitation = as.character(metadata$Limitation)

physeq@sam_data$Outing = as.factor(physeq@sam_data$Outing)
physeq@sam_data$Location = as.factor(physeq@sam_data$Location)

#For the RDA we cannot have NA values, nor should we have any negative values (Based on a priori knowledge of this particular data set)
#Because of this, we are going to fill in NA values with median based on location AND replaces negative values with 0



#This replaces NA with median based on location AND replaces negative values with 0, and then takes sqrt of the final variables

metadata_RDA = metadata %>% 
  dplyr::select(-c('SO4_mg.L_X',"CollectionYear","DATE","CollectionDate","DATE","DO_mg.L_X", 'B_mgL', "TN_TP",#"tot_cya",
                   "NA_mg.L_X","K_mg.L_X","CA_mg.L_X","MG_mg.L_X","CHLA.LC_ug.L_X","PC_ugL","Chl_ugL", 
                   "MCLR_ugL", "MCRR_ugL", "MC_dimethylLR_ugL","MC_980LR_deriv_ugL","Nodularin_811_ugL",
                   "Anatoxin_a_ugL","Euglenaphycin_ugL",'Anabaenopeptins',
                   'Time', 'Station_Id','TSI_Chla',
                   'TSI_TN','TSI_TP','TSI_TN2','TSI_TP2','SPCOND_uS.cm_x','TURB_NTU_X','T.SUS.SD_mg.L_X','NH4_mg.L_X',
                   'OPO4_mg.L_X','TPO4_mg.L_X','ALKALNTY_mg.L_X','TSI','Cyanopeptide','WaterBody','Secchi_Depth'
  )) %>%
  #data.table::setDT() %>%
  group_by(Site) %>%
  mutate_all(funs(ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  mutate_all(funs(ifelse(. < 0, 0, .))) %>% 
  mutate_all(funs(ifelse(is.infinite(.), median(., na.rm = TRUE), .))) %>%
  mutate("DIN" = (NOx_mgL+NH3_mgL)) %>%
  mutate("DIN_DIP" = (NOx_mgL+NH3_mgL)/TRP_mgL) %>%
  mutate_all(funs(ifelse(is.infinite(.), median(., na.rm = TRUE), .))) %>%
  dplyr::select(-c('NOX_mg.L_X','TOTAL_N_mg.L_X',"TN_mgL","TN_TP_X","Salinity_state","tot_cya",
                   "DO.","Sample","TRP_mgL", "TP_mgL", "TN_mgL","NOx_mgL","NO2_mgL","NO3_mgL","NH3_mgL","Depth_m")) %>%
  mutate_at(c("Photic_Depth","Temp","DO_mgL","pH","Conductivity_SPC","Turbidity_FNU","OP_mgL",
              "Al_mgL","Ca_mgL","Cu_mgL","Fe_mgL","K_mgL","Mg_mgL","Mn_mgL","Na_mgL","Zn_mgL",
              "DIN","DIN_DIP"), sqrt) %>%
  as.data.frame()

#This is a secret tool that will help us later
#metadata_RDA1 = metadata_RDA %>% column_to_rownames("SampleID")
#physeq1@sam_data <- sample_data(metadata_RDA1)


metadata_RDA$Site = as.factor(metadata_RDA$Site)
metadata_RDA$Limitation = as.factor(metadata_RDA$Limitation)

#Based on unifrac beta diversity there appears to be a difference between wet/dry season and between certain locations.
#Were going to make a subset of wet, dry, and perhaps each location and outing, but need to confirm with beta

physeq_dry = subset_samples(physeq, `Season` == "Dry")
physeq_wet = subset_samples(physeq, `Season` == "Wet")
physeq_N = subset_samples(physeq, `Limitation` == "N")
physeq_Co = subset_samples(physeq, `Limitation` == "Co")

#Subset limitations sans P-limitation
metadata_L = metadata_RDA %>% filter(Limitation != "P")
physeq_L = subset_samples(physeq, `Limitation` != "P")
physeq_L_BAR = subset_samples(physeq_bar, `Limitation` != "P")
physeq_alpha_L = subset_samples(physeq_alpha, `Limitation` != "P")


#p = metadata %>% dplyr::filter(Cyanopeptide == "Y") %>% dplyr::filter(WaterBody == "Lake") %>% dplyr::select("Outing","Site","MCLR_ugL","MCRR_ugL","MC_dimethylLR_ugL","MC_980LR_deriv_ugL","Nodularin_811_ugL","Anatoxin_a_ugL","Anabaenopeptins","Season")
#write_csv(p, "~/Dropbox (UFL)/Lefler_Dissertation/LakeO_Phylo/toxins.csv")



#####Export for network analysis#####
#This is kind of shit show and far from perfect code, but it works. 
physeq_tot = readRDS("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023/physeq_objects/LakeO_Feb8.rds") 



#Fill in taxa names with lowest classified rank
physeq2 = physeq_tot %>% tax_fix() %>% tax_glom(taxrank="Genus", NArm=F)

c = as.data.frame(physeq2@otu_table)

#Filter out low abundance taxa
physeq_Lake22 = microViz::tax_filter(physeq2,
                                     #min_prevalence = , #must occur in at least TWO samples we have 49 samples, 2/49 = 0.04, not using this since so many diverse environments 
                                     min_total_abundance = 408 #THE ASV MUST OCCUR AT LEAST 100 TIMES IN TOTAL ACROSS ALL SAMPLES
) #This manages to be the top 500 ASVs
physeq_Lake22

gfns = as.data.frame(physeq_Lake22@tax_table) #%>% rownames_to_column() %>% dplyr::select(c("Phylum","Genus")) 

gg = as.data.frame(physeq_Lake22@otu_table) %>% sjmisc::rotate_df()

gg$abundance <- rowSums(gg)
ggg = gg %>% dplyr::select(abundance)

ff = merge(gfns, ggg, by=0)

fff = ff %>% dplyr::select(c(Phylum, Genus, abundance))

write_csv(fff,"~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023/data/ELSA_500_taxonomy.csv")



#Stolen code lol but basically makes your names in "Dolichospermum ASV123"
# generate a vector containing the full taxonomy path for all OTUs
wholetax <- do.call(paste, c(as.data.frame(tax_table(physeq_Lake22))
                             [("Genus")], 
                             sep = "_"))  # to distinguish from "_" within tax ranks

# turn the otu_table into a data.frame
otu_export <- as.data.frame(otu_get(physeq_Lake22))
tmp <- names(otu_export)

# paste wholetax and OTU_ids together
for(i in 1:length(tmp)){
  names(tmp)[i] = paste(wholetax[i], sep = "_")
}

# overwrite old names
names(otu_export) <- names(tmp)

head(otu_export)[4]

thing = otu_export %>% tibble::rownames_to_column(var = "SampleID") %>% arrange(SampleID) 

otu_export = thing %>% column_to_rownames(var="SampleID")

dune.Hellinger_otu_export <- BiodiversityR::disttransform(otu_export, method='hellinger') 

pp = dune.Hellinger_otu_export %>% tibble::rownames_to_column(var = "SampleID")


data_frame_merge_3 <- merge(pp, metadata_RDA,
                            by = 'SampleID', all = TRUE) 

#Makes the row names into what ELSA requires
rownames(data_frame_merge_3) <- paste0("t",data_frame_merge_3$Outing,"r",data_frame_merge_3$Location) 

data_frame_merge_3 = data_frame_merge_3 %>% dplyr::select(-c(Outing,Limitation,Location,Site,Season,DO_mgL,SampleID)) %>% sjmisc::rotate_df()

#If youre missing samples youll need to create a column for the missing samples otherwise ELSA wont run
data3 = data_frame_merge_3 %>%
  mutate(t3r4 = c(0)) %>%
  mutate(t10r1 = c(0)) %>%
  mutate(t10r3 = c(0)) %>%
  rownames_to_column()

#Ensure that you have NO CHARACTERS 

#Write as a tsv
write_tsv(data3, col_names = T, "~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023/ELSA_500_Feb8.txt")





#####BARPLOTS#####

#Phylum
phylum = physeq_tot %>% tax_fix() %>%
  ps_mutate(subject_timepoint = interaction(Outing)) %>%
  ps_arrange(Outing, desc(Outing)) %>%
  comp_barplot(tax_level = "Phylum", n_taxa = 10, facet_by = "Site", nrow = 1, sample_order = "default", label = "subject_timepoint") +
  labs(x = "Outing", 
       y = "Relative Abundance") 

Order = physeq_bar %>% tax_fix() %>%
  ps_mutate(subject_timepoint = interaction(Outing)) %>%
  ps_arrange(Outing, desc(Outing)) %>%
  comp_barplot(tax_level = "Order", n_taxa = 10, facet_by = "Site", nrow = 1, sample_order = "default", label = "subject_timepoint") +
  labs(x = "Outing", 
       y = "Relative Abundance") 

Family = physeq %>% tax_fix() %>%
  ps_mutate(subject_timepoint = interaction(Outing)) %>%
  ps_arrange(Outing, desc(Outing)) %>%
  comp_barplot(tax_level = "Family", n_taxa = 9, facet_by = "Site", nrow = 1, sample_order = "default", label = "subject_timepoint") +
  labs(x = "Outing", 
       y = "Relative Abundance") 

library(patchwork)

(phylum / Order / Family) + plot_annotation(tag_levels = 'A')

#Cyano Genus
Bac_Cyano_Genus = physeq %>% tax_fix() %>%
  ps_mutate(subject_timepoint = interaction(Outing)) %>%
  ps_arrange(Outing, desc(Outing)) %>%
  comp_barplot(tax_level = "Genus", n_taxa = 19, facet_by = "Site", nrow = 1, sample_order = "default", label = "subject_timepoint") +
  labs(x = "Outing", 
       y = "Relative Abundance")


#
#####Season#####

#Faiths PD
dd = metadata %>% filter(WaterBody == "Lake")
FD = picante::pd(physeq_alpha@otu_table, physeq_alpha@phy_tree, include.root=TRUE) %>% tibble::rownames_to_column(var  = "SampleID") %>% full_join(., dd, by = "SampleID" )%>% filter(WaterBody == "Lake")
FD$Outing = as.factor(FD$Outing)

FD_Season = FD %>% 
  ggplot(aes(x = Season, y = PD, fill = Season)) +
  geom_violin(width=1) +
  geom_boxplot(width=0.2, alpha=0.5) +
  #geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme_bw() +
  labs(title = "Faith's PD",
       x = "Season",
       y = "Diversity Measure")

res_aov <- aov(FD$PD ~ FD$Season)
summary(res_aov) #0.84




#Chao
chao = richness(physeq_alpha, index = "chao1", detection = 0) %>% tibble::rownames_to_column(var = "SampleID") %>% full_join(.,dd, by = "SampleID" ) %>% filter(WaterBody == "Lake")
chao$Season = as.factor(chao$Season)


C_Season = chao %>% 
  ggplot(aes(x = Season, y = chao1, fill = Season, group = Season)) +
  geom_violin(width=1) +
  geom_boxplot(width=0.2, alpha=0.5) +
  #geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme_bw() +
  labs(title = "Chao1",
       x = "Season",
       y = "")

res_aov1 <- aov(chao$chao1 ~ chao$Season)
summary(res_aov1) #0.75





#OLD
#z = plot_richness(physeq_alpha, x = "Season",  measures = c("Chao1", "Simpson")) + geom_violin(width=1.4) +
#  geom_boxplot(width=0.2, alpha=0.5) + theme_minimal()


#Barplot, ok so like 60% of this is pico
a = physeq_bar %>% tax_fix() %>%
  ps_select(Site, Season) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Season") %>%
  comp_barplot(tax_level = "Family", n_taxa = 10, bar_width = 0.8) +
  coord_flip()

#Lets remove the Pico
physeq_cyanobacteria2 = physeq_bar %>%
  subset_taxa(Family !="Prochlorococcaceae") 

physeq_cyanobacteria3 = physeq %>%
  subset_taxa(Family =="Prochlorococcaceae")

#And try again
b = physeq_cyanobacteria2 %>% tax_fix() %>%
  ps_select(Site, Season) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Season") %>%
  comp_barplot(tax_level = "Genus", n_taxa = 10, bar_width = 0.8) +
  coord_flip()

#Lets just peek at the peek-o
physeq_cyanobacteria3 %>% tax_fix() %>%
  ps_select(Site, Season) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Season") %>%
  comp_barplot(tax_level = "unique", n_taxa = 19, bar_width = 0.8) +
  coord_flip()

physeq_cyanobacteria2 %>% tax_fix() %>%
  ps_select(Site, Season) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Season") %>%
  comp_barplot(tax_level = "Family", n_taxa = 10, bar_width = 0.8) +
  coord_flip()


#If we do log10 transform, it looks like this. This is more significant than hellinger, althought they are both < 0.05
c = physeq %>% tax_fix() %>%
  tax_transform("log10", rank = "unique") %>% 
  dist_calc(dist = "gunifrac", gunifrac_alpha = 0.5) %>%
  ord_calc(method = "NMDS") %>%
  ord_plot(colour = "Season", size = 2, auto_caption = NA) +
  stat_ellipse(aes(colour = Season), type = "norm") 

x <- physeq %>%  tax_fix() %>%
  tax_transform("log10", rank = "unique") %>% 
  dist_calc(dist = "gunifrac", gunifrac_alpha = 0.5)

# the more permutations you request, the longer it takes
# but also the more stable and precise your p-values become
w <- x %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 9999, # you should use at least 999!
    variables = "Season"
  )

# view the permanova results
perm_get(w) %>% as.data.frame()
#         Df  SumOfSqs         R2        F Pr(>F)
#Season    1 0.03811361 0.06347494 3.592186 0.0069
#Residual 53 0.56233774 0.93652506       NA     NA
#Total    54 0.60045135 1.00000000       NA     NA




#If we do hellinger transform, it looks like this
c = physeq %>% tax_fix() %>%
  tax_transform("hellinger", rank = "unique") %>% 
  dist_calc(dist = "gunifrac", gunifrac_alpha = 0.5) %>%
  ord_calc(method = "NMDS") %>%
  ord_plot(colour = "Season", size = 2, auto_caption = NA) +
  stat_ellipse(aes(colour = Season), type = "norm") +
  theme(legend.position = c(0.15, 0.87),
   legend.background = element_rect(fill = "white", color = "black"))
c

x1 <- physeq %>%  tax_fix() %>%
  tax_transform("hellinger", rank = "unique") %>% 
  dist_calc(dist = "gunifrac", gunifrac_alpha = 0.5)
#> Proportional min_prevalence given: 0.1 --> min 23/222 samples.

# the more permutations you request, the longer it takes
# but also the more stable and precise your p-values become
w1 <- x1 %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 9999, # you should use at least 999!
    variables = "Season"
  )

# view the permanova results
perm_get(w1) %>% as.data.frame()
#         Df  SumOfSqs         R2        F Pr(>F)
#Season    1 0.04053127 0.0414625 2.37908  0.034
#Residual 55 0.93700918 0.9585375      NA     NA
#Total    56 0.97754045 1.0000000      NA     NA



library(patchwork)
qq = (FD_Season | C_Season | c) + plot_layout(guides = 'collect')
(qq) / b  + plot_annotation(tag_levels = 'A')




#https://cran.r-project.org/web/packages/indicspecies/vignettes/IndicatorSpeciesAnalysis.html
#install.packages("indicspecies")
library(indicspecies)
library(tidyverse)

genus <- physeq %>% tax_glom(taxrank="Genus", NArm=T) %>% tax_fix() %>%
  tax_transform("log10", rank = "Genus")
genus_abund = as.data.frame(genus@otu_table) %>% arrange()
metadata_season1 = metadata_magic %>% arrange(SampleID)
genus_abund_mp = indicspecies::multipatt(genus_abund, metadata_season1$Season, func = "r.g", control = how(nperm=9999))
summary(genus_abund_mp)

Group Dry  #sps.  13 
                              stat p.value    
  Microcystaceae Cluster 2   0.470  0.0003 ***
  Aphanizomenon              0.462  0.0007 ***
  Woronichinia               0.444  0.0013 ** 
  Regnicoccus                0.441  0.0012 ** 
  Lacustricoccus             0.439  0.0006 ***
  Prochlorococcaceae_XX      0.421  0.0014 ** 
  Raphidiopsis               0.383  0.0046 ** 
  Synechocystis              0.380  0.0058 ** 
  Metis                      0.345  0.0128 *  
  Pegethrix                  0.316  0.0235 *  
  Cuspidothrix               0.312  0.0216 *  
  Leptolyngbyaceae Cluster 1 0.303  0.0285 *  
  Nodosilinea                0.276  0.0447 *  
  
  Group Wet  #sps.  3 
stat p.value    
  Synechococcaceae Cluster 1 0.664  0.0001 ***
  Apatinema                  0.365  0.0062 ** 
  Neocylindrospermum         0.301  0.0228 *  



#
#####Limitation#####

#Faiths PD
dd_L = metadata_L 
FD_L = picante::pd(physeq_alpha_L@otu_table, physeq_alpha_L@phy_tree, include.root=TRUE) %>% tibble::rownames_to_column(var  = "SampleID") %>% full_join(., dd_L, by = "SampleID" )
FD_L$Outing = as.factor(FD_L$Outing)

FD_L_plot = FD_L %>% 
  ggplot(aes(x = Limitation, y = PD, fill = Limitation)) +
  geom_violin(width=1) +
  geom_boxplot(width=0.2, alpha=0.5) +
  #geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme_bw() +
  labs(title = "Faith's PD",
       x = "Limitation",
       y = "Diversity Measure")

res_aov <- aov(FD_L$PD ~ FD_L$Limitation)
summary(res_aov)




#Chao
chao = richness(physeq_alpha_L, index = "chao1", detection = 0) %>% tibble::rownames_to_column(var = "SampleID") %>% full_join(.,dd_L, by = "SampleID" )
chao$Limitation = as.factor(chao$Limitation)


C_L = chao %>% 
  ggplot(aes(x = Limitation, y = chao1, fill = Limitation, group = Limitation)) +
  geom_violin(width=1) +
  geom_boxplot(width=0.2, alpha=0.5) +
  #geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme_bw() +
  labs(title = "Chao1",
       x = "Limitation",
       y = "")

res_aov1 <- aov(chao$chao1 ~ chao$Limitation)
summary(res_aov1)


#Barplot, ok so like 50% of this is pico
a2 = physeq_L_BAR %>% tax_fix() %>%
  ps_select(Site, Limitation) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Limitation") %>%
  comp_barplot(tax_level = "Family", n_taxa = 10, bar_width = 0.8) +
  coord_flip()

#Lets remove the Pico
physeq_L_cyanobacteria2 = physeq_L_BAR %>%
  subset_taxa(Family !="Prochlorococcaceae") 

physeq_L_cyanobacteria3 = physeq_L_BAR %>%
  subset_taxa(Family =="Prochlorococcaceae")

#And try again
b = physeq_L_cyanobacteria2 %>% tax_fix() %>%
  ps_select(Site, Limitation) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Limitation") %>%
  comp_barplot(tax_level = "Genus", n_taxa = 10, bar_width = 0.8) +
  coord_flip()

b2 = physeq_L_BAR %>% tax_fix() %>%
  ps_select(Site, Limitation) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Limitation") %>%
  comp_barplot(tax_level = "Family", n_taxa = 10, bar_width = 0.8) +
  coord_flip()

(a / a2) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')

#Lets just peek at the peek-o
physeq_L_cyanobacteria3 %>% tax_fix() %>%
  ps_select(Site, Limitation) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Limitation") %>%
  comp_barplot(tax_level = "Genus", n_taxa = 30, bar_width = 0.8) +
  coord_flip()

physeq %>% tax_fix() %>%
  ps_select(Site, Limitation) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Limitation") %>%
  comp_barplot(tax_level = "Genus", n_taxa = 10, bar_width = 0.8) +
  coord_flip()


#If we do log10 transform, it looks like this.
c = physeq_L %>% tax_fix() %>%
  tax_transform("log10", rank = "unique") %>% 
  dist_calc(dist = "gunifrac", gunifrac_alpha = 0.5) %>%
  ord_calc(method = "NMDS") %>%
  ord_plot(colour = "Limitation", size = 2, auto_caption = NA) +
  stat_ellipse(aes(colour = Limitation), type = "norm") 

x <- physeq_L %>%  tax_fix() %>%
  tax_transform("log10", rank = "unique") %>% 
  dist_calc(dist = "gunifrac", gunifrac_alpha = 0.5)
  
w <- x %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 99999, # you should use at least 999!
    variables = c("Limitation")
  )

# view the permanova results
perm_get(w) %>% as.data.frame()
#           Df   SumOfSqs         R2        F Pr(>F)
#Limitation  1 0.01036068 0.01788898 0.9471711 0.4204
#Residual   52 0.56880472 0.98211102        NA     NA
#Total      53 0.57916540 1.00000000        NA     NA

library(patchwork)

qq = (FD_L_plot | C_L) + plot_layout(guides = 'collect')
(qq | c) / b  + plot_annotation(tag_levels = 'A')

(a / a2)+ plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')


#https://cran.r-project.org/web/packages/indicspecies/vignettes/IndicatorSpeciesAnalysis.html
#install.packages("indicspecies")
library(indicspecies)

genus <- physeq_L %>% tax_glom(taxrank="Genus", NArm=F) %>% tax_fix() %>%
  tax_transform("hellinger", rank = "Genus")
genus_abund = as.data.frame(genus@otu_table)
metadata_season1_L = metadata_magic_L %>% arrange(SampleID)
genus_abund_mp = indicspecies::multipatt(genus_abund, metadata_season1_L$Limitation, func = "r.g", control = how(nperm=9999))
summary(genus_abund_mp)

  
#  Group N 
#              stat p.value   
#Synechococcaceae Cluster 1 0.403  0.0069 **
#Planktothrix               0.309  0.0271 * 


#
#####Site#####

#Barplot, ok so like 50% of this is pico
a2 = physeq_bar %>% tax_fix() %>%
  ps_select(Site, Limitation) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Site") %>%
  comp_barplot(tax_level = "Genus", n_taxa = 15, bar_width = 0.8) +
  coord_flip()

#Faiths PD
dd = metadata %>% filter(WaterBody == "Lake")
FD = picante::pd(physeq_alpha@otu_table, physeq_alpha@phy_tree, include.root=TRUE) %>% tibble::rownames_to_column(var  = "SampleID") %>% full_join(., dd, by = "SampleID" )%>% filter(WaterBody == "Lake")
FD$Outing = as.factor(FD$Outing)

FD_Site = FD %>% 
  ggplot(aes(x = Site, y = PD, fill = Site)) +
  geom_violin(width=1) +
  geom_boxplot(width=0.2, alpha=0.5) +
  #geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title = "Faith's PD",
       x = "Season",
       y = "Diversity Measure")

res_aov <- aov(FD$PD ~ FD$Site)
summary(res_aov)
TukeyHSD(res_aov, conf.level=.95) 



#Chao
chao = richness(physeq_alpha, index = "chao1", detection = 0) %>% tibble::rownames_to_column(var = "SampleID") %>% full_join(.,dd, by = "SampleID" ) %>% filter(WaterBody == "Lake")
chao$Season = as.factor(chao$Season)


C_Season = chao %>% 
  ggplot(aes(x = Site, y = chao1, fill = Site, group = Site)) +
  geom_violin(width=1) +
  geom_boxplot(width=0.2, alpha=0.5) +
  #geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title = "Chao1",
       x = "Site",
       y = "")

res_aov1 <- aov(chao$chao1 ~ chao$Site)
summary(res_aov1)
TukeyHSD(res_aov1, conf.level=.95) 




#z = plot_richness(physeq, x = "Site",  measures = c("Chao1", "Simpson")) + geom_violin(width=1.4) +
#  geom_boxplot(width=0.2, alpha=0.5) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


#If we do log10 transform, it looks like this. This is more significant than hellinger, althought they are both < 0.05
c = physeq %>% tax_fix() %>%
  tax_transform("log10", rank = "unique") %>% 
  dist_calc(dist = "gunifrac", gunifrac_alpha = 0.5) %>%
  ord_calc(method = "NMDS") %>%
  ord_plot(colour = "Site", size = 2) +
  stat_ellipse(aes(colour = Site)) 

x <- physeq %>%  tax_fix() %>%
  tax_transform("log10", rank = "unique") %>% 
  dist_calc(dist = "gunifrac", gunifrac_alpha = 0.5)

w <- x %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 9999, # you should use at least 999!
    variables = "Site"
  )

b = dist_get(x)
b

metadata4 = metadata %>% arrange(SampleID) %>% column_to_rownames(var="SampleID")


pairwise_site = pairwiseAdonis::pairwise.adonis(b, metadata4[,"Site"],perm = 99999, p.adjust.m = "holm")
pairwise_site %>% filter(p.adjusted <= 0.05)
#NA




library(patchwork)

(FD_Site | C_Season | c) / a2 + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')




#https://cran.r-project.org/web/packages/indicspecies/vignettes/IndicatorSpeciesAnalysis.html
#install.packages("indicspecies")
library(indicspecies)

genus <- physeq %>% tax_glom(taxrank="Genus", NArm=T) %>% tax_fix() %>%
  tax_transform("log10", rank = "Genus")
genus_abund = as.data.frame(genus@otu_table)
Impact_Site = metadata$Site
genus_abund_mp = multipatt(genus_abund, metadata$Site, duleg = TRUE, func = "r.g", control = how(nperm=9999))
summary(genus_abund_mp)

#Group Clewiston  #sps.  1 
#           stat p.value    
#Apatinema 0.643   1e-04 ***
  
#Group Eagle Bay  #sps.  5 
#                       stat p.value    
#Metis                 0.377  0.0187 *
#Prochlorococcaceae_XX 0.364  0.0302 *
#Snowella              0.360  0.0378 *
#Aphanizomenon         0.358  0.0212 *
#Microcystaceae Family 0.345  0.0256 *
  
#Group Kissimmee Mouth  #sps.  2 
#                   stat p.value    
#Radiocystis       0.427  0.0133 *
#Planktothricoides 0.347  0.0326 *
  
#Group South Lake  #sps.  3 
#                       stat  p.value    
#Prochlorococcus        0.353  0.0294 *
#Prochlorococcaceae_XXX 0.315  0.0379 *
#Parasynechococcus      0.302  0.0443 *



#db-pRDA
physeq %>% tax_fix() %>%
  tax_transform("hellinger", rank = "unique") %>% 
  dist_calc(dist = "gunifrac", gunifrac_alpha = 0.5) %>%
  ord_calc(
    #constraints = c("DIN_DIP"),
    conditions = c("Outing", "Location"),
    method = "CAP", 
    scale_cc = F 
  ) %>%
  ord_plot(
    colour = "Season", size = 2, 
    plot_taxa = 1:8
  ) +
  stat_ellipse(aes(colour = Season))






#####
#####Toxin correlations######
physeq_tox = subset_samples(physeq_bar, `MCLR_ugL` >= "0") 


physeq_tox2 = physeq_tox %>%
  subset_taxa(Family != "Prochlorococcaceae")

physeq_tox2 %>% tax_fix() %>%
  tax_transform("log10", zero_replace = "halfmin") %>% 
  tax_agg("Genus") %>% 
  cor_heatmap(vars = c("MCLR_ugL", "MCRR_ugL","Nodularin_811_ugL"), cor = "kendall", tax_anno = taxAnnotation(
    Prev. = anno_tax_prev(ylim = 0:1)))

#####RDA#####
#This code runs the RDA
#
#Make sure youre using the correct metadata, with NA's filled in

#A quick note
#There is a lot of code just produce one plot, and several objects are made.
#Because of this, each iteration of this code is IDENTICAL
#the final plot is given a specific name


#Transform data and glom as genus
x <- physeq %>% tax_glom(taxrank="Genus", NArm=T) %>% tax_fix() %>%
  tax_transform("log10", rank = "Genus")

#Pull otu table as dataframe
x2 = as.data.frame(x@otu_table)

#Arrange based on identifier
x3 = x2 %>% tibble::rownames_to_column() %>% arrange(rowname)

#Subsample metadata, only needed if subsampling ASVs
metadata4 = metadata_RDA %>% filter(SampleID %in% x3$rowname)#Only needed if subsampling metadata
#Protist_Order = Protist_Order %>% filter(rowname %in% metadata4$SampleID)#Only needed if subsampling metadata

#Arrange based on identifier
metadata4 = metadata4 %>% arrange(SampleID)  %>% dplyr::select(-c("Limitation","Season","Location"))

#Make the identifiers the row names.
#THE ORDER OF THE ROW NAMES MUST BE IDENTICAL BETWEEN THE TWO DATAFRAMES

metadata5 = metadata4 %>% column_to_rownames(var="SampleID") 
x4 = x3 %>% column_to_rownames(var="rowname")

#Remove the junk variables

rm(list = c("x", "x2", "x3", "metadata4"))


###

Ordination.model2 <- rda(x4 ~ ., data = metadata5, scaling="species") 
summary(Ordination.model2)
#Partitioning of variance:
#              Inertia Proportion
#Total          11.182     1.0000
#Conditioned     1.324     0.1184
#Constrained     4.161     0.3721
#Unconstrained   5.697     0.5095

anova.cca(Ordination.model2, step = 1000) #p = 0.21
plot2 <- ordiplot(Ordination.model2, choices=c(1,2))

#These six lines of code make it possible to plot species data in ggplot
sites.long2 <- sites.long(plot2, env.data=metadata5)
head(sites.long2)



species.long2 <- species.long(plot2)
axis.long2 <- axis.long(Ordination.model2, choices=c(1, 2))



#These  4 lines of code add env data into ordination plot.
spec.envfit2 <- envfit(plot2, env=metadata5)
spec.data.envfit2 <- data.frame(r=spec.envfit2$vectors$r, p=spec.envfit2$vectors$pvals, spec.envfit2$vectors$arrows)
spec.data.envfit2$FG <- rownames(spec.data.envfit2)
Env_Arrow <-filter(spec.data.envfit2, p <= .1)
Env_Arrow


spec.envfit <- envfit(plot2, env=x4)
spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)
species.long2 <- species.long(plot2, spec.data=spec.data.envfit)
species.long2

#Removes taxa that don't make up at least 50% of community variance, i.e. drive diversity
species.long3 <- species.long2[species.long2$r >= .65, ]
species.long3


w1 = ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +   
  geom_point(data=sites.long2, 
             aes(x=axis1, y=axis2, colour= Site), 
             size=3, position = position_jitter(.05)) + labs(color = "Site", shape = "Site") +
  geom_segment(data=species.long3, 
               aes(x=0, y=0, xend=axis1*1, yend=axis2*1), 
               colour="blue", size=0.3, arrow=arrow()) +
  geom_text_repel(data=species.long3, 
                  aes(x=axis1*1, y=axis2*1, label=labels),
                  colour="blue") + 
  geom_segment(data = Env_Arrow, mapping = aes(x=0, y=0, xend=RDA1*1, yend=RDA2*1),
               size = .5, color = "red", arrow=arrow(length=unit(2, "cm")*Env_Arrow$r)) +
  geom_text_repel(data = Env_Arrow, aes(x=RDA1, y=RDA2, label= FG), colour = "red")+
  #ggforce::geom_mark_ellipse(data=sites.long2, aes(x=axis1, y=axis2, colour=Limitation)) +  
  theme_minimal() + 
  coord_fixed(ratio=1) 


Ordination.model2 <- rda(x4 ~ ., data = metadata5, scaling="species") 
(fit <- envfit(Ordination.model2, metadata5, perm = 999))

metadata6 = metadata5 %>% dplyr::select(-c("Site","Outing"))


x = bioenv(x5, metadata6, method = "spearman", partial = c("Site", "Outing"),
       metric = "euclidean",
       parallel = 6)
#Photic_Depth DO_mgL OP_mgL Cu_mgL
summary(x)

#####CDMV_RDA#####

#This is just a general genus level RDA
#Transform data and glom as genus
x <- physeq %>% tax_glom(taxrank="Genus", NArm=T) %>% tax_fix() %>%
  tax_transform("log10", rank = "Genus")

#Pull otu table as dataframe
x2 = as.data.frame(x@otu_table)

#Arrange based on identifier
x3 = x2 %>% tibble::rownames_to_column() %>% arrange(rowname)

#Subsample metadata, only needed if subsampling ASVs
metadata4 = metadata_RDA %>% filter(SampleID %in% x3$rowname)#Only needed if subsampling metadata
#Protist_Order = Protist_Order %>% filter(rowname %in% metadata4$SampleID)#Only needed if subsampling metadata

#Arrange based on identifier
metadata4 = metadata4 %>% arrange(SampleID)  %>% dplyr::select(-c("Limitation","Season","Location"))

#Make the identifiers the row names.
#THE ORDER OF THE ROW NAMES MUST BE IDENTICAL BETWEEN THE TWO DATAFRAMES

metadata5 = metadata4 %>% column_to_rownames(var="SampleID") 
x4 = x3 %>% column_to_rownames(var="rowname")

#Remove the junk variables

rm(list = c("x", "x2", "x3", "metadata4"))
#GENUS RDA
x5 = x4 %>% dplyr::select(c("Cyanobium","Dolichospermum","Microcystis", "Cuspidothrix","Vulcanococcus", "Raphidiopsis"))
metadata6  = metadata5 %>% dplyr::select(c("Outing","Site","pH","DO_mgL","Photic_Depth","OP_mgL","OP_mgL","Cu_mgL", "Fe_mgL","K_mgL","DIN", "Zn_mgL","DIN_DIP"))


Ordination.model2 <- rda(x5 ~ . + Condition(Outing+Site), data = metadata6, scaling="species") 
summary(Ordination.model2)
#Partitioning of variance:
#              Inertia Proportion
#Total          1.9446     1.0000
#Conditioned    0.2149     0.1105
#Constrained    0.6684     0.3437
#Unconstrained  1.0612     0.5457

anova.cca(Ordination.model2, step = 1000) 

#Df Variance      F Pr(>F)   
#Model     9  0.64966 2.6067  0.001 **
#Residual 39  1.07999  

plot2 <- ordiplot(Ordination.model2, choices=c(1,2))

#These six lines of code make it possible to plot species data in ggplot
sites.long2 <- sites.long(plot2, env.data=metadata6)
head(sites.long2)

species.long2 <- species.long(plot2)
axis.long2 <- axis.long(Ordination.model2, choices=c(1, 2))



#These  4 lines of code add env data into ordination plot.
spec.envfit2 <- envfit(plot2, env=metadata6)
spec.data.envfit2 <- data.frame(r=spec.envfit2$vectors$r, p=spec.envfit2$vectors$pvals, spec.envfit2$vectors$arrows)
spec.data.envfit2$FG <- rownames(spec.data.envfit2)
Env_Arrow <-filter(spec.data.envfit2, p <= 0.99)
Env_Arrow


spec.envfit <- envfit(plot2, env=x5)
spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)
species.long2 <- species.long(plot2, spec.data=spec.data.envfit)
species.long2


w2 = ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +   
  geom_point(data=sites.long2, 
             aes(x=axis1, y=axis2, colour= Site), 
             size=3, position = position_jitter(.05)) + labs(color = "Site", shape = "Site") +
  geom_segment(data=species.long2, 
               aes(x=0, y=0, xend=axis1*1, yend=axis2*1), 
               colour="blue", size=0.3, arrow=arrow()) +
  geom_text_repel(data=species.long2, 
                  aes(x=axis1*1, y=axis2*1, label=labels),
                  colour="blue") + 
  geom_segment(data = Env_Arrow, 
               aes(x=0, y=0, xend=RDA1*1, yend=RDA2*1),
               size = .5, color = "red", arrow=arrow(length=unit(2, "cm")*Env_Arrow$r)) +
  geom_text_repel(data = Env_Arrow, aes(x=RDA1, y=RDA2, label= FG), colour = "red")+
  #ggforce::geom_mark_ellipse(data=sites.long2, aes(x=axis1, y=axis2, colour=HAT)) +  
  theme_minimal() + 
  coord_fixed(ratio=1) 
w2

(w1 | w2) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')
















#####RDA Season#####

#Transform data and glom as genus
x <- physeq %>% tax_glom(taxrank="Genus", NArm=T) %>% tax_fix() %>%
  tax_transform("log10", rank = "Genus")

#Pull otu table as dataframe
x2 = as.data.frame(x@otu_table)

#Arrange based on identifier
x3 = x2 %>% tibble::rownames_to_column() %>% arrange(rowname)

#Subsample metadata, only needed if subsampling ASVs
metadata4 = metadata_RDA %>% filter(SampleID %in% x3$rowname)#Only needed if subsampling metadata
#Protist_Order = Protist_Order %>% filter(rowname %in% metadata4$SampleID)#Only needed if subsampling metadata

#Arrange based on identifier
metadata4 = metadata4 %>% arrange(SampleID) %>% dplyr::select(-c("Limitation","Location"))

#Make the identifiers the row names.
#THE ORDER OF THE ROW NAMES MUST BE IDENTICAL BETWEEN THE TWO DATAFRAMES

metadata5 = metadata4 %>% column_to_rownames(var="SampleID")
x9 = x3 %>% column_to_rownames(var="rowname")

#Remove the junk variables

rm(list = c("x", "x2", "x3", "metadata4"))



res.man = manova(cbind(Photic_Depth,Temp,DO_mgL,pH,Conductivity_SPC,Salinity_PPT,Turbidity_FNU,
                         OP_mgL,Al_mgL,Ca_mgL,Cu_mgL,Fe_mgL,K_mgL,Mg_mgL,Mn_mgL,Na_mgL,Zn_mgL,DIN,DIN_DIP) ~ Season,
                   data = metadata5)


summary(res.man,test = c("Wilks"))


summary.aov(res.man)

metadata_xyz = metadata5 %>% dplyr::select(c("Outing", "Site", "Season", "Temp", "DO_mgL", "pH", "OP_mgL", 
                                             "Al_mgL", "Cu_mgL", "Fe_mgL", "K_mgL", "Mg_mgL", "Mn_mgL","Na_mgL", "Zn_mgL", "DIN_DIP"))


Ordination.model2 <- rda(x9 ~ .+ Condition(Outing + Site + Season), data = metadata_xyz, scaling="species") 
summary(Ordination.model2)

#Partitioning of variance:
#  Inertia Proportion
#Total          11.182     1.0000
#Conditioned     2.016     0.1803
#Constrained     2.577     0.2305
#Unconstrained   6.589     0.5893

plot2 <- ordiplot(Ordination.model2, choices=c(1,2))
anova.cca(Ordination.model2, step = 1000) #p = 0.341

#These six lines of code make it possible to plot species data in ggplot
sites.long2 <- sites.long(plot2, env.data=metadata_xyz)
head(sites.long2)

species.long2 <- species.long(plot2)
axis.long2 <- axis.long(Ordination.model2, choices=c(1, 2))



#These  4 lines of code add env data into ordination plot.
spec.envfit2 <- envfit(plot2, env=metadata_xyz)
spec.data.envfit2 <- data.frame(r=spec.envfit2$vectors$r, p=spec.envfit2$vectors$pvals, spec.envfit2$vectors$arrows)
spec.data.envfit2$FG <- rownames(spec.data.envfit2)
Env_Arrow <-filter(spec.data.envfit2, p <= 1)
Env_Arrow


spec.envfit <- envfit(plot2, env=x9)
spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)
species.long2 <- species.long(plot2, spec.data=spec.data.envfit)
species.long2


#Removes taxa that don't make up at least 65% of community variance, i.e. drive diversity
species.long3 <- species.long2[species.long2$r >= .65, ] 
species.long3


season = ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +   
  geom_point(data=sites.long2, 
             aes(x=axis1, y=axis2, colour= Site), 
             size=3, position = position_jitter(.05)) + labs(color = "Site", shape = "Site") +
  geom_segment(data=species.long2, 
               aes(x=0, y=0, xend=axis1*1, yend=axis2*1), 
               colour="blue", size=0.3, arrow=arrow()) +
  geom_text_repel(data=species.long2, 
                  aes(x=axis1*1, y=axis2*1, label=labels),
                  colour="blue") + 
  geom_segment(data = Env_Arrow, mapping = aes(x=0, y=0, xend=RDA1*1, yend=RDA2*1),
               size = .5, color = "red", arrow=arrow(length=unit(2, "cm")*Env_Arrow$r)) +
  geom_text_repel(data = Env_Arrow, aes(x=RDA1, y=RDA2, label= FG), colour = "red")+
  #ggforce::geom_mark_ellipse(data=sites.long2, aes(x=axis1, y=axis2, colour=Season)) +  
  theme_minimal() + 
  coord_fixed(ratio=1) 
season


rm(list = c("metadata_xyz","res.man","species.long3", "species.long2", "sites.long2", "spec.data.envfit", "Ordination.model2", "plot2", "axis.long2", "spec.envfit2","spec.data.envfit2","Env_Arrow", "spec.envfit"))





#####RDA Site#####

#Transform data and glom as genus
x <- physeq %>% tax_glom(taxrank="Genus", NArm=T) %>% tax_fix() %>%
  tax_transform("log10", rank = "Genus")

#Pull otu table as dataframe
x2 = as.data.frame(x@otu_table)

#Arrange based on identifier
x3 = x2 %>% tibble::rownames_to_column() %>% arrange(rowname)

#Subsample metadata, only needed if subsampling ASVs
metadata4 = metadata_RDA %>% filter(SampleID %in% x3$rowname)#Only needed if subsampling metadata
#Protist_Order = Protist_Order %>% filter(rowname %in% metadata4$SampleID)#Only needed if subsampling metadata

#Arrange based on identifier
metadata4 = metadata4 %>% arrange(SampleID) %>% dplyr::select(-c("Limitation","Location"))

#Make the identifiers the row names.
#THE ORDER OF THE ROW NAMES MUST BE IDENTICAL BETWEEN THE TWO DATAFRAMES

metadata5 = metadata4 %>% column_to_rownames(var="SampleID")
x9 = x3 %>% column_to_rownames(var="rowname")

#Remove the junk variables

rm(list = c("x", "x2", "x3", "metadata4"))



res.man = manova(cbind(Photic_Depth,Temp,DO_mgL,pH,Conductivity_SPC,Salinity_PPT,Turbidity_FNU,
                       OP_mgL,Al_mgL,Ca_mgL,Cu_mgL,Fe_mgL,K_mgL,Mg_mgL,Mn_mgL,Na_mgL,Zn_mgL,DIN,DIN_DIP) ~ Site,
                 data = metadata5)


summary(res.man,test = c("Wilks"))

summary.aov(res.man)


TukeyHSD(x)

x = aov(DIN_DIP ~ Site,  data = metadata5)
summary.aov(x)
TukeyHSD(x)


metadata_RDA %>% 
  ggplot(aes(x = Site, y = DIN_DIP, fill = Site)) +
  geom_boxplot(width=0.2, ) +
  #geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme_bw()








metadata_xyz = metadata5 %>% dplyr::select(c("Outing", "Site", "Season", "Photic_Depth", "DO_mgL", "pH", "Salinity_PPT", 
                                             "Turbidity_FNU", "Ca_mgL", "Fe_mgL", "K_mgL", "Mg_mgL", "Na_mgL", "DIN" ,"DIN_DIP"))


Ordination.model2 <- rda(x9 ~ .+ Condition(Outing + Site + Season), data = metadata_xyz, scaling="species") 
summary(Ordination.model2)

#Partitioning of variance:
#  Inertia Proportion
#Total          11.182     1.0000
#Conditioned     2.016     0.1803
#Constrained     2.577     0.2305
#Unconstrained   6.589     0.5893

plot2 <- ordiplot(Ordination.model2, choices=c(1,2))
anova.cca(Ordination.model2, step = 1000) #p = 0.341

#These six lines of code make it possible to plot species data in ggplot
sites.long2 <- sites.long(plot2, env.data=metadata_xyz)
head(sites.long2)

species.long2 <- species.long(plot2)
axis.long2 <- axis.long(Ordination.model2, choices=c(1, 2))



#These  4 lines of code add env data into ordination plot.
spec.envfit2 <- envfit(plot2, env=metadata_xyz)
spec.data.envfit2 <- data.frame(r=spec.envfit2$vectors$r, p=spec.envfit2$vectors$pvals, spec.envfit2$vectors$arrows)
spec.data.envfit2$FG <- rownames(spec.data.envfit2)
Env_Arrow <-filter(spec.data.envfit2, p <= 0.1)
Env_Arrow


spec.envfit <- envfit(plot2, env=x9)
spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)
species.long2 <- species.long(plot2, spec.data=spec.data.envfit)
species.long2


#Removes taxa that don't make up at least 65% of community variance, i.e. drive diversity
species.long3 <- species.long2[species.long2$r >= .65, ] 
species.long3


site_RDA = ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +   
  geom_point(data=sites.long2, 
             aes(x=axis1, y=axis2, colour= Site), 
             size=3, position = position_jitter(.05)) + labs(color = "Site", shape = "Site") +
  geom_segment(data=species.long3, 
               aes(x=0, y=0, xend=axis1*1, yend=axis2*1), 
               colour="blue", size=0.3, arrow=arrow()) +
  geom_text_repel(data=species.long3, 
                  aes(x=axis1*1, y=axis2*1, label=labels),
                  colour="blue") + 
  geom_segment(data = spec.data.envfit2, mapping = aes(x=0, y=0, xend=RDA1*1, yend=RDA2*1),
               size = .5, color = "red", arrow=arrow(length=unit(2, "cm")*Env_Arrow$r)) +
  geom_text_repel(data = Env_Arrow, aes(x=RDA1, y=RDA2, label= FG), colour = "red")+
  #ggforce::geom_mark_ellipse(data=sites.long2, aes(x=axis1, y=axis2, colour=Site)) +  
  theme_minimal() + 
  coord_fixed(ratio=1) 
site_RDA


rm(list = c("metadata_xyz","res.man","species.long3", "species.long2", "sites.long2", "spec.data.envfit", "Ordination.model2", "plot2", "axis.long2", "spec.envfit2","spec.data.envfit2","Env_Arrow", "spec.envfit"))






#####RDA for limitation#####
#Transform data and glom as genus
x <- physeq_L %>% tax_glom(taxrank="Genus", NArm=T) %>% tax_fix() %>%
  tax_transform("log10", rank = "Genus")

#Pull otu table as dataframe
x2 = as.data.frame(x@otu_table)

#Arrange based on identifier
x3 = x2 %>% tibble::rownames_to_column() %>% arrange(rowname)

#Subsample metadata, only needed if subsampling ASVs
metadata4 = metadata_L %>% filter(SampleID %in% x3$rowname)#Only needed if subsampling metadata
#Protist_Order = Protist_Order %>% filter(rowname %in% metadata4$SampleID)#Only needed if subsampling metadata

#Arrange based on identifier
metadata4 = metadata4 %>% arrange(SampleID)  %>% dplyr::select(-c("Season","Location"))

#Make the identifiers the row names.
#THE ORDER OF THE ROW NAMES MUST BE IDENTICAL BETWEEN THE TWO DATAFRAMES

metadata5 = metadata4 %>% column_to_rownames(var="SampleID") 
x4 = x3 %>% column_to_rownames(var="rowname")

#Remove the junk variables

rm(list = c("x", "x2", "x3", "metadata4"))

res.man <-  manova(cbind(Photic_Depth,Temp,DO_mgL,pH,Conductivity_SPC,Salinity_PPT,Turbidity_FNU,
                                     OP_mgL,Al_mgL,Ca_mgL,Cu_mgL,Fe_mgL,K_mgL,Mg_mgL,Mn_mgL,Na_mgL,Zn_mgL,DIN,DIN_DIP) ~ Limitation,
                   data = metadata5)


summary(res.man,test = c("Wilks"))


summary.aov(res.man)

metadata_xyz = metadata5 %>% dplyr::select(c("Outing", "Site", "Limitation", "Photic_Depth", "Salinity_PPT", "Turbidity_FNU", 
                                             "Ca_mgL", "Mg_mgL", "Na_mgL", "DIN"))

Ordination.model2 <- rda(x4 ~ . + Condition(Outing + Site + Limitation), data = metadata_xyz, scaling="species") 
summary(Ordination.model2)
#Partitioning of variance:
#              Inertia Proportion
#Total          11.090     1.0000
#Conditioned     1.462     0.1318
#Constrained     4.355     0.3927
#Unconstrained   5.273     0.4755

anova.cca(Ordination.model2, step = 1000) #p = 0.405
plot2 <- ordiplot(Ordination.model2, choices=c(1,2))

#These six lines of code make it possible to plot species data in ggplot
sites.long2 <- sites.long(plot2, env.data=metadata_xyz)
head(sites.long2)



species.long2 <- species.long(plot2)
axis.long2 <- axis.long(Ordination.model2, choices=c(1, 2))



#These  4 lines of code add env data into ordination plot.
spec.envfit2 <- envfit(plot2, env=metadata_xyz)
spec.data.envfit2 <- data.frame(r=spec.envfit2$vectors$r, p=spec.envfit2$vectors$pvals, spec.envfit2$vectors$arrows)
spec.data.envfit2$FG <- rownames(spec.data.envfit2)
Env_Arrow <-filter(spec.data.envfit2, p <= .1)
Env_Arrow


spec.envfit <- envfit(plot2, env=x4)
spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)
species.long2 <- species.long(plot2, spec.data=spec.data.envfit)
species.long2

#Removes taxa that don't make up at least 50% of community variance, i.e. drive diversity
species.long3 <- species.long2[species.long2$r >= .65, ]
species.long3


w_lim = ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +   
  geom_point(data=sites.long2, 
             aes(x=axis1, y=axis2, colour= Site), 
             size=3, position = position_jitter(.05)) + labs(color = "Site", shape = "Site") +
  geom_segment(data=species.long3, 
               aes(x=0, y=0, xend=axis1*1, yend=axis2*1), 
               colour="blue", size=0.3, arrow=arrow()) +
  geom_text_repel(data=species.long3, 
                  aes(x=axis1*1, y=axis2*1, label=labels),
                  colour="blue") + 
  geom_segment(data = Env_Arrow, mapping = aes(x=0, y=0, xend=RDA1*1, yend=RDA2*1),
               size = .5, color = "red", arrow=arrow(length=unit(2, "cm")*Env_Arrow$r)) +
  geom_text_repel(data = Env_Arrow, aes(x=RDA1, y=RDA2, label= FG), colour = "red")+
  ggforce::geom_mark_ellipse(data=sites.long2, aes(x=axis1, y=axis2, colour=Limitation)) +  
  theme_minimal() + 
  coord_fixed(ratio=1) 
w_lim


rm(list = c("metadata_xyz","res.man","species.long3", "species.long2", "sites.long2", "spec.data.envfit", "Ordination.model2", "plot2", "axis.long2", "spec.envfit2","spec.data.envfit2","Env_Arrow", "spec.envfit"))


#####ExportASVs#####


f= physeq_Lake3 %>% tax_fix() %>% subset_taxa(Family != "Prochlorococcaceae") %>% subset_taxa(Family != "Aphanizomenonaceae")  %>% subset_taxa(Family != "Microcystaceae" ) %>%
  tax_select(tax_list = c("Family", "Order","Class"), n_typos = 0, ranks_searched = "Genus")

v = as.data.frame(f@tax_table)

fasta = f@refseq
Genbank_accession = names(fasta)
Sequence = paste(fasta)
df = data.frame(Genbank_accession,Sequence)

iop = as.data.frame(f@tax_table) %>% dplyr::select(Genus) %>% tibble::rownames_to_column("Genbank_accession")
df2 = left_join(iop, df) %>% unite(NAME, c("Genbank_accession", "Genus"))

seqs = Biostrings::DNAStringSet(df2$Sequence)
names(seqs) = paste(df2$NAME)
# Saving the sequences as a fastq file
Biostrings::writeXStringSet(seqs, "~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/CyanoASV.fasta")






#####Exporting for GAMs####


physeq_tot = readRDS("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/physeq_objects/Physeq_fullApr24.rds") 

physeq_rare = phyloseq::rarefy_even_depth(physeq_tot, rngseed = 42069, replace=FALSE)
x = as.data.frame(phyloseq::sample_sums(physeq_rare)) 

physeq_cyanobacteria_rare = physeq_rare %>%
  phyloseq::subset_taxa(Class =="Cyanophyceae") 

metadata <- read.csv("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023/data/LakeO_Bac_Meta_w_DBHYDRO_limitation.csv", header = T, row.names = 1) %>% tibble::rownames_to_column(var = "SampleID") %>% filter(WaterBody == "Lake")

metadata_1 = metadata %>% 
  dplyr::select(-c("Site",'SO4_mg.L_X',"CollectionYear","DATE","CollectionDate","DATE","DO_mg.L_X", 'B_mgL', "TN_TP",#"tot_cya",
                   "NA_mg.L_X","K_mg.L_X","CA_mg.L_X","MG_mg.L_X","CHLA.LC_ug.L_X","PC_ugL","Chl_ugL", 
                   "MCLR_ugL", "MCRR_ugL", "MC_dimethylLR_ugL","MC_980LR_deriv_ugL","Nodularin_811_ugL",
                   "Anatoxin_a_ugL","Euglenaphycin_ugL",'Anabaenopeptins',
                   'Time', 'Station_Id','TSI_Chla',
                   'TSI_TN','TSI_TP','TSI_TN2','TSI_TP2','SPCOND_uS.cm_x','TURB_NTU_X','T.SUS.SD_mg.L_X','NH4_mg.L_X',
                   'OPO4_mg.L_X','TPO4_mg.L_X','ALKALNTY_mg.L_X','TSI','Cyanopeptide','WaterBody','Secchi_Depth',
                   'NOX_mg.L_X','TOTAL_N_mg.L_X',"TN_mgL","TN_TP_X","Salinity_state","tot_cya",
                   "DO.","Sample", "TP_mgL", "TN_mgL" 
  )) %>%
  mutate_all(funs(ifelse(. < 0, 0, .))) %>% 
  mutate("NO2_mgL" = na_if(NO2_mgL, 0)) %>%
  mutate("NO3_mgL" = na_if(NO3_mgL, 0)) %>%
  mutate("NH3_mgL" = na_if(NH3_mgL, 0)) %>%
  mutate("TRP_mgL" = na_if(TRP_mgL, 0)) %>%
  mutate("DIN" = (NOx_mgL+NH3_mgL)) %>%
  mutate("DIN_DIP" = (NOx_mgL+NH3_mgL)/TRP_mgL) %>%
  mutate("NH4_NOX" = (NH3_mgL/NOx_mgL)) %>%
  mutate("depth_thing" = (Photic_Depth/Depth_m)) %>%
  mutate_all(funs(ifelse(is.infinite(.), median(., na.rm = TRUE), .))) %>%
  as.data.frame()


physeq_veg_Genus = physeq_cyanobacteria_rare %>% tax_glom(taxrank="Genus", NArm=T) %>% tax_fix()

GenusNames = tax_name(physeq_veg_Genus, prefix = "", rank = "Genus", sep = "&") %>% taxa_names()

GenusNames = gsub(".*&","",GenusNames)

#replaces OTU identifer with Genus name

taxa_names(physeq_veg_Genus) = GenusNames

Protist_Genus = ranacapa::vegan_otu(physeq_veg_Genus) %>% as.data.frame() %>% tibble::rownames_to_column(var = "SampleID") %>% arrange(SampleID)



ttt = full_join(Protist_Genus, metadata_1, by = "SampleID")
write_csv(ttt, "~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/For_GAM_genus_RARE.csv")





