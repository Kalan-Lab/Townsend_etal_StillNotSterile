# Still not sterile R code 

library(tidyverse)
library(data.table)
library(phyloseq)
library(speedyseq)
library(ggplot2)
library(maditr)
library(microbiome)
library(metagMisc)
library(qiime2R)
library(xlsx)
library(RColorBrewer)
library(microshades)
library(decontam)
library(vegan)
library(DESeq2)
library(ANCOMBC) 
library(Maaslin2)
library(pheatmap)
library(MMUPHin)
library(MBECS)
library(pivottabler)
library(reshape2)

setwd(dir = "https://github.com/Kalan-Lab/Townsend_etal_StillNotSterile")

# Figure 2B - PMAxx optomizatin 
PCR <- read.csv("/AllForearmResultsFIN.csv")
PCR <- as.data.frame(PCR)
PCR2 <- subset.data.frame(PCR, LHKC != "Heated")

ggplot(data=PCR2, aes(x=reorder(LHKC, LHKCorder), fill = factor(LHKC), color = factor(LHKC))) + 
  geom_boxplot(aes(y = MeanPercent), size = 1.3, alpha = 0.7)+
  geom_point(aes(y = MeanPercent),  size = 2.5)+
  #geom_point(aes(y = pCTMean, shape = factor(Replicate)), color = "gray",  size = 2.5)+
  #scale_x_discrete( scales = "free_x") + 
  facet_wrap(~PMAxxConcentration, scales = "free_x", ncol = 4)+
  scale_y_continuous(limits = c(0,1),breaks = c(0, .2, .4, .6, 0.8, 1.0)) + 
  scale_fill_manual(values = c("#EF7234", "#843369","#483470"))+
  scale_color_manual(values = c("#D05011", "#491D3A", "#241A38"))+
  theme_light() +
  ggtitle("Viability PCR POptomizing PMAxx concentration on Human Forearm Samples; 'Percent Live'") +
  ylab("Percent Live") +
  theme(plot.title = element_text(hjust = 0.5))

# Figure 2C - PMAxx optomization body sites 
PCR <- read.csv("/ECT-multipleBodyPMAxx.csv")
PCR <- as.data.frame(PCR)
PCR$Group<- factor(PCR$Group, levels = c("BackLive10", "BackHeatkilled10", "ForearmLive10", "ForearmHeatkilled10", "UmbillicusLive10", "UmbillicusHeatkilled10", "BackLive25", "BackHeatkilled25", "ForearmLive25", "ForearmHeatkilled25", "UmbillicusLive25", "UmbillicusHeatkilled25"))
PCR <- subset.data.frame(PCR, Concentration == "10" )
ggplot(data=PCR, aes(x = Group, y = MeanPercent, color = Group, fill = Group)) + 
  geom_boxplot(alpha = 0.8, size = 1.3)+
  geom_point( size = 2.5)+
  scale_fill_manual(values = c("#483470","#843369","#483470","#843369","#483470", "#843369"))+
  scale_color_manual(values = c( "#241A38","#491D3A", "#241A38",  "#491D3A", "#241A38", "#491D3A"))+
  scale_y_continuous(limits = c(0,1),breaks = c(0, .2, .4, .6, 0.8, 1.0)) + 
  facet_wrap(~BodySite, scales = "free_x", ncol = 4)+
  theme_light() +
  ggtitle("Viability PCR Multiple Body Site Percent live bacteria") +
  ylab("MEAN Perent live") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = -90))


# Figure 3A - qPCR for CHG subjects 
PCR <- read.csv("CHGstudy_qPCR_df.csv")
PCR<- as.data.frame(PCR)

PCR$Timepoint <- factor(PCR$Timepoint, levels = c("Base", "Pre", "OR", "Post", "CV", "Pre2", "OR2", "Post2"))
PCR$TotalBacteria <- gsub(",", "", PCR$TotalBacteria)
PCR$TotalBacteria = as.numeric(as.character(PCR$TotalBacteria))
PCR$ViableBacteria <- gsub(",", "", PCR$ViableBacteria)
PCR$ViableBacteria = as.numeric(as.character(PCR$ViableBacteria))
PCR2<-subset.data.frame(PCR, Timepoint %in% c("Base", "Pre", "OR", "Post", "CV"))
PCR2$MDsite <- paste(PCR2$MoistDry, PCR2$Site)
#PCR_Tbase <- subset.data.frame(PCR2, Timepoint == "Base")
ggplot(PCR2, aes(x= Timepoint))+
  geom_smooth(aes(x= Timepoint, y = log10(V_FC_from_Tbase +1), group = Site),color = "#491D3A", formula = "y~x", fun = "median", stat = 'summary' , method = "loess", position = position_nudge(x=.15 ))+
  geom_smooth(aes(x= Timepoint, y = log10(T_FC_from_Tbase +1), group = Site),color = "#AE8509", formula = "y~x", fun = "median", stat = 'summary' , method = "loess", position = position_nudge(x=-.15 ))+
  geom_boxplot(aes(x= Timepoint, y = log10(T_FC_from_Tbase +1)), fill = "#E9B10C",color = "#AE8509", alpha = 0.8,shape = 25, size =0.7, width = 0.25, position = position_nudge(x=-.15))+
  geom_boxplot(aes(x= Timepoint, y = log10(V_FC_from_Tbase +1)), fill = "#843369",color = "#491D3A",alpha = 0.8, shape =19, size = 0.7, width = 0.25, position = position_nudge(x=.15))+
  #geom_boxplot(data = CFU_Tbase, aes(x= Timepoint, y = log10(T_FC_from_TbaseReplicates +1)), fill = "#E9B10C",color = "#AE8509", alpha = 0.8,shape = 25, size =0.7, width = 0.25, position = position_nudge(x=-.15))+
  geom_hline(yintercept = 0, color = "#C1C1C1")+
  facet_grid(cols = vars(MDsite), scales = "free")+
  scale_y_continuous(limits = c(-4,1))+
  theme_light()+
  ylab("Log10(Fold Change) from total Bacterial at Baseline- w/ replicates")+
  ggtitle("Log10(Fold Change) from total Bacterial at Baseline w/ replicates")

ggplot(PCR2, aes(x= Timepoint))+ # Plot for presentation (drop from viable)
  geom_smooth(aes(x= Timepoint, y = log10(V_FC_from_Vbase +1), group = Site),color = "#491D3A", formula = "y~x", fun = "median", stat = 'summary' , method = "loess")+
  geom_boxplot(aes(x= Timepoint, y = log10(V_FC_from_Vbase +1)), fill = "#843369",color = "#491D3A",alpha = 0.8, shape =19, size = 0.7, width = 0.7)+
  #geom_boxplot(data = CFU_Tbase, aes(x= Timepoint, y = log10(T_FC_from_TbaseReplicates +1)), fill = "#E9B10C",color = "#AE8509", alpha = 0.8,shape = 25, size =0.7, width = 0.25, position = position_nudge(x=-.15))+
  geom_hline(yintercept = 0, color = "#C1C1C1")+
  facet_grid(cols = vars(Site), scales = "free")+
  scale_y_continuous(limits = c(-4,1))+
  theme_light()+
  ylab("Log10(Fold Change) from total Bacterial at Baseline- w/ replicates")+
  ggtitle("Log10(Fold Change) from total Bacterial at Baseline w/ replicates")

# Supplemental Figure 1A - percent live bacteria in surgical subjects 
ggplot(PCR2, aes(x= Timepoint, y = PercentLive))+
  geom_smooth(aes(group = Site),color = "#241A38", formula = "y~x", fun = "median", stat = 'summary' , method = "loess")+
  geom_boxplot(fill = "#483470",alpha = 0.8, shape =19, size = 1, width = 0.6)+
  geom_point(color = "#241A38", shape = 19, size = 1)+
  facet_grid(rows = vars(MoistDry), cols = vars(Site), scales = "free")+
  scale_y_continuous(limits = c(0,1))+
  theme_light()+
  ylab("Percent Live bacteria ")+
  ggtitle("Percent Live bacteria CHG subjects")



### CFU 
CHG_CFU <- read.csv("CHG-CFU-DF.csv")
CHG_CF <- as.data.frame(CHG_CFU)
CHG_CF <- subset.data.frame(CHG_CF, Site != "Leg Rash")
CHG_CF$Timepoint <- factor(CHG_CF$Timepoint, levels = c("Baseline","Pre-OR", "Post-OR", "CV1"))

CFU_complete <- subset.data.frame(CHG_CF, !SubjectID %in% c("CHG-002", "CHG-004","CHG-005","CHG-007","CHG-014","CHG-022","CHG-023","CHG-029","CHG-037","CHG-038") ) # removing subjects that did not endup undergoing surgery 
CFU2 <- subset.data.frame(CFU_complete, Timepoint %in% c("Post-OR", "CV1"))
ggplot(CFU2, aes(x=Timepoint, y=log10(ColonyFormingUnits+1), color = MoistDry, fill = MoistDry))+
  geom_boxplot(position = position_dodge2(width = 0.8), alpha = 0.7, size = 1.4)+
  geom_point(size = 2.5, alpha = 0.9, position = position_dodge(width = 0.75))+
  facet_wrap(~Site, nrow = 1, scales = "free_x")+
  theme_light()+
  scale_fill_manual(values = c("#E9B10C", "#843369"))+
  scale_color_manual(values = c( "#AE8509","#491D3A"))+
  scale_y_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5))+
 # scale_fill_manual(values = c("#F07E45", "#843369","#392F75"))+
 # scale_color_manual(values = c("#D05011", "#491D3A","#1D173A"))+
  theme(text = element_text(size = 14))



### 16S analysis 
ps.CHG.Part1 <- qza_to_phyloseq(features = "table-CHGstudyPart1-dada2.qza",
                              taxonomy = "taxonomy-Silva2022-CHGstudyPart1.qza",
                              metadata = "LK16S009_CHG_Study_Manufest_Server_HealthData.txt")
ps.CHG.Part2 <- qza_to_phyloseq(features = "table-CHGstudyPart2-dada2.qza",
                                taxonomy = "taxonomy-Silva2022-CHGstudyPart2.qza",
                                metadata = "LK16S011_CHG_Study_Manufest_Server_HealthData.txt")

ps.CHG.NotClean <- merge_phyloseq(ps.CHG.Part1, ps.CHG.Part2)
CHG_Study.meta <-  data.frame(sample_data(ps.CHG.NotClean))
CHG_Study.meta$SampleID <- row.names(CHG_Study.meta)
#write.csv(CHG_Study.meta, "CHG_Study_FULL_metadata.csv")

# Negative check 
ps.neg<- subset_samples(ps.CHG.NotClean , PMAxx == "NegEx")
NegREL <- as.data.frame(psmelt(transform(ps.neg, "compositional")))
NegREL2 <- subset.data.frame(NegREL, Abundance > 0)
#write.csv(NegREL2, "/Users/liztown/Documents/KalanLab/CHG-SSI_Study/16S_Seq/CHGstudy16Sanalysis_FINAL/Negatives/NegAbundance.csv")

# cleening up the stuff
ps.Decontam <- subset_taxa(ps.CHG.NotClean, Family !="mitochondria")
ps.Decontam <- subset_taxa(ps.Decontam, Class !="Chloroplast")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !=" ")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !="Cyanobacteria")
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="0319-6G20") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Acidibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Anaeromyxobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Bradyrhizobium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Curtobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Chryseobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Curvibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Dyadobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Enhydrobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Entotheonellaceae") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Exiguobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Ferruginibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Gaiella") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Gracilibacteria") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Herbaspirillum") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="LWQ8") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Marvinbryantia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylobacterium")
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="MB-A2-108") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylotenera") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Motilibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Nakamurella") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="NS9_marine_group") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Ralstonia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Rathayibacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Rhodoferax") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sandaracinobacter") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sediminibacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingomonas")
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingopyxis") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="SWB02") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Tepidiphilus") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Thermomonas") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Undibacterium") 
ps.Decontam <- subset_samples(ps.Decontam, BodySite != "Posterior Elbow")
ps.Decontam  <- subset_samples(ps.Decontam, BodySite != "24 Strain Mock")
ps.CHG <- subset_samples(ps.Decontam, PMAxx!= "NegEx") # removing the negative control
ps.CHG <- phyloseq_filter_prevalence(ps.CHG, prev.trh = 0.02, abund.trh = 100, threshold_condition = "AND") 

ps.CHG# THE CLEANED FILE 
# output summative tables of total and relativce abundance 
CHG.clean <- psmelt(ps.CHG)
#write.csv(CHG.clean, "/Users/liztown/Documents/KalanLab/CHG-SSI_Study/16S_Seq/CHGstudy16Sanalysis_FINAL/CHGsilva_abundanceTotalsOut.csv")
CHG.clean.OTU <- psmelt(ps.CHG@otu_table)
#write.csv(CHG.clean.OTU, "/Users/liztown/Documents/KalanLab/CHG-SSI_Study/16S_Seq/CHGstudy16Sanalysis_FINAL/PICRUST2/CHGforPICRUST_TotalREADS_OTU.csv")

relAlls <- microbiome::transform(ps.CHG, "compositional")
relAllsout <- psmelt(relAlls)
#write.csv(relAllsout, "/Users/liztown/Documents/KalanLab/CHG-SSI_Study/16S_Seq/CHGstudy16Sanalysis_FINAL/CHGsilva_RelativeOUT.csv")
CHG.REL.df<-as.data.frame(relAllsout)

Select.taxa <- subset_taxa(relAlls, Genus %in% c("Pseudomonas", "Staphylococcus", "Escherichia-Shigella", "Bacillus", "Acinetobacter"))
relsout <- psmelt(Select.taxa)
#write.csv(relsout, "/Users/liztown/Documents/KalanLab/CHG-SSI_Study/16S_Seq/CHGstudy16Sanalysis_FINAL/CHGsilva_SELECT_TAXA_RelativeOUT.csv")
Select.taxa.df<-as.data.frame(relsout)

                           
# library size 
df <- as.data.frame(sample_data(ps.CHG)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps.CHG)
df <- df[order(df$LibrarySize),]
med <- median(df$LibrarySize)
med
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=log10(LibrarySize), color=Plate, shape = PMAxx)) + 
  geom_point(size = 3) +
  #scale_y_continuous(breaks = c(0, 1000, 5000, 10000, 20000, 30000))+
  theme_bw()+
  theme(text = element_text(size = 18))+
  ggtitle("Library Size for the CHG")# size for whole data frame

# Supplemental Figure 9 Alpha Abundance 
ps.CHG.gen <- subset_samples(ps.CHG, BodySite != "24 Strain Mock")
ps.CHG.gen@sam_data$Timepoint[ps.CHG.gen@sam_data$Timepoint == "Pre2"] <- "Pre"
ps.CHG.gen@sam_data$Timepoint[ps.CHG.gen@sam_data$Timepoint == "OR2"] <- "OR"
ps.CHG.gen@sam_data$Timepoint[ps.CHG.gen@sam_data$Timepoint == "Post2"] <- "Post"
ps.CHG.gen@sam_data$Timepoint <- factor(ps.CHG.gen@sam_data$Timepoint, c("Base", "Pre","OR", "Post", "CV1"))
tab <-microbiome::alpha(ps.CHG.gen, index = c("observed", "chao1", "diversity_inverse_simpson", "diversity_gini_simpson","diversity_shannon",	"diversity_coverage", "evenness_camargo",	"evenness_pielou",
                                         "evenness_simpson", "evenness_evar", "evenness_bulla", "dominance_dbp", "dominance_dmn", "dominance_absolute", "dominance_relative", "dominance_simpson", 
                                         "dominance_core_abundance", "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_rare_abundance"))
#write.csv(tab, "/Users/liztown/Documents/KalanLab/CHG-SSI_Study/16S_Seq/CHGstudy16Sanalysis_FINAL/CHGStudy_alphaTable.csv")

AlphaCHG <- tab
AlphaCHG$Sample<- row.names(AlphaCHG)
AlphaCHG <- merge(AlphaCHG, CHG.REL.df, by = "Sample", all = T) # the alpha table above but with the metadata added back in 
AlphaCHG$Timepoint[AlphaCHG$Timepoint == "Pre2"] <- "Pre"
AlphaCHG$Timepoint[AlphaCHG$Timepoint == "OR2"] <- "OR"
AlphaCHG$Timepoint[AlphaCHG$Timepoint == "Post2"] <- "Post"
AlphaCHG$Timepoint <- factor(AlphaCHG$Timepoint, levels = c("Base", "Pre","OR", "Post", "CV1"))

ggplot(AlphaCHG, aes(x = Timepoint, y =diversity_shannon, color = MoistDry, fill = MoistDry))+
  #geom_point(size = 1)+
  geom_boxplot(alpha = 0.9,position =position_dodge2(preserve = "single"))+
  stat_summary(fun = "median", linewidth = 1, geom = "line", aes(group = MoistDry, color = MoistDry))+
  facet_grid(rows = vars(rev(PMAxx)), cols = vars(Site), scales = "free_x")+
  scale_fill_manual(values = c("#E9B10C", "#843369"))+
  scale_color_manual(values = c( "#AE8509","#491D3A"))+
  ylim(0,3.5)+
  theme_light()+
  ggtitle("Shannon Alpha Diversity (ASV level) CHG study")


# Figure 3B -rlative abundance plots 
gen <- tax_glom(ps.CHG, "Genus")
relAlls <-microbiome::transform(gen, "compositional")
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "Pre2"] <- "Pre"
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "OR2"] <- "OR"
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "Post2"] <- "Post"
CHG.REL <- psmelt(relAlls)
OtherCHGsilva <- as.data.frame(CHG.REL)
OtherCHGsilva <- subset.data.frame(OtherCHGsilva, MoistDry != "24 Strain Mock")
OtherCHGsilva <- subset.data.frame(OtherCHGsilva, Abundance >= 0.001)
OtherCHGsilva <- subset.data.frame(OtherCHGsilva, Site != "NegEx")
OtherCHGsilva$SiteMD <- paste(OtherCHGsilva$MoistDry, OtherCHGsilva$Site) 
OtherCHGsilva <- OtherCHGsilva %>% mutate(Genus = ifelse(Abundance < 0.05 , "Other", Genus)) # replacing the genus catagory with "Other" if the abundance is less than 1% 
OtherCHGsilva <- OtherCHGsilva %>% mutate(Phylum = ifelse(Abundance < 0.01, "Other", Phylum)) 
table(OtherCHGsilva$Genus, OtherCHGsilva$Phylum)
OtherCHGsilva <- OtherCHGsilva %>% mutate(Genus = case_when(Genus == "uncultured" ~ Family,
                                                            TRUE ~ Genus))
OtherCHGsilva <- OtherCHGsilva %>% mutate(Genus = case_when(Genus == "Caulobacteraceae" ~ "Caulobacteraceae uncultured genus",
                                                            Genus == "Neisseriaceae" ~ "Neisseriaceae uncultured genus",
                                                            TRUE ~ Genus))
table(OtherCHGsilva$Genus, OtherCHGsilva$Phylum)


OtherCHGsilva$Timepoint <- factor(OtherCHGsilva$Timepoint, levels = c("Base", "Pre","OR", "Post", "CV1"))
#write.csv(OtherCHGsilva, "/Users/liztown/Documents/KalanLab/CHG-SSI_Study/16S_Seq/CHGstudy16Sanalysis_FINAL/CHGsilvaGenusWithOTHER.REl.csv")

OtherCHGsilva$Genus<- factor(OtherCHGsilva$Genus, levels = c("Other","uncultured",
                                                             "Actinomyces", "Brachybacterium","Brevibacterium", "Corynebacterium","Cutibacterium","Dermabacter", "Dietzia", "Gardnerella","Gordonia",
                                                             "Kocuria","Kytococcus", "Lawsonella", "Microbacterium", "Micrococcus","Millisia","Mycobacterium", "Ornithinimicrobium"," Pseudoclavibacter",
                                                             "Rhodococcus","Rothia","Streptomyces", #Actinobacteria / actinomycetota (20)
                                                             
                                                             "Abiotrophia", "Aerococcus","Anaerococcus","Anoxybacillus", "Bacillus","Blautia","Brevibacillus", "Clostridium_sensu_stricto_1", "Enterococcus",
                                                             "Eremococcus", "Ezakiella","Facklamia", "Fenollaria", "Finegoldia","Gemella", "Granulicatella","Jeotgalicoccus",  "Lactobacillus", "Lactococcus",
                                                             "Leuconostoc", "Macrococcus","Murdochiella", "Mycoplasma","Negativicoccus","Oribacterium",
                                                             "Peptoniphilus","Peptostreptococcus","Romboutsia", "Staphylococcus", "Streptococcus","Subdoligranulum","Veillonella",  # Firmicutes/ Bacillota (21)
                                                             
                                                             "Bacteroides","Pedobacter","Porphyromonas","Prevotella", # Bacteroidetes / Bacteroidota (4)
                                                             "Campylobacter", #Campilobacteria
                                                             "Deinococcus", #Deinococcus
                                                             "Fusobacterium", #fuso
                                                             
                                                             "Acinetobacter","Aeromonas","Afipia","Amaricoccus","Aureimonas", "Bosea","Brevundimonas","Burkholderia-Caballeronia-Paraburkholderia",
                                                             "Caulobacter","Caulobacteraceae uncultured genus", "Comamonas","Delftia","Denitratisoma","DSSD61", "Enhydrobacter","Enterobacter","Escherichia-Shigella","Haematobacter", "Haemophilus", "Herbaspirillum", 
                                                             "Janthinobacterium","Klebsiella","Luteimonas","Massilia", "Methylobacterium-Methylorubrum","Methylophilus","Neisseria","Neisseriaceae uncultured genus",  "Ottowia","Pantoea", "Paracoccus",
                                                             "Pseudomonas","Raoultella","Roseomonas","Serratia","Sphingopyxis","Stenotrophomonas", "Xanthomonas")) # Proteobacteria /pseudomonadota (30)
sunsetshades <- c("#C1C1C1", #Gray
                  colorRampPalette(c("#F9E39F", "#E9B10C","#EC9220", "#EF7234","#D55D45",
                                     "#F5BAA3","#D77870","#BB4855","#843369","#62386A",
                                     "#DACCE1","#796796","#483470","#2B2358","#15112C"))(60))  # purples          

ggplot(OtherCHGsilva, aes(x= Timepoint, y =Abundance, fill = Genus))+
  geom_bar(stat="identity", position="fill") +
  facet_grid(cols = vars(SiteMD), rows = vars(PMAxx), scales = "free_x")+
  ylab("Relative Abundance")+
  scale_x_discrete(labels = c("Base", "Pre","OR", "Post", "CV1"))+
  theme_light()+
  scale_fill_manual(values = sunsetshades)+
  scale_color_manual(values = sunsetshades)+
  theme(axis.text.x = element_text(size = 7)) +
  ggtitle("Genus present > .5% of reads in a sample")


### DOT PLOTS
# Suplemental figure 6 - dotplot of key taxa for moist and dry sites separatly
gen <- tax_glom(ps.CHG, "Genus")
relAlls <- microbiome::transform(gen, "compositional")
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "Pre2"] <- "Pre"
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "OR2"] <- "OR"
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "Post2"] <- "Post"
relAlls@sam_data$SiteMD <- paste(relAlls@sam_data$MoistDry, relAlls@sam_data$Site) 
relAlls@sam_data$TimSitPxMD <- paste(relAlls@sam_data$MoistDry, relAlls@sam_data$TimeSitePx) 
relAlls.Grouped <- merge_samples2(relAlls, "TimSitPxMD",  fun_otu = mean)
CHG.Grouped <- psmelt(relAlls.Grouped)
CHG.Grouped  <- as.data.frame(CHG.Grouped)
CHG.Grouped <- subset.data.frame(CHG.Grouped , Site != "22 Strain Mock")
CHG.Grouped  <- subset.data.frame(CHG.Grouped , Site != "NegEx")
CHG.Grouped$Phylym<-factor(CHG.Grouped$Phylum, levels = c("Other", "Actinobacteria", "Firmicutes", "Bacteroidota", "Proteobacteria"))
CHG.Grouped$Timepoint <- factor(CHG.Grouped$Timepoint, levels = c("Base", "Pre","OR", "Post","CV1"))
CHG.Grouped$PMAxx <- factor(CHG.Grouped$PMAxx, levels = c("PMAxx", "None"))

CHG.Grouped.AllTop <- subset.data.frame(CHG.Grouped, Genus %in% c("Acinetobacter", "Afipia","Anaerococcus", "Bacillus","Brachybacterium","Brevundimonas", "Corynebacterium", "Cutibacterium", "Dietzia", 
                                                                  "Escherichia-Shigella", "Kocuria", "Finegoldia", "Micrococcus","Staphylococcus", "Streptococcus", "Streptomyces",
                                                                  "Paracoccus", "Pseudomonas"))

ggplot(CHG.Grouped.AllTop, aes(x= Timepoint, y =Genus, Phylum, fill = Phylum, color = Phylum))+
  geom_point(aes(size = Abundance, fill = Phylum, color = Phylum)) +
  facet_grid(rows = vars(PMAxx), cols = vars(SiteMD), scales = "free")+
  ylab("Mean relative Abundance")+
  scale_size(range=c(-1,9), breaks=c(0,.1,.2,.3, 0.4,0.5,0.6))+
  scale_x_discrete(labels = c("Baseline", "Pre-OR","OR", "Post-OR", "Follow-Up"))+
  scale_y_discrete(limits = rev(c("Brachybacterium", "Corynebacterium", "Cutibacterium", "Dietzia","Kocuria","Micrococcus","Streptomyces",
                                  "Anaerococcus", "Bacillus", "Finegoldia","Staphylococcus",  "Streptococcus",
                                  "Acinetobacter", "Afipia","Brevundimonas", "Escherichia-Shigella","Paracoccus", "Pseudomonas")))+
  theme_light()+
  #scale_fill_manual(values = c("#EF7234","#62386A","#BB4855", "#C1C1C1","#2B2358"))+
  #scale_color_manual(values =c("#EF7234","#62386A","#BB4855", "#C1C1C1","#2B2358"))+
  scale_fill_manual(values = c("#EC9220","#843369","#D55D45", "#C1C1C1","#483470"))+
  scale_color_manual(values =c("#EC9220","#843369","#D55D45", "#C1C1C1","#483470"))+
  theme(axis.text.x = element_text(size = 8, angle = 30, hjust =1, vjust = 1), axis.text.y = element_text(size = 8)) +
  ggtitle("Mean Realtive abundance of genera")


# Figure 4C - dotplot for key taxa
gen <- tax_glom(ps.CHG, "Genus")
relAlls <- microbiome::transform(gen, "compositional")
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "Pre2"] <- "Pre"
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "OR2"] <- "OR"
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "Post2"] <- "Post"
relAlls.Grouped <- merge_samples2(relAlls, "TimeSitePx",  fun_otu = mean)
CHG.Grouped <- psmelt(relAlls.Grouped)
CHG.Grouped  <- as.data.frame(CHG.Grouped)
CHG.Grouped <- subset.data.frame(CHG.Grouped , Site != "22 Strain Mock")
CHG.Grouped  <- subset.data.frame(CHG.Grouped , Site != "NegEx")
CHG.Grouped$Phylym<-factor(CHG.Grouped$Phylum, levels = c("Other", "Actinobacteria", "Firmicutes", "Bacteroidota", "Proteobacteria"))
CHG.Grouped$Timepoint <- factor(CHG.Grouped$Timepoint, levels = c("Base", "Pre","OR", "Post","CV1"))
CHG.Grouped$PMAxx <- factor(CHG.Grouped$PMAxx, levels = c("PMAxx", "None"))

CHG.Grouped.AllTop <- subset.data.frame(CHG.Grouped, Genus %in% c("Acinetobacter", "Afipia","Anaerococcus", "Bacillus","Brachybacterium","Brevundimonas", "Corynebacterium", "Cutibacterium", "Dietzia", 
                                                                  "Escherichia-Shigella", "Kocuria", "Finegoldia", "Micrococcus","Staphylococcus", "Streptococcus", "Streptomyces",
                                                                  "Paracoccus", "Pseudomonas"))


ggplot(CHG.Grouped.AllTop, aes(x= Timepoint, y =Genus, Phylum, fill = Phylum, color = Phylum))+
  geom_point(aes(size = Abundance, fill = Phylum, color = Phylum)) +
  facet_grid(rows = vars(PMAxx), cols = vars(Site), scales = "free")+
  ylab("Mean relative Abundance")+
  scale_size(range=c(-1,7.5), breaks=c(0,0.01,.1,.2,.3, 0.4,0.5))+
  scale_x_discrete(labels = c("Base", "Pre","OR", "Post", "CV1"))+
  scale_y_discrete(limits = rev(c("Brachybacterium", "Corynebacterium", "Cutibacterium", "Dietzia","Kocuria","Micrococcus","Streptomyces",
                                  "Anaerococcus", "Bacillus", "Finegoldia","Staphylococcus",  "Streptococcus",
                                  "Acinetobacter", "Afipia","Brevundimonas", "Escherichia-Shigella","Paracoccus", "Pseudomonas")))+
  theme_light()+
  scale_fill_manual(values = c("#EC9220","#BB4855","#62386A", "#C1C1C1","#2B2358"))+
  scale_color_manual(values =c("#EC9220","#BB4855","#62386A", "#C1C1C1","#2B2358"))+
  #scale_fill_manual(values = c("#EC9220","#D55D45", "#843369","#C1C1C1","#483470"))+
  #scale_color_manual(values =c("#EC9220","#D55D45", "#843369", "#C1C1C1","#483470"))+
  theme(axis.text.x = element_text(size = 8, angle = 30, hjust =1, vjust = 1), axis.text.y = element_text(size = 8),legend.position = "none") +
  ggtitle("Mean Realtive abundance of genera")


# Supplemental Figure 5G - Staphylococcus species Reltative abundance plot and dot plot 
relAlls <- microbiome::transform(ps.CHG, "compositional")
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "Pre2"] <- "Pre"
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "OR2"] <- "OR"
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "Post2"] <- "Post"
relAlls@sam_data$SiteMD <- paste(relAlls@sam_data$MoistDry, relAlls@sam_data$Site) 
relAlls@sam_data$TimSitPx2 <- paste(relAlls@sam_data$Timepoint, relAlls@sam_data$Site, relAlls@sam_data$PMAxx) 
relAlls@sam_data$TimSitPxMD <- paste(relAlls@sam_data$MoistDry, relAlls@sam_data$TimeSitePx2) 
relAlls.Grouped <- merge_samples2(relAlls, "TimSitPx2",  fun_otu = mean)
StaphSpp <- subset_taxa(relAlls.Grouped, Genus== "Staphylococcus")
StaphSpp.df <- psmelt(relAlls.Grouped)
StaphSpp.df  <- as.data.frame(StaphSpp.df)
StaphSpp.df <- subset.data.frame(StaphSpp.df , Site != "22 Strain Mock")
StaphSpp.df  <- subset.data.frame(StaphSpp.df  , Site != "NegEx")

StaphBLAST <- read.csv("StaphSpeciesBLAST.csv")
StaphBLAST <-as.data.frame(StaphBLAST)
StaphBLAST$BLAST_Species[StaphBLAST$BLAST_Species == "Staphylococcus arueus"] <- "Staphylococcus aureus"

Staph.Spp.df <- merge(StaphSpp.df, StaphBLAST, by = "OTU")
Staph.Spp.df$Timepoint <- factor(Staph.Spp.df$Timepoint, levels = c("Base", "Pre", "OR", "Post", "CV1"))
Staph.Spp.df$BLAST_Species <- factor(Staph.Spp.df$BLAST_Species, levels = c("Staphylococcus aureus","CoNS Staphylococcus spp.", "Staphylococcus auricularis","Staphylococcus caprae or S. capitis", "Staphylococcus chromogenes ", 
                                                                            "Staphylococcus epidermidis", "Staphylococcus haemolyticus","Staphylococcus hominis", "Staphylococcus lugdunensis","Staphylococcus pasteuri",
                                                                            "Staphylococcus pettenkoferi", "Staphylococcus warneri"))
ggplot(Staph.Spp.df, aes(x= Timepoint, y =BLAST_Species,  fill = BLAST_Species, color = BLAST_Species))+
  geom_point(aes(size = Abundance, fill = BLAST_Species, color = BLAST_Species)) +
  facet_grid(rows = vars(PMAxx), cols = vars(Site), scales = "free")+
  ylab("Mean relative Abundance")+
  scale_size(range=c(-1,9), breaks=c(0,0.01,0.05,.1,.2, 0.3))+
  scale_x_discrete(labels = c("Base", "Pre","OR", "Post", "CV1"))+
  scale_y_discrete(limits = rev(c("Staphylococcus aureus","CoNS Staphylococcus spp.", "Staphylococcus auricularis","Staphylococcus caprae or S. capitis", "Staphylococcus chromogenes ",
                                  "Staphylococcus epidermidis", "Staphylococcus haemolyticus","Staphylococcus hominis", "Staphylococcus lugdunensis","Staphylococcus pasteuri",
                                  "Staphylococcus pettenkoferi", "Staphylococcus warneri")))+
  theme_light()+
  scale_fill_manual(values = c("#E9B10C", "#D55D45", colorRampPalette(c("#BB4855","#A03E5F","#843369","#62386A","#483470"))(10)))+
  scale_color_manual(values = c("#C2940A", "#D55D45", colorRampPalette(c("#BB4855","#A03E5F","#843369","#62386A","#483470"))(10)))+
  theme(axis.text.x = element_text(size = 8, angle = 30, hjust =1, vjust = 1), axis.text.y = element_text(size = 8))+ #,legend.position = "none") +
  ggtitle("Mean Realtive abundance of Staph Species (determined by BLAST)")


### Supplemental Figure 5H - Dotplot for Acinetobacter Bacillus Escherichia and pseudomonas 
PathSpp <- subset_taxa(relAlls.Grouped, Genus %in% c( "Acinetobacter", "Bacillus", "Escherichia-Shigella", "Pseudomonas"))
PathSpp.df <- psmelt(PathSpp)
PathSpp.df  <- as.data.frame(StaphSpp.df)
PathSpp.df <- subset.data.frame(PathSpp.df , Site != "22 Strain Mock")
PathSpp.df  <- subset.data.frame(PathSpp.df  , Site != "NegEx")

PathBLAST <- read.csv("PathBlast.csv")
PathBLAST <-as.data.frame(PathBLAST)

Path.Spp.df <- merge(PathSpp.df, PathBLAST, by = "OTU")
Path.Spp.df$BLAST_Species <- factor(Path.Spp.df$BLAST_Species, levels = c("Acinetobacter spp. ","Acinetobacter johnsonii", "Acinetobacter junii","Acinetobacter lwoffii", "Acinetobacter pittii", 
                                                                          "Acinetobacter radioresistens", "Acinetobacter ursingii","Acinetobacter variabilis", 
                                                                          "Bacillus spp. ","Bacillus cereus","Bacillus coagulans", "Bacillus subtilis",
                                                                          "Escherichia-Shigella spp.","Escherichia coli",
                                                                          "Pseudomonas spp.", "Pseudomonas aeruginosa", "Pseudomonas fluorescens","Pseudomonas lactis","Pseudomonas stutzeri" ))


Path.Spp.df.Top<-Path.Spp.df # reclassifying really lowly abundant species into general species
Path.Spp.df.Top$BLAST_Species[Path.Spp.df.Top$BLAST_Species== "Acinetobacter johnsonii"] <- "Acinetobacter spp. "
Path.Spp.df.Top$BLAST_Species[Path.Spp.df.Top$BLAST_Species== "Acinetobacter junii"] <- "Acinetobacter spp. "
Path.Spp.df.Top$BLAST_Species[Path.Spp.df.Top$BLAST_Species== "Acinetobacter johnsonii"] <- "Acinetobacter spp. "
Path.Spp.df.Top$BLAST_Species[Path.Spp.df.Top$BLAST_Species== "Acinetobacter ursingii"] <- "Acinetobacter spp. "
Path.Spp.df.Top$BLAST_Species[Path.Spp.df.Top$BLAST_Species== "Bacillus coagulans"] <- "Bacillus spp. "
Path.Spp.df.Top$Timepoint <- factor(Path.Spp.df.Top$Timepoint, levels = c("Base", "Pre", "OR", "Post", "CV1"))


ggplot(Path.Spp.df.Top, aes(x= Timepoint, y =BLAST_Species, fill = Genus, color = Genus))+
  geom_point(aes(size = Abundance, fill = BLAST_Species, color = BLAST_Species)) +
  facet_grid(rows = vars(PMAxx), cols = vars(Site), scales = "free")+
  ylab("Mean relative Abundance")+
  scale_size(range=c(-1,9), breaks=c(0, 0.01,0.05,.1))+
  scale_x_discrete(labels = c("Baseline", "Pre-OR","OR", "Post-OR", "Follow-UP"))+
  scale_y_discrete(limits = rev(c("Acinetobacter spp. ","Acinetobacter lwoffii", "Acinetobacter pittii", 
                                  "Acinetobacter radioresistens", "Acinetobacter variabilis", 
                                  "Bacillus spp. ","Bacillus cereus", "Bacillus subtilis",
                                  "Escherichia-Shigella spp.","Escherichia coli",
                                  "Pseudomonas spp.", "Pseudomonas aeruginosa",  "Pseudomonas fluorescens",
                                  "Pseudomonas lactis","Pseudomonas stutzeri" )))+
  theme_light()+
  
  scale_fill_manual(values = c(colorRampPalette(c("#E9B10C","#EC9220","#EF7234"))(5),
                               colorRampPalette(c("#BB4855","#A03E5F","#843369"))(5),
                               colorRampPalette(c("#62386A","#483470","#2B2358"))(5)))+
  scale_color_manual(values = c(colorRampPalette(c("#E9B10C","#EC9220","#EF7234"))(5),
                                colorRampPalette(c("#BB4855","#A03E5F","#843369"))(5),
                                colorRampPalette(c("#62386A","#483470","#2B2358"))(5)))+
  theme(axis.text.x = element_text(size = 8, angle = 30, hjust =1, vjust = 1), axis.text.y = element_text(size = 8)))+ #,legend.position = "none") +
  ggtitle("Mean Realtive abundance of Pathogenic Species (determined by BLAST)")
                                                 
# Supplemental figure 3 B Dot plot for BASELINE- viable microbiome - samples 
#Dot plot:
ps.CHG.BASE<-subset_samples(ps.CHG, Timepoint == "Base")
ps.CHG.gen.base <- tax_glom(ps.CHG.BASE, "Genus") #THIS IS WHAT YOU START WITH THE DOT PLOT
relAlls <- microbiome::transform(ps.CHG.gen.Base, "compositional") #cannot convert into a dataframe
relAlls@sam_data$BodyTimePx <- paste(relAlls@sam_data$BodySite, relAlls@sam_data$Timepoint, relAlls@sam_data$PMAxx)
relAlls.Grouped <- merge_samples2(relAlls, "BodyTimePx",  fun_otu = mean)
CHG.Grouped <- psmelt(relAlls.Grouped)
CHG.Grouped  <- as.data.frame(CHG.Grouped)
CHG.Grouped <- subset.data.frame(CHG.Grouped , Site != "22 Strain Mock")
CHG.Grouped  <- subset.data.frame(CHG.Grouped , Site != "NegEx")
CHG.Grouped$Phylum<-factor(CHG.Grouped$Phylum, levels = c("Actinobacteriota", "Firmicutes","Gemmatimonadota", "Bacteroidota", "Patescibacteria" ,"Proteobacteria"))

CHG.Grouped.AllTop <- subset.data.frame(CHG.Grouped, Genus %in% c("Corynebacterium", "Cutibacterium","Kocuria", "Micrococcus",
                                                                  "Anaerococcus", "Bacillus","Staphylococcus",  "Streptococcus",
                                                                  "Acinetobacter", "Afipia", "Escherichia-Shigella","Pseudomonas"))

CHG.Grouped.AllTop$BodySite <-factor(CHG.Grouped.AllTop$BodySite, levels = c("Forearm", "Armpit","Abdomen", "Knee","Lower Abdomen; Inquinal", "Umbillicus"))


ggplot(CHG.Grouped.AllTop, aes(x= BodySite, y= Genus, fill = Phylum, color = Phylum))+
  geom_point(aes(size = Abundance, fill = Phylum, color = Phylum)) +
  facet_grid(rows = vars(PMAxx), cols = vars(MoistDry), scales = "free")+
  ylab("Mean relative Abundance")+
  xlab("BodySites")+
  scale_x_discrete()+
  scale_y_discrete(limits = rev(c("Corynebacterium", "Cutibacterium", "Kocuria", "Micrococcus",
                                  "Anaerococcus", "Bacillus","Staphylococcus",  "Streptococcus",
                                  "Acinetobacter", "Afipia", "Escherichia-Shigella","Pseudomonas")))+
  scale_size(range=c(-1,12), breaks=c(0,.01, .1,.2,.3, 0.4,0.5,0.6))+
  scale_fill_manual(values = c("#EC9220","#843369","#D55D45", "#C1C1C1","#483470"))+
  scale_color_manual(values =c("#EC9220","#843369","#D55D45", "#C1C1C1","#483470"))+
  theme_light()+
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) +
  ggtitle("Mean Relative abundance of genera")



# Beta diversity plots 
Sunset<-c("#C1C1C1", #Gray
          colorRampPalette(c("#F8DD8C", "#E9B10C","#EC9220", "#EF7234","#D55D45", "#BB4855", "#A03E5F",
                             "#843369", "#62386A","#483470","#2B2358","#15112C"))(76))#yellows #Actinobacteria / actinomycetota (7)

ps.CHG.gen <- tax_glom(ps.CHG, "Genus")
ps.CHG.gen <- subset_samples(ps.CHG.gen, BodySite != "24 Strain Mock")
ps.CHG.gen@sam_data$Timepoint[ps.CHG.gen@sam_data$Timepoint == "Pre2"] <- "Pre"
ps.CHG.gen@sam_data$Timepoint[ps.CHG.gen@sam_data$Timepoint == "OR2"] <- "OR"
ps.CHG.gen@sam_data$Timepoint[ps.CHG.gen@sam_data$Timepoint == "Post2"] <- "Post"
min(sample_sums(ps.CHG.gen))# minimum sample read is 27
median(sample_sums(ps.CHG.gen)) # 5081
max(sample_sums(ps.CHG.gen)) #2567160
table(sample_sums(ps.CHG.gen))

rarecurve(t(otu_table(ps.CHG)), step=100) 
rarecurve(t(otu_table(ps.CHG.gen)), step=100, ylim =c(0,50), xlim=c(0,900))  ## considering using a read cut off of 5000 for beta diversity metrics 

ps.rarefiedsmall = rarefy_even_depth(ps.CHG.gen, rngseed=1, sample.size=500, replace=F)

# Supplemental Figure 2 Total microb iome vs, viable microbiome (PMAxx) over time 
ps.rar.Base<- subset_samples(ps.rarefiedsmall, Timepoint == "Base")
ps.rar.Pre<- subset_samples(ps.rarefiedsmall, Timepoint == "Pre")
ps.rar.OR<- subset_samples(ps.rarefiedsmall, Timepoint == "OR")
ps.rar.Post<- subset_samples(ps.rarefiedsmall, Timepoint == "Post")
ps.rar.CV<- subset_samples(ps.rarefiedsmall, Timepoint == "CV1")

GP = ps.rar.Base
GP.ord <- ordinate(GP, "NMDS",  "bray")
plot_ordination(GP, GP.ord, type="samples", color="PMAxx")+
  geom_point(size=1) + ggtitle("Bray Curtis")+
  stat_ellipse(type = "t", level = 0.9)+
  geom_point(aes( color =PMAxx), size = 4.5)+
  scale_color_manual(values = c( "#E9B10C", "#843369"))+
  scale_y_continuous(limits = c(-1.25, 1.25))+
  scale_x_continuous(limits = c(-1.25, 1.25))+
  ggtitle("Baseline")+
  theme_bw()
GP = ps.rar.Pre
GP.ord <- ordinate(GP, "NMDS",  "bray") 
plot_ordination(GP, GP.ord, type="samples", color="PMAxx")+
  geom_point(size=1) + ggtitle("Bray Curtis")+
  stat_ellipse(type = "t", level = 0.9)+
  geom_point(aes( color =PMAxx), size = 4.5)+
  scale_color_manual(values = c( "#E9B10C", "#843369"))+
  ggtitle("Pre")+
  scale_y_continuous(limits = c(-1.25, 1.25))+
  scale_x_continuous(limits = c(-1.25, 1.25))+
  theme_bw()
GP = ps.rar.OR
GP.ord <- ordinate(GP, "NMDS",  "bray") 
plot_ordination(GP, GP.ord, type="samples", color="PMAxx")+ 
  geom_point(size=1) + ggtitle("Bray Curtis")+
  stat_ellipse(type = "t", level = 0.9)+
  geom_point(aes( color =PMAxx), size = 4.5)+
  scale_color_manual(values = c( "#E9B10C", "#843369"))+
  scale_y_continuous(limits = c(-1.25, 1.25))+
  scale_x_continuous(limits = c(-1.25, 1.25))+
  ggtitle("OR")+
  theme_bw()
GP = ps.rar.Post
GP.ord <- ordinate(GP, "NMDS",  "bray") 
plot_ordination(GP, GP.ord, type="samples", color="PMAxx")+ + 
  geom_line(aes(group = SpecificSubjectSample), color = "#CCCCCC")+
  geom_point(size=1) + ggtitle("Bray Curtis")+
  stat_ellipse(type = "t", level = 0.9)+
  geom_point(aes( color =PMAxx), size = 4.5)+
  scale_color_manual(values = c( "#E9B10C", "#843369"))+
  scale_y_continuous(limits = c(-1.25, 1.25))+
  scale_x_continuous(limits = c(-1.25, 1.25))+
  ggtitle("Post")+
  theme_bw()
GP = ps.rar.CV
GP.ord <- ordinate(GP, "NMDS",  "bray") 
plot_ordination(GP, GP.ord, type="samples", color="PMAxx")+
  geom_line(aes(group = SpecificSubjectSample), color = "#CCCCCC")+
  geom_point(size=1) + ggtitle("Bray Curtis")+
  stat_ellipse(type = "t", level = 0.9)+
  geom_point(aes( color =PMAxx), size = 4.5)+
  scale_color_manual(values = c( "#E9B10C", "#843369"))+
  scale_y_continuous(limits = c(-1.25, 1.25))+
  scale_x_continuous(limits = c(-1.25, 1.25))+
  ggtitle("CV")+
  theme_bw()

# Supplemental Table 4 - Total vs. viable microbiome PERMANOVAs
sampledf <- data.frame(sample_data(ps.rarefiedsmall) )
CHGsilva_bray <- phyloseq::distance(ps.rarefiedsmall, method = "bray")
Bray_clust <-hclust(CHGsilva_bray)
adonis2(CHGsilva_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) #

sampledf <- data.frame(sample_data(ps.rar.Base) )
CHGsilva_bray <- phyloseq::distance(ps.rar.Base, method = "bray")
Bray_clust <-hclust(CHGsilva_bray)
adonis2(CHGsilva_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

sampledf <- data.frame(sample_data(ps.rar.Pre) )
CHGsilva_bray <- phyloseq::distance(ps.rar.Pre, method = "bray")
Bray_clust <-hclust(CHGsilva_bray)
adonis2(CHGsilva_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

sampledf <- data.frame(sample_data(ps.rar.OR) )
CHGsilva_bray <- phyloseq::distance(ps.rar.OR, method = "bray")
Bray_clust <-hclust(CHGsilva_bray)
adonis2(CHGsilva_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

sampledf <- data.frame(sample_data(ps.rar.Post) )
CHGsilva_bray <- phyloseq::distance(ps.rar.Post, method = "bray")
Bray_clust <-hclust(CHGsilva_bray)
adonis2(CHGsilva_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 

sampledf <- data.frame(sample_data(ps.rar.CV) )
CHGsilva_bray <- phyloseq::distance(ps.rar.CV, method = "bray")
Bray_clust <-hclust(CHGsilva_bray)
adonis2(CHGsilva_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 



# Supplemental Figure 3A and 3G Subject Features Associations with Baseline Microbiome 
ps.CHG.gen.Base<-subset_samples(ps.CHG.gen, Timepoint== "Base")
ps.rarefied.Base = rarefy_even_depth(ps.CHG.gen.Base, rngseed=1, sample.size=750, replace=F)
ps.rare.B.Viable = subset_samples(ps.rarefied.Base, PMAxx =="PMAxx")
GP = ps.rare.B.Viable 
GP.ord <- ordinate(GP, "NMDS",  "bray") 
plot_ordination(GP, GP.ord, type="samples", color="BodySite")+
  geom_point(size=0.1) + 
  stat_ellipse(type = "t", level = 0.8)+
  geom_point(aes(color =BodySite, fill = BodySite, shape = Site), size = 5.5, alpha = 0.7)+
  scale_shape_manual(values = c(21,23, 24,25,22))+
  scale_color_manual(values = c("#EF7234","#BB4855", "#E9B10C","#D55E45", "#483470","#843369"))+
  scale_fill_manual(values = c("#EF7234","#BB4855", "#E9B10C","#D55E45", "#483470","#843369"))+
  theme_bw()+
  ggtitle("Bray Curtis NMDS; VIABLE baseline (all)")
plot_ordination(GP, GP.ord, type="samples", color="Gender")+# shape = "EdgeCenter") + 
  geom_point(size=0.1) + 
  stat_ellipse(type = "t", level = 0.8)+
  geom_point(aes(color =Gender, fill = Gender, shape = Gender), size = 5.5, alpha = 0.7)+
  scale_shape_manual(values = c(21,23, 24,25,22))+
  scale_color_manual(values = c("#BB4855", "#483470"))+
  scale_fill_manual(values = c("#BB4855", "#483470"))+
  theme_bw()+
  ggtitle("Bray Curtis NMDS; VIABLE surgery site (all)")

# Supplemental Table 5 - PERMANOVAS for subject factors associations with baseline microbial communities 
ps.rare.B.Viable = subset_samples(ps.rarefied.Base, PMAxx =="PMAxx")
sampledf <- data.frame(sample_data(ps.rare.B.Viable) )
CHGsilva_bray <- phyloseq::distance(ps.rare.B.Viable, method = "bray")

Bray_clust <-hclust(CHGsilva_bray)

vegan::adonis2(CHGsilva_bray ~ Age, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ MoistDry, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Site, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ PMAxx, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ BodySite, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ BMI, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ SmokingYN, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ SignificantAlcohol, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ PostQ1.SkinTypeODC, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ PostQ2.HowFrequentlyDoYouBathe, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ FrequencyOfBaths.NOperWeek, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ PostQ3.HoursSinceLastShower, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Acne, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ AtopicDermatitis, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ HxMRSA, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ PostHW6.UrbanSuburbanRural, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ FarmAnimalExposure, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ AnyAnimalExposure, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Cats, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Dogs, by= "margin", data = sampledf, permutations = 9999)


adonis2(CHGsilva_bray ~ Gender+BodySite + Site, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite,  by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ BodySite  + Gender, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+ BodySite+ MoistDry, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + Age, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + BMI, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + SmokingYN, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + SignificantAlcohol, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + PostQ1.SkinTypeODC, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + PostQ2.HowFrequentlyDoYouBathe, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + FrequencyOfBaths.NOperWeek, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + PostQ3.HoursSinceLastShower, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + Acne, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + AtopicDermatitis, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + HxMRSA, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + PostHW6.UrbanSuburbanRural, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + FarmAnimalExposure, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + Cats, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + Dogs, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+BodySite  + AnyAnimalExposure, by= "terms", data = sampledf, permutations = 9999)


# Supplemental table 6 Subject factors associated with surgery day microbiome - whats driving differences in Viable communities the day of surgery 
ps.rar.V <- subset_samples(ps.rarefiedsmall, PMAxx == "PMAxx")
ps.rar.SurgDay <- subset_samples(ps.rar.V, BC.SurgD == "Surgery Day")
ps.rar.SurgDay.S <- subset_samples(ps.rar.SurgDay, Site != "Control")
ps.rar.SurgDay.C <- subset_samples(ps.rar.SurgDay, Site == "Control")

sampledf <- data.frame(sample_data(ps.rar.SurgDay.S ))
CHGsilva_bray <- phyloseq::distance(ps.rar.SurgDay.S, method = "bray")
Bray_clust <-hclust(CHGsilva_bray)

adonis2(CHGsilva_bray ~ BodySite, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ MoistDry, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ SurgeryType , by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ SurgicalSiteInfection, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ AntibioticProphylaxisYN, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ HairRemoval, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ Antiseptic, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ Open.Lap, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~LengthOfSurgery.IncisionToClose , by= "margin", data = sampledf, permutations = 9999) # 

adonis2(CHGsilva_bray ~ Gender+ BodySite +HairRemoval, by= "terms", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ Gender+ BodySite +SurgeryType, by= "terms", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ Gender+ BodySite +SurgicalSiteInfection, by= "terms", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ Gender+ BodySite +AntibioticProphylaxisYN, by= "terms", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ Gender+ BodySite +Open.Lap, by= "terms", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ Gender+ BodySite +AntibioticProphylaxisYN+Open.Lap, by= "terms", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ Gender+ BodySite +LengthOfSurgery.IncisionToClose, by= "terms", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ Gender+ BodySite:MoistDry, by= "terms", data = sampledf, permutations = 9999) # 

sampledf <- data.frame(sample_data(ps.rar.SurgDay.C ))
CHGsilva_bray <- phyloseq::distance(ps.rar.SurgDay.C, method = "bray")
Bray_clust <-hclust(CHGsilva_bray)

adonis2(CHGsilva_bray ~ BodySite, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ SurgeryType , by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ SurgicalSiteInfection, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ AntibioticProphylaxisYN, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ HairRemoval, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ Antiseptic, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ Open.Lap, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~LengthOfSurgery.IncisionToClose , by= "margin", data = sampledf, permutations = 9999) # 



# Surgical Site Microbial communities 
ps.rar.Surg<- subset_samples(ps.rarefiedsmall, Site == "SurgeSite")
ps.rar.SurgV<- subset_samples(ps.rar.Surg, PMAxx == "PMAxx")
ps.rar.SurgT<- subset_samples(ps.rar.Surg, PMAxx != "PMAxx")
ps.rar.SurgV.M <- subset_samples(ps.rar.SurgV, MoistDry == "Moist")
ps.rar.SurgV.D <- subset_samples(ps.rar.SurgV, MoistDry == "Dry")

# Figure 4A - viable surgical site microbiome over time 
GP = ps.rar.SurgV
GP.ord <- ordinate(GP, "NMDS",  "bray")
plot_ordination(GP, GP.ord, type="samples", color="Timepoint")+# shape = "EdgeCenter") + 
  geom_point(size=0.1) + 
  stat_ellipse(type = "t", level = 0.9)+
  geom_point(aes(color =Timepoint, fill = Timepoint, shape = BodySite), size = 5.5, alpha = 0.7)+
  scale_shape_manual(values = c(24,25,21,22))+
  scale_color_manual(values = c("#483470","#843369","#EF7234", "#E9B10C", "#D55D45"))+
  scale_fill_manual(values = c("#483470","#843369","#EF7234", "#E9B10C", "#D55D45"))+
  theme_bw()+
  ggtitle("Bray Curtis NMDS; VIABLE surgery site (all)")

### Figure 4B  vector plot for the viable surgery site bray curtis NMDS
GP = ps.rar.SurgV
GP.ord <- ordinate(GP, "NMDS",  "bray")
BrayDF <- as.data.frame(GP.ord$points)
BrayPositions <- subset.data.frame(BrayDF, select = c("MDS1", "MDS2")) %>% as_tibble(rownames = "Sample")
SurgV_OUT <-as.data.frame(psmelt(ps.rar.SurgV)) # this is the CSV with the full metadata and the count of the abundance for each otu in each sample
OTUtaxa <- subset.data.frame(SurgV_OUT , select = c("OTU", "Phylum", "Genus"))
Brayshaired<- merge(SurgV_OUT, BrayPositions, by = "Sample") # merging the two dataframes to get corelations 
cor_x <- Brayshaired %>%
  nest(data = -OTU) %>%
  mutate(cor_x = map(data, ~ cor.test(.x$Abundance, .x$MDS1, method = "spearman",exact = FALSE) %>% tidy())) %>%
  select(OTU, cor_x) %>%
  unnest(cor_x) %>%
  select(OTU, estimate, p.value)
cor_y  <- Brayshaired %>%
  nest(data = -OTU) %>%
  mutate(cor_y = map(data, ~ cor.test(.x$Abundance, .x$MDS2, method = "spearman",exact = FALSE) %>% tidy())) %>%
  select(OTU, cor_y) %>%
  unnest(cor_y) %>%
  select(OTU, estimate, p.value)
correlations <- inner_join(cor_x, cor_y, by = "OTU") 
filtered <- correlations %>%
  filter(p.value.x < 0.1 | p.value.y < 0.1)
corTaxa <- merge(filtered, OTUtaxa, by = "OTU") 
corTaxa <- unique(corTaxa)
corTaxa # THIS IS THE Data fram with the list of the significant OTU
OTUtaxaAnnotations <- subset.data.frame(corTaxa, select = c("OTU", "Phylum", "Genus"))
high_corr <- filter(corTaxa, abs(estimate.x >0.001) | abs(estimate.y>0.001))
high_corr%>% 
  ggplot(aes(x=0, xend=estimate.x, y=0, yend =estimate.y)) +
  geom_segment() 
BrayPositions %>%
  ggplot(aes(x=MDS1, y = MDS2)) +
  geom_point()+
  geom_segment(data = high_corr, aes(x=0, xend=estimate.x, y=0, yend =estimate.y))+
  geom_text_repel(data = high_corr, aes(x=estimate.x, y =estimate.y, label = OTU), min.segment.length = 0.01, max.overlaps = 15, size = 2) +
  theme_bw()
BrayPositions %>%
  ggplot(aes(x=MDS1, y = MDS2)) +
  geom_point()+
  geom_segment(data = high_corr, aes(x=0, xend=estimate.x, y=0, yend =estimate.y, color = Genus), size= 1.1)+
  geom_text_repel(data = high_corr, aes(x=estimate.x, y =estimate.y, label = Genus), min.segment.length = 0.01, max.overlaps = 15, size = 3) +
  theme_bw()
BrayPositions %>%
  ggplot(aes(x =NMDS1, y = NMDS2)) +
  geom_segment(data = high_corr, aes(x=0, xend=estimate.x, y=0, yend =estimate.y), size= 1.1)+
  geom_text_repel(data = high_corr, aes(x=estimate.x, y =estimate.y, label = Genus), min.segment.length = 0.01, max.overlaps = 15, size = 3) +
  scale_y_continuous(limits = c(-0.6,.6))+
  scale_x_continuous(limits = c(-0.6,.6))+
  theme_bw()


# Supplemental Table 7 (code part 1/2)- surgical site samples 
sampledf <- data.frame(sample_data(ps.rar.SurgV) )
CHGsilva_bray <- phyloseq::distance(ps.rar.SurgV, method = "bray")
adonis2(CHGsilva_bray ~ Timepoint, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ BodySite, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ MoistDry, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ BC.SurgD, by= "margin", data = sampledf, permutations = 9999)
adonis2(CHGsilva_bray ~ AntibioticProphylaxisYN, by= "margin", data = sampledf, permutations = 9999)

adonis2(CHGsilva_bray ~ Gender+BodySite +AntibioticProphylaxisYN+ Timepoint, by= "term", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~Gender+BodySite +AntibioticProphylaxisYN+ BC.SurgD, by= "term", data = sampledf, permutations = 9999)

adonis2(CHGsilva_bray ~ Gender+BodySite + Timepoint, by= "term", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~Gender+BodySite + BC.SurgD, by= "term", data = sampledf, permutations = 9999)

ps.rar.surgV.BFU <- subset_samples(ps.rar.SurgV, Timepoint %in% c("Base", "CV1"))
CHGsilva_bray <- phyloseq::distance(ps.rar.surgV.BFU , method = "bray")
sampledf <- data.frame(sample_data(ps.rar.surgV.BFU ) )
adonis2(CHGsilva_bray ~ Timepoint, by= "margin", data = sampledf, permutations = 9999) # difference between baseline and followup
adonis2(CHGsilva_bray ~ Gender+BodySite + Timepoint, by= "term", data = sampledf, permutations = 9999) 

sampledf <- data.frame(sample_data(ps.rar.SurgT) )
CHGsilva_bray <- phyloseq::distance(ps.rar.SurgT, method = "bray")
adonis2(CHGsilva_bray ~ Timepoint, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ BodySite, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ MoistDry, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ BC.SurgD, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ AntibioticProphylaxisYN, by= "margin", data = sampledf, permutations = 9999) # 

adonis2(CHGsilva_bray ~ Gender+BodySite +AntibioticProphylaxisYN+ Timepoint, by= "term", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~Gender+BodySite +AntibioticProphylaxisYN+ BC.SurgD, by= "term", data = sampledf, permutations = 9999)

adonis2(CHGsilva_bray ~ Gender+BodySite + Timepoint, by= "term", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~Gender+BodySite + BC.SurgD, by= "term", data = sampledf, permutations = 9999)

ps.rar.surgT.BFU <- subset_samples(ps.rar.SurgT, Timepoint %in% c("Base", "CV1"))
CHGsilva_bray <- phyloseq::distance(ps.rar.surgT.BFU , method = "bray")
sampledf <- data.frame(sample_data(ps.rar.surgT.BFU ) )
adonis2(CHGsilva_bray ~ Timepoint, by= "margin", data = sampledf, permutations = 9999) # difference between baseline and followup
adonis2(CHGsilva_bray ~ Gender+BodySite + Timepoint, by= "term", data = sampledf, permutations = 9999) 



# Supplemental Figure 8 A-B - Control Site
ps.rar.Con<- subset_samples(ps.rarefiedsmall, Site == "Control")
ps.rar.ConV<- subset_samples(ps.rar.Con, PMAxx == "PMAxx")
ps.rar.ConT<- subset_samples(ps.rar.Con, PMAxx != "PMAxx")

GP = ps.rar.ConV
GP.ord <- ordinate(GP, "NMDS",  "bray")
plot_ordination(GP, GP.ord, type="samples", color="Timepoint")+# shape = "EdgeCenter") + 
  geom_point(size=0.1) + 
  stat_ellipse(type = "t", level = 0.8)+
  geom_point(aes(color =Timepoint, fill = Timepoint, shape = BodySite), size = 5.5, alpha = 0.7)+
  scale_shape_manual(values = c(21,22))+
  scale_color_manual(values = c("#483470","#843369","#EF7234", "#E9B10C", "#D55D45"))+
  scale_fill_manual(values = c("#483470","#843369","#EF7234", "#E9B10C", "#D55D45"))+
  theme_bw()+
  ggtitle("Bray Curtis NMDS; Viable DNA CONTROL site (all)")

plot_ordination(GP, GP.ord, type="samples", color="BodySite")+# shape = "EdgeCenter") + 
  geom_point(size=0.1) + 
  stat_ellipse(type = "t", level = 0.85)+
  geom_point(aes(color =BodySite, fill = BodySite, shape = BodySite), size = 5.5, alpha = 0.7)+
  scale_shape_manual(values = c(21,22))+
  scale_color_manual(values = c("#483470", "#D55D45"))+
  scale_fill_manual(values = c("#483470", "#D55D45"))+
  theme_bw()+
  ggtitle("Bray Curtis NMDS; Viable DNA CONTROL site (all)")

# supplemental Table 7 - part 2/2 - Control site microbiome
sampledf <- data.frame(sample_data(ps.rar.ConV) )
CHGsilva_bray <- phyloseq::distance(ps.rar.ConV, method = "bray")
adonis2(CHGsilva_bray ~ Timepoint, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ BodySite, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ MoistDry, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ BC.SurgD, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ AntibioticProphylaxisYN, by= "margin", data = sampledf, permutations = 9999) # 

adonis2(CHGsilva_bray ~ Gender+BodySite +AntibioticProphylaxisYN+ Timepoint, by= "term", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~Gender+BodySite +AntibioticProphylaxisYN+ BC.SurgD, by= "term", data = sampledf, permutations = 9999)

adonis2(CHGsilva_bray ~ Gender+BodySite + Timepoint, by= "term", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~Gender+BodySite + BC.SurgD, by= "term", data = sampledf, permutations = 9999)

ps.rar.conV.BFU <- subset_samples(ps.rar.ConV, Timepoint %in% c("Base", "CV1"))
CHGsilva_bray <- phyloseq::distance(ps.rar.conV.BFU , method = "bray")
sampledf <- data.frame(sample_data(ps.rar.conV.BFU ) )
adonis2(CHGsilva_bray ~ Timepoint, by= "margin", data = sampledf, permutations = 9999) # difference between baseline and followup
adonis2(CHGsilva_bray ~ Gender+BodySite + Timepoint, by= "term", data = sampledf, permutations = 9999) 

sampledf <- data.frame(sample_data(ps.rar.ConT) )
CHGsilva_bray <- phyloseq::distance(ps.rar.ConT, method = "bray")
adonis2(CHGsilva_bray ~ Timepoint, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ BodySite, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ BC.SurgD, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~ AntibioticProphylaxisYN, by= "margin", data = sampledf, permutations = 9999) # 

adonis2(CHGsilva_bray ~ Gender+BodySite +AntibioticProphylaxisYN+ Timepoint, by= "term", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~Gender+BodySite +AntibioticProphylaxisYN+ BC.SurgD, by= "term", data = sampledf, permutations = 9999)

adonis2(CHGsilva_bray ~ Gender+BodySite + Timepoint, by= "term", data = sampledf, permutations = 9999) # 
adonis2(CHGsilva_bray ~Gender+BodySite + BC.SurgD, by= "term", data = sampledf, permutations = 9999)

ps.rar.conT.BFU <- subset_samples(ps.rar.ConT, Timepoint %in% c("Base", "CV1"))
CHGsilva_bray <- phyloseq::distance(ps.rar.conT.BFU , method = "bray")
sampledf <- data.frame(sample_data(ps.rar.conT.BFU ) )
adonis2(CHGsilva_bray ~ Timepoint, by= "margin", data = sampledf, permutations = 9999) # difference between baseline and followup
adonis2(CHGsilva_bray ~ Gender+BodySite + Timepoint, by= "term", data = sampledf, permutations = 9999) 



# Supplemental Figure 7 A -  Change in Bray Curtis beta diversity VIABLE ONLY at each timepoint from baseline 
ps.rare.S.viable <- subset_samples(ps.rarefiedsmall, PMAxx == "PMAxx")
delta_bray <- phyloseq::distance(ps.rare.S.viable, method = "bray")
delta_bray_Mat <- as.matrix(delta_bray)
sampledf <- data.frame(sample_data(ps.rare.S.viable))
sampledf$SampleID <- row.names(sampledf)

delta_BrayD<- as.data.frame(as.table(delta_bray_Mat)) 
#delta_BrayD <- delta_BrayD %>% rename("Var1"= "SampleID", 
                                     # "Var2"="SampleB" ,
                                     # "Freq" ="BrayCurtisDistance")
delta_BrayD <- delta_BrayD %>% rename("SampleID" = "Var1", 
                                     "SampleB" = "Var2",
                                     "BrayCurtisDistance" = "Freq")
delta_BrayDF <-merge(delta_BrayD, sampledf, by = "SampleID", all = T)
delta_BrayDF <- subset.data.frame(delta_BrayDF, Timepoint %in% c("Base", "Pre", "OR", "Post","CV1"))
delta_BrayDF$SubjectSCpma <- paste(delta_BrayDF$Subject, delta_BrayDF$Site, delta_BrayDF$PMAxx)

Subject.001 <- subset.data.frame(delta_BrayDF, SubjectSCpma == "CHG-001 SurgeSite PMAxx" )
Base.samp<- unique(droplevels(Subject.001$SampleID[Subject.001$Timepoint =="Base"]))
Subject.001.toBase<-subset.data.frame(Subject.001, SampleB %in% Base.samp)
delta_Bray_toTVBaseline<-Subject.001.toBase
SubSCPx.lst <- unique(delta_BrayDF$SubjectSCpma)
for (SubSCPx in SubSCPx.lst){
  Subject.X <- subset.data.frame(delta_BrayDF, SubjectSCpma == SubSCPx)
  Base.samp<- unique(droplevels(Subject.X$SampleID[Subject.X$Timepoint =="Base"]))
  Subject.X.toBase<-subset.data.frame(Subject.X, SampleB %in% Base.samp)
  delta_Bray_toTVBaseline<-rbind(delta_Bray_toTVBaseline,Subject.X.toBase)
}
delta_Bray_toViabBaseline<- delta_Bray_toTVBaseline
delta_Bray_toViabBaseline <- subset.data.frame(delta_Bray_toViabBaseline, !SampleID %in% c("LK16S011-043","LK16S011-063" )) # removing CHG-016 lower abdomen inquinal samples 
delta_Bray_toViabBaseline <- subset.data.frame(delta_Bray_toViabBaseline, !SampleB %in% c("LK16S011-043","LK16S011-063" ))

delta_Bray_toViabBaseline$Timepoint <- factor(delta_Bray_toViabBaseline$Timepoint, levels = c("Base", "Pre", "OR", "Post","CV1"))
delta_Bray_toViabBaseline$SiteMD <- paste(delta_Bray_toViabBaseline$MoistDry, delta_Bray_toViabBaseline$Site) 

#write.csv(delta_Bray_toViabBaseline, "./BrayCHANGE_fromRespectiveViableCSS_BASELINE.csv")

ggplot(delta_Bray_toViabBaseline , aes(x = Timepoint, y = BrayCurtisDistance, fill = Site, color = Site))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.15))+
  stat_summary(fun = "median", linewidth = 1, geom = "line", aes(group = Site))+
  facet_grid(cols = vars(Site), scales = "free_x")+
  #geom_hline(yintercept=0.625, color = "gray")+
  ylab("Degree of Change in Bray Curtis From Baseline")+
  scale_x_discrete(labels = c("Baseline","Pre-OR", "OR", "Post-OR","Follow-Up"))+
  theme_light()+
  scale_color_manual(values = c("#AE8509","#86323B"))+
  scale_fill_manual(values = c("#E9B10C", "#BB4855"))+
  #theme(axis.text.x = element_text(size = 7)) +
  ggtitle("Change in bray curtis community from subject V baseline sample")


# Supplemental Figure 7B CV to BASELINE CHange in beta diversity from SELF vs others 
### looking to see if someone is more like someone else's baseline (so taking the MINIMUM distance between one subjects Cv and anyone elses baseline sample)
ps.rare.S.viable <- subset_samples(ps.rarefiedsmall, PMAxx == "PMAxx")
ps.rare.S.viable.BFU <- subset_samples(ps.rare.S.viable, Timepoint %in% c("Base", "CV1"))
delta_bray <- phyloseq::distance(ps.rare.S.viable.BFU, method = "bray")
delta_bray_Mat <- as.matrix(delta_bray)
sampledf <- data.frame(sample_data(ps.rare.S.viable.BFU))
sampledf$SampleID <- row.names(sampledf)

delta_BrayD<- as.data.frame(as.table(delta_bray_Mat)) 
#delta_BrayD <- delta_BrayD %>% rename("Var1"= "SampleID", 
# "Var2"="SampleB" ,
# "Freq" ="BrayCurtisDistance")
delta_BrayD <- delta_BrayD %>% rename("SampleCV" = "Var1", 
                                      "SampleBase" = "Var2",
                                      "BrayCurtisDistance" = "Freq")
delta_BrayDF <-merge(delta_BrayD, sampledf, by.x = "SampleCV", by.y = "SampleID", all = T)
delta_BrayDF$SubjectSCpma <- paste(delta_BrayDF$Subject, delta_BrayDF$Site, delta_BrayDF$PMAxx)
delta_BrayDF$SurgConPma <- paste(delta_BrayDF$Site, delta_BrayDF$PMAxx)

Subject.001 <- subset.data.frame(delta_BrayDF, SubjectSCpma == "CHG-001 SurgeSite PMAxx" ) # to self baseline 
Base.samp<- unique(droplevels(Subject.001$SampleCV[Subject.001$Timepoint =="Base"]))
Subject.001.toBase<-subset.data.frame(Subject.001, SampleBase %in% Base.samp)
Subject.001.toBase<-subset.data.frame(Subject.001.toBase, Timepoint == "CV1")
delta_Bray_toTVBaseline<-Subject.001.toBase
SubSCPx.lst <- unique(delta_BrayDF$SubjectSCpma)
for (SubSCPx in SubSCPx.lst){
  Subject.X <- subset.data.frame(delta_BrayDF, SubjectSCpma == SubSCPx)
  Base.samp<- unique(droplevels(Subject.X$SampleCV[Subject.X$Timepoint =="Base"]))
  Subject.X.toBase<-subset.data.frame(Subject.X, SampleBase %in% Base.samp)
  Subject.X.toBase<-subset.data.frame(Subject.X.toBase, Timepoint == "CV1")
  delta_Bray_toTVBaseline<-rbind(delta_Bray_toTVBaseline,Subject.X.toBase)
}
delta_Bray_toViabBaseline<- delta_Bray_toTVBaseline
delta_Bray_toViabBaseline <- subset.data.frame(delta_Bray_toViabBaseline, !SampleCV %in% c("LK16S011-043","LK16S011-063" )) # removing CHG-016 lower abdomen inquinal samples 
delta_Bray_toViabBaseline <- subset.data.frame(delta_Bray_toViabBaseline, !SampleBase %in% c("LK16S011-043","LK16S011-063" ))

delta_Bray_toViabBaseline$Timepoint <- factor(delta_Bray_toViabBaseline$Timepoint, levels = c("Base", "CV1"))
delta_Bray_toViabBaseline$SiteMD <- paste(delta_Bray_toViabBaseline$MoistDry, delta_Bray_toViabBaseline$Site) 
SubSelf.PMA.CV.toBase.SELF<-delta_Bray_toViabBaseline

SampB.Base.DF<- select(sampledf, c(SampleID, Subject, Timepoint,Site)) # combining it all 
SampB.Base.DF <- SampB.Base.DF %>% rename("SampleBase" = "SampleID", 
                                          "SubjectBase" = "Subject",
                                          "TimepointB" = "Timepoint",
                                          "SiteB"= "Site")
Viable.CV.toBase.SELF.DF<-merge(SubSelf.PMA.CV.toBase.SELF, SampB.Base.DF , by = "SampleBase", all = F)
Viable.CV.toBase.SELF.DF<-Viable.CV.toBase.SELF.DF[Viable.CV.toBase.SELF.DF$Site==Viable.CV.toBase.SELF.DF$SiteB, ]
Viable.CV.toBase.SELF.DF$SelfElse <- with(Viable.CV.toBase.SELF.DF, ifelse(Subject == SubjectBase, "to Self Baseline", "to Someone Elses Baseline"))


SS.Samples <- subset.data.frame(delta_BrayDF, Site =="SurgeSite")
All.Base.SS.samples <- unique(droplevels(SS.Samples$SampleCV[SS.Samples$Timepoint =="Base"])) # b/c right now the metafile is associated w/ the sample in the SampleCV col, and this col currently includes all samples

Subject.001 <- subset.data.frame(SS.Samples, SubjectSCpma == "CHG-001 SurgeSite PMAxx" ) # to SS EVERYONEELSE baseline 
Base.samp<- unique(droplevels(Subject.001$SampleCV[Subject.001$Timepoint =="Base"]))
Subject.001.toBa<-subset.data.frame(Subject.001, SampleBase %in% All.Base.SS.samples)
Subject.001.toBase<-subset.data.frame(Subject.001.toBa, !SampleBase %in% Base.samp)
Subject.001.toBase<-subset.data.frame(Subject.001.toBase, Timepoint == "CV1")
delta_Bray_toTVBaseline<-Subject.001.toBase
SubSCPx.lst <- unique(SS.Samples$SubjectSCpma)
for (SubSCPx in SubSCPx.lst){
  Subject.X <- subset.data.frame(SS.Samples, SubjectSCpma == SubSCPx)
  Base.samp<- unique(droplevels(Subject.X$SampleCV[Subject.X$Timepoint =="Base"]))
  Subject.X.toBa<-subset.data.frame(Subject.X, SampleBase %in% All.Base.SS.samples)
  Subject.X.toBase<-subset.data.frame(Subject.X.toBa, !SampleBase %in% Base.samp)
  Subject.X.toBase<-subset.data.frame(Subject.X.toBase, Timepoint == "CV1")
  delta_Bray_toTVBaseline<-rbind(delta_Bray_toTVBaseline,Subject.X.toBase)
}
delta_Bray_toViabBaseline<- delta_Bray_toTVBaseline
delta_Bray_toViabBaseline <- subset.data.frame(delta_Bray_toViabBaseline, !SampleCV %in% c("LK16S011-043","LK16S011-063" )) # removing CHG-016 lower abdomen inquinal samples 
delta_Bray_toViabBaseline <- subset.data.frame(delta_Bray_toViabBaseline, !SampleBase %in% c("LK16S011-043","LK16S011-063" ))
delta_Bray_toViabBaseline.SSe <- delta_Bray_toViabBaseline

SampB.Base.DF<- select(sampledf, c(SampleID, Subject, Timepoint,Site)) # combining it all 
SampB.Base.DF <- SampB.Base.DF %>% rename("SampleBase" = "SampleID", 
                                          "SubjectBase" = "Subject",
                                          "TimepointB" = "Timepoint",
                                          "SiteB"= "Site")
Viable.CV.toBase.SS.ELSE.DF<-merge(delta_Bray_toViabBaseline.SSe, SampB.Base.DF , by = "SampleBase", all = F)
Viable.CV.toBase.SS.ELSE.DF<-Viable.CV.toBase.SS.ELSE.DF[Viable.CV.toBase.SS.ELSE.DF$Site==Viable.CV.toBase.SS.ELSE.DF$SiteB, ]
Viable.CV.toBase.SS.ELSE.DF.min<-Viable.CV.toBase.SS.ELSE.DF %>% group_by(SampleCV) %>% slice(which.min(BrayCurtisDistance))
Viable.CV.toBase.SS.ELSE.DF.min$Timepoint <- factor(Viable.CV.toBase.SS.ELSE.DF.min$Timepoint, levels = c("Base", "CV1"))
Viable.CV.toBase.SS.ELSE.DF.min$SiteMD <- paste(Viable.CV.toBase.SS.ELSE.DF.min$MoistDry, Viable.CV.toBase.SS.ELSE.DF.min$Site) 
Viable.CV.toBase.SS.ELSE.DF.min$SelfElse <- with(Viable.CV.toBase.SS.ELSE.DF.min, ifelse(Subject == SubjectBase, "to Self Baseline", "to Someone Elses Baseline"))


C.Samples <- subset.data.frame(delta_BrayDF, Site =="Control")
All.Base.C.samples <- unique(droplevels(C.Samples$SampleCV[C.Samples$Timepoint =="Base"])) # b/c right now the metafile is associated w/ the sample in the SampleCV col, and this col currently includes all samples

Subject.001 <- subset.data.frame(C.Samples, SubjectSCpma == "CHG-001 Control PMAxx" ) # to SS EVERYONEELSE baseline 
Base.samp<- unique(droplevels(Subject.001$SampleCV[Subject.001$Timepoint =="Base"]))
Subject.001.toBase<-subset.data.frame(Subject.001, SampleBase %in% All.Base.C.samples)
Subject.001.toBase<-subset.data.frame(Subject.001.toBase, !SampleBase %in% Base.samp)
Subject.001.toBase<-subset.data.frame(Subject.001.toBase, Timepoint == "CV1")
delta_Bray_toTVBaseline<-Subject.001.toBase
SubSCPx.lst <- unique(C.Samples$SubjectSCpma)
for (SubSCPx in SubSCPx.lst){
  Subject.X <- subset.data.frame(C.Samples, SubjectSCpma == SubSCPx)
  Base.samp<- unique(droplevels(Subject.X$SampleCV[Subject.X$Timepoint =="Base"]))
  Subject.X.toBase<-subset.data.frame(Subject.X, SampleBase %in% All.Base.C.samples)
  Subject.X.toBase<-subset.data.frame(Subject.X.toBase, !SampleBase %in% Base.samp)
  Subject.X.toBase<-subset.data.frame(Subject.X.toBase, Timepoint == "CV1")
  delta_Bray_toTVBaseline<-rbind(delta_Bray_toTVBaseline,Subject.X.toBase)
}
delta_Bray_toViabBaseline<- delta_Bray_toTVBaseline
delta_Bray_toViabBaseline <- subset.data.frame(delta_Bray_toViabBaseline, !SampleCV %in% c("LK16S011-043","LK16S011-063" )) # removing CHG-016 lower abdomen inquinal samples 
delta_Bray_toViabBaseline <- subset.data.frame(delta_Bray_toViabBaseline, !SampleBase %in% c("LK16S011-043","LK16S011-063" ))
delta_Bray_toViabBaseline.Ce <- delta_Bray_toViabBaseline

SampB.Base.DF<- select(sampledf, c(SampleID, Subject, Timepoint,Site)) # combining it all 
SampB.Base.DF <- SampB.Base.DF %>% rename("SampleBase" = "SampleID", 
                                          "SubjectBase" = "Subject",
                                          "TimepointB" = "Timepoint",
                                          "SiteB"= "Site")
Viable.CV.toBase.C.ELSE.DF<-merge(delta_Bray_toViabBaseline.Ce, SampB.Base.DF , by = "SampleBase", all = F)
Viable.CV.toBase.C.ELSE.DF<-Viable.CV.toBase.C.ELSE.DF[Viable.CV.toBase.C.ELSE.DF$Site==Viable.CV.toBase.C.ELSE.DF$SiteB, ]
Viable.CV.toBase.C.ELSE.DF.min<-Viable.CV.toBase.C.ELSE.DF %>% group_by(SampleCV) %>% slice(which.min(BrayCurtisDistance))
Viable.CV.toBase.C.ELSE.DF.min$Timepoint <- factor(Viable.CV.toBase.C.ELSE.DF.min$Timepoint, levels = c("Base", "CV1"))
Viable.CV.toBase.C.ELSE.DF.min$SiteMD <- paste(Viable.CV.toBase.C.ELSE.DF.min$MoistDry, Viable.CV.toBase.C.ELSE.DF.min$Site) 
Viable.CV.toBase.C.ELSE.DF.min$SelfElse <- with(Viable.CV.toBase.C.ELSE.DF.min, ifelse(Subject == SubjectBase, "to Self Baseline", "to Someone Elses Baseline"))

Viable.CV.toBase.SELF.ELSE.DF<-rbind(Viable.CV.toBase.SELF.DF,Viable.CV.toBase.C.ELSE.DF.min)
Viable.CV.toBase.SELF.ELSE.DF<-rbind(Viable.CV.toBase.SELF.ELSE.DF,Viable.CV.toBase.SS.ELSE.DF.min)

ggplot(Viable.CV.toBase.SELF.ELSE.DF, aes(x = SelfElse, y = BrayCurtisDistance))+
  geom_boxplot(aes(color = SelfElse, fill = SelfElse), linewidth = 1)+
  geom_point(aes(color = SelfElse), position=position_jitter(width=0.15))+
  facet_grid(cols =  vars(Site), scales = "free_x")+
  theme_light()+
  scale_color_manual(values = c("#241A38", "#AE8509"))+
  scale_fill_manual(values = c("#483470", "#E9B10C"))+
  #theme(axis.text.x = element_text(size = 7)) +
  ggtitle("Change in bray curtis community from subject V baseline sample")

#write.csv(./FollowUP_to_SELF_ELSE_baseline.csv" )



## MAASLIN 
### maaslin for JUST the PMAxx treated samples only  (viable microbiome) 
input_genus <- read.csv("./CHGsilva_RelativeOUT.genus.csv")
input_genus_data <- as.data.frame(input_genus )
rownames(input_genus_data) <-input_genus_data[,1]

input_meta_file <-CHG_Study.meta
input_meta_file <-subset.data.frame(input_meta_file, Timepoint != "NegEx")
input_meta_file$Timepoint[input_meta_file$Timepoint == "Pre2"] <- "Pre"
input_meta_file$Timepoint[input_meta_file$Timepoint == "OR2"] <- "OR"
input_meta_file$Timepoint[input_meta_file$Timepoint == "Post2"] <- "Post"

Surg.Day.meta<- subset.data.frame(input_meta_file, BC.SurgD == "Surgery Day")

PMAxx_meta_file <- subset.data.frame(input_meta_file, PMAxx == "PMAxx")
SSpma_MF <- subset.data.frame(PMAxx_meta_file, Site =="SurgeSite")
V.MoistS_MF <- subset.data.frame(SSpma_MF, MoistDry == "Moist" )
V.DryS_MF <- subset.data.frame(SSpma_MF, MoistDry == "Dry")
SSpma_MF$Timepoint <- factor(SSpma_MF$Timepoint, levels = c("Base", "Pre","OR", "Post", "CV1"))

V.MoistS_MF$Timepoint <- factor(V.MoistS_MF$Timepoint, levels = c("Base", "Pre","OR", "Post", "CV1"))
V.DryS_MF$Timepoint <- factor(V.DryS_MF$Timepoint, levels = c("Base", "Pre","OR", "Post", "CV1"))

# Supplemental Figure 2 difference between PMAxx and none (controled for by subject)
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN/PMAxxV.None/CHG_PMAxxNone_Genus_maslin", 
  fixed_effects = c("PMAxx"),
  reference = c("PMAxx", "None"),
  random_effects = c("SpecificSubjectSample")) # controling for the specific sample it origionated from 

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Surg.Day.meta, 
  output = "./MAASLIN/PMAxxV.None/SurgeryDay_CHG_PMAxxNone_Genus_maslin", 
  min_abundance = 0.01,
  fixed_effects = c("PMAxx"),
  reference = c("PMAxx", "None"),
  random_effects = c("SpecificSubjectSample")) # controling for the specific sample it origionated from 

# Supplemental Figure 3 - Subject features associated with Baseline communities 
PMAxx_Base_meta_file <- subset.data.frame(PMAxx_meta_file, Timepoint == "Base")
fit_data = Maaslin2( #ask about this one forgot where site is and can't find on the graph
  input_data = input_genus_data, 
  input_metadata = PMAxx_Base_meta_file, 
  output = "./MAASLIN/BaselineSubjectFeatures/Gender_AdjBodySite_test", 
  fixed_effects = c("Gender"),
  reference = c("Gender", "M"),
  random_effects = c("BodySite")) # controling for the specific sample it origionated from 
fit_data = Maaslin2( #ask about this one forgot where site is and can't find on the graph
  input_data = input_genus_data, 
  input_metadata = PMAxx_Base_meta_file, 
  output = "./MAASLIN/BaselineSubjectFeatures/Gender_AdjMoistDry", 
  fixed_effects = c("Gender"),
  reference = c("Gender", "M"),
  random_effects = c("MoistDry"))

#Supplemental Figure 4 - antibiotic prophylaxis 
SSpma_MF <- subset.data.frame(PMAxx_meta_file, Site =="SurgeSite")
V.SS_MF.SD <- subset.data.frame(SSpma_MF, BC.SurgD == "Surgery Day")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = V.SS_MF.SD , 
  output = "./MAASLIN/SurgeryDaySurgeryFeatures/AntibioticProphylaxis_adjBodySite+Gender", 
  min_abundance = 0.01,
  fixed_effects = c("Antibiotics"),
  reference = c("Antibiotics", "None"))
  random_effects = c("Bodysite", "Gender"))

#Figure 4 and Supplemental Table 8: Surgical Site - Viable Microbiome 
# Evaluating differences from BASELINE in viable microbiome samples 
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = SSpma_MF, 
  output = "./MAASLIN/SurgViable/CHG_V.ALL_overtime_ControledForSubject+Gender+BodySiteABX_Genus_maslin", 
  fixed_effects = c("Timepoint"),
  reference = c("Timepoint", "Base"),
  random_effects = c("Subject", "BodySite", "Gender", "AntibioticProphylaxisYN"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = V.DryS_MF, 
  output = "./MAASLIN/SurgViable/CHG_V.Dry_overtime_ControledForSubjectBSgenABX_Genus_maslin", 
  fixed_effects = c("Timepoint"),
  reference = c("Timepoint", "Base"),
  random_effects = c("Subject","BodySite", "Gender", "AntibioticProphylaxisYN"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = V.MoistS_MF, 
  output = "./MAASLIN/SurgViable/CHG_V.Moist_overtime_ControledForSubjectBSgenABx_Genus_maslin", 
  fixed_effects = c("Timepoint"),
  reference = c("Timepoint", "Base"),
  random_effects = c("Subject","BodySite", "Gender", "AntibioticProphylaxisYN"))

#Supplemental Figure 8 and supplemental Table 8 Control site - viable microbiome 
## evaluating control site only (PMAxx treated)
Cpma_MF <- subset.data.frame(PMAxx_meta_file, Site =="Control")
V.MoistC_MF <- subset.data.frame(Cpma_MF, MoistDry == "Moist")
V.DryC_MF <- subset.data.frame(Cpma_MF, MoistDry == "Dry")

V.MoistC_MF$Timepoint <- factor(V.MoistC_MF$Timepoint, levels = c("Base", "Pre","OR", "Post", "CV1"))
V.DryC_MF$Timepoint <- factor(V.DryC_MF$Timepoint, levels = c("Base", "Pre","OR", "Post", "CV1"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Cpma_MF, 
  output = "./MAASLIN/ContViable/CHG_V.ALL_C_overtime_ControledForSubject+bodySite+GenderABX_Genus_maslin", 
  fixed_effects = c("Timepoint"),
  reference = c("Timepoint", "Base"),
  random_effects = c("Subject", "BodySite","Gender","AntibioticProphylaxisYN"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = V.DryC_MF, 
  output = "./MAASLIN/ContViable/CHG_V.Dry_C_overtime_ControledForSubjectBSgenABX_Genus_maslin", 
  fixed_effects = c("Timepoint"),
  reference = c("Timepoint", "Base"),
  random_effects = c("Subject","BodySite","Gender","AntibioticProphylaxisYN"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = V.MoistC_MF, 
  output = "./MAASLIN/ContViable/CHG_V.Moist_C_overtime_ControledForSubjectBSgenABX_Genus_maslin", 
  fixed_effects = c("Timepoint"),
  reference = c("Timepoint", "Base"),
  random_effects = c("Subject","BodySite","Gender","AntibioticProphylaxisYN"))

## Total microbiome (no PMAxx)
None_meta_file <- subset.data.frame(input_meta_file, PMAxx != "PMAxx")
SSn_MF <- subset.data.frame(None_meta_file, Site =="SurgeSite")
T.MoistS_MF <- subset.data.frame(SSn_MF, MoistDry == "Moist" )
T.DryS_MF <- subset.data.frame(SSn_MF, MoistDry == "Dry")

T.MoistS_MF$Timepoint <- factor(T.MoistS_MF$Timepoint, levels = c("Base", "Pre","OR", "Post", "CV1"))
T.DryS_MF$Timepoint <- factor(T.DryS_MF$Timepoint, levels = c("Base", "Pre","OR", "Post", "CV1"))

#Supplemental Table 8 - Surgical Site TOTAL MICROBIOME 
# Evaluating differences from BASELINE in TOTAL microbiome samples 
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = SSn_MF, 
  output = "/Users/liztown/Documents/KalanLab/CHG-SSI_Study/16S_Seq/CHGstudy16Sanalysis_FINAL/MAASLIN/SurgTotal/CHG_T.ALL_overtime_ControledForSubjectBodySiteGenderABX_Genus_maslin", 
  min_abundance = 0.01,
  fixed_effects = c("Timepoint"),
  reference = c("Timepoint", "Base"),
  random_effects = c("Subject","BodySite", "Gender", "AntibioticProphylaxisYN"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = T.MoistS_MF, 
  output = "/Users/liztown/Documents/KalanLab/CHG-SSI_Study/16S_Seq/CHGstudy16Sanalysis_FINAL/MAASLIN/SurgTotal/CHG_T.Moist_overtime_ControledForSubjectBSgenABX_Genus_maslin", 
  min_abundance = 0.01,
  fixed_effects = c("Timepoint"),
  reference = c("Timepoint", "Base"),
  random_effects = c("Subject","BodySite", "Gender", "AntibioticProphylaxisYN"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = T.DryS_MF, 
  output = "/Users/liztown/Documents/KalanLab/CHG-SSI_Study/16S_Seq/CHGstudy16Sanalysis_FINAL/MAASLIN/SurgTotal/CHG_T.Dry_overtime_ControledForSubjectBSgenABX_Genus_maslin", 
  min_abundance = 0.01,
  fixed_effects = c("Timepoint"),
  reference = c("Timepoint", "Base"),
  random_effects = c("Subject","BodySite", "Gender", "AntibioticProphylaxisYN"))

#Supplemental Table 8 - Control site TOTAL MICROBIOME 
## evaluating control site only TOTAL 
CNone_MF <- subset.data.frame(None_meta_file, Site =="Control")
T.MoistC_MF <- subset.data.frame(CNone_MF, MoistDry == "Moist")
T.DryC_MF <- subset.data.frame(CNone_MF, MoistDry == "Dry")

T.MoistC_MF$Timepoint <- factor(T.MoistC_MF$Timepoint, levels = c("Base", "Pre","OR", "Post", "CV1"))
T.DryC_MF$Timepoint <- factor(T.DryC_MF$Timepoint, levels = c("Base", "Pre","OR", "Post", "CV1"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = CNone_MF, 
  output = "/Users/liztown/Documents/KalanLab/CHG-SSI_Study/16S_Seq/CHGstudy16Sanalysis_FINAL/MAASLIN/ContTotal/CHG_T.All_C_overtime_ControledForSubjectBodySiteGenderABX_Genus_maslin", 
  min_abundance = 0.01,
  fixed_effects = c("Timepoint"),
  reference = c("Timepoint", "Base"),
  random_effects = c("Subject","BodySite", "Gender","AntibioticProphylaxisYN"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = T.DryC_MF, 
  output = "/Users/liztown/Documents/KalanLab/CHG-SSI_Study/16S_Seq/CHGstudy16Sanalysis_FINAL/MAASLIN/ContTotal/CHG_T.Dry_C_overtime_ControledForSubjectBSgenABX_Genus_maslin", 
  min_abundance = 0.01,
  fixed_effects = c("Timepoint"),
  reference = c("Timepoint", "Base"),
  random_effects = c("Subject", "BodySite","Gender","AntibioticProphylaxisYN"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = T.MoistC_MF, 
  output = "/Users/liztown/Documents/KalanLab/CHG-SSI_Study/16S_Seq/CHGstudy16Sanalysis_FINAL/MAASLIN/ContTotal/CHG_T.Moist_C_overtime_ControledForSubjectBSGenABX_Genus_maslin", 
  min_abundance = 0.01,
  fixed_effects = c("Timepoint"),
  reference = c("Timepoint", "Base"),
  random_effects = c("Subject", "BodySite","Gender","AntibioticProphylaxisYN"))

# Additional MAASLIN plots 
# Subject surgery features associated with Microbial community on day of surgery (viable microbiome only)
SSpma_MF <- subset.data.frame(PMAxx_meta_file, Site =="SurgeSite")
V.SS_MF.SD <- subset.data.frame(SSpma_MF, BC.SurgD == "Surgery Day")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = V.SS_MF.SD , 
  output = "./MAASLIN/SurgeryDaySurgeryFeatures/AntibioticProphylaxis", 
  min_abundance = 0.01,
  fixed_effects = c("Antibiotics"),
  reference = c("Antibiotics", "None"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = V.SS_MF.SD , 
  output = "./MAASLIN/SurgeryDaySurgeryFeatures/AntibioticProphylaxis_adjMoistDry", 
  min_abundance = 0.01,
  fixed_effects = c("Antibiotics"),
  reference = c("Antibiotics", "None"),
  random_effects = c("MoistDry"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = V.SS_MF.SD , 
  output = "./MAASLIN/SurgeryDaySurgeryFeatures/AntibioticProphylaxis_adjMoistDry+Gender+Subject", 
  min_abundance = 0.01,
  fixed_effects = c("Antibiotics"),
  reference = c("Antibiotics", "None"),
  random_effects = c("MoistDry", "Gender", "Subject"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = V.SS_MF.SD , 
  output = "./MAASLIN/SurgeryDaySurgeryFeatures/OpenLap.vLap_adjBodySite+Gender+Subject", 
  min_abundance = 0.01,
  fixed_effects = c("Open.Lap"),
  reference = c("Open.Lap", "Laproscopic"),
  random_effects = c("BodySite", "Gender", "Subject"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = V.SS_MF.SD , 
  output = "./MAASLIN/SurgeryDaySurgeryFeatures/OpenLap.vOpen_adjMoistDry+Gender+Subject", 
  min_abundance = 0.01,
  fixed_effects = c("Open.Lap"),
  reference = c("Open.Lap", "Open"),
  random_effects = c("MoistDry", "Gender", "Subject"))

## MAASLIN for STAPHYLOCCAL ASVs
### maaslin for JUST the PMAxx treated samples only 
relAlls <- microbiome::transform(ps.CHG, "compositional")
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "Pre2"] <- "Pre"
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "OR2"] <- "OR"
relAlls@sam_data$Timepoint[relAlls@sam_data$Timepoint == "Post2"] <- "Post"
StaphSpp2.df <- psmelt(relAlls)
StaphSpp2.df  <- as.data.frame(StaphSpp2.df)
StaphSpp2.df <- subset.data.frame(StaphSpp2.df , Site != "22 Strain Mock")
StaphSpp2.df  <- subset.data.frame(StaphSpp2.df  , Site != "NegEx")


StaphBLAST <- read.csv("./StaphSpeciesBLAST.csv")
StaphBLAST <-as.data.frame(StaphBLAST)
StaphBLAST$BLAST_Species[StaphBLAST$BLAST_Species == "Staphylococcus arueus"] <- "Staphylococcus aureus"

Staph.Spp2.df <- merge(StaphSpp2.df, StaphBLAST, by = "OTU", all = T)
Staph.Spp2.df$BLAST_Species <- factor(Staph.Spp2.df$BLAST_Species, levels = c("Staphylococcus aureus","CoNS Staphylococcus spp.", "Staphylococcus auricularis","Staphylococcus caprae or S. capitis", "Staphylococcus chromogenes ", 
                                                                            "Staphylococcus epidermidis", "Staphylococcus haemolyticus","Staphylococcus hominis", "Staphylococcus lugdunensis","Staphylococcus pasteuri",
                                                                            "Staphylococcus pettenkoferi", "Staphylococcus warneri", "NA"))

Staph.Spp2.df.Top<-Staph.Spp2.df # reclassifying really lowly abundant species to CoNS staph
Staph.Spp2.df.Top$BLAST_Species[Staph.Spp2.df.Top$BLAST_Species== "Staphylococcus auricularis"] <- "CoNS Staphylococcus spp."
Staph.Spp2.df.Top$BLAST_Species[Staph.Spp2.df.Top$BLAST_Species== "Staphylococcus chromogenes "] <- "CoNS Staphylococcus spp."

#write.csv(Staph.Spp2.df.Top, "./CHG_StaphBlast_Merged.csv")


input_Staph <- read.csv("./CHG_StaphBlast_Merged_pivot.csv")
input_Staph_data <- as.data.frame(input_Staph )
rownames(input_Staph_data) <-input_Staph_data[,1]

input_meta_file <-as.data.frame(CHG_Study.meta)
input_meta_file <-subset.data.frame(input_meta_file, Timepoint != "NegEx")
input_meta_file$Timepoint[input_meta_file$Timepoint == "Pre2"] <- "Pre"
input_meta_file$Timepoint[input_meta_file$Timepoint == "OR2"] <- "OR"
input_meta_file$Timepoint[input_meta_file$Timepoint == "Post2"] <- "Post"

PMAxx_meta_file <- subset.data.frame(input_meta_file, PMAxx == "PMAxx")
SSpma_MF <- subset.data.frame(PMAxx_meta_file, Site =="SurgeSite")
SSpma_MF$Timepoint <- factor(SSpma_MF$Timepoint, levels = c("Base", "Pre","OR", "Post", "CV1"))

Cpma_MF <- subset.data.frame(PMAxx_meta_file, Site =="Control")
Cpma_MF$Timepoint <- factor(Cpma_MF$Timepoint, levels = c("Base", "Pre","OR", "Post", "CV1"))

fit_data = Maaslin2(
  input_data = input_Staph_data, 
  input_metadata = SSpma_MF, 
  output = "./MAASLIN/SurgViable/Staph_time_MAASLIN", 
  min_abundance = 0.01,
  fixed_effects = c("Timepoint"),
  reference = c("Timepoint", "Base"),
  random_effects = c("Subject","BodySite", "Gender", "AntibioticProphylaxisYN"))

fit_data = Maaslin2(
  input_data = input_Staph_data, 
  input_metadata = Cpma_MF, 
  output = "./MAASLIN/ContViable/Staph_time_MAASLIN", 
  min_abundance = 0.01,
  fixed_effects = c("Timepoint"),
  reference = c("Timepoint", "Base"),
  random_effects = c("Subject","BodySite", "Gender", "AntibioticProphylaxisYN"))


# Supplemental Figure 10 - CARD on healthy subject skin metagenomes - prevelence of antiseptic genes 
trim <- function (x) gsub("^\\s+|\\s+$", "", x) # remove leading and trailing white spaces from cells 

meta.df<- as.data.frame(read_csv("./LKMB002_Metadata_UPDATED.csv"))

CardResults.df <- as.data.frame(read_csv("./Concatenated_Results_KMA.csv"))
CardResults.df3 <- CardResults.df %>% separate_rows(DrugClass, sep = ";")
CardResults.df3$DrugClass <- trim(CardResults.df3$DrugClass)  
CardResults.df3 <- subset.data.frame(CardResults.df3, DrugClass != "glycopeptide antibiotic")
CardResults.df3<- CardResults.df3 %>% mutate(AMRclass = case_when(DrugClass == "aminocoumarin antibiotic" ~ "Aminocoumarin",
                                                                  DrugClass == "aminoglycoside antibiotic" ~ "Aminoglycoside",
                                                                  DrugClass == "antibacterial free fatty acids" ~ "Antibacterial Free Fatty Acids",
                                                                  DrugClass == "bicyclomycin-like antibiotic" ~ "Bicyclomycin-like",
                                                                  DrugClass == "carbapenem" ~ "Beta-Lactam",
                                                                  DrugClass == "cephalosporin" ~ "Beta-Lactam",
                                                                  DrugClass == "cephamycin" ~ "Beta-Lactam",
                                                                  DrugClass == "diaminopyrimidine antibiotic" ~ "Diaminopyrimidine",
                                                                  DrugClass == "disinfecting agents and antiseptics" ~ "Disinfecting Agents and Antiseptics",
                                                                  DrugClass == "fluoroquinolone antibiotic" ~ "Fluoroquinolone",
                                                                  DrugClass == "glycopeptide antibiotic" ~ "Glycopeptide",
                                                                  DrugClass == "glycylcycline" ~ "Glycylcycline",
                                                                  DrugClass == "isoniazid-like antibiotic" ~ "Isoniazid-like",
                                                                  DrugClass == "lincosamide antibiotic" ~ "Lincosamide",
                                                                  DrugClass == "macrolide antibiotic" ~ "Macrolide",
                                                                  DrugClass == "monobactam" ~ "Beta-Lactam",
                                                                  DrugClass == "nitrofuran antibiotic" ~ "Nitrofuran",
                                                                  DrugClass == "nitroimidazole antibiotic" ~ "Nitroimidazole",
                                                                  DrugClass == "nucleoside antibiotic" ~ "Nucleoside Antibiotics",
                                                                  DrugClass == "oxazolidinone antibiotic" ~ "Oxazolidinone",
                                                                  DrugClass == "penam" ~ "Beta-Lactam",
                                                                  DrugClass == "penem" ~ "Beta-Lactam",
                                                                  DrugClass == "peptide antibiotic" ~ "Peptide Antibiotics",
                                                                  DrugClass == "phenicol antibiotic" ~ "Phenicol",
                                                                  DrugClass == "phosphonic acid antibiotic" ~ "Phosphonic Acid",
                                                                  DrugClass == "pleuromutilin antibiotic" ~ "Pleuromutilin",
                                                                  DrugClass == "rifamycin antibiotic" ~ "Rifamycin",
                                                                  DrugClass == "streptogramin antibiotic" ~ "Streptogramin",
                                                                  DrugClass == "sulfonamide antibiotic" ~ "Sulfonamide",
                                                                  DrugClass == "tetracycline antibiotic" ~ "Tetracycline",))
CHG.res <- subset.data.frame(CardResults.df3, AMRclass == "Disinfecting Agents and Antiseptics")

CHG.res$Hit <- ifelse((CHG.res$PercentCoverage >=70) & (CHG.res$AverageMAPQ_CompletelyMappedReads >= 50), 1, 0) # 1 == hit, 0 = no hit
CHG.res <- subset.data.frame(CHG.res, Hit == 1) # keeping only hits with 75% coverage
CHG.res <- merge(meta.df, CHG.res, by = "SampleID", all = T)

CHG.res.FIN <-distinct(CHG.res, SampleID, AROTerm, .keep_all= TRUE) # removing any duplicate hits (ie a sample with multiple efflux pumps hits or aminoglycoside hits - so there is only one aminoglycoside count "presnet" for each sample)
CHG.res.FIN  <- CHG.res.FIN  %>% drop_na(AROTerm)
CHG.res.FIN  <-subset.data.frame(CHG.res.FIN , SkinSite != "Neg")
CHG.res.FIN  <-subset.data.frame(CHG.res.FIN , SkinSite != "Mock")
CHG.res.FIN  <-subset.data.frame(CHG.res.FIN , AROTerm != "NA")
CHG.res.FIN  <-subset.data.frame(CHG.res.FIN , SubjectID != "NA")
CHG.res.BodySite.freq <- CHG.res.FIN  %>% group_by(SkinSite) %>% 
  count(AROTerm) %>% 
  mutate(freq = n / 34 * 100) # divided by the 34 subjects
CHG.res.BodySite.freq
CHG.res.Subject.freq <- CHG.res.FIN %>% group_by(SubjectID) %>% 
  count(AROTerm) %>% 
  mutate(freq = n / 8 * 100) # divided by the 8 skin sites for each subject
CHG.res.Subject.freq 
CHG.res.Collective.freq <- CHG.res.FIN %>% group_by(AROTerm) %>% 
  count(AROTerm) %>% 
  mutate(freq = n / 272 * 100) # divided by the 8 skin sites for each subject
CHG.res.Collective.freq

AROm <- na.omit(subset.data.frame(CHG.res, select = c(AROTerm, AMRclass, Resistomes_and_Variants_Observed_Pathogens,ResistanceMechanism )))
AROm <-distinct(AROm, AROTerm, AMRclass,Resistomes_and_Variants_Observed_Pathogens, ResistanceMechanism, .keep_all= TRUE)
rownames(AROm) <-AROm$AROTerm
ARO.matrix<- AROm[-1]

BodySite.Matrix <- acast(CHG.res.BodySite.freq, AROTerm ~ SkinSite, value.var="freq")
BodySite.Matrix[is.na(BodySite.Matrix)] <- 0

pheatmap(BodySite.Matrix, 
         color = c("white","#FFF38F", "#FDEA45" ,"#91B959" , "#248A8D","#344571","#440054"),
         breaks = c(0,1,10, 20, 40,60, 80, 100),
         annotation_row =  ARO.matrix,
         #annotation_colors = ExtColor,
         border_color = F,
         fontsize_row = 5, fontsize_col = 10)

sunsetshades <- c("#C1C1C1", #Gray
                  colorRampPalette(c("#F9E39F", "#E9B10C","#EC9220", "#EF7234","#D55D45",
                                     "#F5BAA3","#D77870","#BB4855","#843369","#62386A",
                                     "#DACCE1","#796796","#483470","#2B2358","#15112C"))(76))  # purples          

pheatmap(BodySite.Matrix, 
         color = c("white","#F8DD8C", "#E9B10C" ,"#EF7234" , "#BB4855","#843369","#261F4D"),
         breaks = c(0,1,10, 20, 40,60, 80, 100),
         annotation_row =  ARO.matrix,
         #annotation_colors = ExtColor,
         border_color = F,
         fontsize_row = 5, fontsize_col = 10)

Subject.Matrix <- acast(CHG.res.Subject.freq , AROTerm ~ SubjectID, value.var="freq")
Subject.Matrix[is.na(Subject.Matrix)] <- 0

pheatmap(Subject.Matrix, 
         color = c("white","#FFF38F", "#FDEA45" ,"#91B959" , "#248A8D","#344571","#440054"),
         breaks = c(0,1,10, 20, 40,60, 80, 100),
         border_color = F,
         annotation_row =  ARO.matrix,
         #annotation_colors = ExtColor,
         fontsize_row = 5, fontsize_col = 5)
