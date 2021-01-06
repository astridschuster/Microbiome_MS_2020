##########################
#Microbiome data analysis for Schuster and Strehlow et al., 2020 mSphare 
##########################

library(dplyr)
library(gridExtra)
library(ggnet)


# Check ampvis2 version
if(packageVersion("ampvis2") < "2.4.0") {
  stop(paste("The accepted ampvis2 version is 2.4.0. Your version is ", packageVersion("ampvis2"), sep=""))
}

#tell r where to look for files

#setwd("~/Desktop/my_RScripts16S")

#working directory and specific libraries for your laptop computer

setwd("C:/Users/username/Google Drive/my_RScripts16S")
library(ampvis2, lib="/Users/strehlow/Documents/R_packages")

Sys.setenv(LANG = "en")

# Primer 1
P1_metadata <- readxl::read_excel(path="metadata.xls") %>% 
  subset(Libtype == "abV4-C")
P1_otutable <- read.delim(file = "otutable.txt", 
                          check.names = FALSE)
P1_d <- amp_load(otutable = P1_otutable,
                 metadata = P1_metadata,
                 fasta = "otus.fa")
P1_d_s <- amp_subset_samples(P1_d, !(Samplename %in% c("PC","NC", "NEG", "PC1", "PC2", "NC1", "NC2"))) #remove sequencing control samples
P1_d_s <- amp_subset_samples(P1_d_s, (Season %in% c("Summer","Spring"))) #filter for seasonal data

#as of 13/3, also includes suberities, tethya and halichondria as a reference

#Filter for only anoxic samples (9 samples and 3368 OTUs including Filter)
#P1_d_a <- amp_subset_samples(P1_d_s, (Oxy_cat %in% c("Anoxic")))

# The SILVA taxonomy we use have spaces, which makes various errors. Run this code to fix it - remove unknown etc.
P1_d_s$tax <- mutate_all(P1_d_s$tax, function(x) {
  x <- as.character(x)
  x[grepl("uncultured|unknown|incertae sedis", tolower(x))] <- ""
  x <- gsub(" +", "_", x)
  return(x)
})
rownames(P1_d_s$tax) <- P1_d_s$tax$OTU


# Export the fasta file(s) and OTU table with a nice name

# Primer 1
amp_export_fasta(data = P1_d_s, 
                 filename = paste0("[", format(Sys.Date(), "%Y-%m-%d"), "] DNASense ", projectID, " V4 OTUs.fa"))

amp_export_otutable(data = P1_d_s, sep = ",",
                    filename = paste0("[", format(Sys.Date(), "%Y-%m-%d"), "] DNASense ", projectID, " V4 OTU-table"),
                    id = "Samplename", normalise = TRUE)

# Export a corrected OTUtable for the app - remember to rename to 'otutable.txt' and upload for use as otutable in ampvis2
amp_export_otutable(data = P1_d_s, sep = "\t",
                    filename = "otutable_V4",
                    id = "SeqID", normalise = FALSE)

#note: all data for field anoxia study = P1_d_s

#updated data to only include H. stillefera, eurypon sp2, water and sediment
focusTaxa <- amp_subset_samples(P1_d_s, (Species_name %in% c("Hymeraphia_stellifera", 'Eurypon_sp2', 'Sediment', 'Water')))

#45 samples and 257 OTUs have been filtered 
#Before: 100 samples and 4775 OTUs
#After: 55 samples and 4518 OTUs

#heatmap alldata
amp_heatmap(P1_d_s, 
            group_by = "Species_name", 
            #order_x_by = "cluster",
            tax_aggregate = 'Phylum', 
            tax_show = 25,
            normalise = T,
            tax_add = "OTU", 
            tax_class = "p__Proteobacteria",
            plot_colorscale = "sqrt", 
            color_vector = c("White", "Red"),
            plot_na= T,
            plot_values = T,
            tax_empty = "best",
            min_abundance = 0.1,
            max_abundance = 20,
            plot_values_size = 3,
            facet_by = "Oxy_cat"
) +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "right")
####
#all data by class Fig. S7
####
amp_heatmap(P1_d_s, 
            group_by = "Species_name", 
            #order_x_by = "cluster",
            tax_aggregate = c("Class"), 
            tax_show = 30,
            normalise = T,
            #tax_add = "OTU", 
            #tax_class = "p__Proteobacteria",
            plot_colorscale = "sqrt", 
            color_vector = c("White", "Red"),
            plot_na= T,
            plot_values = T,
            tax_empty = "best",
            min_abundance = 0.1,
            max_abundance = 20,
            plot_values_size = 3,
            facet_by = "Oxy_cat"
) +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "right")
#size: 1300 X 700

######
#used for figure7 - gather top otus from all species
######

#From looking at heatmaps of individual sponge species, not all high abundance OTUS were
#represented just from the top 25 most abundant overall, according to ampvis
All_top <- amp_subset_taxa(P1_d_s, tax_vector = c('OTU_1', 'OTU_2', 'OTU_3', 'OTU_7',
                                                  'OTU_6', 'OTU_1075', 'OTU_17', 'OTU_10', 'OTU_31', 'OTU_39', 'OTU_4', 'OTU_9', 'OTU_19',
                                                  'OTU_5', 'OTU_21', 'OTU_145', 'OTU_8', 'OTU_28', 'OTU_56', 'OTU_53','OTU_389', 'OTU_27',
                                                  'OTU_47', 'OTU_66', 'OTU_271', 'OTU_113'), normalise = TRUE)

#Figure7-  simplified, with OTUs in front
amp_heatmap(All_top, 
            group_by = "Species_name", 
            #order_x_by = "cluster",
            #tax_aggregate = c("Genus"), 
            tax_show = 26,
            normalise = F,
            tax_add = "OTU", 
            tax_class = "p__Proteobacteria",
            plot_colorscale = "sqrt", 
            color_vector = c("White", "Red"),
            plot_na= T,
            plot_values = T,
            tax_empty = "best",
            min_abundance = 0.1,
            max_abundance = 20,
            plot_values_size = 4,
            facet_by = "Oxy_cat"
) +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 11),
        legend.position = "right")
#size:1300x700

#extra info for OTUs from Figure 7 to add to final figure post review
amp_heatmap(All_top, 
            group_by = "Species_name", 
            #order_x_by = "cluster",
            #tax_aggregate = c("Genus"), 
            tax_show = 26,
            normalise = F,
            tax_add =  "Order", "Genus",
            tax_class = "p__Proteobacteria",
            plot_colorscale = "sqrt", 
            color_vector = c("White", "Red"),
            plot_na= T,
            plot_values = T,
            tax_empty = "best",
            min_abundance = 0.1,
            max_abundance = 20,
            plot_values_size = 4,
            facet_by = "Oxy_cat"
) +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 11),
        legend.position = "right")
#size:1300x700
#size:1300x700

######
#shannon index- focus taxa 
######
stats <- focusTaxa %>%
  amp_alphadiv(measure = c("observed", "shannon"), rarefy = 16000) %>%
  select(SeqID, 
         Samplename, 
         ConcExt, 
         ConcLib, 
         RawReads, 
         ObservedOTUs, 
         Shannon) %>%
  mutate(ConcLib = round(ConcLib, 1),
         ConcExt = round(ConcExt, 1),
         Shannon = round(Shannon, 1)) %>%
  mutate(ConcLib = replace(ConcLib, ConcLib < 0.2, "< 0.2"),
         ConcExt = replace(ConcExt, ConcExt < 2, "< 2")) %>%
  arrange(as.numeric(substring(SeqID,10)))


#table of information/ diversity statistics 
colnames(stats) <- c("Sequencing\nID", 
                     "Sample_nName", 
                     "Extraction\nConc. [ng/uL]", 
                     "Library\nConc. [ng/uL]", 
                     "Reads", 
                     "Observed_nOTUs", 
                     "Shannon_nIndex"
)


tt2 <- ttheme_default(core=list(fg_params=list(hjust=1, 
                                               x = 0.95, 
                                               fontsize = 6)),
                      colhead=list(fg_params=list(fontsize = 2)))

grid.table(stats[01:89,], rows= NULL, theme = tt2)

#only ten rows to fit screen better: grid.table(stats[01:10,], rows= NULL, theme = tt2)


#boxplot of shannon index 
ggplot(stats, aes(y=Shannon_nIndex)) + geom_boxplot()

#shannon index  ~3 

ggplot(stats, aes(y=Observed_nOTUs)) + geom_boxplot()

#shannon index  ~600 

######
#Shannon index - all data
######

#all data 
stats <- P1_d_s %>%
  amp_alphadiv(measure = c("observed", "shannon"), rarefy = 16000) %>%
  select(SeqID, 
         Samplename, 
         ConcExt, 
         ConcLib, 
         RawReads, 
         ObservedOTUs, 
         Shannon) %>%
  mutate(ConcLib = round(ConcLib, 1),
         ConcExt = round(ConcExt, 1),
         Shannon = round(Shannon, 1)) %>%
  mutate(ConcLib = replace(ConcLib, ConcLib < 0.2, "< 0.2"),
         ConcExt = replace(ConcExt, ConcExt < 2, "< 2")) %>%
  arrange(as.numeric(substring(SeqID,10)))

#The following sample(s) have not been rarefied (less than 16000 reads): MQ191007-263

#table of information/ diversity statistics 
colnames(stats) <- c("Sequencing\nID", 
                     "Sample_nName", 
                     "Extraction\nConc. [ng/uL]", 
                     "Library\nConc. [ng/uL]", 
                     "Reads", 
                     "Observed_nOTUs", 
                     "Shannon_nIndex"
)


tt2 <- ttheme_default(core=list(fg_params=list(hjust=1, 
                                               x = 0.95, 
                                               fontsize = 6)),
                      colhead=list(fg_params=list(fontsize = 2)))

grid.table(stats[01:89,], rows= NULL, theme = tt2)

#only ten rows to fit screen better: grid.table(stats[01:10,], rows= NULL, theme = tt2)


#boxplot of shannon index 
ggplot(stats, aes(y=Shannon_nIndex)) + geom_boxplot()

#shannon index still approx 3 

ggplot(stats, aes(y=Observed_nOTUs)) + geom_boxplot()

#observed OTUs approx 550 per sample 

#####
#Figure3:
#####
#Figure 3A PCA
amp_ordinate(data = focusTaxa, 
             type = "PCA", 
             transform = "hellinger", 
             #sample_label_by = "Samplename",
             sample_shape_by = "Oxy_cat",
             sample_color_by = "Species_name",
             #species_plot = TRUE,
             #species_nlabels = 5,
             #filter_species = 2.0,
             sample_colorframe = "TRUE", #or TRUE for the same as sample_color_by
             #sample_colorframe_label = "Groups",
             sample_point_size = 2.5,
             sample_label_size = 1) +
  theme_classic() +
  theme(legend.position = "right") 
#scale_color_discrete(name = "")

#export size: 948X656

#Figure 3B CCA based on oxygen condition
result <- amp_ordinate(focusTaxa, 
                       type = "CCA", 
                       constrain = "Oxy_cat",
                       transform = "hellinger", 
                       #sample_label_by = "Season",
                       sample_shape_by = "Species_name",
                       sample_color_by = "Oxy_cat",
                       species_plot = TRUE,
                       #species_nlabels = 5,
                       sample_point_size = 4,
                       sample_colorframe = "TRUE", #or TRUE for the same as sample_color_by
                       sample_colorframe_label = "Oxy_cat",
                       #sample_plotly = 'all',
                       #sample_colorframe = "AnotherGroup" #or TRUE for the same as sample_color_by
                       sample_label_size = 3) +
  theme_classic() +
  theme(legend.position = "right") 
#scale_color_discrete(name = "")

plot(result)

#export size: 948X656

#Figure 3C heatmap
#Drivers = OTUs within the red shape in CCA, i.e. driving separation
#indentifed by hand using DNAsense shiny app https://dnasense.shinyapps.io/dnasense/
drivers <- c('OTU_951', 'OTU_120', 'OTU_484', 'OTU_838', 'OTU_2072', 'OTU_125', 'OTU_209', 'OTU_124', 'OTU_176', 'OTU_345', 'OTU_422', 'OTU_68', 'OTU_927', 'OTU_1150', 'OTU_587', 'OTU_189', 'OTU_1042', 'OTU_397', 'OTU_394', 'OTU_337', 'OTU_205', 'OTU_810', 'OTU_926', 'OTU_676', 'OTU_337', 'OTU_541', 'OTU_1127', 'OTU_144', 'OTU_885', 'OTU_148', 'OTU_247', 'OTU_2427', 'OTU_4541', 'OTU_3760', 'OTU_682', 'OTU_921', 'OTU_516', 'OTU_718', 'OTU_211', 'OTU_732', 'OTU_3811', 'OTU_1255', 'OTU_1072', 'OTU_732', 'OTU_3811', 'OTU_1255', 'OTU_697', 'OTU_1394', 'OTU_1842', 'OTU_790', 'OTU_488', 'OTU_699', 'OTU_1794', 'OTU_872', 'OTU_1132', 'OTU_260', 'OTU_1672', 'OTU_547', 'OTU_408', 'OTU_3710', 'OTU_4180', 'OTU_1352', 'OTU_208', 'OTU_1205', 'OTU_424', 'OTU_584', 'OTU_1455', 'OTU_1460', 'OTU_847', 'OTU_813', 'OTU_311', 'OTU_1122', 'OTU_1216', 'OTU_296', 'OTU_2293', 'OTU_4180', 'OTU_395', 'OTU_531', 'OTU_517', 'OTU_530', 'OTU_442', 'OTU_359', 'OTU_1833', 'OTU_881', 'OTU_1060', 'OTU_1220', 'OTU_607', 'OTU_767', 'OTU_2675', 'OTU_606', 'OTU_298', 'OTU_1559', 'OTU_717', 'OTU_858', 'OTU_2652', 'OTU_968', 'OTU_4497', 'OTU_206', 'OTU_1388', 'OTU_819', 'OTU_736', 'OTU_692', 'OTU_1665', 'OTU_1360', 'OTU_2607', 'OTU_1094', 'OTU_505', 'OTU_242', 'OTU_201', 'OTU_2299', 'OTU_938', 'OTU_1155', 'OTU_1756', 'OTU_944', 'OTU_1133', 'OTU_996', 'OTU_9362', 'OTU_1605', 'OTU_994', 'OTU_379', 'OTU_1492', 'OTU_623', 'OTU_521', 'OTU_942', 'OTU_214', 'OTU_694', 'OTU_545', 'OTU_1509', 'OTU_778', 'OTU_755', 'OTU_551',  'OTU_4796', 'OTU_535', 'OTU_647', 'OTU_367', 'OTU_199', 'OTU_645', 'OTU_326', 'OTU_635', 'OTU_115', 'OTU_122', 'OTU_1422', 'OTU_818', 'OTU_297', 'OTU_1896', 'OTU_261', 'OTU_174', 'OTU_248', 'OTU_1054', 'OTU_100', 'OTU_1474', 'OTU_483', 'OTU_413', 'OTU_142', 'OTU_283', 'OTU_323', 'OTU_910', 'OTU_401', 'OTU_640', 'OTU_2428', 'OTU_1105', 'OTU_770', 'OTU_640', 'OTU_910', 'OTU_191', 'OTU_30', 'OTU_105', 'OTU_432', 'OTU_186', 'OTU_230', 'OTU_1779', 'OTU_1571', 'OTU_168', 'OTU_12', 'OTU_2672', 'OTU_207', 'OTU_2775', 'OTU_69', 'OTU_827', 'OTU_324', 'OTU_171', 'OTU_63', 'OTU_4101', 'OTU_2927', 'OTU_55', 'OTU_1717', 'OTU_67')

#removed any possible duplicates 
udrivers <- unique(drivers)

select_otus <- amp_subset_taxa(focusTaxa, tax_vector = udrivers, normalise = TRUE)

amp_heatmap(select_otus, 
            group_by = "Species_name",
            order_x_by = "cluster",
            tax_aggregate = c("Order"), 
            tax_show = 25,
            tax_add = "OTU", 
            normalise = FALSE,
            tax_class = "p__Proteobacteria",
            plot_colorscale = "sqrt", 
            color_vector = c("White", "Red"),
            plot_na= T,
            plot_values = T,
            tax_empty = "best",
            min_abundance = 0.1,
            max_abundance = 20,
            plot_values_size = 3,
            facet_by = "Oxy_cat"
) +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "right")

#size:1000, 700

#rerun to get additional taxanomic information post review 
amp_heatmap(select_otus, 
            group_by = "Species_name",
            #order_x_by = "cluster",
            #tax_aggregate = c("Order"), 
            tax_show = 25,
            tax_add =  "OTU", "Phylum", 
            normalise = FALSE,
            tax_class = "p__Proteobacteria",
            plot_colorscale = "sqrt", 
            color_vector = c("White", "Red"),
            plot_na= T,
            plot_values = T,
            tax_empty = "best",
            min_abundance = 0.1,
            max_abundance = 20,
            plot_values_size = 3,
            facet_by = "Oxy_cat"
) +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "right")

######
#Fig s5- CCA all data (no outgroup)
######
no_outgroup <- amp_subset_samples(P1_d_s, !(Oxy_cat %in% c('Outgroup')))

result <- amp_ordinate(no_outgroup, 
                       type = "CCA", 
                       constrain = "Oxy_cat",
                       transform = "hellinger", 
                       sample_label_by = "Samplename",
                       #sample_shape_by = "Species_name",
                       sample_color_by = "Oxy_cat",
                       species_plot = TRUE,
                       #species_nlabels = 5,
                       sample_point_size = 4,
                       sample_colorframe = "TRUE", #or TRUE for the same as sample_color_by
                       sample_colorframe_label = "Oxy_cat",
                       #sample_plotly = 'all',
                       #sample_colorframe = "AnotherGroup" #or TRUE for the same as sample_color_by
                       sample_label_size = 3) +
  theme_classic() +
  theme(legend.position = "right") 
#scale_color_discrete(name = "")

plot(result)

#export size: 948X656

#######
#PERMANOVA (table 1)
#######

otus = t(focusTaxa$abund)
meta = focusTaxa$metadata
taxtab = focusTaxa$tax

library(vegan)

#two factor permanova, with sample type and oxygen condition
permanova = adonis(otus~Oxy_cat*Species_name, data = meta, distance = "bray", permutations = 10000)

permanova$aov.tab

####
#Figure4 -top 6 most abundant OTUs in focus taxa - determined through explorations of subset, e.g. heatmaps above
####

top7 <- c('OTU_1', 'OTU_2', 'OTU_3', 'OTU_7', 'OTU_1075', 'OTU_17', 'OTU_6')

key_otus <- amp_subset_taxa(focusTaxa, tax_vector =top7, normalise = T)

amp_heatmap(key_otus, 
            group_by = "Species_name",
            order_x_by = "cluster",
            #tax_aggregate = c("Genus"), 
            tax_show = 9,
            tax_add = "OTU",
            normalise = F,
            tax_class = "p__Proteobacteria",
            plot_colorscale = "sqrt", 
            color_vector = c("White", "Red"),
            plot_na= T,
            plot_values = T,
            tax_empty = "best",
            min_abundance = 0.1,
            max_abundance = 20,
            plot_values_size = 3,
            #plot_legendbreaks = c(1, 5, 10),
            facet_by = "Oxy_cat"
) +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "right")

######################
#Comparative statistics of top 7 OTUS (individually) across oxygen conditions. Table S1.
#######################

######
#OTU1 statistical analysis
######
#all data for OTU1
#Eurypon
logitTransform <- function(p) { log(p/(1-p)) }
asinTransform <- function(p) { asin(sqrt(p)) }

OTU1 <- amp_subset_taxa(focusTaxa, tax_vector = c('OTU_1'), normalise = TRUE)

otus = t(OTU1$abund)
meta = OTU1$metadata
taxtab = OTU1$tax

#test for normalily
shapiro.test(otus)

#subsample only eurpon for normality
OTU1_norm <- amp_subset_samples(OTU1, (Species_name %in% c('Eurypon_sp2')))


abund = data.frame(t(OTU1_norm$abund))
meta = OTU1_norm$metadata
taxtab = OTU1_norm$tax

meta$abund <- abund$OTU_1

hist(abund$OTU_1)

#still not normal, so used Kruskal-wallis test for OTU1, eurypon

kruskal.test(meta$abund~meta$Oxy_cat)

#OTU1 - hymer - is abundance different from zero in anoxia?

#hymeraphia stillefera
OTU1_hymer <-amp_subset_samples(OTU1, (Species_name %in% c('Hymeraphia_stellifera')))


#OTU 1 in hymeraphia should still be significant?y different from zero though - try t test on just the anoxic data
OTU1_a <-amp_subset_samples(OTU1_hymer, (Oxy_cat %in% c('Anoxic')))
#two samples, one with a relative abundance of 7.64 and one with 0.0122

t.test(c(7.64, 0.0122), mu=0)

#######
#OTU2
#######

OTU2 <- amp_subset_taxa(focusTaxa, tax_vector = c('OTU_2'), normalise = TRUE)

#only eurypon
OTU2_norm <- amp_subset_samples(OTU2, (Species_name %in% c('Eurypon_sp2')))

abund = data.frame(t(OTU2_norm$abund))
meta = OTU2_norm$metadata
taxtab = OTU2_norm$tax

meta$abund <- abund$OTU_2

hist(abund$OTU_2)
#not normal

kruskal.test(meta$abund~meta$Oxy_cat)

#OTU2 - H. stellifera

OTU2_hymer <-amp_subset_samples(OTU2, (Species_name %in% c('Hymeraphia_stellifera')))

OTU2_hymer_a <-amp_subset_samples(OTU2_hymer, (Oxy_cat %in% c('Anoxic')))
#two samples with abundances of 5.75 and 0.0109

t.test(c(5.75, 0.0109), mu=0)

####
#OTU3
####

OTU3 <- amp_subset_taxa(focusTaxa, tax_vector = c('OTU_3'), normalise = TRUE)
OTU3_eury <-amp_subset_samples(OTU3, (Species_name %in% c('Eurypon_sp2')))
OTU3_hymer <-amp_subset_samples(OTU3, (Species_name %in% c('Hymeraphia_stellifera')))
OTU3_water <-amp_subset_samples(OTU3, (Species_name %in% c('Water')))
OTU3_sed <-amp_subset_samples(OTU3, (Species_name %in% c('Sediment')))

#split by species, check if each species has abundaces greater than zeros
#looks like Hs is the most important

otus = data.frame(t(OTU3_eury$abund))

t.test(otus, mu=0)

otus = data.frame(t(OTU3_water$abund))

t.test(otus, mu=0)


#OTU3 comparisons

otus = data.frame(t(OTU3$abund))
meta = OTU3$metadata
taxtab = OTU3$tax

meta$abund <- otus$OTU_3

kruskal.test(meta$abund~meta$Species_name)

pairwise.wilcox.test(meta$abund, meta$Species_name, p.adjust.method = "BH")


#######
#OTU7
#######
OTU7 <- amp_subset_taxa(focusTaxa, tax_vector = c('OTU_7'), normalise = TRUE)

OTU7_eury <-amp_subset_samples(OTU7, (Species_name %in% c('Eurypon_sp2')))

otus = data.frame(t(OTU7_eury$abund))
meta = OTU7_eury$metadata
taxtab = OTU7_eury$tax

meta$abund <- otus$OTU_7

#not normal for eurypon

kruskal.test(meta$abund~meta$Oxy_cat)

#OTU7 H- Stellifera

OTU7_hymer <-amp_subset_samples(OTU7, (Species_name %in% c('Hymeraphia_stellifera')))

otus = data.frame(t(OTU7_hymer$abund))
meta = OTU7_hymer$metadata
taxtab = OTU7_hymer$tax

shapiro.test(otus$OTU_7)
#normal :)

hist(otus$OTU_7)

meta$abund <- otus$OTU_7

anova = aov(meta$abund~meta$Oxy_cat)

summary(anova)

TukeyHSD(anova)

####
#OTU6
####

OTU6 <- amp_subset_taxa(focusTaxa, tax_vector = c('OTU_6'), normalise = TRUE)
#zero in all environmental samples, nt normal otherwise

otus = t(OTU6$abund)
meta = OTU6$metadata
taxtab = OTU6 $tax

OTU6_e <-amp_subset_samples(OTU6, (Species_name %in% c('Eurypon_sp2')))

#t-test done for eurypon as it is only in anoxia

otus = data.frame(t(OTU6_e$abund))
meta = data.frame(OTU6_e$metadata)
taxtab = OTU6_e$tax

#normal

meta$abund <- otus$OTU_6

anova = aov(meta$abund~meta$Oxy_cat)

summary(anova)

TukeyHSD(anova)

#t-test for Hymeraphia

OTU6_eury <-amp_subset_samples(OTU6, (Species_name %in% c('Hymeraphia_stellifera')))

OTU6_eury_a <-amp_subset_samples(OTU6_eury, (Oxy_cat %in% c('Anoxic')))

t.test(c(0.549,0.002), mu=0)

####
#OTU1075
####
OTU1075 <- amp_subset_taxa(focusTaxa, tax_vector = c('OTU_1075'), normalise = TRUE)

OTU1075_e <-amp_subset_samples(OTU1075, (Species_name %in% c('Eurypon_sp2')))

otus = data.frame(t(OTU1075_e$abund))
meta = data.frame(OTU1075_e$metadata)
taxtab = OTU1075_e$tax

meta$abund <- otus$OTU_1075

hist(meta$abund)
shapiro.test(meta$abund)
#yes, normal

anova = aov(abund~Oxy_cat, data = meta)

summary(anova)


TukeyHSD(anova)


OTU1075_h <-amp_subset_samples(OTU1075, (Species_name %in% c('Hymeraphia_stellifera')))

OTU1075_ha <-amp_subset_samples(OTU1075_h, (Oxy_cat %in% c('Anoxic')))

#t.test for hymeraphia 

t.test(c(1.89,0), mu=0)

####
#OTU17
####
OTU17 <- amp_subset_taxa(focusTaxa, tax_vector = c('OTU_17'), normalise = TRUE)

OTU17_s <-amp_subset_samples(OTU17, (Species_name %in% c('Hymeraphia_stellifera', 'Eurypon_sp2')))

otus = data.frame(t(OTU17_s$abund))
meta = data.frame(OTU17_s$metadata)
taxtab = OTU17_s$tax

meta$abund <- otus$OTU_17

hist(meta$abund)
shapiro.test(meta$abund)

#not nornmal with both sponge species 

kruskal.test(meta$abund~meta$Species_name)

#	Kruskal-Wallis rank sum test

#data:  meta$abund by meta$Species_name
#Kruskal-Wallis chi-squared = 19.705, df = 1, p-value = 0.000009035

wilcox.test(meta$abund~meta$Species_name)

#separated both species

#not normal

OTU17_s <-amp_subset_samples(OTU17, (Species_name %in% c('Hymeraphia_stellifera')))

otus = data.frame(t(OTU17_s$abund))
meta = data.frame(OTU17_s$metadata)
taxtab = OTU17_s$tax

meta$abund <- otus$OTU_17

hist(meta$abund)
shapiro.test(meta$abund)

#h. stellifera normal

anova = aov(abund~Oxy_cat, data = meta)
summary(anova)

#no significant differences for h. stellifera OTU17

OTU17_s <-amp_subset_samples(OTU17, (Species_name %in% c('Eurypon_sp2')))

otus = data.frame(t(OTU17_s$abund))
meta = data.frame(OTU17_s$metadata)
taxtab = OTU17_s$tax

meta$abund <- otus$OTU_17

hist(meta$abund)
shapiro.test(meta$abund)

#eurypon = not normal
kruskal.test(meta$abund~meta$Oxy_cat)

#######################################
#Physical data analysis
#######################################

#CTD data anaylsis 

Sys.setenv(LANG = "en")

library(dplyr)
library(tidyverse)
library(magrittr)

CTD <- read.csv('CTD_LH_ALL_DATA.csv', sep=';')

head(CTD)

ggplot(CTD, aes(x = Depth, y = Oxygen_uM)) +
  geom_point()

#sediment on probe ruined labhra upcast data on 27-7-2019
#filter out these data 

CTDF <-filter(CTD, !(Date %in% c("27-07-2019") & Cast %in% c('Up')))

#plot of filtered data for oxygen concentration
#summer 2019, hypoxic conditions at 27m
O2norm_plot <-ggplot(CTDF, aes(x = Oxygen_uM, y = Depth, colour=Location)) +
  geom_point()+
  #geom_smooth()+
  scale_y_reverse()+
  geom_hline(yintercept = 27)+
  xlab("Oxygen (uM)") + ylab("Depth (m)")
#horizontal line is at 27m, where we collected the sponges

O2norm_plot

#plotting the same for trace sensor
O2trace_plot <- ggplot(CTDF, aes(x = Oxygen_nM, y = Depth, colour=Location)) +
  geom_point()+
  #geom_smooth()+
  scale_y_reverse()+
  ylim(45, 30)+
  xlim(0, 200)+
  #scale_y_reverse()+
  xlab("Oxygen concentration (nM)") + ylab("Depth (m)")

O2trace_plot

#geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs") for both plots

#true anoxia as low as 35 m during the summer of 2019, but there was also lost of mixing, so daily changes were present and sulphide probably could not accumulate
ggplot(CTDF, aes(x = Cond, y = Depth, colour=Location)) +
  geom_point()+
  #geom_smooth()+
  scale_y_reverse()+
  geom_hline(yintercept = 27)+
  xlab("Conductivity") + ylab("Depth (m)")

#Conductivity plot
#one west cliff sample had a lot less conductivity, exluded for now. Might have had sediment on it
Cond_plot <- ggplot(CTDF, aes(x = Cond, y = Depth, colour=Location)) +
  geom_point()+
  geom_smooth()+
  scale_y_reverse()+
  xlim(39, 45)+
  geom_hline(yintercept = 27)+
  xlab("Conductivity") + ylab("Depth (m)")

Cond_plot

#Chlorophyll plot
#upcast data bad for chlorophyll at west cliff on 25-07-2019, so filtered it out
CTDF2 <- filter(CTDF, !(Date %in% c("25-07-2019") & Cast %in% c('Up')))
CTDF2 <- filter(CTDF2, !(Location %in% c("Labhra") & Cast %in% c('Up')))

#there are strange dots that seem out of order in the plot, but it may be an artifact of 
#the point of inflection in the cast

Chla_plot <- ggplot(CTDF2, aes(x = Chl_A_cal, y = Depth, colour=Location)) +
  #geom_smooth()+
  geom_point()+
  #geom_line()+
  scale_y_reverse()+
  #xlim(39, 45)+
  geom_hline(yintercept = 27)+
  xlab("Chlorophyl A") + ylab("Depth (m)")

Chla_plot

#definite chloropyll maximum at 10m, not sure what the strange deep maximum is though... '

#PAR (light plot)
PAR_plot <- ggplot(CTDF2, aes(x = PAR, y = Depth, colour=Location)) +
  #geom_smooth()+
  geom_point()+
  #geom_line()+
  scale_y_reverse()+
  #xlim(39, 45)+
  geom_hline(yintercept = 27)+
  xlab("PAR") + ylab("Depth (m)")

PAR_plot

#Turbidity plot
Turb_plot <- ggplot(CTDF, aes(x = Turb, y = Depth, colour=Location)) +
  #geom_smooth()+
  geom_point()+
  #geom_line()+
  scale_y_reverse()+
  xlim(0, 25)+
  geom_hline(yintercept = 27)+
  xlab("Turbidity") + ylab("Depth (m)")

Turb_plot

#turbidity relatively uniform, purhaps it increases around the thermocline though

###Temperature plot
Temp_plot <- ggplot(CTDF, aes(x = Temp, y = Depth, colour=Location)) +
  #geom_smooth()+
  geom_point()+
  #geom_line()+
  scale_y_reverse()+
  #xlim(0, 25)+
  geom_hline(yintercept = 27)+
  xlab("Temperature (C)") + ylab("Depth (m)")

Temp_plot

#Yes, turbidity has a local peak right around the thermocline :)

#Salinity
Salin_plot <- ggplot(CTDF, aes(x = SALIN, y = Depth, colour=Location)) +
  #geom_smooth()+
  geom_point()+
  #geom_line()+
  scale_y_reverse()+
  xlim(34.4, 35)+
  geom_hline(yintercept = 27)+
  xlab("Salinity (PSU)") + ylab("Depth (m)")

Salin_plot

#salinity was relatively constant, increaded by aproximately 2 units

#combined all desired plots into one figure using the ggarrange 

#7 graphs of interest 

#install.packages('ggpubr')
library(ggpubr)

O2trace_plot <- O2trace_plot + xlab("Trace oxygen (nM)")
Turb_plot <-Turb_plot + xlab("Turbidity (FTU)")

Figure_CTD <- ggarrange(O2norm_plot, O2trace_plot, Temp_plot, Turb_plot,
                        Salin_plot,
                        labels = c('A', 'B', 'C', 'D', 'E'),
                        common.legend = TRUE,
                        ncol = 5, nrow = 1, 
                        legend = 'bottom')
###FigureS1
Figure_CTD

#maybe move chlorophyl A and PAR into graph with Lisa's data?
Pig <- read.csv('Lisa_pigments.csv', sep=';')


Pig_plot <- ggplot(Pig, aes(x = Concentration, y = Depth, colour=Location, shape=Pigment)) +
  #geom_smooth()+
  geom_point(size=3)+
  #geom_point(shape=17, size=3, aes(x= Pheophytin_a, y = Depth))+
  #geom_line()+
  scale_y_reverse()+
  #xlim(34.4, 35)+
  geom_hline(yintercept = 27)+
  #scale_fill_manual(name = "Location", labels = c("West", 'Labhra'), values = c('#56B4E9', '#F8766d'))+
  #scale_shape_manual(name = "Pigment", labels = c("Chlorophyll a", 'Pheophytin_a'), values = c(1, 2))+
  xlab("Pigment concentration (ug/L)") + ylab("Depth (m)")

Pig_plot 

#pigment and light figure


#fix axes and add unit labels 
PAR_plot<-PAR_plot  + ylim(45, 0) + xlab("PAR (?E)") 
Chla_plot<-Chla_plot  + ylim(45, 0) + xlab("Chlorophyll a CTD (?g/L)") 
Pig_plot<-Pig_plot + ylim(45, 0) + xlab("Pigment concentration (?g/L)") 

Figure_Pig <- ggarrange(PAR_plot, Chla_plot, Pig_plot,
                        labels = c('A', 'B', 'C'),
                        common.legend = FALSE,
                        legend = 'right', ncol = 3, nrow = 1) 
###
#FigureS2
Figure_Pig 

#saved as image with aspect ratio 948, 656

#CTD pigment values differ by a factor of 10 from filter values

####profile data from Rob 
Prof <- read.csv('Key_oxygen_profiles.csv', sep=';')

head(Prof)

#rotated graph instead in order to get the lines to connect. So one could do that above. 

Profo2_plot <- ggplot(Prof, aes(x = ?..Depth, y = O2, colour=Condition)) +
  #geom_smooth(method = 'gam')+
  geom_point(size=3)+
  geom_line(size=2)+
  scale_x_reverse()+
  #xlim(0, 25)+
  #geom_hline(yintercept = 27)+
  coord_flip()+
  xlab("Depth (m)") + ylab("Oxygen (mg/L)")

Profo2_plot 

#temperature

Proftemp_plot <- ggplot(Prof, aes(x = ?..Depth, y = Temp, colour=Condition)) +
  #geom_smooth(method = 'gam')+
  geom_point(size=3)+
  geom_line(size=2)+
  scale_x_reverse()+
  #xlim(0, 25)+
  #geom_hline(yintercept = 27)+
  coord_flip()+
  xlab("Depth (m)") + ylab("Temperature (C)")

Proftemp_plot

#Profiles from rob figure
Figure_Prof <- ggarrange(Profo2_plot, Proftemp_plot,
                         labels = c('A', 'B'),
                         common.legend = TRUE,
                         legend = 'bottom', ncol = 2, nrow = 1) 
###for Figure1
Figure_Prof