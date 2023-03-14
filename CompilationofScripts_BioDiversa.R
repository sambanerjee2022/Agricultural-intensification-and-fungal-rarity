########### Script Compilation  #################

#####Libraries Used. 
library(tidyverse)
library(microbiome)
library(phyloseq)
library(metagMisc)
library(MiscMetabar)
library(EnhancedVolcano)
library(DESeq2)



#Part 0 Match sequencing samples names
rm(list=ls())
# OTU_pre contains irrelevant OTU data
OTU_pre <- read.csv("DATA/otu_table.csv", row.names = 1)
# change site names and delete unuseful data
Primer_match <- read.csv("TABLES/Primer_Sample_assignment.csv", na.strings = c("NA","","-"))
Primer_match <- na.omit(Primer_match)
row.names(Primer_match) <- Primer_match$Sample.ID
Primer_match$sequencingName <- paste("s",Primer_match$new.Plate,"_",Primer_match$F.Primer,".",Primer_match$R.Primer, sep = "")
Name_match <- match(colnames(OTU_pre),Primer_match$sequencingName)
OTU <- OTU_pre[,!is.na(Name_match)]
colnames(OTU) <- row.names(Primer_match) [Name_match[!is.na(Name_match)]]

# delete singletons and all-zero
OTU <- OTU[rowSums(OTU)>1,]
OTU <- as.matrix(OTU)
OTU <- OTU[,order(colnames(OTU),decreasing = F)]
OTU_pre <- OTU

## import taxonomy, and delete singleton or all-zero OTUs
otu_taxa <- read.csv("DATA/otu_taxa.csv",row.names = 1, stringsAsFactors = F)
otu_taxa <- otu_taxa[which(rownames(otu_taxa) %in% rownames(OTU)),]
otu_taxa_pre <- otu_taxa

## import design file
design <- read.csv("DATA/design.csv",row.names = 1, stringsAsFactors = F)
row.names(design) <- design$Sample_ID

## correct typos
design$Treatment[design$Treatment == "High "]<- "High"
design$Treatment[design$Treatment == "Low "]<- "Low"

design$Treatment[design$Treatment == "High" | design$Treatment == "Low" ] <- "Arable"

design[design$Treatment == "Organic"| design$Treatment == "Conventional","Treatment"] <- "Arable"
design_pre <- design
## create phyloseq object
physeq_table<- phyloseq(otu_table(OTU, taxa_are_rows=T),
                        tax_table(as.matrix(otu_taxa)),
                        sample_data(design))


########## Alluvial Plots ########################
#### AQ ####
AQ=subset_samples(physeq_table, Country=="France")
AQ
AQ<-filter_taxa(AQ, function(x) sum(x) > 0.000000000, prune=TRUE)
AQ

top20otus=names(sort(taxa_sums(AQ), TRUE)[1:10])
top20otus
taxtab20=cbind(tax_table(AQ), OTU10=NA)
taxtab20[top20otus, "OTU10"] <- as(tax_table(AQ)[top20otus, "Trophic.Mode"], "character")
tax_table(AQ)<-tax_table(taxtab20)

TEST<-as.data.frame(tax_table(merged20))
TESTmelt<-psmelt(merged20)
write.csv(TESTmelt, file="AQ_TEST.csv")

merged<-merge_samples(AQ, "Treatment")
sample_data(merged)$Treatments<-levels(sample_data(AQ)$Treatment)
merged=transform_sample_counts(merged,function(x)100*x/sum(x))

merged20 = prune_taxa(top20otus, merged)

TEST<-as.data.frame(tax_table(merged20))
TESTmelt<-psmelt(merged20)
write.csv(TESTmelt, file="FR_TEST.csv")



#### CY ####
CY=subset_samples(physeq_table, Country=="Germany")
CY
CY<-filter_taxa(CY, function(x) sum(x) > 0.000000000, prune=TRUE)
CY

top20otus=names(sort(taxa_sums(CY), TRUE)[1:10])
top20otus
taxtab20=cbind(tax_table(CY), OTU10=NA)
taxtab20[top20otus, "OTU10"] <- as(tax_table(CY)[top20otus, "Class"], "character")
tax_table(CY)<-tax_table(taxtab20)

merged<-merge_samples(CY, "Treatment")
sample_data(merged)$Treatments<-levels(sample_data(CY)$Treatment)
merged=transform_sample_counts(merged,function(x)100*x/sum(x))

merged20 = prune_taxa(top20otus, merged)

TESTmelt<-psmelt(merged20)
write.csv(TESTmelt, file="GERM_alluvial.csv")

#### CO ####
CO=subset_samples(physeq_table, Country=="Spain")
CO
CO<-filter_taxa(CO, function(x) sum(x) > 0.000000000, prune=TRUE)
CO

top20otus=names(sort(taxa_sums(CO), decreasing = TRUE)[1:10])
top20otus
taxtab20=cbind(tax_table(CO), OTU10=NA)
taxtab20[top20otus, "OTU10"] <- as(tax_table(CO)[top20otus, "Class"], "character")
tax_table(CO)<-tax_table(taxtab20)

merged<-merge_samples(CO, "Treatment")
sample_data(merged)$Treatments<-levels(sample_data(CO)$Treatment)
merged=transform_sample_counts(merged,function(x)100*x/sum(x))

merged20 = prune_taxa(top20otus, merged)

TESTmelt<-psmelt(merged20)
write.csv(TESTmelt, file="Spain_alluvial.csv")

###Switzserland
CO1=subset_samples(physeq_table, Country=="Switzerland")
CO1
CO1<-filter_taxa(CO1, function(x) sum(x) > 0.000000000, prune=TRUE)
CO1

top20otus=names(sort(taxa_sums(CO1), decreasing = TRUE)[1:10])
top20otus
taxtab20=cbind(tax_table(CO1), OTU10=NA)
taxtab20[top20otus, "OTU10"] <- as(tax_table(CO1)[top20otus, "Class"], "character")
tax_table(CO1)<-tax_table(taxtab20)

merged<-merge_samples(CO1, "Treatment")
sample_data(merged)$Treatments<-levels(sample_data(CO1)$Treatment)
merged=transform_sample_counts(merged,function(x)100*x/sum(x))

merged20 = prune_taxa(top20otus, merged)

TESTmelt<-psmelt(merged20)
write.csv(TESTmelt, file="SWZ_alluvial.csv")

###SWEDEN
CO3=subset_samples(physeq_table, Country=="Sweden")
CO3
CO3<-filter_taxa(CO1, function(x) sum(x) > 0.000000000, prune=TRUE)
CO3

top20otus=names(sort(taxa_sums(CO3), decreasing = TRUE)[1:10])
top20otus
taxtab20=cbind(tax_table(CO3), OTU10=NA)
taxtab20[top20otus, "OTU10"] <- as(tax_table(CO3)[top20otus, "Class"], "character")
tax_table(CO3)<-tax_table(taxtab20)

merged<-merge_samples(CO3, "Treatment")
sample_data(merged)$Treatments<-levels(sample_data(CO3)$Treatment)
merged=transform_sample_counts(merged,function(x)100*x/sum(x))

merged20 = prune_taxa(top20otus, merged)

TESTmelt<-psmelt(merged20)
write.csv(TESTmelt, file="SWEDEN_alluvial.csv")


####### Top country tables were later merged and input in https://www.rawgraphs.io/learning/how-to-make-an-alluvial-diagram ##








## Filter to according sites. 
Physeq_2orless <- filter_taxa(physeq_table,
                              function(x){sum(x > 0) <= 2}, ### not seen more than 0 time and in 10 samples or less. 
                              prune = TRUE)


### Getting Soil Properties. 
soilProp <- readxl::read_excel("SoilProperties_and_qPCR_Subset_dataset.xlsx",sheet = 1)
weatherProp <- readxl::read_excel("SoilProperties_and_qPCR_Subset_dataset.xlsx",sheet = 3)
### Combining Tables 
combined_soilProp <-bact.pcoa %>% select(Sample_ID,Axis.1,Axis.2) %>% left_join(.,soilProp,by="Sample_ID") %>% left_join(.,weatherProp,by="Sample_ID") 
combined_weatherProp <-bact.pcoa %>% select(Sample_ID,Axis.1,Axis.2)  
### Removing NAs 
soilProp_noNA2 <-combined_soilProp %>% select(Sample_ID,Country,Axis.1,Axis.2,Sand:avg_precip_longterm) %>% drop_na(.)

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

library(Hmisc)

soilProp_weather<- select(soilProp_noNA,c(avg_temp_longterm_1,avg_temp_longterm_2,avg_precip_longterm))

res2<-rcorr(as.matrix(soilProp_noNA),type = "spearman")

testcorr <- flattenCorrMatrix(res2$r, res2$P)

filtered_testcorr_axis1<-testcorr %>% filter(row=="Axis.1") %>% arrange(desc(cor))# %>% filter(p<0.05)

write_csv(filtered_testcorr_axis1,"OverallSamples2sites_spearmancorr_soilprop_axis.1_cor_nofilter.csv")

columnsoil <- filtered_testcorr_axis1$column[2:12]


#### Selecting positive variables for scatterplot. 

selected_values <- soilProp_noNA %>% select(all_of(columnsoil)) %>% select(-c("Calcium","WHC"))


#### Scaling of values. 
forscalinf_PCA <- selected_values %>% select(basal_respiration:AMP)

scaled_values_functions <- scale(forscalinf_PCA,
                                 center = TRUE,scale = TRUE)

chem_pca <- vegan::capscale(scaled_values_functions ~ 1, method = "euc") ### PCA plot

library(ecodist)

varespec.euc <- vegan::vegdist(scaled_values_functions, method = "euc") # dissimilarity matrix using bray-curtis distance indices on the varespec dataset native to vegan

pcoaVS <- pco(varespec.euc, #negvals = "zero",
              dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999


plot(pcoaVS$vectors[,1], pcoaVS$vectors[,2], #type = "n", 
     xlab = "PCoA1", ylab = "PCoA2",
     axes = TRUE, main = "PCoA (Euclidean distance) on Soil Properties")
text(pcoaVS$vectors[,1], pcoaVS$vectors[,2], labels(varespec.euc), 
     cex = 0.9, xpd = TRUE)
pcoaVS$values # eigenvalue for each component. This is a measure of the variance explained by each dimension
pcoaVS$vectors # eigenvectors. Each column contains the scores for that dimension.


library(ape)
library(vegan)

mite.D <- vegdist(scaled_values_functions, "euc")
res <- pcoa(mite.D)
res$values
biplot(res)
biplot(res,scaled_values_functions)

PCoA_soilAxis<- res[["vectors"]]
autoplot()

biplot(chem_pca)

vegan::scores(res, display = "sites", choices = c(1:8)) %>%
  head %>% round(2)

summary(chem_pca)$cont %>% 
  as.data.frame %>%       
  round(2)


screeplot(chem_pca, type="lines")


PC1<- vegan::scores(chem_pca, display = "sites", choices = c(1:6))


pruebita <-data.frame(Sample_ID=soilProp$Sample_ID,PCoA.1=soilProp_noNA$Axis.1,PCoA.2=soilProp_noNA$Axis.2,
                      PCoA_soilAxis)


test_long <- pruebita %>% pivot_longer(cols = Axis.1:Axis.8) 

soil_long <- data.frame(PCoA.1 =soilProp_noNA$Axis.1,PCoA.2=soilProp_noNA$Axis.2,forscalinf_PCA,weatherZscore=weatherZ$zscoreWeather,zscoreSoil=scalednewZ$zscore) %>% select(-c(C_tot_mg))%>% 
  pivot_longer(cols = basal_respiration:AMP) 

soil_long$name <- factor(soil_long$name,levels = c("MWD","basal_respiration",
                                                   "C_tot_.","C_microbial","AMP","DON","N_microbial","P_microbial"))

ggplot(soil_long, aes(value,PCoA.1)) +  
  geom_point(size=5, color="black",# alpha = .85
             shape=21,fill="chartreuse4")+
  geom_smooth(method = lm, se = FALSE, colour="black", size=1.5, formula = y ~ x)+
  scale_color_manual(values = "#FDC086")+
  facet_wrap(~name, scales = "free_x" , 
             strip.position="bottom",ncol = 4)+
  ggpubr::stat_cor(label.x = 0, method = "spearman",label.sep = "\n", size = 5,cor.coef.name = "rho",
  ) +
  theme_bw(24)+
  labs(x="Soil Properties PCoA axis",y="PCoA 1 of Rare taxa")+
  theme(
    #legend.position = c(0.90, 0.10),
    strip.placement = "outside",
    strip.background = element_blank(),
    #strip.text.x = element_blank(),
    #strip.text.y = element_blank(),
    axis.line = element_line(colour = "black"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    #axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank()
  )

summary(lm(Axis.1~DON,selected_values))

scalednewZ<- data.frame(scaled_values_functions)%>% select(-c("C_tot_mg")) %>% 
  mutate(zscore = rowSums(.))

weatherZ<- soilProp_weather  %>% scale(.,
                                       center = TRUE,scale = TRUE) %>% data.frame(.) %>%  mutate(zscoreWeather = rowSums(.))


soil_newz <- data.frame(PCoA.1 =soilProp_noNA$Axis.1,PCoA.2=soilProp_noNA$Axis.2,zscore_soil=scalednewZ$zscore,)

ggplot(soil_newz, aes(zscore_soil,PCoA.1)) +  
  geom_point(size=5, color="black",# alpha = .85
             shape=21,fill="chartreuse4")+
  geom_smooth(method = lm, se = FALSE, colour="black", size=1.5, formula = y ~ x)+
  # scale_color_manual(values = "#FDC086")+
  # facet_wrap(~name, scales = "free_x" , 
  #            strip.position="bottom",ncol = 4)+
  ggpubr::stat_cor(label.x = 0, method = "spearman",label.sep = "\n", size = 10,cor.coef.name = "rho",
  ) +
  theme_bw(24)+
  labs(x="Soil Properties Z-score",y="PCoA 1 of Rare taxa")+
  theme(
    #legend.position = c(0.90, 0.10),
    strip.placement = "outside",
    strip.background = element_blank(),
    #strip.text.x = element_blank(),
    #strip.text.y = element_blank(),
    #axis.line = element_line(colour = "black"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    #panel.border = element_blank(), 
    panel.background = element_blank(),
    #axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank()
  )

summary(lm(PCoA.1~zscore,soil_newz))

stable_C <- nlme::lme(PCoA.1 ~ basal_respiration, random= ~1 |weatherZscore,data=soil_tablemaster)
nlme::lme(PCoA.1 ~ N_microbial, random= ~1 |weatherZscore,data=soil_tablemaster)


soil_tablemaster <- data.frame(Sample_ID=soilProp_noNA2$Sample_ID,Country=soilProp_noNA2$Country,PCoA.1 =soilProp_noNA$Axis.1,PCoA.2=soilProp_noNA$Axis.2,forscalinf_PCA,weatherZscore=weatherZ$zscoreWeather,zscoreSoil=scalednewZ$zscore) %>% select(-c(C_tot_mg))
write_csv(soil_tablemaster,"lme_masterTable_weather_asZscore.csv")

####### Number of OTUS #########################
Arable2<- psadd::subset_samples_no_zero(Physeq_2orless,Treatment=="Arable")
Grassland2<- psadd::subset_samples_no_zero(Physeq_2orless,Treatment=="Grassland")

plot_bar(Physeq_2orless, x="Treatment", fill="Treatment")



ntaxa(physeq = Physeq_2orless)
ntaxa(Grassland2)


venn2sites<-ps_venn(Physeq_2orless,"Treatment",
                    fill = c("#FDC086","#7FC97F"))

df <- data.frame(Treatment=c("Arable", "Grassland"),
                 Total_OTU=c(558, 769))

Barplot2sitesorless<-ggplot(df,aes(x=Treatment,y=Total_OTU,fill=Treatment))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#FDC086","#7FC97F"))+
  theme_bw(24)+labs(x="",y="Total number of OTUs")+theme(legend.position = "none")

#### Composition
no_unidentified_order_depleted <-subset_taxa(Physeq_2orless,Order != "o__")
ps1.com.fam_depleted <- microbiomeutilities::aggregate_top_taxa2(no_unidentified_order_depleted, 
                                                                 "Order", top = 19)
#ps2_depleted<- transform_sample_counts(ps1.com.fam_depleted, function(x) (x / sum(x))*100)

plot_composition(ps1.com.fam_depleted,
                 #sample.sort = "scientific_name",
                 x.label = "env_material",
                 average_by = "Treatment",
                 otu.sort = "abundance")+
  scale_fill_tableau(palette = "Tableau 20")+ 
  theme_bw(15)+
  theme(legend.position = "right") + labs(fill="Order", y="Relative abundance (%)", x="",
                                          title = "Depleted specialist top orders",
                                          subtitle = "(Unidentified removed)")


############ Total Composition ########################
Spain10<- psadd::subset_samples_no_zero(Physeq_2orless,Country=="Spain" )
France10<- psadd::subset_samples_no_zero(Physeq_2orless,Country=="France" )
Germany10<- psadd::subset_samples_no_zero(Physeq_2orless,Country=="Germany" )
Sweden10<- psadd::subset_samples_no_zero(Physeq_2orless,Country=="Sweden" )
Switzerland10<- psadd::subset_samples_no_zero(Physeq_2orless,Country=="Switzerland" )
###### Creating absolute abundance table 
Spain10Top <- microbiomeutilities::aggregate_top_taxa2(Spain10,"Order", top = 10)
France10Top <- microbiomeutilities::aggregate_top_taxa2(France10,"Order", top = 10)
Germany10Top <- microbiomeutilities::aggregate_top_taxa2(Germany10,"Order", top = 10)
Sweden10Top <- microbiomeutilities::aggregate_top_taxa2(Sweden10,"Order", top = 10)
Switz10Top <- microbiomeutilities::aggregate_top_taxa2(Switzerland10,"Order", top = 10)

ps_venn(Spain10,"Treatment")

no_unidentified_order_depleted <-subset_taxa(Physeq_2orless,Order != "o__")
only2<- microbiomeutilities::aggregate_top_taxa2(no_unidentified_order_depleted,"Order", top = 20

SP10<- Spain10Top@otu_table %>% as.data.frame() %>% rownames_to_column(.,"OTU")
FR10 <- France10Top@otu_table %>% as.data.frame() %>% rownames_to_column(.,"OTU")
GERM10 <-Germany10Top@otu_table %>% as.data.frame() %>% rownames_to_column(.,"OTU")
Sweden10 <- Sweden10Top@otu_table %>% as.data.frame() %>% rownames_to_column(.,"OTU")
SWITZ10 <- Switz10Top@otu_table %>% as.data.frame() %>% rownames_to_column(.,"OTU")
                                                 
SP_FR<- full_join(SP10,FR10, "OTU")


SP_FR_GERM <- full_join(SP_FR,GERM10, "OTU")

SP_FR_GERM_SWEDEN <- full_join(SP_FR_GERM,Sweden10, "OTU")


SP_FR_GERM_SWEDEN_SWITZ <- full_join(SP_FR_GERM_SWEDEN,SWITZ10, "OTU")



Checking <- sjmisc::rotate_df(SP_FR_GERM_SWEDEN_SWITZ) %>% janitor::row_to_names(row_number = 1)




Checking2 <- cbind(cbind(Checking,
                         Treatment=design$Treatment,
                         Country=design$Country,
                         Sample_ID=design$Sample_ID
))

Checking2<- Checking2 %>% 
  mutate_at(vars(1:21), as.numeric)

Checking2$OtherNEW <- Checking2$o__ + Checking2$Other

Checking2[is.na(Checking2)] <- 0

Checking3<-Checking2 %>% select(.,-c(`o__`,`Other`))

Checking3<-Checking3 %>% relocate(OtherNEW, .before = Treatment)  


TopOrdersMean <- Checking3 %>%
  group_by(Treatment,Country) %>%
  summarise_at(vars(1:20), list(mean))



Checking4 <- pivot_longer(Checking3,1:20 )

pivot_longer(new_df,cols = 2:12)

# Grouped

TopOrdersMean %>% pivot_longer(3:22 ) %>% 
  ggplot(aes(fill=name, y=value, x=Treatment)) + 
  geom_bar(stat="identity")+
  #facet_wrap(~Country,nrow = 5,scales = "free_y" , 
  #strip.position="right")+
  labs(fill="Order",y="Absolute abundace")+
  theme_bw(24)+
  theme(
    #legend.position = c(0.90, 0.10),
    strip.placement = "outside",
    strip.background = element_blank(),
    #strip.text.x = element_blank(),
    #strip.text.y = element_blank(),
    axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank()
  )+scale_fill_tableau(palette = "Tableau 20")


ggplot(Checking4, aes(fill=name, y=mean, x=Treatment)) + 
  geom_bar(stat="identity")+
  #facet_wrap(~Country,nrow = 5,scales = "free_y" , 
  # strip.position="right")+
  labs(fill="Order",y="Absolute abundace")+
  theme_bw(24)+
  theme(
    #legend.position = c(0.90, 0.10),
    strip.placement = "outside",
    strip.background = element_blank(),
    #strip.text.x = element_blank(),
    #strip.text.y = element_blank(),
    axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank()
  )+scale_fill_tableau(palette = "Tableau 20")



Prueba2only<-only2@otu_table %>% as.data.frame() %>% rownames_to_column(.,"OTU") %>% 
  sjmisc::rotate_df(.) %>% janitor::row_to_names(row_number = 1)

Prueb2only2 <- cbind(cbind(Prueba2only,
                           Treatment=design$Treatment,
                           Country=design$Country,
                           Sample_ID=design$Sample_ID
))

Prueb2only2<- Prueb2only2 %>% 
  mutate_at(vars(1:21), as.numeric)

Prueb2only2$OtherNEW <- Prueb2only2$o__ + Prueb2only2$Other

Checking2[is.na(Checking2)] <- 0

Prueba3only2<-Prueb2only2 %>% select(.,-c(`o__`,`Other`)) %>% relocate(OtherNEW, .before = Treatment) 

TopOrdersMeanIverall <- Prueb2only2 %>%
  group_by(Treatment,Country) %>%
  summarise_at(vars(1:21), list(mean))


Treatment_Top<- TopOrdersMeanIverall %>% pivot_longer(3:23 ) %>% 
  ggplot(aes(fill=name, y=value, x=Treatment, color=name)) + 
  geom_bar(stat="identity")+
  #facet_wrap(~Country,nrow = 3,ncol = 2,scales = "free_y" , 
  #strip.position="right")+
  labs(fill="Order",y="Absolute abundace")+
  theme_bw(24)+
  theme(
    legend.position = "right",
    #legend.position = c(0.90, 0.15),
    strip.placement = "outside",
    strip.background = element_blank(),
    #strip.text.x = element_blank(),
    #strip.text.y = element_blank(),
    axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank()
  )+scale_fill_manual(values=c("#4E79A7","#A0CBE8","#F28E2B","#FFBE7D","#59A14F",
                                        "#8CD17D", "#B6992D","#F1CE63", "#499894", "#86BCB6",
                                        "#E15759","#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
                                        "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6","black"))+
                                          scale_color_manual(values=c("#4E79A7","#A0CBE8","#F28E2B","#FFBE7D","#59A14F","#8CD17D", 
                                                                               "#B6992D","#F1CE63", "#499894", "#86BCB6", "#E15759",
                                                                               "#FF9D9A", "#79706E", "#BAB0AC", "#D37295","#FABFD2", 
                                                                               "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6","black"))+
                                                                                 guides(color="none")

ggpubr::ggarrange(TestaPLot,Treatment_Top,
                  #Countries_Top, 
                  common.legend = TRUE, ncol = 2,legend="none",labels = c("A", "B"))





# Extract the legend. Returns a gtable
leg <- ggpubr::get_legend(Treatment_Top)

# Convert to a ggplot and print
ggpubr::as_ggplot(leg)



TestaPLot<- Barplot2sitesorless+
  theme_bw(24)+
  theme(
    legend.position = "none",
    #legend.position = c(0.90, 0.15),
    strip.placement = "outside",
    strip.background = element_blank(),
    #strip.text.x = element_blank(),
    #strip.text.y = element_blank(),
    axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank()
  )+ annotation_custom(grob = venn2sites, xmin = 0.65, xmax = 1.30, 
                       ymin = 600, ymax = 775)


PredicBarplot<-ggplot(predictors,
                      aes(x= reorder(`...1`, predictors$`% increase in mean square error`),
                          y=predictors$`% increase in mean square error`, label= predictors$`...1`)
)+
  geom_bar(fill="orange",color="black", stat = "identity",)+
  #scale_fill_manual(values = c("#FDC086","#7FC97F"))+
  geom_text(#angle = 90, 
    hjust = -.05, size = 3,  fontface="bold")+
  theme_bw(24)+labs(x="",y="")+theme(legend.position = "none")+coord_flip()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank()
  )


scatter_pred<- ggplot(predictors,
                      aes(label= predictors$`...1`,x= predictors$`% increase in mean square error`,
                          y=predictors$`Increase in node purity`
                      ))+
  geom_point(fill="orange", color="black", shape=21, size = 10,stroke = 2
  )+
  #scale_fill_manual(values = c("#FDC086","#7FC97F"))+
  theme_bw(24)+labs(x="% increase in Mean Square Error",y="Increase in node purity")+
  theme(legend.position = "none")+geom_text(hjust=0, 
                                            vjust=-2.2,
                                            nudge_x = 0.01, fontface="bold")+
  theme(#panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), #axis.line = element_line(colour = "black"
  )
annotation_custom(grob = ggplotGrob(PredicBarplot), xmin = 0, xmax =6, 
                  ymin = 0, ymax = 0.4)


ggdraw() +
  draw_plot(scatter_pred) +
  draw_plot(PredicBarplot, x = .1, y = .5, width = .3, height = .3)



##### Differential Abundance Analysis
#ASV_modified <- microbiomeutilities::format_to_besthit(ASV_clean) # Taxa to best hit

ASV_clean_plus_1 <- transform_sample_counts(Physeq_2orless, 
                                            function(x) x+1 ) # Plus One Addition


diagdds = phyloseq_to_deseq2(ASV_clean_plus_1, ~ Treatment)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")


###Comparisons

ArablevGrassland  <- results(diagdds, cooksCutoff = FALSE,
                             contrast = c('Treatment','Arable','Grassland'))
ArablevGrasslandVolcano <-EnhancedVolcano(ArablevGrassland ,
                                          lab = rownames(ArablevGrassland),
                                          #selectLab = top50_sigOE_genes,
                                          title = "OTUs present in 2 sites or less",
                                          subtitle = "Arable vs Grassland",
                                          x = 'log2FoldChange',
                                          y = 'pvalue',
                                          drawConnectors = TRUE,
                                          widthConnectors = 0.75
)


###### Niche Breadth
BacteriaNiche2<- as.data.frame(Physeq_2orless@otu_table) %>% 
  rownames_to_column(.,"OTU")

Bacteria_RareRanks2 <- as.data.frame(Physeq_2orless@tax_table) %>% 
  rownames_to_column(.,"OTU")




Prueba <- levins.Bn(BacteriaNiche2,2,
                    design$Treatment
)



TenSites_result_Levin <- levins.Bn(BacteriaNiche,2,
                                   TenSites_SampleInfo$Treatment
)



ggplot(TenSites_result_Levin,aes(Bn,))+geom_histogram(binwidth = 0.010)+ 
  geom_vline(xintercept = 0.55, linetype="dotted",color = "red", size=1.5)+
  #geom_vline(xintercept = 1.0, linetype="dotted", color = "blue", size=1.5)+
  theme_bw(14)

ggplot(TestLevin,aes(x = reorder(OTU,Bn),y=Bn))+geom_bar(stat = "identity")



result_Levin <- TenSites_result_Levin %>% 
  rownames_to_column(.,"OTU")

TestLevin <- left_join(result_Levin,Bacteria_RareRanks,"OTU")

### Getting Specialist

TenSitesSpecialist <- TestLevin %>% filter(.,Bn <= 0.55)

ggplot(TestLevin,aes(x = reorder(OTU,Bn),y=Bn))+geom_bar(stat = "identity")+
  geom_hline(yintercept = 0.55, linetype="dotted",color = "red", size=1.5)+
  theme_bw(14)

ggplot(TestLevin, aes(x = reorder(OTU,Bn),y=Bn, 
                      fill=ifelse(OTU==c("OTU95","OTU1060","OTU1814",
                                         "OTU1437","OTU170","OTU596",
                                         "OTU2368","OTU716","OTU766","OTU224"), "A","B")) + 
         geom_bar(stat = "identity",show.legend=FALSE) +
         scale_fill_manual(values=c(A="red", B="gray"))
       
       TestLevin  
       
       TestLevin %>% mutate(
         highlight = case_when(OTU == "OTU95" ~ "H",
                               OTU == "OTU1060"  ~ "H",
                               OTU== "OTU1814" ~ "H",
                               OTU=="OTU1437"~"H",
                               OTU=="OTU170"~ "H",
                               OTU=="OTU596"~ "H",
                               OTU== "OTU2368"~ "H",
                               OTU=="OTU716"~ "H",
                               OTU=="OTU766"~ "H",
                               OTU=="OTU224"~ "H"
         )) %>% ggplot(aes(x = reorder(OTU,Bn),y=Bn))+
         geom_bar(aes(fill = highlight),stat = "identity",show.legend=FALSE) +
         geom_hline(yintercept = 0.55, linetype="dotted",
                    color = "red", size=1.5)+
         theme_bw(24)+
         theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())+
         labs(x="OTUs")+ylim(0,0.55)
       
       
       
       TestLevin%>% 
         filter(Bn <=0.55)%>% 
         mutate(
           highlight = case_when(OTU == "OTU95" ~ "H",
                                 OTU == "OTU1060"  ~ "H",
                                 OTU== "OTU1814" ~ "H",
                                 OTU=="OTU1437"~"H",
                                 OTU=="OTU170"~ "H",
                                 OTU=="OTU596"~ "H",
                                 OTU== "OTU2368"~ "H",
                                 OTU=="OTU716"~ "H",
                                 OTU=="OTU766"~ "H",
                                 OTU=="OTU224"~ "H"
           ))  %>% ggplot(aes(x = reorder(OTU,Bn),y=Bn))+
         #geom_bar(aes(fill ="black"),stat = "identity",show.legend=FALSE,size=2)+
         geom_col(aes(fill=highlight),show.legend=FALSE,width = 0.5,
                  #size=30,position = "dodge"
         )+
         # geom_hline(yintercept = 0.55, linetype="dotted",
         # color = "red", size=1.5)+
         theme_bw(24)+
         theme(
           axis.ticks.y=element_blank(),
           axis.text.y = element_text(color="red",face = "bold",
                                      size=10,angle = 360))+labs(x="OTUs")+
         scale_x_discrete(breaks=c("OTU95","OTU1060","OTU1814",
                                   "OTU1437","OTU170","OTU596",
                                   "OTU2368","OTU716","OTU766","OTU224"),
                          labels=c("OTU95","OTU1060","OTU1814",
                                   "OTU1437","OTU170","OTU596",
                                   "OTU2368","OTU716","OTU766","OTU224"),
                          guide = guide_axis(n.dodge = 2))+
         scale_fill_manual(values = c("H"="red", "NA"="darkgray"))+
coord_flip()
       
       
       
       
       
       
TenSitesSpecialist
       

####### MI Scatter plot.
Test5 <- read_csv("2sites_FinalTable_Matrix.csv")

ggplot(Test5, aes(MI,Axis.1)) +  
  geom_point(size=5, color="black"#,# alpha = .85
             ,shape=21,fill="orange")+ 
  geom_smooth(method = lm, se = FALSE, colour="black", size=1.5, formula = y ~ x)+
  #scale_color_manual(values = "#FDC086")+
  #facet_wrap(~name, scales = "free_y" , 
  #strip.position="left")+
  ggpubr::stat_regline_equation(label.x =0, label.y = 0.40,
  )+
  ggpubr::stat_cor(label.x = 0, label.y = 0.45, method = "spearman",label.sep = "\n", size = 12
  ) +
  theme_bw(24)+
  labs(x="Management Intensity",y="PCoA 1 of Rare taxa")+
  theme(
    #legend.position = c(0.90, 0.10),
    strip.placement = "outside",
    strip.background = element_blank(),
    #strip.text.x = element_blank(),
    #strip.text.y = element_blank(),
    #axis.line = element_line(colour = "black"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    #panel.border = element_blank(), 
    panel.background = element_blank(),
    #axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank()
  )





