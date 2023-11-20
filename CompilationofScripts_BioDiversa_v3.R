#This script contains the necessary codes for analyses included in the Fungal Rarity manuscript


#********************************************Loading packages and getting files ready**********************************************************
pacman::p_load(rmarkdown,knitr,reticulate, igraph,edgeR, phyloseq, reshape2, RColorBrewer,ggplot2, RVAideMemoire, vegan, gridExtra, cowplot, stringi, vegan, plyr, microbiome, igraph, minerva, indicspecies, cowplot, DescTools, qgraph, PMCMR, maps, PBSmapping, raster, sp, venn, corrplot, betapart)

knitr::opts_chunk$set(echo = TRUE, warning = FALSE)
#opts_chunk$set(dev = c('pdf'), pdf.options(encoding = "ISOLatin9.enc"),fig.path = "images/") 


#Match sequencing samples names
rm(list=ls())
OTU_pre <- read.csv("otu_table.csv", row.names = 1)
Primer_match <- read.csv("Primer_Sample_assignment.csv", na.strings = c("NA","","-"))
Primer_match <- na.omit(Primer_match)
row.names(Primer_match) <- Primer_match$Sample.ID
Primer_match$sequencingName <- paste("s",Primer_match$new.Plate,"_",Primer_match$F.Primer,".",Primer_match$R.Primer, sep = "")
Name_match <- match(colnames(OTU_pre),Primer_match$sequencingName)
OTU <- OTU_pre[,!is.na(Name_match)]
colnames(OTU) <- row.names(Primer_match) [Name_match[!is.na(Name_match)]]

#Getting data ready for analyses
OTU <- OTU[rowSums(OTU)>1,]
OTU <- as.matrix(OTU)
OTU <- OTU[,order(colnames(OTU),decreasing = F)]
OTU_pre <- OTU

otu_taxa <- read.csv("otu_taxa.csv",row.names = 1, stringsAsFactors = F)
otu_taxa <- otu_taxa[which(rownames(otu_taxa) %in% rownames(OTU)),]

otu_taxa_pre <- otu_taxa
design <- read.csv("design.csv",row.names = 1, stringsAsFactors = F)
row.names(design) <- design$Sample_ID
design$Treatment[design$Treatment == "High "]<- "High"
design$Treatment[design$Treatment == "Low "]<- "Low"
design$Treatment[design$Treatment == "High" | design$Treatment == "Low" ] <- "Arable"

design[design$Treatment == "Organic"| design$Treatment == "Conventional","Treatment"] <- "Arable"
design_pre <- design

## create phyloseq object
physeq_table<- phyloseq(otu_table(OTU, taxa_are_rows=T),
                        tax_table(as.matrix(otu_taxa)),
                        sample_data(design))
#Mapping samples 
world = map_data('world')
world <- world %>%
  # the PBSmapping package expects a data frame with these names
  rename(c(long="X", lat="Y", group="PID", order="POS"))%>%
  clipPolys(xlim=c(-10, 25), ylim=c(35, 70), keepExtra=TRUE) %>%
  mutate(highlight=region %in% c("Spain" ,"France","Switzerland", "Germany","Sweden"))

## Map of arable sites
ggplot() +
  geom_polygon(data = world, aes(x=X, y=Y, group=PID, fill=highlight),color='black') +
  coord_fixed(ratio = 1.3) +
  scale_fill_manual(values=c('floralwhite', 'gray'))+
  geom_point(data=design[design$Treatment == "Arable",],
             aes (x = Long, y = Lat),
             colour = 'red',
             size = 1)+
  geom_rect(aes(xmin = -10, xmax = 25, ymin = 35, ymax = 70), colour="black", fill=NA)+
  xlim(-10, 25)+
  ylim(35,70)+
  theme_classic()

## map of grassland sites
ggplot() +
  geom_polygon(data = world, aes(x=X, y=Y, group=PID, fill=highlight),color='black') +
  coord_fixed(ratio = 1.3) +
  scale_fill_manual(values=c('floralwhite', 'gray'))+
  geom_point(data=design[design$Treatment == "Grassland",],
             aes (x = Long, y = Lat),
             colour = 'green',
             size = 1)+
  theme_classic()


#********************************************Alpha Diversity**********************************************************
#rarefaction curve
#create color vector
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,50), col=col_vector[1:50])
rowSums(t(OTU))[order(rowSums(t(OTU)), decreasing = T)]
S <- specnumber(t(OTU))
raremax <- min(rowSums(t(OTU)))
Srare <- rarefy(t(OTU), raremax)
forrare <- as.data.frame(t(OTU))

## plotting the rarefaction curve, colored by country
rare_country <- rare_color <- design$Country[match(row.names(forrare), row.names(design))]
rare_color[which(rare_country == "Spain")] <- col_vector[14]
rare_color[which(rare_country == "France")] <- col_vector[18]
rare_color[which(rare_country == "Switzerland")] <- col_vector[20]
rare_color[which(rare_country == "Germany")] <- col_vector[22]
rare_color[which(rare_country == "Sweden")] <- col_vector[26]

pdf("country_rarefaction.pdf")
rarecurve(forrare, step = 20, label = F, col = rare_color, legend=T)
dev.off()

#plotting rarecurve  by land use type
rare_landuse <- rare_color <- design$Treatment[match(row.names(forrare), row.names(design))]
rare_color[which(rare_landuse == "Arable")] <- "red"
rare_color[which(rare_landuse == "Grassland")] <- "blue"
pdf("landuse_rarefaction.pdf")
rarecurve(forrare, step = 20, label = F, col = rare_color, legend=T)
dev.off()

#rare plot for each country
forrare_SP <- forrare[which(rare_country == "Spain"),]
rare_landuse <- rare_color <- design$Treatment[match(row.names(forrare_SP), row.names(design))]
rare_color[which(rare_landuse == "Arable")] <- "red"
rare_color[which(rare_landuse == "Grassland")] <- "blue"

rarecurve(forrare_SP, step = 20, label = F, col = rare_color)

forrare_FR <- forrare[which(rare_country == "France"),]
rare_landuse <- rare_color <- design$Treatment[match(row.names(forrare_FR), row.names(design))]
rare_color[which(rare_landuse == "Arable")] <- "red"
rare_color[which(rare_landuse == "Grassland")] <- "blue"
rarecurve(forrare_FR, step = 20, label = F, col = rare_color)

forrare_CH <- forrare[which(rare_country == "Switzerland"),]
rare_landuse <- rare_color <- design$Treatment[match(row.names(forrare_CH), row.names(design))]
rare_color[which(rare_landuse == "Arable")] <- "red"
rare_color[which(rare_landuse == "Grassland")] <- "blue"
rarecurve(forrare_CH, step = 20, label = F, col = rare_color)

forrare_GE <- forrare[which(rare_country == "Germany"),]
rare_landuse <- rare_color <- design$Treatment[match(row.names(forrare_GE), row.names(design))]
rare_color[which(rare_landuse == "Arable")] <- "red"
rare_color[which(rare_landuse == "Grassland")] <- "blue"
rarecurve(forrare_GE, step = 20, label = F, col = rare_color)

forrare_SW <- forrare[which(rare_country == "Sweden"),]
rare_landuse <- rare_color <- design$Treatment[match(row.names(forrare_SW), row.names(design))]
rare_color[which(rare_landuse == "Arable")] <- "red"
rare_color[which(rare_landuse == "Grassland")] <- "blue"
rarecurve(forrare_SW, step = 20, label = F, col = rare_color)

plot.new()
legend("bottomleft",legend=c("Spain","France","Switzerland","Germany","Sweden"),col= col_vector[c(14,18,20,22,26)],
       lty=1)

plot.new()
legend("bottomleft",legend=c("Arable land", "Grassland"),col= c("red", "blue"),
       lty=1)

#choose threshold and rarefy
set.seed(1111)
rarefied_OTU_df <- rarefy_even_depth(physeq_table, sample.size =  2000)
rarefied_OTU <- as.data.frame(rarefied_OTU_df@otu_table@.Data)

thrown_samples <- names(colSums(OTU)[colSums(OTU)<2000])
write.csv(
  data.frame(sample.name = thrown_samples,
             country = design$Country[match(thrown_samples, row.names(design))],
             land.use = design$Treatment[match(thrown_samples, row.names(design))],
             sample.depth = unname(colSums(OTU)[colSums(OTU)<2000])),
  "thrown samples.csv")

OTU <- rarefied_OTU
otu_taxa <- otu_taxa[which(rownames(otu_taxa) %in% rownames(OTU)),]
design <- design[which(rownames(design) %in% colnames(OTU)),]
#write.csv(rarefied_OTU,"OTU_rarefied_3thrown.csv", row.names = T)

## create phyloseq object of rarefied OTU
physeq_table<- phyloseq(otu_table(OTU, taxa_are_rows=T),
                        tax_table(as.matrix(otu_taxa)),
                        sample_data(design))

#Land Use and Interaction with Country Identity
alpha_est <- data.frame(
  Richness = estimate_richness(physeq_table , measures = "Observed")[,1],
  Chao1 = estimate_richness(physeq_table , measures = "Chao1")[,1],
  ACE = estimate_richness(physeq_table , measures = "ACE")[,1],
  Shannon = estimate_richness(physeq_table , measures = "Shannon")[,1],
  InvSimpson = estimate_richness(physeq_table , measures = "InvSimpson")[,1],
  Evenness = evenness(physeq_table,'pielou')
)

# simplify treatment to only arable and grassland
alpha_est <- merge(alpha_est, design, by = "row.names")
alpha_est$Treatment[alpha_est$Treatment == "Organic"| alpha_est$Treatment == "Conventional"] <- "Arable"
alpha_est_melt <- melt(alpha_est, measure.vars = c("Treatment","Country"))
alpha_est_melt$value <- factor(alpha_est_melt$value,levels = unique(alpha_est_melt$value))

alpha_est_melt$value <- as.character(alpha_est_melt$value)

alpha_est_melt$value <- factor(alpha_est_melt$value, levels = c("Arable","Grassland","Spain","France","Switzerland","Germany","Sweden"))


rich<-boxplot(alpha_est_melt[,"Richness"] ~ alpha_est_melt$value,
              ylab = "Richness",
              main = "Overall",
              cex.axis = 0.8,
              col = c("red", "blue"))

boxplot(alpha_est_melt[,"InvSimpson"] ~ alpha_est_melt$value,
        ylab = "InvSimpson",
        main = "Overall",
        cex.axis = 0.8,
        col = c("red", "blue"))
boxplot(alpha_est_melt[,"Shannon"] ~ alpha_est_melt$value,
        ylab = "Shannon",
        main = "Overall",
        cex.axis = 0.8,
        col = c("red", "blue"))
boxplot(alpha_est_melt[,"ACE"] ~ alpha_est_melt$value,
        ylab = "ACE",
        main = "Overall",
        cex.axis = 0.8,
        col = c("red", "blue"))
boxplot(alpha_est_melt[alpha_est_melt$variable == "Treatment",][,"Richness"] ~ as.character(alpha_est_melt[alpha_est_melt$variable == "Treatment",]$value),
        ylab = "Richness",
        main = "Overall",
        cex.axis = 0.8,
        col = c("red", "blue"))
boxplot(alpha_est_melt[alpha_est_melt$variable == "Treatment",][,"pielou"] ~ as.character(alpha_est_melt[alpha_est_melt$variable == "Treatment",]$value),
        ylab = "pielou",
        main = "Overall",
        cex.axis = 0.8,
        col = c("red", "blue"))

write.csv(alpha_est, "alpha_est.csv")

#plot arable vs. grassland for each country

CountryNames <- c("Spain","France","Switzerland","Germany","Sweden")
for (i in c(2,4,5,6,7))
{
  boxplot(alpha_est[,i] ~ factor(alpha_est$Treatment) + factor(alpha_est$Country, levels = CountryNames) , col = rep(col_vector[c(3,1)],5))
}

## loss of OTUs in arable land
(mean(alpha_est[alpha_est$Country == "Spain" & alpha_est$Treatment == "Grassland",2])-
    mean(alpha_est[alpha_est$Country == "Spain" & alpha_est$Treatment == "Arable",2]))
mean(alpha_est[alpha_est$Country == "Spain" & alpha_est$Treatment == "Grassland",2])

(mean(alpha_est[alpha_est$Country == "France" & alpha_est$Treatment == "Grassland",2])-
    mean(alpha_est[alpha_est$Country == "France" & alpha_est$Treatment == "Arable",2]))
mean(alpha_est[alpha_est$Country == "France" & alpha_est$Treatment == "Grassland",2])

(mean(alpha_est[alpha_est$Country == "Switzerland" & alpha_est$Treatment == "Grassland",2])-
    mean(alpha_est[alpha_est$Country == "Switzerland" & alpha_est$Treatment == "Arable",2]))
mean(alpha_est[alpha_est$Country == "Switzerland" & alpha_est$Treatment == "Grassland",2])

(mean(alpha_est[alpha_est$Country == "Germany" & alpha_est$Treatment == "Grassland",2])-
    mean(alpha_est[alpha_est$Country == "Germany" & alpha_est$Treatment == "Arable",2]))
mean(alpha_est[alpha_est$Country == "Germany" & alpha_est$Treatment == "Grassland",2])

(mean(alpha_est[alpha_est$Country == "Sweden" & alpha_est$Treatment == "Grassland",2])-
    mean(alpha_est[alpha_est$Country == "Sweden" & alpha_est$Treatment == "Arable",2]))
mean(alpha_est[alpha_est$Country == "Sweden" & alpha_est$Treatment == "Grassland",2])

(mean(alpha_est[alpha_est$Treatment == "Grassland",2])-
    mean(alpha_est[alpha_est$Treatment == "Arable",2]))/
  mean(alpha_est[alpha_est$Treatment == "Grassland",2])


#Geographical Locations
#plot alpha diversity along latitude
row.names(alpha_est)<- alpha_est$Row.names
alpha_est$shape <- 3
alpha_est$shape[alpha_est$Treatment == "Grassland"] <- 16

col_vector[1:5] <-  col_vector[c(14,18,20,22,26)]

alpha_est$color <- col_vector[1]
alpha_est$color[alpha_est$Country == "Germany"] <- col_vector[2]
alpha_est$color[alpha_est$Country == "Spain"] <- col_vector[3]
alpha_est$color[alpha_est$Country == "Sweden"] <- col_vector[4]
alpha_est$color[alpha_est$Country == "Switzerland"] <- col_vector[5]

for (i in c(2,4,5,6))
{
  fit.arable <- lm(alpha_est[row.names(design[design$Treatment == "Arable",]),i] ~ alpha_est[row.names(design[design$Treatment == "Arable",]),]$Lat)
  print(paste("Adj R2 = ",signif(summary(fit.arable)$adj.r.squared, 5), "\nIntercept =",signif(fit.arable$coef[[1]],5 ),                                   " \nSlope =",signif(fit.arable$coef[[2]], 5),                                    " \nP =",signif(summary(fit.arable)$coef[2,4], 5)))
  
  fit.grassland <- lm(alpha_est[row.names(design[design$Treatment == "Grassland",]),i] ~ alpha_est[row.names(design[design$Treatment == "Grassland",]),]$Lat)
  print(paste("Adj R2 = ",signif(summary(fit.grassland)$adj.r.squared, 5), "\nIntercept =",signif(fit.grassland$coef[[1]],5 ),                                   " \nSlope =",signif(fit.grassland$coef[[2]], 5),                                    " \nP =",signif(summary(fit.grassland)$coef[2,4], 5)))
  print(cor.test(alpha_est[row.names(design[design$Treatment == "Grassland",]),i] , alpha_est[row.names(design[design$Treatment == "Grassland",]),]$Lat)$estimate)
  
  
  
  
  print(ggplot()+  geom_point(aes(x = alpha_est$Lat, y = alpha_est[,i]),shape = alpha_est$shape, color = alpha_est$color))
  
  fit.arable <- lm(alpha_est[row.names(design[design$Treatment == "Arable",]),i] ~ poly(alpha_est[row.names(design[design$Treatment == "Arable",]),]$Lat,2,raw = T))
  print(paste("Adj R2 = ",signif(summary(fit.arable)$adj.r.squared, 5),                                    " \nP =",signif(summary(fit.arable)$coef[2,4], 5)))
  
  fit.grassland <- lm(alpha_est[row.names(design[design$Treatment == "Grassland",]),i] ~ poly(alpha_est[row.names(design[design$Treatment == "Grassland",]),]$Lat,2,raw = T))
  print(paste("Adj R2 = ",signif(summary(fit.grassland)$adj.r.squared, 5),                                    " \nP =",signif(summary(fit.grassland)$coef[2,4], 5)))
  
  print(ggplot()+  geom_point(aes(x = alpha_est$Lat, y = alpha_est[,i], colour = alpha_est$Treatment))
        + stat_smooth(aes(x = alpha_est$Lat, y = alpha_est[,i], colour = alpha_est$Treatment),method="lm", formula = y ~ poly(x, 2, raw =TRUE), se = F))
}


#plotting without Switzerland
i = 2
fit.arable <- lm(alpha_est[row.names(design[design$Treatment == "Arable" & design$Country != "Switzerland",]),i] ~ poly(alpha_est[row.names(design[design$Treatment == "Arable"& design$Country != "Switzerland",]),]$Lat,2,raw = T))
print(paste("Adj R2 = ",signif(summary(fit.arable)$adj.r.squared, 5),                                    " \nP =",signif(summary(fit.arable)$coef[2,4], 5)))

fit.grassland <- lm(alpha_est[row.names(design[design$Treatment == "Grassland"& design$Country != "Switzerland",]),i] ~ poly(alpha_est[row.names(design[design$Treatment == "Grassland"& design$Country != "Switzerland",]),]$Lat,2,raw = T))

print(ggplot()+  geom_point(aes(x = alpha_est[row.names(design[design$Country != "Switzerland",]),]$Lat, y = alpha_est[row.names(design[ design$Country != "Switzerland",]),i], colour = alpha_est[row.names(design[design$Country != "Switzerland",]),]$Treatment))
      + stat_smooth(aes(x = alpha_est[row.names(design[design$Country != "Switzerland",]),]$Lat, y = alpha_est[row.names(design[ design$Country != "Switzerland",]),i], colour = alpha_est[row.names(design[design$Country != "Switzerland",]),]$Treatment),method="lm", formula = y ~ poly(x, 2, raw =TRUE), se = F))


#plot alpha diversity along longitude
for (i in c(2,4,5,6))
{
  fit.arable <- lm(alpha_est[row.names(design[design$Treatment == "Arable",]),i] ~ alpha_est[row.names(design[design$Treatment == "Arable",]),]$Long)
  print(paste("Adj R2 = ",signif(summary(fit.arable)$adj.r.squared, 5), "\nIntercept =",signif(fit.arable$coef[[1]],5 ),                                   " \nSlope =",signif(fit.arable$coef[[2]], 5),                                    " \nP =",signif(summary(fit.arable)$coef[2,4], 5)))
  
  fit.grassland <- lm(alpha_est[row.names(design[design$Treatment == "Grassland",]),i] ~ alpha_est[row.names(design[design$Treatment == "Grassland",]),]$Long)
  print(paste("Adj R2 = ",signif(summary(fit.grassland)$adj.r.squared, 5), "\nIntercept =",signif(fit.grassland$coef[[1]],5 ),                                   " \nSlope =",signif(fit.grassland$coef[[2]], 5),                                    " \nP =",signif(summary(fit.grassland)$coef[2,4], 5)))
  
  print(ggplot()+  geom_point(aes(x = alpha_est$Long, y = alpha_est[,i], colour = alpha_est$Treatment)))
}


#********************************************Soil Properties**********************************************************
SoilProperties <- c("C_tot_.","N_tot_.","P_tot_mg","P_Olsen","pH","Organic_C_.","C_microbial","N_microbial","basal_respiration","Sand","Silt","Clay","moisture","CEC")

for (j in SoilProperties)
{
  for (i in c(2,4,5,6))
  {
    print(j)
    print(colnames(alpha_est)[i])
    ## linear fit
    fit.arable <- lm(alpha_est[row.names(design[design$Treatment == "Arable",]),i] ~ alpha_est[row.names(design[design$Treatment == "Arable",]),j])
    print(paste("Adj R2 = ",signif(summary(fit.arable)$adj.r.squared, 5), "\nIntercept =",signif(fit.arable$coef[[1]],5 ),                                   " \nSlope =",signif(fit.arable$coef[[2]], 5),                                    " \nP =",signif(summary(fit.arable)$coef[2,4], 5)))
    
    fit.grassland <- lm(alpha_est[row.names(design[design$Treatment == "Grassland",]),i] ~ alpha_est[row.names(design[design$Treatment == "Grassland",]),j])
    print(paste("Adj R2 = ",signif(summary(fit.grassland)$adj.r.squared, 5), "\nIntercept =",signif(fit.grassland$coef[[1]],5 ),                                   " \nSlope =",signif(fit.grassland$coef[[2]], 5),                                    " \nP =",signif(summary(fit.grassland)$coef[2,4], 5)))
    
    #print(cor.test(alpha_est[row.names(design[design$Treatment == "Arable",]),i] , alpha_est[row.names(design[design$Treatment == "Arable",]),j], method = "pearson")$estimate)
    
    #print(cor.test(alpha_est[row.names(design[design$Treatment == "Grassland",]),i] , alpha_est[row.names(design[design$Treatment == "Grassland",]),j], method = "pearson")$estimate)
    
    print(ggplot()+  geom_point(aes(x = alpha_est[,j], y = alpha_est[,i], colour = alpha_est$Treatment))
          + geom_smooth(aes(x = alpha_est[,j], y = alpha_est[,i], colour = alpha_est$Treatment),method = "lm",se = F))
    
    
    ## second order fit
    fit.arable <- lm(alpha_est[row.names(design[design$Treatment == "Arable",]),i] ~ poly(alpha_est[row.names(design[design$Treatment == "Arable",]),j],2,raw = T))
    print(paste("Adj R2 = ",signif(summary(fit.arable)$adj.r.squared, 5),                                    " \nP =",signif(summary(fit.arable)$coef[2,4], 5)))
    
    fit.grassland <- lm(alpha_est[row.names(design[design$Treatment == "Grassland",]),i] ~ poly(alpha_est[row.names(design[design$Treatment == "Grassland",]),j],2,raw = T))
    print(paste("Adj R2 = ",signif(summary(fit.grassland)$adj.r.squared, 5),                                    " \nP =",signif(summary(fit.grassland)$coef[2,4], 5)))
    
    print(ggplot()+  geom_point(aes(x = alpha_est[,j], y = alpha_est[,i], colour = alpha_est$Treatment))
          + stat_smooth(aes(x = alpha_est[,j], y = alpha_est[,i], colour = alpha_est$Treatment),method="lm", formula = y ~ poly(x, 2, raw =TRUE), se = F))
    
  }
}

#plot correlation
SoilProperties <- c("C_tot_.","N_tot_.","P_tot_mg","P_Olsen","pH","Organic_C_.","C_microbial","N_microbial","basal_respiration","Sand","Silt","Clay","moisture","CEC","MIN")

Arable <- row.names(design[design$Treatment=="Arable",])
Grassland <- row.names(design[design$Treatment=="Grassland",])

corr_div_soil <- data.frame(matrix(ncol = 6, nrow = 15))
row.names(corr_div_soil) <- SoilProperties
colnames(corr_div_soil) <- c("Richness_arable","Shannon-Wiener_arable","ACE_arable","Richness_grassland","Shannon-Wiener_grassland","ACE_grassland")

DiversityIndex <- c("Richness","Shannon","ACE")
NnEnzyme <- read.csv("N and Enzymes.csv", row.names = 1)
design$MIN <- NnEnzyme[row.names(design),"MIN"]
soildf <- design[,SoilProperties]

for (i in 1:15)
  for (j in 1:3)
  {
    corr_div_soil[i,j] <- cor.test(~alpha_est[Arable,DiversityIndex[j]] +
                                     soildf[Arable,i], method = "spearman"
    )$estimate
  }

for (i in 1:15)
  for (j in 4:6)
  {
    corr_div_soil[i,j] <- cor.test(~alpha_est[Grassland,DiversityIndex[j-3]] +
                                     soildf[Grassland,i], method = "spearman"
    )$estimate
  }
row.names(corr_div_soil) <- c("Total C","Total N","Total P","Olsen P","pH","Organic C","Microbial C","Microbial N","Basal Respiration", "Sand","Silt","Clay","Moisture","CEC","MIN")
corrplot(t(as.matrix(corr_div_soil)), is.corr=FALSE, tl.srt=41, method = "circle", tl.cex =0.9)

## correlation among bioclimatic variables
colnames(soildf) <- row.names(corr_div_soil)
corrplot(as.matrix(cor(na.omit(soildf), method = "spearman")), method = "number", type = "lower", diag = F)

for (i in 1:15)
  for (j in 1:3)
  {
    corr_div_soil[i,j] <- cor.test(~alpha_est[Arable,DiversityIndex[j]] +
                                     soildf[Arable,i], method = "spearman"
    )$p.value
  }

for (i in 1:15)
  for (j in 4:6)
  {
    corr_div_soil[i,j] <- cor.test(~alpha_est[Grassland,DiversityIndex[j-3]] +
                                     soildf[Grassland,i], method = "spearman"
    )$p.value
  }
corrplot(t(as.matrix(corr_div_soil)), is.corr=FALSE, tl.srt=41, method = "number", tl.cex =0.9)

plot(x = alpha_est$C_tot_., y = alpha_est$Richness)
cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot(data = alpha_est)+
  geom_point(aes(x = C_tot_., y = Richness, color = Treatment)) + 
  scale_color_manual(values = cols[c(3,1)])+
  xlab("Total Carbon Content")+
  geom_smooth(aes(x = C_tot_., y = Richness, color = Treatment), method = "lm", se = F)


#PCA of soil properties

col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pca1<-princomp(na.omit(soildf[,c(1:15)]),cor=T, scale = T)
summary(pca1)
LoadingComp1n2 <- as.data.frame(pca1$loadings[1:15,1:2])
ScoreComp1n2 <- as.data.frame(pca1$scores[,1:2])
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

alpha_est <- merge(alpha_est[,-1],-ScoreComp1n2,by = "row.names")
row.names(alpha_est) <- alpha_est$Row.names
for (j in c("Comp.1","Comp.2"))
{
  for (i in c(2,4,5,6))
  {
    print(j)
    print(colnames(alpha_est)[i])
    fit.arable <- lm(alpha_est[row.names(design[design$Treatment == "Arable",]),i] ~ alpha_est[row.names(design[design$Treatment == "Arable",]),j])
    print(paste("Adj R2 = ",signif(summary(fit.arable)$adj.r.squared, 5), "\nIntercept =",signif(fit.arable$coef[[1]],5 ),                                   " \nSlope =",signif(fit.arable$coef[[2]], 5),                                    " \nP =",signif(summary(fit.arable)$coef[2,4], 5)))
    
    fit.grassland <- lm(alpha_est[row.names(design[design$Treatment == "Grassland",]),i] ~ alpha_est[row.names(design[design$Treatment == "Grassland",]),j])
    print(paste("Adj R2 = ",signif(summary(fit.grassland)$adj.r.squared, 5), "\nIntercept =",signif(fit.grassland$coef[[1]],5 ),                                   " \nSlope =",signif(fit.grassland$coef[[2]], 5),                                    " \nP =",signif(summary(fit.grassland)$coef[2,4], 5)))
    print(ggplot()+  geom_point(aes(x = alpha_est[,j], y = alpha_est[,i], colour = alpha_est$Treatment))+ geom_smooth(aes(x = alpha_est[,j], y = alpha_est[,i],color = alpha_est$Treatment),method = "lm",se = F)+ 
            scale_color_manual(values = col_vector[c(3,1)]) )
  }
}

soildf$Country <- design[row.names(soildf),"Country"]
soildf$Treatment <- design[row.names(soildf),"Treatment"]
soildf <- na.omit(soildf)
library(ggfortify)
autoplot(pca1, data = soildf, colour = "Country", shape = "Treatment", loadings = T, loadings.label = T, loadings.label.size = 4)


#********************************************Climate data**********************************************************
r <- raster::getData("worldclim",var="bio",res=10)
coords <- data.frame(design$Long,design$Lat)
points <- SpatialPoints(coords, proj4string =r@crs)
climatedf <- as.data.frame(extract(r,points))
row.names(climatedf)<- row.names(design)

#plot correlation
Arable <- row.names(design[design$Treatment=="Arable",])
Grassland <- row.names(design[design$Treatment=="Grassland",])

corr_div_clim <- data.frame(matrix(ncol = 6, nrow = 19))
row.names(corr_div_clim) <- paste("bio",c(1:19),sep = "")
colnames(corr_div_clim) <- c("Richness_arable","Shannon-Wiener_arable","ACE_arable","Richness_grassland","Shannon-Wiener_grassland","ACE_grassland")

DiversityIndex <- c("Richness","Shannon","ACE")


for (i in 1:19)
  for (j in 1:3)
  {
    corr_div_clim[i,j] <- cor.test(~alpha_est[Arable,DiversityIndex[j]] +
                                     climatedf[Arable,i], method = "spearman"
    )$estimate
  }

for (i in 1:19)
  for (j in 4:6)
  {
    corr_div_clim[i,j] <- cor.test(~alpha_est[Grassland,DiversityIndex[j-3]] +
                                     climatedf[Grassland,i], method = "spearman"
    )$estimate
  }
corrplot(t(as.matrix(corr_div_clim)), is.corr=FALSE, tl.srt=45, method = "number")

## correlation among bioclimatic variables
corrplot(as.matrix(cor(climatedf, method = "spearman")), method = "number", type = "lower", diag = F)

write.csv(climatedf, "climatedf.csv")

## replace bio with real names
row.names(corr_div_clim)<- c("Annual Mean Temperature","Mean Diurnal Range","Isothermality","Temperature Seasonality","Max Temperature of Warmest Month","Min Temperature of Coldest Month","Temperature Annual Range","Mean Temperature of Wettest Quarter","Mean Temperature of Driest Quarter","Mean Temperature of Warmest Quarter","Mean Temperature of Coldest Quarter","Annual Precipitation","Precipitation of Wettest Month","Precipitation of Driest Month","Precipitation Seasonality","Precipitation of Wettest Quarter","Precipitation of Driest Quarter","Precipitation of Warmest Quarter","Precipitation of Coldest Quarter")
corrplot(t(as.matrix(corr_div_clim)), is.corr=FALSE, tl.srt=45, method = "circle", tl.cex =0.7)

pca2<-princomp(climatedf[,c(1:9,12:15,18,19)],cor=T, scale = T)
summary(pca2)
LoadingComp1n2 <- as.data.frame(pca2$loadings[1:15,1:2])
row.names(LoadingComp1n2) <- c("Annual Mean Temperature","Mean Diurnal Range (Mean of monthly (max temp - min temp))","Isothermality (BIO2/BIO7) (* 100)","Temperature Seasonality (standard deviation *100)","Max Temperature of Warmest Month","Min Temperature of Coldest Month","Temperature Annual Range (BIO5-BIO6)","Mean Temperature of Wettest Quarter","Mean Temperature of Driest Quarter","Mean Temperature of Warmest Quarter","Mean Temperature of Coldest Quarter","Annual Precipitation","Precipitation of Wettest Month","Precipitation of Driest Month","Precipitation Seasonality (Coefficient of Variation)","Precipitation of Wettest Quarter","Precipitation of Driest Quarter","Precipitation of Warmest Quarter","Precipitation of Coldest Quarter")[c(1:9,12:15,18,19)]
write.csv(LoadingComp1n2, "LoadingComp1n2.csv")
ScoreComp1n2 <- as.data.frame(pca2$scores[,1:2])


alpha_est <- merge(alpha_est[,-1],ScoreComp1n2,by = "row.names")
row.names(alpha_est) <- alpha_est$Row.names
for (j in c("Comp.1.y","Comp.2.y"))
{
  for (i in c(2,4,5,6))
  {
    print(j)
    print(colnames(alpha_est)[i])
    fit.arable <- lm(alpha_est[row.names(design[design$Treatment == "Arable",]),i] ~ alpha_est[row.names(design[design$Treatment == "Arable",]),j])
    print(paste("Adj R2 = ",signif(summary(fit.arable)$adj.r.squared, 5), "\nIntercept =",signif(fit.arable$coef[[1]],5 ),                                   " \nSlope =",signif(fit.arable$coef[[2]], 5),                                    " \nP =",signif(summary(fit.arable)$coef[2,4], 5)))
    
    fit.grassland <- lm(alpha_est[row.names(design[design$Treatment == "Grassland",]),i] ~ alpha_est[row.names(design[design$Treatment == "Grassland",]),j])
    print(paste("Adj R2 = ",signif(summary(fit.grassland)$adj.r.squared, 5), "\nIntercept =",signif(fit.grassland$coef[[1]],5 ),                                   " \nSlope =",signif(fit.grassland$coef[[2]], 5),                                    " \nP =",signif(summary(fit.grassland)$coef[2,4], 5)))
    print(ggplot()+  geom_point(aes(x = alpha_est[,j], y = alpha_est[,i], colour = alpha_est$Treatment))+ geom_smooth(aes(x = alpha_est[,j], y = alpha_est[,i],color = alpha_est$Treatment),method = "lm",se = F))
  }
}

ggplot()+
  geom_bar(fill = "black",aes(y = -LoadingComp1n2$Comp.1, x = factor(row.names(LoadingComp1n2), levels = rev(row.names(LoadingComp1n2)))), stat = "identity")+ coord_flip()+ xlab("")+ylab("Loadings for Principal Component 1")

climatedf$Country <- design[row.names(climatedf),"Country"]
climatedf$Treatment <- design[row.names(climatedf),"Treatment"]
#soildf <- na.omit(soildf)
#library(ggfortify)
autoplot(pca2, data = climatedf, colour = "Country", shape = "Treatment", loadings = T, loadings.label = T, loadings.label.size = 4)

#plot relation between mean temperature, mean precipitation and alpha diversity
row.names(alpha_est)<- alpha_est$Row.names
alpha_est <- merge(alpha_est[,-1], climatedf[,c("bio1","bio12")],by="row.names")
row.names(alpha_est)<- alpha_est$Row.names
for (j in c("bio1","bio12"))
{
  for (i in c(2,4,5,6))
  {
    print(j)
    print(colnames(alpha_est)[i])
    fit.arable <- lm(alpha_est[row.names(design[design$Treatment == "Arable",]),i] ~ alpha_est[row.names(design[design$Treatment == "Arable",]),j])
    print(paste("Adj R2 = ",signif(summary(fit.arable)$adj.r.squared, 5), "\nIntercept =",signif(fit.arable$coef[[1]],5 ),                                   " \nSlope =",signif(fit.arable$coef[[2]], 5),                                    " \nP =",signif(summary(fit.arable)$coef[2,4], 5)))
    
    fit.grassland <- lm(alpha_est[row.names(design[design$Treatment == "Grassland",]),i] ~ alpha_est[row.names(design[design$Treatment == "Grassland",]),j])
    print(paste("Adj R2 = ",signif(summary(fit.grassland)$adj.r.squared, 5), "\nIntercept =",signif(fit.grassland$coef[[1]],5 ),                                   " \nSlope =",signif(fit.grassland$coef[[2]], 5),                                    " \nP =",signif(summary(fit.grassland)$coef[2,4], 5)))
    print(ggplot()+  geom_point(aes(x = alpha_est[,j], y = alpha_est[,i], colour = alpha_est$Treatment)))
  }
}

#********************************************Community composition**********************************************************
#Class level - relative abundance arable vs. fungi in all and five countries
otu_taxa <- otu_taxa_pre
OTU<- OTU_pre
design <- design_pre
rm(list=ls()[-match(c("OTU","otu_taxa","design"),ls())])

OTU_RA <- t(t(OTU)/colSums(OTU)) * 100
OTU_norm <- sqrt(OTU)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
## plot frequency distribution of relative abundance
ggplot(as.data.frame(rowMeans(OTU_RA)), aes(x=as.data.frame(rowMeans(OTU_RA)))) + geom_density()

## interdependence of number of sites occupied and mean abundance
presence <- c(rep(0,length(OTU_RA[,1])))
for (i in 1:length(OTU_RA[,1]))
{
  presence[i] <- sum(OTU_RA[i,] != 0)
}
Arable <- row.names(design[design$Treatment == "Arable",])
Grassland <- row.names(design[design$Treatment == "Grassland",])
#presencedf <- data.frame(presence = presence, log_abundance = log(rowMeans(OTU_RA)))
#row.names(presencedf) <- row.names(OTU_RA)


#ggplot(data = presencedf, aes(x = presence, y = log(rowMeans(OTU_RA))))+
#  geom_point(aes(color = ))
plot(presence, log(rowSums(OTU_norm)/presence))

## arable

OTU_RA_arable <- OTU_RA[,Arable]
OTU_norm_arable <- OTU_norm[,Arable]

presence_arable <- c(rep(0,length(OTU_RA[,1])))
for (i in 1:length(OTU_RA[,1]))
{
  presence_arable[i] <- sum(OTU_RA[i,Arable] != 0)
}

plot(presence_arable, log(rowSums(OTU_norm_arable)/presence_arable))

## grassland

OTU_RA_grassland <- OTU_RA[,Grassland]
OTU_norm_grassland <- OTU_norm[,Grassland]

presence_grassland <- c(rep(0,length(OTU_RA_grassland[,1])))
for (i in 1:length(OTU_RA[,1]))
{
  presence_grassland[i] <- sum(OTU_RA[i,Grassland] != 0)
}

plot(presence_grassland, log(rowSums(OTU_RA_grassland)/presence_grassland))

ggplot()+
  geom_point(aes(x = presence_grassland, y = log(rowMeans(OTU_RA_grassland))), color = "blue")+
  geom_point(aes(x = presence_arable, y = log(rowMeans(OTU_RA_arable))), color = "red")

## number of arbla sites are too much bigger, randomly select 61
set.seed(2018)
s <- sample(c(1:156),61)
OTU_RA_arable_61 <- OTU_RA_arable[,s]
OTU_norm_arable_61 <- OTU_norm_arable[,s]
presence_arable <- c(rep(0,length(OTU_RA_arable_61[,1])))
for (i in 1:length(OTU_RA_arable_61[,1]))
{
  presence_arable[i] <- sum(OTU_RA_arable_61[i,] != 0)
}

## take away 0
OTU_norm_grassland <- OTU_norm_grassland[presence_grassland!=0,]
presence_grassland <- presence_grassland[presence_grassland!=0]
OTU_norm_arable_61 <- OTU_norm_arable_61[presence_arable!=0,]
presence_arable <- presence_arable[presence_arable!=0]


ggplot()+
  geom_point(aes(x = presence_grassland, y = log(rowSums(OTU_norm_grassland)/61)), color = "black", size =2,alpha = 0.5)+
  geom_point(aes(x = presence_arable, y = log(rowSums(OTU_norm_arable_61)/61)), color = col_vector[3], size =2,alpha = 0.5)+
  ylab("log (Normalized OTU Abundance)")+
  xlab("Number of Sites Occupied")

## check significance
df1 <- data.frame(presence = presence_grassland, abundance = rowSums(OTU_norm_grassland)/61,type = "grassland")
df2 <- data.frame(presence = presence_arable, abundance = rowSums(OTU_norm_arable_61)/61, type = "arable")
df <- rbind(df1, df2)



#  stat_smooth(aes(x = presence_grassland, y = log10(rowSums(OTU_norm_grassland)/61)), color = col_vector[1],method = "lm", formula = y ~ poly(x, 2, raw =TRUE),  se = F)
#  geom_smooth(aes(x = presence_grassland, y = log(rowMeans(OTU_RA_grassland))),method = "lm", se = F)+
#  geom_smooth(aes(x = presence_arable, y = log(rowMeans(OTU_RA_arable_61))),method = "lm", se = F)

# get names of fungal classes
names_fun <- names(sort(table(otu_taxa[,"Class"]), decr=T))

## Preparation of matrix with relative abundance by CLASS
y <- NULL
otunames <- rownames(OTU_RA)
for (i in names_fun){
  x <- array(colSums(OTU_RA[rownames(otu_taxa)[which(otu_taxa$Class == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(names_fun)
colnames(y) <- paste(colnames(OTU_RA))
CLASS_mat_fun <- as.data.frame(y)

CLASS_mat_fun$mean <- rowMeans(CLASS_mat_fun)
CLASS_mat_fun <- CLASS_mat_fun[order(CLASS_mat_fun$mean,decreasing = TRUE),]

plot(x= log(rowMeans(CLASS_mat_fun[,Grassland])), y =(rowMeans(CLASS_mat_fun[,Arable])-rowMeans(CLASS_mat_fun[,Grassland]))/rowMeans(CLASS_mat_fun[,Grassland]))

## select classes with MEAN abundances top 9
abundant_class_fun <- rownames(CLASS_mat_fun)[1:10]
CLASS_mat_fun <- CLASS_mat_fun[abundant_class_fun,]
others <- 100-colSums(CLASS_mat_fun)
CLASS_mat_fun["c__",]<- CLASS_mat_fun["c__",] + others
row.names(CLASS_mat_fun)[which(abundant_class_fun == "c__")] <- "others"
CLASS_mat_fun <- CLASS_mat_fun[c(c(1:10)[-which(abundant_class_fun == "c__")],which(abundant_class_fun == "c__")),]

OTU_special_present <- as.data.frame(OTU_RA)
OTU_special_present$present_grassland <- rowMeans(OTU_RA[,Grassland])
OTU_special_present$present_arable <- rowMeans(OTU_RA[,Arable])
OTU_special_present$special <- OTU_special_present$present_grassland * OTU_special_present$present_arable
OTU_special_present <- OTU_special_present[OTU_special_present$special == 0,]

length(OTU_special_present$present_grassland[OTU_special_present$present_grassland!=0])
length(OTU_special_present$present_arable[OTU_special_present$present_arable!=0])

OTU_special_arable <- data.frame(Class = otu_taxa[row.names(OTU_special_present[OTU_special_present$present_arable!=0,]),]$Class, RA = OTU_special_present$present_arable[OTU_special_present$present_arable!=0])

OTU_special_grassland <- data.frame(Class = otu_taxa[row.names(OTU_special_present[OTU_special_present$present_grassland!=0,]),]$Class, RA = OTU_special_present$present_grassland[OTU_special_present$present_grassland!=0])

cols <- c("grey50",col_vector[c(1:7,9,22)])

class_special_arable <- aggregate(OTU_special_arable$RA, by = list(as.character(OTU_special_arable$Class)), sum)

class_special_grassland <- aggregate(OTU_special_grassland$RA, by = list(as.character(OTU_special_grassland$Class)), sum)

class_special_arable_top10 <- class_special_arable[order(class_special_arable$x, decreasing = T),][1:10,]
class_special_arable_top10$type <- "Arable"
class_special_arable_top10$Group.1 <- factor(class_special_arable_top10$Group.1, level = c("c__",rev(class_special_arable_top10$Group.1[-1])) )
#class_special_arable_top10$color <- cols

ggplot(class_special_arable_top10, aes(fill=Group.1, y=x, x=type)) + 
  geom_bar( stat="identity") + 
  scale_fill_manual(values = cols)

class_special_grassland_top10 <- class_special_grassland[order(class_special_grassland$x, decreasing = T),][1:10,]
class_special_grassland_top10$type <- "Grassland"

class_special_grassland_top10$Group.1 <- factor(class_special_grassland_top10$Group.1, level = c("c__",rev(class_special_grassland_top10$Group.1[-2])) )
#class_special_grassland_top10$color <- cols[c(2,1,3:10)]



ggplot(class_special_grassland_top10, aes(fill=Group.1, y=x, x=type)) + 
  geom_bar( stat="identity") + 
  scale_fill_manual(values = cols)


#check significant differences
```{r}
fun_mat <- as.data.frame(t(CLASS_mat_fun))
fun_mat$Country <- design$Country[match(row.names(fun_mat),row.names(design))]
fun_mat$Treatment <- design$Treatment[match(row.names(fun_mat),row.names(design))]
fun_mat[,"Country"] <- as.factor(fun_mat[,"Country"] )
fun_mat[,"Treatment"] <- as.factor(fun_mat[,"Treatment"] )
fun_mat$Treatment[fun_mat$Treatment == "Conventional" | fun_mat$Treatment == "Organic"] <- "Arable"

for (i in 1:9)
{
  print(colnames(fun_mat)[i])
  
  print(kruskal.test(fun_mat[,i] ~ fun_mat[,"Country"] ))
  print(kruskal.test(fun_mat[,i] ~ fun_mat[,"Treatment"] ))
  print(posthoc.kruskal.nemenyi.test(fun_mat[,i],fun_mat[,"Country"]))
  print(posthoc.kruskal.nemenyi.test(fun_mat[,i],fun_mat[,"Treatment"]))
  #print(PostHocTest(aov1, method = "duncan") )
}  


#calculate mean and s.e. for arable vs. grassland and for each country
# overall
CLASS_mat_fun$overall_mean <- rowMeans(CLASS_mat_fun)
CLASS_mat_fun$overall_se <- apply(CLASS_mat_fun,1, se)
CLASS_mat_fun$overall_begin <- CLASS_mat_fun$overall_mean - CLASS_mat_fun$overall_se
CLASS_mat_fun$overall_end <- CLASS_mat_fun$overall_mean + CLASS_mat_fun$overall_se

## simplify treatment to arable and grassland
design$Treatment[design$Treatment == "Conventional" | design$Treatment == "Organic"] <- "Arable"

## Arable
CLASS_mat_fun$Arable_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable",])])
CLASS_mat_fun$Arable_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable",])],1,se)
CLASS_mat_fun$Arable_begin <- CLASS_mat_fun$Arable_mean - CLASS_mat_fun$Arable_se
CLASS_mat_fun$Arable_end <- CLASS_mat_fun$Arable_mean + CLASS_mat_fun$Arable_se


## Grassland
CLASS_mat_fun$Grassland_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland",])])
CLASS_mat_fun$Grassland_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland",])],1,se)
CLASS_mat_fun$Grassland_begin <- CLASS_mat_fun$Grassland_mean - CLASS_mat_fun$Grassland_se
CLASS_mat_fun$Grassland_end <- CLASS_mat_fun$Grassland_mean + CLASS_mat_fun$Grassland_se


## sp
CLASS_mat_fun$sp_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Country == "Spain",])])
CLASS_mat_fun$sp_se <- apply(CLASS_mat_fun[,row.names(design[design$Country == "Spain",])],1, se)
CLASS_mat_fun$sp_begin <- CLASS_mat_fun$sp_mean - CLASS_mat_fun$sp_se
CLASS_mat_fun$sp_end <- CLASS_mat_fun$sp_mean + CLASS_mat_fun$sp_se


## fr
CLASS_mat_fun$fr_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Country == "France",])])
CLASS_mat_fun$fr_se <- apply(CLASS_mat_fun[,row.names(design[design$Country == "France",])],1,se)
CLASS_mat_fun$fr_begin <- CLASS_mat_fun$fr_mean - CLASS_mat_fun$fr_se
CLASS_mat_fun$fr_end <- CLASS_mat_fun$fr_mean + CLASS_mat_fun$fr_se


## ch
CLASS_mat_fun$ch_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Country == "Switzerland",])])
CLASS_mat_fun$ch_se <- apply(CLASS_mat_fun[,row.names(design[design$Country == "Switzerland",])],1,se)
CLASS_mat_fun$ch_begin <- CLASS_mat_fun$ch_mean - CLASS_mat_fun$ch_se
CLASS_mat_fun$ch_end <- CLASS_mat_fun$ch_mean + CLASS_mat_fun$ch_se


## ge
CLASS_mat_fun$ge_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Country == "Germany",])])
CLASS_mat_fun$ge_se <- apply(CLASS_mat_fun[,row.names(design[design$Country == "Germany",])],1,se)
CLASS_mat_fun$ge_begin <- CLASS_mat_fun$ge_mean - CLASS_mat_fun$ge_se
CLASS_mat_fun$ge_end <- CLASS_mat_fun$ge_mean + CLASS_mat_fun$ge_se

## sw
CLASS_mat_fun$sw_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Country == "Sweden",])])
CLASS_mat_fun$sw_se <- apply(CLASS_mat_fun[,row.names(design[design$Country == "Sweden",])],1, se)
CLASS_mat_fun$sw_begin <- CLASS_mat_fun$sw_mean - CLASS_mat_fun$sw_se
CLASS_mat_fun$sw_end <- CLASS_mat_fun$sw_mean + CLASS_mat_fun$sw_se


## data transformation
forplotting <- data.frame(Var1 = rep(row.names(CLASS_mat_fun),times=8),
                          Var2 = rep(c("overall","sp","fr","ch","ge","sw","Arable","Grassland"),each=length(CLASS_mat_fun$ge_begin)),
                          mean=c(CLASS_mat_fun$overall_mean,CLASS_mat_fun$sp_mean,CLASS_mat_fun$fr_mean,CLASS_mat_fun$ch_mean,CLASS_mat_fun$ge_mean,CLASS_mat_fun$sw_mean,CLASS_mat_fun$Arable_mean,CLASS_mat_fun$Grassland_mean),
                          begin=c(CLASS_mat_fun$overall_begin, CLASS_mat_fun$sp_begin,CLASS_mat_fun$fr_begin,CLASS_mat_fun$ch_begin,CLASS_mat_fun$ge_begin,CLASS_mat_fun$sw_begin,CLASS_mat_fun$Arable_begin,CLASS_mat_fun$Grassland_begin),
                          end=c(CLASS_mat_fun$overall_end, CLASS_mat_fun$sp_end,CLASS_mat_fun$fr_end,CLASS_mat_fun$ch_end,CLASS_mat_fun$ge_end,CLASS_mat_fun$sw_end,CLASS_mat_fun$Arable_end,CLASS_mat_fun$Grassland_end))

mean_fun <- forplotting$mean
l <- length(CLASS_mat_fun$overall_mean)
for (i in 2:l)
{
  for(j in 0:7)
  {
    forplotting[i+j*l,"begin"]= forplotting[i+j*l,"begin"] + forplotting[i-1+j*l,"mean"]
    forplotting[i+j*l,"end"]=forplotting[i+j*l,"end"] + forplotting[i-1+j*l,"mean"]
    forplotting[i+j*l,"mean"]=forplotting[i+j*l,"mean"] + forplotting[i-1+j*l,"mean"]
  }
}

forplotting <- cbind(forplotting,real_mean = mean_fun)
forplotting$Var1 <- factor(forplotting$Var1, levels = rev(forplotting$Var1[1:length(CLASS_mat_fun$ge_begin)]))
forplotting$Var2<-factor(forplotting$Var2, levels = unique(forplotting$Var2))


##### Plot bar chart
dodge <- position_dodge(width = 0.9)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

cols <- c("grey50",col_vector[c(1:7,9,22)])
limits <- aes(ymax = begin , ymin = end)
ggplot(forplotting, aes(x = Var2, y = real_mean, fill = Var1)) + geom_bar(stat="identity", color = "black" ) + 
  scale_fill_manual(values = cols) +
  geom_errorbar(limits, colour = "black", width = 1,  position = dodge) +
  ylab("Relative Abundance")+ 
  labs(fill="Class")+
  xlab("")+
  theme_bw() +
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank()) 

#class level - for each country plot arable vs. grassland
## Arable Spain
CLASS_mat_fun$Arable_Spain_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable" & design$Country == "Spain",])])
CLASS_mat_fun$Arable_Spain_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable"& design$Country == "Spain",])],1,se)
CLASS_mat_fun$Arable_Spain_begin <- CLASS_mat_fun$Arable_Spain_mean - CLASS_mat_fun$Arable_Spain_se
CLASS_mat_fun$Arable_Spain_end <- CLASS_mat_fun$Arable_Spain_mean + CLASS_mat_fun$Arable_Spain_se


## Grassland Spain
CLASS_mat_fun$Grassland_Spain_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland" & design$Country == "Spain",])])
CLASS_mat_fun$Grassland_Spain_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland"& design$Country == "Spain",])],1,se)
CLASS_mat_fun$Grassland_Spain_begin <- CLASS_mat_fun$Grassland_Spain_mean - CLASS_mat_fun$Grassland_Spain_se
CLASS_mat_fun$Grassland_Spain_end <- CLASS_mat_fun$Grassland_Spain_mean + CLASS_mat_fun$Grassland_Spain_se

## Arable France
CLASS_mat_fun$Arable_France_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable" & design$Country == "France",])])
CLASS_mat_fun$Arable_France_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable"& design$Country == "France",])],1,se)
CLASS_mat_fun$Arable_France_begin <- CLASS_mat_fun$Arable_France_mean - CLASS_mat_fun$Arable_France_se
CLASS_mat_fun$Arable_France_end <- CLASS_mat_fun$Arable_France_mean + CLASS_mat_fun$Arable_France_se

## Grassland France
CLASS_mat_fun$Grassland_France_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland" & design$Country == "France",])])
CLASS_mat_fun$Grassland_France_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland"& design$Country == "France",])],1,se)
CLASS_mat_fun$Grassland_France_begin <- CLASS_mat_fun$Grassland_France_mean - CLASS_mat_fun$Grassland_France_se
CLASS_mat_fun$Grassland_France_end <- CLASS_mat_fun$Grassland_France_mean + CLASS_mat_fun$Grassland_France_se

## Arable Switzerland
CLASS_mat_fun$Arable_Switzerland_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable" & design$Country == "Switzerland",])])
CLASS_mat_fun$Arable_Switzerland_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable"& design$Country == "Switzerland",])],1,se)
CLASS_mat_fun$Arable_Switzerland_begin <- CLASS_mat_fun$Arable_Switzerland_mean - CLASS_mat_fun$Arable_Switzerland_se
CLASS_mat_fun$Arable_Switzerland_end <- CLASS_mat_fun$Arable_Switzerland_mean + CLASS_mat_fun$Arable_Switzerland_se


## Grassland Switzerland
CLASS_mat_fun$Grassland_Switzerland_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland" & design$Country == "Switzerland",])])
CLASS_mat_fun$Grassland_Switzerland_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland"& design$Country == "Switzerland",])],1,se)
CLASS_mat_fun$Grassland_Switzerland_begin <- CLASS_mat_fun$Grassland_Switzerland_mean - CLASS_mat_fun$Grassland_Switzerland_se
CLASS_mat_fun$Grassland_Switzerland_end <- CLASS_mat_fun$Grassland_Switzerland_mean + CLASS_mat_fun$Grassland_Switzerland_se

## Arable Germany
CLASS_mat_fun$Arable_Germany_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable" & design$Country == "Germany",])])
CLASS_mat_fun$Arable_Germany_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable"& design$Country == "Germany",])],1,se)
CLASS_mat_fun$Arable_Germany_begin <- CLASS_mat_fun$Arable_Germany_mean - CLASS_mat_fun$Arable_Germany_se
CLASS_mat_fun$Arable_Germany_end <- CLASS_mat_fun$Arable_Germany_mean + CLASS_mat_fun$Arable_Germany_se


## Grassland Germany
CLASS_mat_fun$Grassland_Germany_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland" & design$Country == "Germany",])])
CLASS_mat_fun$Grassland_Germany_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland"& design$Country == "Germany",])],1,se)
CLASS_mat_fun$Grassland_Germany_begin <- CLASS_mat_fun$Grassland_Germany_mean - CLASS_mat_fun$Grassland_Germany_se
CLASS_mat_fun$Grassland_Germany_end <- CLASS_mat_fun$Grassland_Germany_mean + CLASS_mat_fun$Grassland_Germany_se

## Arable Sweden
CLASS_mat_fun$Arable_Sweden_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable" & design$Country == "Sweden",])])
CLASS_mat_fun$Arable_Sweden_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable"& design$Country == "Sweden",])],1,se)
CLASS_mat_fun$Arable_Sweden_begin <- CLASS_mat_fun$Arable_Sweden_mean - CLASS_mat_fun$Arable_Sweden_se
CLASS_mat_fun$Arable_Sweden_end <- CLASS_mat_fun$Arable_Sweden_mean + CLASS_mat_fun$Arable_Sweden_se


## Grassland Sweden
CLASS_mat_fun$Grassland_Sweden_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland" & design$Country == "Sweden",])])
CLASS_mat_fun$Grassland_Sweden_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland"& design$Country == "Sweden",])],1,se)
CLASS_mat_fun$Grassland_Sweden_begin <- CLASS_mat_fun$Grassland_Sweden_mean - CLASS_mat_fun$Grassland_Sweden_se
CLASS_mat_fun$Grassland_Sweden_end <- CLASS_mat_fun$Grassland_Sweden_mean + CLASS_mat_fun$Grassland_Sweden_se

## data transformation
forplotting <- data.frame(Var1 = rep(row.names(CLASS_mat_fun),times=10),
                          Var2 = rep(c("Arable_Spain","Grassland_Spain" ,"Arable_France","Grassland_France", "Arable_Switzerland","Grassland_Switzerland", "Arable_Germany","Grassland_Germany", "Arable_Sweden","Grassland_Sweden"  ),each=length(CLASS_mat_fun$ge_begin)),
                          mean=unlist(c(CLASS_mat_fun[,c(62:71)*4+3]),use.names = F),
                          begin=unlist(c(CLASS_mat_fun[,c(62:71)*4+5]),use.names = F),
                          end=unlist(c(CLASS_mat_fun[,c(62:71)*4+6])),use.names = F)

mean_fun <- forplotting$mean
l <- length(CLASS_mat_fun$overall_mean)
for (i in 2:l)
{
  for(j in 0:9)
  {
    forplotting[i+j*l,"begin"]= forplotting[i+j*l,"begin"] + forplotting[i-1+j*l,"mean"]
    forplotting[i+j*l,"end"]=forplotting[i+j*l,"end"] + forplotting[i-1+j*l,"mean"]
    forplotting[i+j*l,"mean"]=forplotting[i+j*l,"mean"] + forplotting[i-1+j*l,"mean"]
  }
}

forplotting <- cbind(forplotting,real_mean = mean_fun)
forplotting$Var1 <- factor(forplotting$Var1, levels = rev(forplotting$Var1[1:length(CLASS_mat_fun$ge_begin)]))
forplotting$Var2<-factor(forplotting$Var2, levels = unique(forplotting$Var2))
 
##### Plot barchart
dodge <- position_dodge(width = 0.9)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

cols <- c("grey50",col_vector[c(1:7,9,22)])
limits <- aes(ymax = begin , ymin = end)
ggplot(forplotting, aes(x = Var2, y = real_mean, fill = Var1)) + geom_bar(stat="identity", color = "black") + 
  scale_fill_manual(values = cols) +
  geom_errorbar(limits, colour = "black", width = 1,  position = dodge) +
  ylab("Relative Abundance")+ 
  labs(fill="Class")+
  xlab("")+
  theme_bw() +
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank()) 

#Phylum level - relative abundance arable vs. fungi in all and five countries
# get names of fungal phyla
names_fun <- names(sort(table(otu_taxa[,"Phylum"]), decr=T))

## Preparation of matrix with relative abundance by Phylum
y <- NULL
otunames <- rownames(OTU_RA)
for (i in names_fun){
  x <- array(colSums(OTU_RA[rownames(otu_taxa)[which(otu_taxa$Phylum == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(names_fun)
colnames(y) <- paste(colnames(OTU_RA))
Phylum_mat_fun <- as.data.frame(y)

Phylum_mat_fun$mean <- rowMeans(Phylum_mat_fun)
Phylum_mat_fun <- Phylum_mat_fun[order(Phylum_mat_fun$mean,decreasing = TRUE),]

## select Phylumes with MEAN abundances top 9
abundant_Phylum_fun <- rownames(Phylum_mat_fun)[1:10]
Phylum_mat_fun <- Phylum_mat_fun[abundant_Phylum_fun,]
others <- 100-colSums(Phylum_mat_fun)
Phylum_mat_fun["p__",]<- Phylum_mat_fun["p__",] + others
row.names(Phylum_mat_fun)[which(abundant_Phylum_fun == "p__")] <- "others"
Phylum_mat_fun <- Phylum_mat_fun[c(c(1:10)[-which(abundant_Phylum_fun == "p__")],which(abundant_Phylum_fun == "p__")),]

#calculate mean and s.e. for arable vs. grassland and for each country
# overall
Phylum_mat_fun$overall_mean <- rowMeans(Phylum_mat_fun)
Phylum_mat_fun$overall_se <- apply(Phylum_mat_fun,1, se)
Phylum_mat_fun$overall_begin <- Phylum_mat_fun$overall_mean - Phylum_mat_fun$overall_se
Phylum_mat_fun$overall_end <- Phylum_mat_fun$overall_mean + Phylum_mat_fun$overall_se

## simplify treatment to arable and grassland
design$Treatment[design$Treatment == "Conventional" | design$Treatment == "Organic"] <- "Arable"

## Arable
Phylum_mat_fun$Arable_mean <- rowMeans(Phylum_mat_fun[,row.names(design[design$Treatment == "Arable",])])
Phylum_mat_fun$Arable_se <- apply(Phylum_mat_fun[,row.names(design[design$Treatment == "Arable",])],1,se)
Phylum_mat_fun$Arable_begin <- Phylum_mat_fun$Arable_mean - Phylum_mat_fun$Arable_se
Phylum_mat_fun$Arable_end <- Phylum_mat_fun$Arable_mean + Phylum_mat_fun$Arable_se


## Grassland
Phylum_mat_fun$Grassland_mean <- rowMeans(Phylum_mat_fun[,row.names(design[design$Treatment == "Grassland",])])
Phylum_mat_fun$Grassland_se <- apply(Phylum_mat_fun[,row.names(design[design$Treatment == "Grassland",])],1,se)
Phylum_mat_fun$Grassland_begin <- Phylum_mat_fun$Grassland_mean - Phylum_mat_fun$Grassland_se
Phylum_mat_fun$Grassland_end <- Phylum_mat_fun$Grassland_mean + Phylum_mat_fun$Grassland_se


## sp
Phylum_mat_fun$sp_mean <- rowMeans(Phylum_mat_fun[,row.names(design[design$Country == "Spain",])])
Phylum_mat_fun$sp_se <- apply(Phylum_mat_fun[,row.names(design[design$Country == "Spain",])],1, se)
Phylum_mat_fun$sp_begin <- Phylum_mat_fun$sp_mean - Phylum_mat_fun$sp_se
Phylum_mat_fun$sp_end <- Phylum_mat_fun$sp_mean + Phylum_mat_fun$sp_se


## fr
Phylum_mat_fun$fr_mean <- rowMeans(Phylum_mat_fun[,row.names(design[design$Country == "France",])])
Phylum_mat_fun$fr_se <- apply(Phylum_mat_fun[,row.names(design[design$Country == "France",])],1,se)
Phylum_mat_fun$fr_begin <- Phylum_mat_fun$fr_mean - Phylum_mat_fun$fr_se
Phylum_mat_fun$fr_end <- Phylum_mat_fun$fr_mean + Phylum_mat_fun$fr_se


## ch
Phylum_mat_fun$ch_mean <- rowMeans(Phylum_mat_fun[,row.names(design[design$Country == "Switzerland",])])
Phylum_mat_fun$ch_se <- apply(Phylum_mat_fun[,row.names(design[design$Country == "Switzerland",])],1,se)
Phylum_mat_fun$ch_begin <- Phylum_mat_fun$ch_mean - Phylum_mat_fun$ch_se
Phylum_mat_fun$ch_end <- Phylum_mat_fun$ch_mean + Phylum_mat_fun$ch_se


## ge
Phylum_mat_fun$ge_mean <- rowMeans(Phylum_mat_fun[,row.names(design[design$Country == "Germany",])])
Phylum_mat_fun$ge_se <- apply(Phylum_mat_fun[,row.names(design[design$Country == "Germany",])],1,se)
Phylum_mat_fun$ge_begin <- Phylum_mat_fun$ge_mean - Phylum_mat_fun$ge_se
Phylum_mat_fun$ge_end <- Phylum_mat_fun$ge_mean + Phylum_mat_fun$ge_se

## sw
Phylum_mat_fun$sw_mean <- rowMeans(Phylum_mat_fun[,row.names(design[design$Country == "Sweden",])])
Phylum_mat_fun$sw_se <- apply(Phylum_mat_fun[,row.names(design[design$Country == "Sweden",])],1, se)
Phylum_mat_fun$sw_begin <- Phylum_mat_fun$sw_mean - Phylum_mat_fun$sw_se
Phylum_mat_fun$sw_end <- Phylum_mat_fun$sw_mean + Phylum_mat_fun$sw_se


## data transformation
forplotting <- data.frame(Var1 = rep(row.names(Phylum_mat_fun),times=8),
                          Var2 = rep(c("overall","sp","fr","ch","ge","sw","Arable","Grassland"),each=length(Phylum_mat_fun$ge_begin)),
                          mean=c(Phylum_mat_fun$overall_mean,Phylum_mat_fun$sp_mean,Phylum_mat_fun$fr_mean,Phylum_mat_fun$ch_mean,Phylum_mat_fun$ge_mean,Phylum_mat_fun$sw_mean,Phylum_mat_fun$Arable_mean,Phylum_mat_fun$Grassland_mean),
                          begin=c(Phylum_mat_fun$overall_begin, Phylum_mat_fun$sp_begin,Phylum_mat_fun$fr_begin,Phylum_mat_fun$ch_begin,Phylum_mat_fun$ge_begin,Phylum_mat_fun$sw_begin,Phylum_mat_fun$Arable_begin,Phylum_mat_fun$Grassland_begin),
                          end=c(Phylum_mat_fun$overall_end, Phylum_mat_fun$sp_end,Phylum_mat_fun$fr_end,Phylum_mat_fun$ch_end,Phylum_mat_fun$ge_end,Phylum_mat_fun$sw_end,Phylum_mat_fun$Arable_end,Phylum_mat_fun$Grassland_end))

mean_fun <- forplotting$mean
l <- length(Phylum_mat_fun$overall_mean)
for (i in 2:l)
{
  for(j in 0:7)
  {
    forplotting[i+j*l,"begin"]= forplotting[i+j*l,"begin"] + forplotting[i-1+j*l,"mean"]
    forplotting[i+j*l,"end"]=forplotting[i+j*l,"end"] + forplotting[i-1+j*l,"mean"]
    forplotting[i+j*l,"mean"]=forplotting[i+j*l,"mean"] + forplotting[i-1+j*l,"mean"]
  }
}

forplotting <- cbind(forplotting,real_mean = mean_fun)
forplotting$Var1 <- factor(forplotting$Var1, levels = rev(forplotting$Var1[1:length(Phylum_mat_fun$ge_begin)]))
forplotting$Var2<-factor(forplotting$Var2, levels = unique(forplotting$Var2))


##### Plot barchart
dodge <- position_dodge(width = 0.9)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

cols <- c("grey50",col_vector[c(1:7,9,22)])
limits <- aes(ymax = begin , ymin = end)
ggplot(forplotting, aes(x = Var2, y = real_mean, fill = Var1)) + geom_bar(stat="identity", color = "black") + 
  scale_fill_manual(values = cols) +
  geom_errorbar(limits, colour = "black", width = 1,  position = dodge) +
  ylab("Relative Abundance")+ 
  labs(fill="Phylum")+
  xlab("")+
  theme_bw() +
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank()) 

#********************************************Beta Diversity**********************************************************
rm(list=ls()[-match(c("OTU","otu_taxa","design"),ls())])
otu_norm <- sqrt(OTU)
physeq_norm <- phyloseq(otu_table(otu_norm, taxa_are_rows=T),
                        tax_table(as.matrix(otu_taxa)),
                        sample_data(design))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


#PCoA for all
col_vector[1:5] <-  col_vector[c(14,18,20,22,26)]
pcoa_norm <- ordinate(physeq_norm,"PCoA","bray")
pcoa_all <- plot_ordination(physeq_norm, pcoa_norm, type="sites", color="Country", shape="Treatment")
pcoa_all <- pcoa_all+
  geom_point(size=1.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
  scale_shape_manual(name="Treatment",values=c(15,16))+
  scale_color_manual(name="Country",values= col_vector[c(2,4,1,5,3)])+
  xlab(paste("PCo 1", paste("(",round(pcoa_norm$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
  ylab(paste("PCo 2", paste("(",round(pcoa_norm$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  guides(shape=guide_legend(nrow=2,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("PCoA analysis")
pcoa_all

pcoa_all <- plot_ordination(physeq_norm, pcoa_norm, type="sites", color="Treatment")
pcoa_all <- pcoa_all+
  geom_point(size=1.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
  scale_color_manual(name="Treatment",values= col_vector[1:2])+
  xlab(paste("PCo 1", paste("(",round(pcoa_norm$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
  ylab(paste("PCo 2", paste("(",round(pcoa_norm$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  guides(shape=guide_legend(nrow=2,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("PCoA analysis")
pcoa_all



pcoa_all <- plot_ordination(physeq_norm, pcoa_norm, type="sites", color="Lat")
pcoa_all <- pcoa_all+
  geom_point(size=1.5)+
  scale_color_continuous(low="blue", high="goldenrod2")+
  theme_bw()+
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
  xlab(paste("PCo 1", paste("(",round(pcoa_norm$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
  ylab(paste("PCo 2", paste("(",round(pcoa_norm$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(shape=guide_legend(nrow=2,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("PCoA analysis")
pcoa_all


#PERMANOVA for all
## Perform PERMANVOA testing for depth and type effects on soil community
dis <- phyloseq::distance(physeq_norm, method = "bray")
#dis <- vegdist(t(otu_table(physeq_norm)),method="bray")
sampledf <- data.frame(sample_data(physeq_norm))
paov <- adonis(dis ~ Treatment * Country, data=sampledf, permutations=999)
paov

## Perform pairwise comparision
pairwise.perm.manova(dis, sampledf$Treatment, nperm=999)
pairwise.perm.manova(dis, sampledf$Country, nperm=999)

## Perform BETADISP test of multivariate dispersions
bdi <- betadisper(dis,sampledf$Treatment,type = "centroid")
permutest(bdi,pairwise=T,permutations=how(nperm=999))

bdi <- betadisper(dis,sampledf$Country,type = "centroid")
permutest(bdi,pairwise=T,permutations=how(nperm=999))

cap_Treatment <- ordinate(physeq_norm,"CAP","bray", ~Treatment)
source("variance_functions.R")
var_tab_Treatment <- variability_table(cap_Treatment)

set.seed(2222)  
permutest(cap_Treatment, permutations=how(nperm=9999))

cap_Treatment_plot <- plot_ordination(physeq_norm,cap_Treatment,color="Treatment",shape="Treatment")+
  geom_point(size=2.5)+
  scale_shape_manual(values=c(15,16,17,18))+
  scale_color_manual(values=col_vector[1:4])+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle(bquote(paste(.(round(var_tab_Treatment["constrained","proportion"],2)*100),"% of total variation (CI=",
                       .(round(cca_ci(cap_Treatment)[1],2)*100),"%, ",
                       .(round(cca_ci(cap_Treatment)[2],2)*100),"%); ", "p < 0.001" )))+
  stat_ellipse(aes(group = Treatment),type = "t")+
  theme(plot.title = element_text(size = 10, face = "bold"))

## Perform CAP analysis constrained to cropping system on bulk soil bacteria community
cap_Country <- ordinate(physeq_norm,"CAP","bray", ~Country)

## Calculate proporition of variation table
source("variance_functions.R")
var_tab_Country <- variability_table(cap_Country)

set.seed(3333)  
permutest(cap_Country, permutations=how(nperm=9999))

cap_Country_plot <- plot_ordination(physeq_norm,cap_Country,color="Country",shape="Country")+
  geom_point(size=2.5)+
  scale_shape_manual(values=c(15,16,17,18,19))+
  scale_color_manual(values=col_vector[1:5])+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle(bquote(paste(.(round(var_tab_Country["constrained","proportion"],2)*100),"% of total variation (CI=",
                       .(round(cca_ci(cap_Country)[1],2)*100),"%, ",
                       .(round(cca_ci(cap_Country)[2],2)*100),"%); ", "p < 0.001" )))+
  stat_ellipse(aes(group = Country),Country = "t")+
  theme(plot.title = element_text(size = 10, face = "bold"))

grid.arrange(cap_Treatment_plot, cap_Country_plot, nrow=2, ncol=1)

#Mantel test geographical location and community composition
all_site <- row.names(design)
arable_site <- row.names(design[design$Treatment == "Arable",])
grassland_site <- row.names(design[design$Treatment == "Grassland",])
commdist<-vegdist(t(otu_norm[,all_site]), method="bray")
geodist<-vegdist(design[all_site,c("Lat","Long")], method="euclidean")
mantel(commdist, geodist, method="pearson", permutations=999)
plot(geodist,commdist,pch=16,cex=0.5,col="black",bty="l",xlab="Spatial distance",ylab="Species composition dissimilarity")
commdist<-vegdist(t(otu_norm[,arable_site]), method="bray")
geodist<-vegdist(design[arable_site,c("Lat","Long")], method="euclidean")
mantel(commdist, geodist, method="pearson", permutations=999)
plot(geodist,commdist,pch=16,cex=0.5,col="black",bty="l",xlab="Spatial distance",ylab="Species composition dissimilarity")
commdist<-vegdist(t(otu_norm[,grassland_site]), method="bray")
geodist<-vegdist(design[grassland_site,c("Lat","Long")], method="euclidean")
mantel(commdist, geodist, method="pearson", permutations=999)
plot(geodist,commdist,pch=16,cex=0.5,col="black",bty="l",xlab="Spatial distance",ylab="Species composition dissimilarity")

r <- raster::getData("worldclim",var="bio",res=10)
coords <- data.frame(design$Long,design$Lat)
points <- SpatialPoints(coords, proj4string =r@crs)
climatedf <- as.data.frame(extract(r,points))
row.names(climatedf)<- row.names(design)

## climate and community composition
commdist<-vegdist(t(otu_norm[,all_site]), method="bray")
climdist<-vegdist(climatedf[all_site,c(1,3,8,12)], method="euclidean")
mantel(commdist, climdist, method="pearson", permutations=999)
plot(climdist,commdist,pch=16,cex=0.5,col="black",bty="l",xlab="Climate distance",ylab="Species composition dissimilarity")

## arable
commdist<-vegdist(t(otu_norm[,arable_site]), method="bray")
climdist<-vegdist(climatedf[arable_site,c(1,3,8,12)], method="euclidean")
mantel(commdist, climdist, method="pearson", permutations=999)
plot(climdist,commdist,pch=16,cex=0.5,col="black",bty="l",xlab="Climate distance",ylab="Species composition dissimilarity")

## grassland
commdist<-vegdist(t(otu_norm[,grassland_site]), method="bray")
climdist<-vegdist(climatedf[grassland_site,c(1,3,8,12)], method="euclidean")
mantel(commdist, climdist, method="pearson", permutations=999)
plot(climdist,commdist,pch=16,cex=0.5,col="black",bty="l",xlab="Climate distance",ylab="Species composition dissimilarity")

#Mantel test of soil properties and community composition 
SoilProperties <- c("C_tot_.","N_tot_.","P_tot_mg","P_Olsen","pH","Organic_C_.","C_microbial","N_microbial","basal_respiration","Sand","Silt","Clay","moisture","CEC")

soildf <- design[,SoilProperties]
all_site <- row.names(na.omit(soildf[all_site,]))
arable_site <- row.names(na.omit(soildf[arable_site,]))
grassland_site <- row.names(na.omit(soildf[grassland_site,]))

commdist<-vegdist(t(otu_norm[,all_site]), method="bray")
soildist<-vegdist(soildf[all_site,], method="euclidean")
mantel(commdist, soildist, method="pearson", permutations=999)
plot(soildist,commdist,pch=16,cex=0.5,col="black",bty="l",xlab="Soil properties distance",ylab="Species composition dissimilarity")

## arable
commdist<-vegdist(t(otu_norm[,arable_site]), method="bray")
soildist<-vegdist(soildf[arable_site,], method="euclidean")
mantel(commdist, soildist, method="pearson", permutations=999)
plot(soildist,commdist,pch=16,cex=0.5,col="black",bty="l",xlab="Soil properties distance",ylab="Species composition dissimilarity")

## grassland
commdist<-vegdist(t(otu_norm[,grassland_site]), method="bray")
soildist<-vegdist(soildf[grassland_site,], method="euclidean")
mantel(commdist, soildist, method="pearson", permutations=999)
plot(soildist,commdist,pch=16,cex=0.5,col="black",bty="l",xlab="Soil properties distance",ylab="Species composition dissimilarity")

#Soil properties and geographical locations
## all
geodist <- vegdist(design[all_site,c("Lat","Long")], method="euclidean")
soildist<-vegdist(soildf[all_site,], method="euclidean")
mantel(geodist, soildist, method="pearson", permutations=999)
plot(geodist,soildist,pch=16,cex=0.5,col="black",bty="l",xlab="Geographical distance",ylab="Soil properties dissimilarity")

## arable
geodist <- vegdist(design[arable_site,c("Lat","Long")], method="euclidean")
soildist<-vegdist(soildf[arable_site,], method="euclidean")
mantel(geodist, soildist, method="pearson", permutations=999)
plot(geodist,soildist,pch=16,cex=0.5,col="black",bty="l",xlab="Geographical distance",ylab="Soil properties dissimilarity")

## grassland
geodist <- vegdist(design[grassland_site,c("Lat","Long")], method="euclidean")
soildist<-vegdist(soildf[grassland_site,], method="euclidean")
mantel(geodist, soildist, method="pearson", permutations=999)
plot(geodist,soildist,pch=16,cex=0.5,col="black",bty="l",xlab="Geographical distance",ylab="Soil properties dissimilarity")

#Soil properties and climate
```{r}
## all
soildist<-vegdist(soildf[all_site,], method="euclidean")
climdist<-vegdist(climatedf[all_site,c(1,3,8,12)], method="euclidean")
mantel(climdist, soildist, method="pearson", permutations=999)
plot(climdist, soildist,pch=16,cex=0.5,col="black",bty="l",xlab="Climate distance",ylab="Soil properties dissimilarity")

## arable
soildist<-vegdist(soildf[arable_site,], method="euclidean")
climdist<-vegdist(climatedf[arable_site,c(1,3,8,12)], method="euclidean")
mantel(climdist, soildist, method="pearson", permutations=999)
plot(climdist, soildist,pch=16,cex=0.5,col="black",bty="l",xlab="Climate distance",ylab="Soil properties dissimilarity")

## grassland
soildist<-vegdist(soildf[grassland_site,], method="euclidean")
climdist<-vegdist(climatedf[grassland_site,c(1,3,8,12)], method="euclidean")
mantel(climdist, soildist, method="pearson", permutations=999)
plot(climdist, soildist,pch=16,cex=0.5,col="black",bty="l",xlab="Climate distance",ylab="Soil properties dissimilarity")

#geographical location and climate
```{r}
## all
geodist <- vegdist(design[all_site,c("Lat","Long")], method="euclidean")
climdist<-vegdist(climatedf[all_site,c(1,3,8,12)], method="euclidean")
mantel(geodist, climdist, method="pearson", permutations=999)
plot(geodist, climdist,pch=16,cex=0.5,col="black",bty="l",xlab="Climate distance",ylab="Soil properties dissimilarity")

## arable
geodist <- vegdist(design[arable_site,c("Lat","Long")], method="euclidean")
climdist<-vegdist(climatedf[arable_site,c(1,3,8,12)], method="euclidean")
mantel(geodist, climdist, method="pearson", permutations=999)
plot(geodist, climdist,pch=16,cex=0.5,col="black",bty="l",xlab="Climate distance",ylab="Soil properties dissimilarity")

## grassland
geodist <- vegdist(design[grassland_site,c("Lat","Long")], method="euclidean")
climdist<-vegdist(climatedf[grassland_site,c(1,3,8,12)], method="euclidean")
mantel(geodist, climdist, method="pearson", permutations=999)
plot(geodist, climdist,pch=16,cex=0.5,col="black",bty="l",xlab="Climate distance",ylab="Soil properties dissimilarity")

#Partial Mantel test  
## All
geodist <- vegdist(design[all_site,c("Lat","Long")], method="euclidean")
climdist<-vegdist(climatedf[all_site,c(1,2,4,7,8,12,15,19)], method="euclidean")
soildist<-vegdist(soildf[all_site,], method="euclidean")
commdist<-vegdist(t(otu_norm[,all_site]), method="bray")

mantel.partial(commdist,climdist, geodist, method="pearson", permutations=999)
mantel.partial(commdist,soildist, geodist, method="pearson", permutations=999)

## Arable
geodist <- vegdist(design[arable_site,c("Lat","Long")], method="euclidean")
climdist<-vegdist(climatedf[arable_site,c(1,2,4,7,8,12,15,19)], method="euclidean")
soildist<-vegdist(soildf[arable_site,], method="euclidean")
commdist<-vegdist(t(otu_norm[,arable_site]), method="bray")

mantel.partial(commdist,climdist, geodist, method="pearson", permutations=999)
mantel.partial(commdist,soildist, geodist, method="pearson", permutations=999)

## Grassland
geodist <- vegdist(design[grassland_site,c("Lat","Long")], method="euclidean")
climdist<-vegdist(climatedf[grassland_site,c(1,2,4,7,8,12,15,19)], method="euclidean")
soildist<-vegdist(soildf[grassland_site,], method="euclidean")
commdist<-vegdist(t(otu_norm[,grassland_site]), method="bray")

mantel.partial(commdist,climdist, geodist, method="pearson", permutations=999)
mantel.partial(commdist,soildist, geodist, method="pearson", permutations=999)

envdf <- cbind(climatedf[all_site,c(1,2,4,7,8,12,15,19)],soildf[all_site,])

## all
geodist <- vegdist(design[all_site,c("Lat","Long")], method="euclidean")
envdist<-vegdist(envdf[all_site,], method="euclidean")
commdist<-vegdist(t(otu_norm[,all_site]), method="bray")
mantel.partial(commdist,envdist, geodist, method="pearson", permutations=999)
mantel.partial(commdist, geodist,envdist, method="pearson", permutations=999)
## arable
geodist <- vegdist(design[arable_site,c("Lat","Long")], method="euclidean")
envdist<-vegdist(envdf[arable_site,], method="euclidean")
commdist<-vegdist(t(otu_norm[,arable_site]), method="bray")
mantel.partial(commdist,envdist, geodist, method="pearson", permutations=999)
mantel.partial(commdist, geodist,envdist, method="pearson", permutations=999)

## grassland
geodist <- vegdist(design[grassland_site,c("Lat","Long")], method="euclidean")
envdist<-vegdist(envdf[grassland_site,], method="euclidean")
commdist<-vegdist(t(otu_norm[,grassland_site]), method="bray")
mantel.partial(commdist,envdist, geodist, method="pearson", permutations=999)
mantel.partial(commdist, geodist,envdist, method="pearson", permutations=999)

alldf <- cbind(envdf[all_site,],design[all_site,c("Lat","Long")])

## all
alldist <- vegdist(alldf[all_site,], method="euclidean")
commdist<-vegdist(t(otu_norm[,all_site]), method="bray")
mantel(commdist,alldist, method="pearson", permutations=999)

## arable
alldist <- vegdist(alldf[arable_site,], method="euclidean")
commdist<-vegdist(t(otu_norm[,arable_site]), method="bray")
mantel(commdist,alldist, method="pearson", permutations=999)

## grassland
alldist <- vegdist(alldf[grassland_site,], method="euclidean")
commdist<-vegdist(t(otu_norm[,grassland_site]), method="bray")
mantel(commdist,alldist, method="pearson", permutations=999)

geoclimdf <- cbind(design[all_site,c("Lat","Long")],climatedf[all_site,c(1,3,8,12)])

## all
geoclimdist <- vegdist(geoclimdf[all_site,], method="euclidean")
commdist<-vegdist(t(otu_norm[,all_site]), method="bray")
mantel(commdist,geoclimdist, method="pearson", permutations=999)

## arable
geoclimdist <- vegdist(geoclimdf[arable_site,], method="euclidean")
commdist<-vegdist(t(otu_norm[,arable_site]), method="bray")
mantel(commdist,geoclimdist, method="pearson", permutations=999)

## grassland
geoclimdist <- vegdist(geoclimdf[grassland_site,], method="euclidean")
commdist<-vegdist(t(otu_norm[,grassland_site]), method="bray")
mantel(commdist,geoclimdist, method="pearson", permutations=999)


#***********************************Venn Diagram*****************************************************************************
rm(list=ls()[-match(c("OTU","otu_taxa","design"),ls())])
physeq_table<- phyloseq(otu_table(OTU, taxa_are_rows=T),
                        tax_table(as.matrix(otu_taxa)),
                        sample_data(design))

sp.arable <- rowSums(OTU[,row.names(design[design$Treatment == "Arable" & design$Country == "Spain",])])
sp.arable <- names(sp.arable[sp.arable != 0])

sp.grassland <- rowSums(OTU[,row.names(design[design$Treatment == "Grassland" & design$Country == "Spain",])])
sp.grassland <- names(sp.grassland[sp.grassland != 0])

fr.arable <- rowSums(OTU[,row.names(design[design$Treatment == "Arable" & design$Country == "France",])])
fr.arable <- names(fr.arable[fr.arable != 0])

fr.grassland <- rowSums(OTU[,row.names(design[design$Treatment == "Grassland" & design$Country == "France",])])
fr.grassland <- names(fr.grassland[fr.grassland != 0])

ch.arable <- rowSums(OTU[,row.names(design[design$Treatment == "Arable" & design$Country == "Switzerland",])])
ch.arable <- names(ch.arable[ch.arable != 0])

ch.grassland <- rowSums(OTU[,row.names(design[design$Treatment == "Grassland" & design$Country == "Switzerland",])])
ch.grassland <- names(ch.grassland[ch.grassland != 0])

ge.arable <- rowSums(OTU[,row.names(design[design$Treatment == "Arable" & design$Country == "Germany",])])
ge.arable <- names(ge.arable[ge.arable != 0])

ge.grassland <- rowSums(OTU[,row.names(design[design$Treatment == "Grassland" & design$Country == "Germany",])])
ge.grassland <- names(ge.grassland[ge.grassland != 0])

sw.arable <- rowSums(OTU[,row.names(design[design$Treatment == "Arable" & design$Country == "Sweden",])])
sw.arable <- names(sw.arable[sw.arable != 0])

sw.grassland <- rowSums(OTU[,row.names(design[design$Treatment == "Grassland" & design$Country == "Sweden",])])
sw.grassland <- names(sw.grassland[sw.grassland != 0])

venn(list(Spain=sp.arable,France=fr.arable,Switzerland=ch.arable,Germany=ge.arable,Sweden=sw.arable), zcolor = "style")
venn(list(Spain=sp.grassland,France=fr.grassland,Switzerland=ch.grassland,Germany=ge.grassland,Sweden=sw.grassland), zcolor = "style")


#***********************************relative abundance*****************************************************************************
otu_taxa <- otu_taxa_pre
OTU<- OTU_pre
design <- design_pre
rm(list=ls()[-match(c("OTU","otu_taxa","design"),ls())])

OTU_RA <- t(t(OTU)/colSums(OTU)) * 100
OTU_norm <- sqrt(OTU)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

presence <- c(rep(0,length(OTU_RA[,1])))
for (i in 1:length(OTU_RA[,1]))
{
  presence[i] <- sum(OTU_RA[i,] != 0)
}
Arable <- row.names(design[design$Treatment == "Arable",])
Grassland <- row.names(design[design$Treatment == "Grassland",])
plot(presence, log(rowSums(OTU_norm)/presence))

## arable
OTU_RA_arable <- OTU_RA[,Arable]
OTU_norm_arable <- OTU_norm[,Arable]

presence_arable <- c(rep(0,length(OTU_RA[,1])))
for (i in 1:length(OTU_RA[,1]))
{
  presence_arable[i] <- sum(OTU_RA[i,Arable] != 0)
}

plot(presence_arable, log(rowSums(OTU_norm_arable)/presence_arable))

## grassland
OTU_RA_grassland <- OTU_RA[,Grassland]
OTU_norm_grassland <- OTU_norm[,Grassland]

presence_grassland <- c(rep(0,length(OTU_RA_grassland[,1])))
for (i in 1:length(OTU_RA[,1]))
{
  presence_grassland[i] <- sum(OTU_RA[i,Grassland] != 0)
}

plot(presence_grassland, log(rowSums(OTU_RA_grassland)/presence_grassland))

ggplot()+
  geom_point(aes(x = presence_grassland, y = log(rowMeans(OTU_RA_grassland))), color = "blue")+
  geom_point(aes(x = presence_arable, y = log(rowMeans(OTU_RA_arable))), color = "red")

## randomly selecting 61 arable sites
set.seed(4444)
s <- sample(c(1:156),61)
OTU_RA_arable_61 <- OTU_RA_arable[,s]
OTU_norm_arable_61 <- OTU_norm_arable[,s]
presence_arable <- c(rep(0,length(OTU_RA_arable_61[,1])))
for (i in 1:length(OTU_RA_arable_61[,1]))
{
  presence_arable[i] <- sum(OTU_RA_arable_61[i,] != 0)
}

OTU_norm_grassland <- OTU_norm_grassland[presence_grassland!=0,]
presence_grassland <- presence_grassland[presence_grassland!=0]
OTU_norm_arable_61 <- OTU_norm_arable_61[presence_arable!=0,]
presence_arable <- presence_arable[presence_arable!=0]


ggplot()+
  geom_point(aes(x = presence_grassland, y = log(rowSums(OTU_norm_grassland)/61)), color = "black", size =2,alpha = 0.5)+
  geom_point(aes(x = presence_arable, y = log(rowSums(OTU_norm_arable_61)/61)), color = col_vector[3], size =2,alpha = 0.5)+
  ylab("log (Normalized OTU Abundance)")+
  xlab("Number of Sites Occupied")

## check significance
df1 <- data.frame(presence = presence_grassland, abundance = rowSums(OTU_norm_grassland)/61,type = "grassland")
df2 <- data.frame(presence = presence_arable, abundance = rowSums(OTU_norm_arable_61)/61, type = "arable")
df <- rbind(df1, df2)


names_fun <- names(sort(table(otu_taxa[,"Class"]), decr=T))

## Preparation of matrix with relative abundance by CLASS
y <- NULL
otunames <- rownames(OTU_RA)
for (i in names_fun){
  x <- array(colSums(OTU_RA[rownames(otu_taxa)[which(otu_taxa$Class == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(names_fun)
colnames(y) <- paste(colnames(OTU_RA))
CLASS_mat_fun <- as.data.frame(y)

CLASS_mat_fun$mean <- rowMeans(CLASS_mat_fun)
CLASS_mat_fun <- CLASS_mat_fun[order(CLASS_mat_fun$mean,decreasing = TRUE),]

plot(x= log(rowMeans(CLASS_mat_fun[,Grassland])), y =(rowMeans(CLASS_mat_fun[,Arable])-rowMeans(CLASS_mat_fun[,Grassland]))/rowMeans(CLASS_mat_fun[,Grassland]))

## select classes with MEAN abundances top 9
abundant_class_fun <- rownames(CLASS_mat_fun)[1:10]
CLASS_mat_fun <- CLASS_mat_fun[abundant_class_fun,]
others <- 100-colSums(CLASS_mat_fun)
CLASS_mat_fun["c__",]<- CLASS_mat_fun["c__",] + others
row.names(CLASS_mat_fun)[which(abundant_class_fun == "c__")] <- "others"
CLASS_mat_fun <- CLASS_mat_fun[c(c(1:10)[-which(abundant_class_fun == "c__")],which(abundant_class_fun == "c__")),]

OTU_special_present <- as.data.frame(OTU_RA)
OTU_special_present$present_grassland <- rowMeans(OTU_RA[,Grassland])
OTU_special_present$present_arable <- rowMeans(OTU_RA[,Arable])
OTU_special_present$special <- OTU_special_present$present_grassland * OTU_special_present$present_arable
OTU_special_present <- OTU_special_present[OTU_special_present$special == 0,]

length(OTU_special_present$present_grassland[OTU_special_present$present_grassland!=0])
length(OTU_special_present$present_arable[OTU_special_present$present_arable!=0])

OTU_special_arable <- data.frame(Class = otu_taxa[row.names(OTU_special_present[OTU_special_present$present_arable!=0,]),]$Class, RA = OTU_special_present$present_arable[OTU_special_present$present_arable!=0])

OTU_special_grassland <- data.frame(Class = otu_taxa[row.names(OTU_special_present[OTU_special_present$present_grassland!=0,]),]$Class, RA = OTU_special_present$present_grassland[OTU_special_present$present_grassland!=0])



cols <- c("grey50",col_vector[c(1:7,9,22)])

class_special_arable <- aggregate(OTU_special_arable$RA, by = list(as.character(OTU_special_arable$Class)), sum)

class_special_grassland <- aggregate(OTU_special_grassland$RA, by = list(as.character(OTU_special_grassland$Class)), sum)

class_special_arable_top10 <- class_special_arable[order(class_special_arable$x, decreasing = T),][1:10,]
class_special_arable_top10$type <- "Arable"
class_special_arable_top10$Group.1 <- factor(class_special_arable_top10$Group.1, level = c("c__",rev(class_special_arable_top10$Group.1[-1])) )
#class_special_arable_top10$color <- cols

ggplot(class_special_arable_top10, aes(fill=Group.1, y=x, x=type)) + 
  geom_bar( stat="identity") + 
  scale_fill_manual(values = cols)

class_special_grassland_top10 <- class_special_grassland[order(class_special_grassland$x, decreasing = T),][1:10,]
class_special_grassland_top10$type <- "Grassland"

class_special_grassland_top10$Group.1 <- factor(class_special_grassland_top10$Group.1, level = c("c__",rev(class_special_grassland_top10$Group.1[-2])) )
#class_special_grassland_top10$color <- cols[c(2,1,3:10)]



ggplot(class_special_grassland_top10, aes(fill=Group.1, y=x, x=type)) + 
  geom_bar( stat="identity") + 
  scale_fill_manual(values = cols)

#check special OTUs at order level
OTU_special_arable <- data.frame(Order = otu_taxa[row.names(OTU_special_present[OTU_special_present$present_arable!=0,]),]$Order, RA = OTU_special_present$present_arable[OTU_special_present$present_arable!=0])
OTU_special_grassland <- data.frame(Order = otu_taxa[row.names(OTU_special_present[OTU_special_present$present_grassland!=0,]),]$Order, RA = OTU_special_present$present_grassland[OTU_special_present$present_grassland!=0])



cols <- c("grey50",col_vector[c(1:7,9,22)])

Order_special_arable <- aggregate(OTU_special_arable$RA, by = list(as.character(OTU_special_arable$Order)), sum)

Order_special_grassland <- aggregate(OTU_special_grassland$RA, by = list(as.character(OTU_special_grassland$Order)), sum)

Order_special_arable_top10 <- Order_special_arable[order(Order_special_arable$x, decreasing = T),][1:10,]
Order_special_arable_top10$type <- "Arable"
Order_special_arable_top10$Group.1 <- factor(Order_special_arable_top10$Group.1, level = c("o__",rev(Order_special_arable_top10$Group.1[-1])) )
#Order_special_arable_top10$color <- cols

ggplot(Order_special_arable_top10, aes(fill=Group.1, y=x, x=type)) + 
  geom_bar( stat="identity") + 
  scale_fill_manual(values = cols)

Order_special_grassland_top10 <- Order_special_grassland[order(Order_special_grassland$x, decreasing = T),][1:10,]
Order_special_grassland_top10$type <- "Grassland"

Order_special_grassland_top10$Group.1 <- factor(Order_special_grassland_top10$Group.1, level = c("o__",rev(Order_special_grassland_top10$Group.1[-1])) )

ggplot(Order_special_grassland_top10, aes(fill=Group.1, y=x, x=type)) + 
  geom_bar( stat="identity") + 
  scale_fill_manual(values = cols)


#check significant differences
fun_mat <- as.data.frame(t(CLASS_mat_fun))
fun_mat$Country <- design$Country[match(row.names(fun_mat),row.names(design))]
fun_mat$Treatment <- design$Treatment[match(row.names(fun_mat),row.names(design))]
fun_mat[,"Country"] <- as.factor(fun_mat[,"Country"] )
fun_mat[,"Treatment"] <- as.factor(fun_mat[,"Treatment"] )
fun_mat$Treatment[fun_mat$Treatment == "Conventional" | fun_mat$Treatment == "Organic"] <- "Arable"

for (i in 1:9)
{
  print(colnames(fun_mat)[i])
  
  print(kruskal.test(fun_mat[,i] ~ fun_mat[,"Country"] ))
  print(kruskal.test(fun_mat[,i] ~ fun_mat[,"Treatment"] ))
  print(posthoc.kruskal.nemenyi.test(fun_mat[,i],fun_mat[,"Country"]))
  print(posthoc.kruskal.nemenyi.test(fun_mat[,i],fun_mat[,"Treatment"]))
  #print(PostHocTest(aov1, method = "duncan") )
}  

#calculate mean and s.e. for arable vs. grassland and for each country

# overall
CLASS_mat_fun$overall_mean <- rowMeans(CLASS_mat_fun)
CLASS_mat_fun$overall_se <- apply(CLASS_mat_fun,1, se)
CLASS_mat_fun$overall_begin <- CLASS_mat_fun$overall_mean - CLASS_mat_fun$overall_se
CLASS_mat_fun$overall_end <- CLASS_mat_fun$overall_mean + CLASS_mat_fun$overall_se

## simplify treatment to arable and grassland
design$Treatment[design$Treatment == "Conventional" | design$Treatment == "Organic"] <- "Arable"

## Arable
CLASS_mat_fun$Arable_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable",])])
CLASS_mat_fun$Arable_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable",])],1,se)
CLASS_mat_fun$Arable_begin <- CLASS_mat_fun$Arable_mean - CLASS_mat_fun$Arable_se
CLASS_mat_fun$Arable_end <- CLASS_mat_fun$Arable_mean + CLASS_mat_fun$Arable_se


## Grassland
CLASS_mat_fun$Grassland_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland",])])
CLASS_mat_fun$Grassland_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland",])],1,se)
CLASS_mat_fun$Grassland_begin <- CLASS_mat_fun$Grassland_mean - CLASS_mat_fun$Grassland_se
CLASS_mat_fun$Grassland_end <- CLASS_mat_fun$Grassland_mean + CLASS_mat_fun$Grassland_se


## sp
CLASS_mat_fun$sp_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Country == "Spain",])])
CLASS_mat_fun$sp_se <- apply(CLASS_mat_fun[,row.names(design[design$Country == "Spain",])],1, se)
CLASS_mat_fun$sp_begin <- CLASS_mat_fun$sp_mean - CLASS_mat_fun$sp_se
CLASS_mat_fun$sp_end <- CLASS_mat_fun$sp_mean + CLASS_mat_fun$sp_se


## fr
CLASS_mat_fun$fr_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Country == "France",])])
CLASS_mat_fun$fr_se <- apply(CLASS_mat_fun[,row.names(design[design$Country == "France",])],1,se)
CLASS_mat_fun$fr_begin <- CLASS_mat_fun$fr_mean - CLASS_mat_fun$fr_se
CLASS_mat_fun$fr_end <- CLASS_mat_fun$fr_mean + CLASS_mat_fun$fr_se


## ch
CLASS_mat_fun$ch_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Country == "Switzerland",])])
CLASS_mat_fun$ch_se <- apply(CLASS_mat_fun[,row.names(design[design$Country == "Switzerland",])],1,se)
CLASS_mat_fun$ch_begin <- CLASS_mat_fun$ch_mean - CLASS_mat_fun$ch_se
CLASS_mat_fun$ch_end <- CLASS_mat_fun$ch_mean + CLASS_mat_fun$ch_se


## ge
CLASS_mat_fun$ge_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Country == "Germany",])])
CLASS_mat_fun$ge_se <- apply(CLASS_mat_fun[,row.names(design[design$Country == "Germany",])],1,se)
CLASS_mat_fun$ge_begin <- CLASS_mat_fun$ge_mean - CLASS_mat_fun$ge_se
CLASS_mat_fun$ge_end <- CLASS_mat_fun$ge_mean + CLASS_mat_fun$ge_se

## sw
CLASS_mat_fun$sw_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Country == "Sweden",])])
CLASS_mat_fun$sw_se <- apply(CLASS_mat_fun[,row.names(design[design$Country == "Sweden",])],1, se)
CLASS_mat_fun$sw_begin <- CLASS_mat_fun$sw_mean - CLASS_mat_fun$sw_se
CLASS_mat_fun$sw_end <- CLASS_mat_fun$sw_mean + CLASS_mat_fun$sw_se


## data transformation
forplotting <- data.frame(Var1 = rep(row.names(CLASS_mat_fun),times=8),
                          Var2 = rep(c("overall","sp","fr","ch","ge","sw","Arable","Grassland"),each=length(CLASS_mat_fun$ge_begin)),
                          mean=c(CLASS_mat_fun$overall_mean,CLASS_mat_fun$sp_mean,CLASS_mat_fun$fr_mean,CLASS_mat_fun$ch_mean,CLASS_mat_fun$ge_mean,CLASS_mat_fun$sw_mean,CLASS_mat_fun$Arable_mean,CLASS_mat_fun$Grassland_mean),
                          begin=c(CLASS_mat_fun$overall_begin, CLASS_mat_fun$sp_begin,CLASS_mat_fun$fr_begin,CLASS_mat_fun$ch_begin,CLASS_mat_fun$ge_begin,CLASS_mat_fun$sw_begin,CLASS_mat_fun$Arable_begin,CLASS_mat_fun$Grassland_begin),
                          end=c(CLASS_mat_fun$overall_end, CLASS_mat_fun$sp_end,CLASS_mat_fun$fr_end,CLASS_mat_fun$ch_end,CLASS_mat_fun$ge_end,CLASS_mat_fun$sw_end,CLASS_mat_fun$Arable_end,CLASS_mat_fun$Grassland_end))

mean_fun <- forplotting$mean
l <- length(CLASS_mat_fun$overall_mean)
for (i in 2:l)
{
  for(j in 0:7)
  {
    forplotting[i+j*l,"begin"]= forplotting[i+j*l,"begin"] + forplotting[i-1+j*l,"mean"]
    forplotting[i+j*l,"end"]=forplotting[i+j*l,"end"] + forplotting[i-1+j*l,"mean"]
    forplotting[i+j*l,"mean"]=forplotting[i+j*l,"mean"] + forplotting[i-1+j*l,"mean"]
  }
}

forplotting <- cbind(forplotting,real_mean = mean_fun)
forplotting$Var1 <- factor(forplotting$Var1, levels = rev(forplotting$Var1[1:length(CLASS_mat_fun$ge_begin)]))
forplotting$Var2<-factor(forplotting$Var2, levels = unique(forplotting$Var2))

##### Plot barchart
dodge <- position_dodge(width = 0.9)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

cols <- c("grey50",col_vector[c(1:7,9,22)])
limits <- aes(ymax = begin , ymin = end)
ggplot(forplotting, aes(x = Var2, y = real_mean, fill = Var1)) + geom_bar(stat="identity", color = "black" ) + 
  scale_fill_manual(values = cols) +
  geom_errorbar(limits, colour = "black", width = 1,  position = dodge) +
  ylab("Relative Abundance")+ 
  labs(fill="Class")+
  xlab("")+
  theme_bw() +
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank()) 


#*********************************************Functional Groups****************************************************************
combined_functions <- paste(as.data.frame(table(otu_taxa$Guild))[,1], sep = '-', collapse = '-')
combined_functions <- unique(unlist(strsplit(combined_functions, "-")))[-1]

functions_df <- data.frame(matrix(ncol = length(combined_functions), nrow = length(row.names(design))))
colnames(functions_df) <- combined_functions

## OTU_functions -- OTU table with only OTUs that are assigned to guilds
OTU_functions <- t(OTU[row.names(otu_taxa[otu_taxa$Guild != "-",]),])
row.names(functions_df)<- row.names(OTU_functions)

## combine all functions that exist in one site
for (i in 1:204)
{
  ## OTUs that exist in this site 
  OTUs<- names(which(OTU_functions[i,]!= 0 ))
  ## combine functions of all OTUs that exist in this site
  functions <- otu_taxa[OTUs,"Guild"]
  functions_df$functions[i] <- as.character(paste(functions, sep = '-', collapse = '-'))
}

for (i in 1:33)
  for (j in 1:204)
  {
    if (grepl(combined_functions[i], functions_df[j,34],fix = T))  functions_df[j,i] = 1
    else functions_df[j,i] = 0 
    
  }

functions_df$funcrich <- rowSums(functions_df[,1:33])
functions_df["arable",] <- NA
functions_df["grassland",] <- NA

Arable <- row.names(design[design$Treatment=="Arable",])
Grassland <- row.names(design[design$Treatment=="Grassland",])

functions_df["arable",1:33] <- colSums( functions_df[Arable,1:33]) * 100/length(Arable)
functions_df["grassland",1:33] <- colSums( functions_df[Grassland,1:33]) * 100/length(Grassland)

x <- as.data.frame(t(functions_df[c("arable","grassland"),1:33]))

row.names(x[x$arable>0 & x$grassland > 0,])
row.names(x[x$arable>0 & x$grassland == 0,])
row.names(x[x$arable == 0 & x$grassland > 0,])

##plot functions
plotfunctions <- data.frame(functions = combined_functions, RA_arable = NA, RA_grassland = NA)
OTU_RA <- t(t(OTU)/colSums(OTU)) * 100
RA_landuse <- data.frame(RA_grassland = rowMeans(OTU_RA[,row.names(design[design$Treatment == "Grassland",])]), RA_arable = rowMeans(OTU_RA[,row.names(design[design$Treatment == "Arable",])]))
for (i in 1:33)
{
  k = grepl(combined_functions[i], otu_taxa[row.names(RA_landuse),]$Guild)
  plotfunctions[i,"RA_arable"] = sum(RA_landuse[k,"RA_arable"])
  plotfunctions[i,"RA_grassland"] = sum(RA_landuse[k,"RA_grassland"])
}

## grab RA of functions for each site
functions_persite <- data.frame(matrix(, nrow=33, ncol=length(row.names(design))+1))
colnames(functions_persite)[2:(length(row.names(design))+1)] <- row.names(design)
colnames(functions_persite)[1] <- "functions"
functions_persite$functions <- plotfunctions$functions

for (i in 1:33)
{
  k = grepl(combined_functions[i], otu_taxa[row.names(RA_landuse),]$Guild)
  for (j in 1:length(row.names(design)))
    functions_persite[i,j+1] <- sum(OTU_RA[k,colnames(functions_persite)[j+1]])
  
}
CLASS_mat_fun <- functions_persite[order(rowMeans(functions_persite[,2:34]), decreasing = T),]
CLASS_mat_fun <-CLASS_mat_fun[-match("Wood Saprotrop",CLASS_mat_fun$functions),]

CLASS_mat_fun <- CLASS_mat_fun[1:10,]

row.names(CLASS_mat_fun) <- CLASS_mat_fun$functions
CLASS_mat_fun <- CLASS_mat_fun[,-1]

## Arable Spain
CLASS_mat_fun$Arable_Spain_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable" & design$Country == "Spain",])])
CLASS_mat_fun$Arable_Spain_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable"& design$Country == "Spain",])],1,se)
CLASS_mat_fun$Arable_Spain_begin <- CLASS_mat_fun$Arable_Spain_mean - CLASS_mat_fun$Arable_Spain_se
CLASS_mat_fun$Arable_Spain_end <- CLASS_mat_fun$Arable_Spain_mean + CLASS_mat_fun$Arable_Spain_se


## Grassland Spain
CLASS_mat_fun$Grassland_Spain_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland" & design$Country == "Spain",])])
CLASS_mat_fun$Grassland_Spain_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland"& design$Country == "Spain",])],1,se)
CLASS_mat_fun$Grassland_Spain_begin <- CLASS_mat_fun$Grassland_Spain_mean - CLASS_mat_fun$Grassland_Spain_se
CLASS_mat_fun$Grassland_Spain_end <- CLASS_mat_fun$Grassland_Spain_mean + CLASS_mat_fun$Grassland_Spain_se

## Arable France
CLASS_mat_fun$Arable_France_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable" & design$Country == "France",])])
CLASS_mat_fun$Arable_France_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable"& design$Country == "France",])],1,se)
CLASS_mat_fun$Arable_France_begin <- CLASS_mat_fun$Arable_France_mean - CLASS_mat_fun$Arable_France_se
CLASS_mat_fun$Arable_France_end <- CLASS_mat_fun$Arable_France_mean + CLASS_mat_fun$Arable_France_se

## Grassland France
CLASS_mat_fun$Grassland_France_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland" & design$Country == "France",])])
CLASS_mat_fun$Grassland_France_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland"& design$Country == "France",])],1,se)
CLASS_mat_fun$Grassland_France_begin <- CLASS_mat_fun$Grassland_France_mean - CLASS_mat_fun$Grassland_France_se
CLASS_mat_fun$Grassland_France_end <- CLASS_mat_fun$Grassland_France_mean + CLASS_mat_fun$Grassland_France_se

## Arable Switzerland
CLASS_mat_fun$Arable_Switzerland_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable" & design$Country == "Switzerland",])])
CLASS_mat_fun$Arable_Switzerland_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable"& design$Country == "Switzerland",])],1,se)
CLASS_mat_fun$Arable_Switzerland_begin <- CLASS_mat_fun$Arable_Switzerland_mean - CLASS_mat_fun$Arable_Switzerland_se
CLASS_mat_fun$Arable_Switzerland_end <- CLASS_mat_fun$Arable_Switzerland_mean + CLASS_mat_fun$Arable_Switzerland_se


## Grassland Switzerland
CLASS_mat_fun$Grassland_Switzerland_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland" & design$Country == "Switzerland",])])
CLASS_mat_fun$Grassland_Switzerland_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland"& design$Country == "Switzerland",])],1,se)
CLASS_mat_fun$Grassland_Switzerland_begin <- CLASS_mat_fun$Grassland_Switzerland_mean - CLASS_mat_fun$Grassland_Switzerland_se
CLASS_mat_fun$Grassland_Switzerland_end <- CLASS_mat_fun$Grassland_Switzerland_mean + CLASS_mat_fun$Grassland_Switzerland_se

## Arable Germany
CLASS_mat_fun$Arable_Germany_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable" & design$Country == "Germany",])])
CLASS_mat_fun$Arable_Germany_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable"& design$Country == "Germany",])],1,se)
CLASS_mat_fun$Arable_Germany_begin <- CLASS_mat_fun$Arable_Germany_mean - CLASS_mat_fun$Arable_Germany_se
CLASS_mat_fun$Arable_Germany_end <- CLASS_mat_fun$Arable_Germany_mean + CLASS_mat_fun$Arable_Germany_se


## Grassland Germany
CLASS_mat_fun$Grassland_Germany_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland" & design$Country == "Germany",])])
CLASS_mat_fun$Grassland_Germany_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland"& design$Country == "Germany",])],1,se)
CLASS_mat_fun$Grassland_Germany_begin <- CLASS_mat_fun$Grassland_Germany_mean - CLASS_mat_fun$Grassland_Germany_se
CLASS_mat_fun$Grassland_Germany_end <- CLASS_mat_fun$Grassland_Germany_mean + CLASS_mat_fun$Grassland_Germany_se

## Arable Sweden
CLASS_mat_fun$Arable_Sweden_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable" & design$Country == "Sweden",])])
CLASS_mat_fun$Arable_Sweden_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Arable"& design$Country == "Sweden",])],1,se)
CLASS_mat_fun$Arable_Sweden_begin <- CLASS_mat_fun$Arable_Sweden_mean - CLASS_mat_fun$Arable_Sweden_se
CLASS_mat_fun$Arable_Sweden_end <- CLASS_mat_fun$Arable_Sweden_mean + CLASS_mat_fun$Arable_Sweden_se


## Grassland Sweden
CLASS_mat_fun$Grassland_Sweden_mean <- rowMeans(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland" & design$Country == "Sweden",])])
CLASS_mat_fun$Grassland_Sweden_se <- apply(CLASS_mat_fun[,row.names(design[design$Treatment == "Grassland"& design$Country == "Sweden",])],1,se)
CLASS_mat_fun$Grassland_Sweden_begin <- CLASS_mat_fun$Grassland_Sweden_mean - CLASS_mat_fun$Grassland_Sweden_se
CLASS_mat_fun$Grassland_Sweden_end <- CLASS_mat_fun$Grassland_Sweden_mean + CLASS_mat_fun$Grassland_Sweden_se

## data transformation
forplotting <- data.frame(Var1 = rep(row.names(CLASS_mat_fun),times=10),
                          Var2 = rep(c("Arable_Spain","Grassland_Spain" ,"Arable_France","Grassland_France", "Arable_Switzerland","Grassland_Switzerland", "Arable_Germany","Grassland_Germany", "Arable_Sweden","Grassland_Sweden"  ),each=10),
                          mean=unlist(c(CLASS_mat_fun[,c(51:60)*4+1]),use.names = F),
                          begin=unlist(c(CLASS_mat_fun[,c(51:60)*4+3]),use.names = F),
                          end=unlist(c(CLASS_mat_fun[,c(51:60)*4+4])),use.names = F)

mean_fun <- forplotting$mean
l <- 10
for (i in 2:l)
{
  for(j in 0:9)
  {
    forplotting[i+j*l,"begin"]= forplotting[i+j*l,"begin"] + forplotting[i-1+j*l,"mean"]
    forplotting[i+j*l,"end"]=forplotting[i+j*l,"end"] + forplotting[i-1+j*l,"mean"]
    forplotting[i+j*l,"mean"]=forplotting[i+j*l,"mean"] + forplotting[i-1+j*l,"mean"]
  }
}

forplotting <- cbind(forplotting,real_mean = mean_fun)
forplotting$Var1 <- factor(forplotting$Var1, levels = rev(forplotting$Var1[1:10]))
forplotting$Var2<-factor(forplotting$Var2, levels = unique(forplotting$Var2))

##### Plot barchart
dodge <- position_dodge(width = 0.9)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

cols <- c("grey50",col_vector[c(1:7,9,22)])
limits <- aes(ymax = begin , ymin = end)
ggplot(forplotting, aes(x = Var2, y = real_mean, fill = Var1)) + geom_bar(stat="identity", color = "black" ) + 
  scale_fill_manual(values = cols) +
  geom_errorbar(limits, colour = "black", width = 1,  position = dodge) +
  ylab("Relative Abundance")+ 
  labs(fill="Class")+
  xlab("")+
  theme_bw() +
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank()) 

k = grepl("polen", otu_taxa[row.names(RA_landuse),]$Guild)
plotfunctions[31,"RA_arable"] = sum(RA_landuse[k,"RA_arable"])
plotfunctions[31,"RA_grassland"] = sum(RA_landuse[k,"RA_grassland"])

plotfunctions <- plotfunctions[-c(1,4,26,31,33),]

plotfunctions_1 <- plotfunctions[which (plotfunctions$functions %in%                                           row.names(x[x$arable>80,])),]
plotfunctions_1 <- plotfunctions_1[order(plotfunctions_1$RA_arable, decreasing = T),]
plotfunctions_1_melt <- melt(plotfunctions_1)
plotfunctions_1_melt$functions <- factor(plotfunctions_1_melt$functions, level = plotfunctions_1$functions)

ggplot(data = plotfunctions_1_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = functions, y = value, fill = variable))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plotfunctions_2 <- plotfunctions[which (plotfunctions$functions %in%                                           row.names(x[x$arable<80,])),]
plotfunctions_2 <- plotfunctions_2[order(plotfunctions_2$RA_arable, decreasing = T),]
plotfunctions_2_melt <- melt(plotfunctions_2)
plotfunctions_2_melt$functions <- factor(plotfunctions_2_melt$functions, level = plotfunctions_2$functions)

ggplot(data = plotfunctions_2_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = functions, y = value, fill = variable))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## plot commonness of functions
toplotcommonness <- t(functions_df[c(205:206),1:33])
toplotcommonness_melt <- melt(toplotcommonness)
ggplot(data = toplotcommonness_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toplotcommonness_1 <- toplotcommonness[as.character(plotfunctions_1$functions),]
toplotcommonness_1_melt <- melt(toplotcommonness_1)
toplotcommonness_1_melt$Var1 <- factor(toplotcommonness_1_melt$Var1, level = as.character(plotfunctions_1$functions))
ggplot(data = toplotcommonness_1_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toplotcommonness_2 <- toplotcommonness[as.character(plotfunctions_2$functions),]
toplotcommonness_2_melt <- melt(toplotcommonness_2)
toplotcommonness_2_melt$Var1 <- factor(toplotcommonness_2_melt$Var1, level = as.character(plotfunctions_2$functions))
ggplot(data = toplotcommonness_2_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Spain
plotfunctions <- data.frame(functions = combined_functions, RA_arable = NA, RA_grassland = NA)
OTU_RA <- t(t(OTU)/colSums(OTU)) * 100
RA_landuse <- data.frame(RA_grassland = rowMeans(OTU_RA[,row.names(design[design$Treatment == "Grassland" & design$Country == "Spain",])]), RA_arable = rowMeans(OTU_RA[,row.names(design[design$Treatment == "Arable" & design$Country == "Spain",])]))
for (i in 1:33)
{
  k = grepl(combined_functions[i], otu_taxa[row.names(RA_landuse),]$Guild)
  plotfunctions[i,"RA_arable"] = sum(RA_landuse[k,"RA_arable"])
  plotfunctions[i,"RA_grassland"] = sum(RA_landuse[k,"RA_grassland"])
}

## something wrong with Plant pathogenic (?) on polen, manually change
k = grepl("polen", otu_taxa[row.names(RA_landuse),]$Guild)
plotfunctions[31,"RA_arable"] = sum(RA_landuse[k,"RA_arable"])
plotfunctions[31,"RA_grassland"] = sum(RA_landuse[k,"RA_grassland"])

plotfunctions <- plotfunctions[-c(1,4,26,31,33),]

plotfunctions_1 <- plotfunctions[which (plotfunctions$functions %in%                                           row.names(x[x$arable>80,])),]
plotfunctions_1 <- plotfunctions_1[order(plotfunctions_1$RA_arable, decreasing = T),]
plotfunctions_1_melt <- melt(plotfunctions_1)
plotfunctions_1_melt$functions <- factor(plotfunctions_1_melt$functions, level = plotfunctions_1$functions)

ggplot(data = plotfunctions_1_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = functions, y = value, fill = variable))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plotfunctions_2 <- plotfunctions[which (plotfunctions$functions %in%                                           row.names(x[x$arable<80,])),]
plotfunctions_2 <- plotfunctions_2[order(plotfunctions_2$RA_arable, decreasing = T),]
plotfunctions_2_melt <- melt(plotfunctions_2)
plotfunctions_2_melt$functions <- factor(plotfunctions_2_melt$functions, level = plotfunctions_2$functions)

ggplot(data = plotfunctions_2_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = functions, y = value, fill = variable))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## plot commonness of functions
toplotcommonness <- t(functions_df[c(205:206),1:33])
toplotcommonness_melt <- melt(toplotcommonness)
ggplot(data = toplotcommonness_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toplotcommonness_1 <- toplotcommonness[as.character(plotfunctions_1$functions),]
toplotcommonness_1_melt <- melt(toplotcommonness_1)
toplotcommonness_1_melt$Var1 <- factor(toplotcommonness_1_melt$Var1, level = as.character(plotfunctions_1$functions))
ggplot(data = toplotcommonness_1_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toplotcommonness_2 <- toplotcommonness[as.character(plotfunctions_2$functions),]
toplotcommonness_2_melt <- melt(toplotcommonness_2)
toplotcommonness_2_melt$Var1 <- factor(toplotcommonness_2_melt$Var1, level = as.character(plotfunctions_2$functions))
ggplot(data = toplotcommonness_2_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#France
plotfunctions <- data.frame(functions = combined_functions, RA_arable = NA, RA_grassland = NA)
OTU_RA <- t(t(OTU)/colSums(OTU)) * 100
RA_landuse <- data.frame(RA_grassland = rowMeans(OTU_RA[,row.names(design[design$Treatment == "Grassland" & design$Country == "France",])]), RA_arable = rowMeans(OTU_RA[,row.names(design[design$Treatment == "Arable" & design$Country == "France",])]))
for (i in 1:33)
{
  k = grepl(combined_functions[i], otu_taxa[row.names(RA_landuse),]$Guild)
  plotfunctions[i,"RA_arable"] = sum(RA_landuse[k,"RA_arable"])
  plotfunctions[i,"RA_grassland"] = sum(RA_landuse[k,"RA_grassland"])
}

## something wrong with Plant pathogenic (?) on polen, manually change
k = grepl("polen", otu_taxa[row.names(RA_landuse),]$Guild)
plotfunctions[31,"RA_arable"] = sum(RA_landuse[k,"RA_arable"])
plotfunctions[31,"RA_grassland"] = sum(RA_landuse[k,"RA_grassland"])

plotfunctions <- plotfunctions[-c(1,4,26,31,33),]

plotfunctions_1 <- plotfunctions[which (plotfunctions$functions %in%                                           row.names(x[x$arable>80,])),]
plotfunctions_1 <- plotfunctions_1[order(plotfunctions_1$RA_arable, decreasing = T),]
plotfunctions_1_melt <- melt(plotfunctions_1)
plotfunctions_1_melt$functions <- factor(plotfunctions_1_melt$functions, level = plotfunctions_1$functions)

ggplot(data = plotfunctions_1_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = functions, y = value, fill = variable))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plotfunctions_2 <- plotfunctions[which (plotfunctions$functions %in%                                           row.names(x[x$arable<80,])),]
plotfunctions_2 <- plotfunctions_2[order(plotfunctions_2$RA_arable, decreasing = T),]
plotfunctions_2_melt <- melt(plotfunctions_2)
plotfunctions_2_melt$functions <- factor(plotfunctions_2_melt$functions, level = plotfunctions_2$functions)

ggplot(data = plotfunctions_2_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = functions, y = value, fill = variable))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## plot commonness of functions
toplotcommonness <- t(functions_df[c(205:206),1:33])
toplotcommonness_melt <- melt(toplotcommonness)
ggplot(data = toplotcommonness_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toplotcommonness_1 <- toplotcommonness[as.character(plotfunctions_1$functions),]
toplotcommonness_1_melt <- melt(toplotcommonness_1)
toplotcommonness_1_melt$Var1 <- factor(toplotcommonness_1_melt$Var1, level = as.character(plotfunctions_1$functions))
ggplot(data = toplotcommonness_1_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toplotcommonness_2 <- toplotcommonness[as.character(plotfunctions_2$functions),]
toplotcommonness_2_melt <- melt(toplotcommonness_2)
toplotcommonness_2_melt$Var1 <- factor(toplotcommonness_2_melt$Var1, level = as.character(plotfunctions_2$functions))
ggplot(data = toplotcommonness_2_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Switzerland
plotfunctions <- data.frame(functions = combined_functions, RA_arable = NA, RA_grassland = NA)
OTU_RA <- t(t(OTU)/colSums(OTU)) * 100
RA_landuse <- data.frame(RA_grassland = rowMeans(OTU_RA[,row.names(design[design$Treatment == "Grassland" & design$Country == "Switzerland",])]), RA_arable = rowMeans(OTU_RA[,row.names(design[design$Treatment == "Arable" & design$Country == "Switzerland",])]))
for (i in 1:33)
{
  k = grepl(combined_functions[i], otu_taxa[row.names(RA_landuse),]$Guild)
  plotfunctions[i,"RA_arable"] = sum(RA_landuse[k,"RA_arable"])
  plotfunctions[i,"RA_grassland"] = sum(RA_landuse[k,"RA_grassland"])
}

## something wrong with Plant pathogenic (?) on polen, manually change
k = grepl("polen", otu_taxa[row.names(RA_landuse),]$Guild)
plotfunctions[31,"RA_arable"] = sum(RA_landuse[k,"RA_arable"])
plotfunctions[31,"RA_grassland"] = sum(RA_landuse[k,"RA_grassland"])

plotfunctions <- plotfunctions[-c(1,4,26,31,33),]

plotfunctions_1 <- plotfunctions[which (plotfunctions$functions %in%                                           row.names(x[x$arable>80,])),]
plotfunctions_1 <- plotfunctions_1[order(plotfunctions_1$RA_arable, decreasing = T),]
plotfunctions_1_melt <- melt(plotfunctions_1)
plotfunctions_1_melt$functions <- factor(plotfunctions_1_melt$functions, level = plotfunctions_1$functions)

ggplot(data = plotfunctions_1_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = functions, y = value, fill = variable))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plotfunctions_2 <- plotfunctions[which (plotfunctions$functions %in%                                           row.names(x[x$arable<80,])),]
plotfunctions_2 <- plotfunctions_2[order(plotfunctions_2$RA_arable, decreasing = T),]
plotfunctions_2_melt <- melt(plotfunctions_2)
plotfunctions_2_melt$functions <- factor(plotfunctions_2_melt$functions, level = plotfunctions_2$functions)

ggplot(data = plotfunctions_2_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = functions, y = value, fill = variable))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## plot commonness of functions
toplotcommonness <- t(functions_df[c(205:206),1:33])
toplotcommonness_melt <- melt(toplotcommonness)
ggplot(data = toplotcommonness_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toplotcommonness_1 <- toplotcommonness[as.character(plotfunctions_1$functions),]
toplotcommonness_1_melt <- melt(toplotcommonness_1)
toplotcommonness_1_melt$Var1 <- factor(toplotcommonness_1_melt$Var1, level = as.character(plotfunctions_1$functions))
ggplot(data = toplotcommonness_1_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toplotcommonness_2 <- toplotcommonness[as.character(plotfunctions_2$functions),]
toplotcommonness_2_melt <- melt(toplotcommonness_2)
toplotcommonness_2_melt$Var1 <- factor(toplotcommonness_2_melt$Var1, level = as.character(plotfunctions_2$functions))
ggplot(data = toplotcommonness_2_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Germany
plotfunctions <- data.frame(functions = combined_functions, RA_arable = NA, RA_grassland = NA)
OTU_RA <- t(t(OTU)/colSums(OTU)) * 100
RA_landuse <- data.frame(RA_grassland = rowMeans(OTU_RA[,row.names(design[design$Treatment == "Grassland" & design$Country == "Germany",])]), RA_arable = rowMeans(OTU_RA[,row.names(design[design$Treatment == "Arable" & design$Country == "Germany",])]))
for (i in 1:33)
{
  k = grepl(combined_functions[i], otu_taxa[row.names(RA_landuse),]$Guild)
  plotfunctions[i,"RA_arable"] = sum(RA_landuse[k,"RA_arable"])
  plotfunctions[i,"RA_grassland"] = sum(RA_landuse[k,"RA_grassland"])
}

## something wrong with Plant pathogenic (?) on polen, manually change
k = grepl("polen", otu_taxa[row.names(RA_landuse),]$Guild)
plotfunctions[31,"RA_arable"] = sum(RA_landuse[k,"RA_arable"])
plotfunctions[31,"RA_grassland"] = sum(RA_landuse[k,"RA_grassland"])

plotfunctions <- plotfunctions[-c(1,4,26,31,33),]

plotfunctions_1 <- plotfunctions[which (plotfunctions$functions %in%                                           row.names(x[x$arable>80,])),]
plotfunctions_1 <- plotfunctions_1[order(plotfunctions_1$RA_arable, decreasing = T),]
plotfunctions_1_melt <- melt(plotfunctions_1)
plotfunctions_1_melt$functions <- factor(plotfunctions_1_melt$functions, level = plotfunctions_1$functions)

ggplot(data = plotfunctions_1_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = functions, y = value, fill = variable))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plotfunctions_2 <- plotfunctions[which (plotfunctions$functions %in%                                           row.names(x[x$arable<80,])),]
plotfunctions_2 <- plotfunctions_2[order(plotfunctions_2$RA_arable, decreasing = T),]
plotfunctions_2_melt <- melt(plotfunctions_2)
plotfunctions_2_melt$functions <- factor(plotfunctions_2_melt$functions, level = plotfunctions_2$functions)

ggplot(data = plotfunctions_2_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = functions, y = value, fill = variable))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## plot commonness of functions
toplotcommonness <- t(functions_df[c(205:206),1:33])
toplotcommonness_melt <- melt(toplotcommonness)
ggplot(data = toplotcommonness_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toplotcommonness_1 <- toplotcommonness[as.character(plotfunctions_1$functions),]
toplotcommonness_1_melt <- melt(toplotcommonness_1)
toplotcommonness_1_melt$Var1 <- factor(toplotcommonness_1_melt$Var1, level = as.character(plotfunctions_1$functions))
ggplot(data = toplotcommonness_1_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toplotcommonness_2 <- toplotcommonness[as.character(plotfunctions_2$functions),]
toplotcommonness_2_melt <- melt(toplotcommonness_2)
toplotcommonness_2_melt$Var1 <- factor(toplotcommonness_2_melt$Var1, level = as.character(plotfunctions_2$functions))
ggplot(data = toplotcommonness_2_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Sweden
plotfunctions <- data.frame(functions = combined_functions, RA_arable = NA, RA_grassland = NA)
OTU_RA <- t(t(OTU)/colSums(OTU)) * 100
RA_landuse <- data.frame(RA_grassland = rowMeans(OTU_RA[,row.names(design[design$Treatment == "Grassland" & design$Country == "Sweden",])]), RA_arable = rowMeans(OTU_RA[,row.names(design[design$Treatment == "Arable" & design$Country == "Sweden",])]))
for (i in 1:33)
{
  k = grepl(combined_functions[i], otu_taxa[row.names(RA_landuse),]$Guild)
  plotfunctions[i,"RA_arable"] = sum(RA_landuse[k,"RA_arable"])
  plotfunctions[i,"RA_grassland"] = sum(RA_landuse[k,"RA_grassland"])
}

## something wrong with Plant pathogenic (?) on polen, manually change
k = grepl("polen", otu_taxa[row.names(RA_landuse),]$Guild)
plotfunctions[31,"RA_arable"] = sum(RA_landuse[k,"RA_arable"])
plotfunctions[31,"RA_grassland"] = sum(RA_landuse[k,"RA_grassland"])

plotfunctions <- plotfunctions[-c(1,4,26,31,33),]

plotfunctions_1 <- plotfunctions[which (plotfunctions$functions %in%                                           row.names(x[x$arable>80,])),]
plotfunctions_1 <- plotfunctions_1[order(plotfunctions_1$RA_arable, decreasing = T),]
plotfunctions_1_melt <- melt(plotfunctions_1)
plotfunctions_1_melt$functions <- factor(plotfunctions_1_melt$functions, level = plotfunctions_1$functions)

ggplot(data = plotfunctions_1_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = functions, y = value, fill = variable))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plotfunctions_2 <- plotfunctions[which (plotfunctions$functions %in%                                           row.names(x[x$arable<80,])),]
plotfunctions_2 <- plotfunctions_2[order(plotfunctions_2$RA_arable, decreasing = T),]
plotfunctions_2_melt <- melt(plotfunctions_2)
plotfunctions_2_melt$functions <- factor(plotfunctions_2_melt$functions, level = plotfunctions_2$functions)

ggplot(data = plotfunctions_2_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = functions, y = value, fill = variable))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## plot commonness of functions
toplotcommonness <- t(functions_df[c(205:206),1:33])
toplotcommonness_melt <- melt(toplotcommonness)
ggplot(data = toplotcommonness_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toplotcommonness_1 <- toplotcommonness[as.character(plotfunctions_1$functions),]
toplotcommonness_1_melt <- melt(toplotcommonness_1)
toplotcommonness_1_melt$Var1 <- factor(toplotcommonness_1_melt$Var1, level = as.character(plotfunctions_1$functions))
ggplot(data = toplotcommonness_1_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toplotcommonness_2 <- toplotcommonness[as.character(plotfunctions_2$functions),]
toplotcommonness_2_melt <- melt(toplotcommonness_2)
toplotcommonness_2_melt$Var1 <- factor(toplotcommonness_2_melt$Var1, level = as.character(plotfunctions_2$functions))
ggplot(data = toplotcommonness_2_melt)+
  geom_bar(stat="identity",position='dodge',aes(x = Var1, y = value, fill = Var2))+
  scale_fill_manual(values = col_vector[c(3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#functional richness and taxonomic richness
alpha_est <- read.csv("alpha_est.csv")
design$funcrich <- functions_df[row.names(design),"funcrich"]
row.names(alpha_est) <- alpha_est$Row.names
design$richness <- alpha_est[row.names(design),"Richness"]

#Taxonomic richness
fit.arable <- lm(alpha_est[row.names(design[design$Treatment == "Arable",]),"Richness"] ~ alpha_est[row.names(design[design$Treatment == "Arable",]),]$Lat)
print(paste("Adj R2 = ",signif(summary(fit.arable)$adj.r.squared, 5), "\nIntercept =",signif(fit.arable$coef[[1]],5 ),                                   " \nSlope =",signif(fit.arable$coef[[2]], 5),                                    " \nP =",signif(summary(fit.arable)$coef[2,4], 5)))
print(AICc(fit.arable))

fit.grassland <- lm(alpha_est[row.names(design[design$Treatment == "Grassland",]),"Richness"] ~ alpha_est[row.names(design[design$Treatment == "Grassland",]),]$Lat)
print(paste("Adj R2 = ",signif(summary(fit.grassland)$adj.r.squared, 5), "\nIntercept =",signif(fit.grassland$coef[[1]],5 ),                                   " \nSlope =",signif(fit.grassland$coef[[2]], 5),                                    " \nP =",signif(summary(fit.grassland)$coef[2,4], 5)))
print(AICc(fit.grassland))

fit.arable <- lm(alpha_est[row.names(design[design$Treatment == "Arable",]),"Richness"] ~ poly(alpha_est[row.names(design[design$Treatment == "Arable",]),]$Lat,2,raw = T))
print(paste("Adj R2 = ",signif(summary(fit.arable)$adj.r.squared, 5),                                    " \nP =",signif(summary(fit.arable)$coef[2,4], 5)))
print(AICc(fit.arable))

fit.grassland <- lm(alpha_est[row.names(design[design$Treatment == "Grassland",]),"Richness"] ~ poly(alpha_est[row.names(design[design$Treatment == "Grassland",]),]$Lat,2,raw = T))
print(paste("Adj R2 = ",signif(summary(fit.grassland)$adj.r.squared, 5),                                    " \nP =",signif(summary(fit.grassland)$coef[2,4], 5)))

print(AICc(fit.grassland))


ggplot()+  geom_point(aes(x = alpha_est$Lat, y = alpha_est[,"Richness"], colour = alpha_est$Treatment))+
  stat_smooth(aes(x = alpha_est$Lat, y = alpha_est[,"Richness"], colour = alpha_est$Treatment),method="lm", formula = y ~ poly(x, 2, raw =TRUE), se = F)+
  stat_smooth(aes(x = alpha_est$Lat, y = alpha_est[,"Richness"], colour = alpha_est$Treatment), method="lm", se = F, linetype = "dashed")+
  scale_color_manual(values= col_vector[c(3,1)])

mean(alpha_est[alpha_est$Treatment == "Grassland","Richness"])/mean(alpha_est[alpha_est$Treatment == "Arable","Richness"])

fit.arable <- lm(alpha_est[row.names(design[design$Treatment == "Arable",]),"Richness"] ~ alpha_est[row.names(design[design$Treatment == "Arable",]),]$Long)
print(paste("Adj R2 = ",signif(summary(fit.arable)$adj.r.squared, 5), "\nIntercept =",signif(fit.arable$coef[[1]],5 ),                                   " \nSlope =",signif(fit.arable$coef[[2]], 5),                                    " \nP =",signif(summary(fit.arable)$coef[2,4], 5)))
print(AICc(fit.arable))

fit.grassland <- lm(alpha_est[row.names(design[design$Treatment == "Grassland",]),"Richness"] ~ alpha_est[row.names(design[design$Treatment == "Grassland",]),]$Long)
print(paste("Adj R2 = ",signif(summary(fit.grassland)$adj.r.squared, 5), "\nIntercept =",signif(fit.grassland$coef[[1]],5 ),                                   " \nSlope =",signif(fit.grassland$coef[[2]], 5),                                    " \nP =",signif(summary(fit.grassland)$coef[2,4], 5)))
print(AICc(fit.grassland))

fit.arable <- lm(alpha_est[row.names(design[design$Treatment == "Arable",]),"Richness"] ~ poly(alpha_est[row.names(design[design$Treatment == "Arable",]),]$Long,2,raw = T))
print(paste("Adj R2 = ",signif(summary(fit.arable)$adj.r.squared, 5),                                    " \nP =",signif(summary(fit.arable)$coef[2,4], 5)))
print(AICc(fit.arable))

fit.grassland <- lm(alpha_est[row.names(design[design$Treatment == "Grassland",]),"Richness"] ~ poly(alpha_est[row.names(design[design$Treatment == "Grassland",]),]$Long,2,raw = T))
print(paste("Adj R2 = ",signif(summary(fit.grassland)$adj.r.squared, 5),                                    " \nP =",signif(summary(fit.grassland)$coef[2,4], 5)))

print(AICc(fit.grassland))


ggplot()+  geom_point(aes(x = alpha_est$Long, y = alpha_est[,"Richness"], colour = alpha_est$Treatment))+
  stat_smooth(aes(x = alpha_est$Long, y = alpha_est[,"Richness"], colour = alpha_est$Treatment),method="lm", formula = y ~ poly(x, 2, raw =TRUE), se = F)+
  stat_smooth(aes(x = alpha_est$Long, y = alpha_est[,"Richness"], colour = alpha_est$Treatment), method="lm", se = F, linetype = "dashed")+
  scale_color_manual(values= col_vector[c(3,1)])

mean(alpha_est[alpha_est$Treatment == "Grassland","Richness"])/mean(alpha_est[alpha_est$Treatment == "Arable","Richness"])

#plot richness, shannon and ACE
boxplot(alpha_est[,"Richness"] ~ factor(alpha_est$Treatment) + factor(alpha_est$Country, levels = CountryNames) , col = rep(col_vector[c(3,1)],5))

kruskal.test(alpha_est[,"Richness"] ~ alpha_est[,"Country"] )
kruskal.test(alpha_est[,"Richness"] ~ alpha_est[,"Treatment"] )

boxplot(alpha_est[,"Shannon"] ~ factor(alpha_est$Treatment) + factor(alpha_est$Country, levels = CountryNames) , col = rep(col_vector[c(3,1)],5))

kruskal.test(alpha_est[,"Shannon"] ~ alpha_est[,"Country"] )
kruskal.test(alpha_est[,"Shannon"] ~ alpha_est[,"Treatment"] )


boxplot(alpha_est[,"ACE"] ~ factor(alpha_est$Treatment) + factor(alpha_est$Country, levels = CountryNames) , col = rep(col_vector[c(3,1)],5))

kruskal.test(alpha_est[,"ACE"] ~ alpha_est[,"Country"] )
kruskal.test(alpha_est[,"ACE"] ~ alpha_est[,"Treatment"] )

#matrix of all
diversityindex <- c("PD","richness","funcrich")
soilProperties <- c("C_tot_.","N_tot_.","P_tot_mg","pH","Sand","Silt","Clay","moisture","CEC")
Climate <- read.csv("climatedf.csv", row.names =1)

variablesSEM <- cbind(design[,c(diversityindex, soilProperties)], Climate[row.names(design),c("bio1","bio3","bio8","bio12")])
corrplot(as.matrix(cor(na.omit(variablesSEM), method = "spearman")), method = "number", type = "lower", diag = F)
write.csv(variablesSEM, "variablesSEM.csv")


#corrplot(as.matrix(cor(na.omit(variablesSEM[Arable,]), method = "spearman")), method = "number", type = "lower", diag = F)
corrplot(as.matrix(cor(na.omit(variablesSEM[Grassland,]), method = "spearman")), method = "number", type = "lower", diag = F)



#land use-specific soil fungi
OTU_special <- as.data.frame(OTU)
OTU_special$special_grassland <- rowMeans(OTU[,Grassland])
OTU_special$special_arable <- rowMeans(OTU[,Arable])
OTU_special$special <- OTU_special$special_grassland * OTU_special$special_arable
OTU_special <- OTU_special[OTU_special$special == 0,]


length(OTU_special$special_grassland[OTU_special$special_grassland!=0])
length(OTU_special$special_arable[OTU_special$special_arable!=0])

OTU_special_arable <- data.frame(Class = otu_taxa[row.names(OTU_special[OTU_special$special_arable!=0,]),]$Class, abundance = OTU_special$special_arable[OTU_special$special_arable!=0])

OTU_special_grassland <- data.frame(Class = otu_taxa[row.names(OTU_special[OTU_special$special_grassland!=0,]),]$Class, abundance = OTU_special$special_grassland[OTU_special$special_grassland!=0])

cols <- c("grey50",col_vector[c(1:7,9,22)])

class_special_arable <- aggregate(OTU_special_arable$abundance, by = list(as.character(OTU_special_arable$Class)), sum)

class_special_grassland <- aggregate(OTU_special_grassland$abundance, by = list(as.character(OTU_special_grassland$Class)), sum)

class_special_arable_top10 <- class_special_arable[order(class_special_arable$x, decreasing = T),][1:10,]
class_special_arable_top10$type <- "Arable"
class_special_arable_top10$Group.1 <- factor(class_special_arable_top10$Group.1, level = c("c__",rev(class_special_arable_top10$Group.1[-1])) )
#class_special_arable_top10$color <- cols

ggplot(class_special_arable_top10, aes(fill=Group.1, y=x, x=type)) + 
  geom_bar( stat="identity") + 
  scale_fill_manual(values = cols)

class_special_grassland_top10 <- class_special_grassland[order(class_special_grassland$x, decreasing = T),][1:10,]
class_special_grassland_top10$type <- "Grassland"

class_special_grassland_top10$Group.1 <- factor(class_special_grassland_top10$Group.1, level = c("c__",rev(class_special_grassland_top10$Group.1[-2])) )
#class_special_grassland_top10$color <- cols[c(2,1,3:10)]



ggplot(class_special_grassland_top10, aes(fill=Group.1, y=x, x=type)) + 
  geom_bar( stat="identity") + 
  scale_fill_manual(values = cols)

## plot special taxa for both land use types
class_special_top10 <- merge(class_special_grassland_top10, class_special_arable_top10, by = "Group.1", all = T)
class_special_top10$type.x <- "Grassland"
class_special_top10$type.y <- "Arable"
class_special_top10$x.x [is.na(class_special_top10$x.x)] <- 0
class_special_top10$x.y [is.na(class_special_top10$x.y)] <- 0
class_special_top10$mean <- rowMeans(class_special_top10[,c("x.x","x.y")])
class_special_top10 <- class_special_top10[order(class_special_top10$mean, decreasing = T),]
class_special_top10$Group.1 <- factor(as.character(class_special_top10$Group.1), c("c__",rev(as.character(class_special_top10$Group.1)[-1]) )     )


cols <- c("grey50",col_vector[c(1:7,9,22)])

ggplot() + 
  geom_bar(data = class_special_top10, aes(fill=Group.1, y=x.x, x=type.x), stat="identity") + 
  geom_bar(data = class_special_top10, aes(fill=Group.1, y=x.y, x=type.y), stat="identity") + 
  scale_fill_manual(values = c("grey50",col_vector[c(10,11,9,3,12,2,13,14,22,5,6,4)]))



## functional composition
## select 56 samples from arable land
shared_n <- rep(0,100)
grassland_n <- rep(0,100)
arable_n <- rep(0,100)


for (i in 1:100)
{
  set.seed(i)
  s <- sample(148,56)
  OTU_special_56 <- OTU[,c(Grassland, Arable[s])]
  OTU_special_56$special_grassland <- rowMeans(OTU[,Grassland])
  OTU_special_56$special_arable <- rowMeans(OTU[,Arable[s]])
  OTU_special_56$special <- OTU_special_56$special_grassland * OTU_special_56$special_arable
  shared_n[i] <- length(OTU_special_56$special[OTU_special_56$special != 0])
  OTU_special_56 <- OTU_special_56[OTU_special_56$special == 0,]
  grassland_n[i] <- length(OTU_special_56$special_grassland[OTU_special_56$special_grassland!=0])
  arable_n[i] <- length(OTU_special_56$special_arable[OTU_special_56$special_arable!=0])
  
  
}

mean(shared_n)
mean(arable_n)
mean(grassland_n)

mean(shared_n) / (mean(shared_n) + mean(arable_n) + mean(grassland_n))
mean(arable_n) / (mean(shared_n) + mean(arable_n) + mean(grassland_n))
mean(grassland_n) / (mean(shared_n) + mean(arable_n) + mean(grassland_n))


##Rank Abundance
# get names of fungal classes
names_fun <- names(sort(table(otu_taxa[,"Class"]), decr=T))

## Preparation of matrix with relative abundance by CLASS
y <- NULL
otunames <- rownames(OTU_RA)
for (i in names_fun){
  x <- array(colSums(OTU_RA[rownames(otu_taxa)[which(otu_taxa$Class == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(names_fun)
colnames(y) <- paste(colnames(OTU_RA))
CLASS_mat_fun <- as.data.frame(y)

CLASS_mat_fun$mean <- rowMeans(CLASS_mat_fun)
CLASS_mat_fun <- CLASS_mat_fun[order(CLASS_mat_fun$mean,decreasing = TRUE),]

plot(x= log(rowMeans(CLASS_mat_fun[,Grassland])), y =(rowMeans(CLASS_mat_fun[,Arable])-rowMeans(CLASS_mat_fun[,Grassland]))/rowMeans(CLASS_mat_fun[,Grassland]))

## select classes with MEAN abundances top 9
abundant_class_fun <- rownames(CLASS_mat_fun)[1:10]
CLASS_mat_fun <- CLASS_mat_fun[abundant_class_fun,]
others <- 100-colSums(CLASS_mat_fun)
CLASS_mat_fun["c__",]<- CLASS_mat_fun["c__",] + others
row.names(CLASS_mat_fun)[which(abundant_class_fun == "c__")] <- "others"
CLASS_mat_fun <- CLASS_mat_fun[c(c(1:10)[-which(abundant_class_fun == "c__")],which(abundant_class_fun == "c__")),]

OTU_special_present <- as.data.frame(OTU_RA)
OTU_special_present$present_grassland <- rowMeans(OTU_RA[,Grassland])
OTU_special_present$present_arable <- rowMeans(OTU_RA[,Arable])
OTU_special_present$special <- OTU_special_present$present_grassland * OTU_special_present$present_arable
OTU_special_present <- OTU_special_present[OTU_special_present$special == 0,]

length(OTU_special_present$present_grassland[OTU_special_present$present_grassland!=0])
length(OTU_special_present$present_arable[OTU_special_present$present_arable!=0])

OTU_special_arable <- data.frame(Class = otu_taxa[row.names(OTU_special_present[OTU_special_present$present_arable!=0,]),]$Class, RA = OTU_special_present$present_arable[OTU_special_present$present_arable!=0])

OTU_special_grassland <- data.frame(Class = otu_taxa[row.names(OTU_special_present[OTU_special_present$present_grassland!=0,]),]$Class, RA = OTU_special_present$present_grassland[OTU_special_present$present_grassland!=0])



cols <- c("grey50",col_vector[c(1:7,9,22)])

class_special_arable <- aggregate(OTU_special_arable$RA, by = list(as.character(OTU_special_arable$Class)), sum)

class_special_grassland <- aggregate(OTU_special_grassland$RA, by = list(as.character(OTU_special_grassland$Class)), sum)

class_special_arable_top10 <- class_special_arable[order(class_special_arable$x, decreasing = T),][1:10,]
class_special_arable_top10$type <- "Arable"
class_special_arable_top10$Group.1 <- factor(class_special_arable_top10$Group.1, level = c("c__",rev(class_special_arable_top10$Group.1[-1])) )
#class_special_arable_top10$color <- cols

ggplot(class_special_arable_top10, aes(fill=Group.1, y=x, x=type)) + 
  geom_bar( stat="identity") + 
  scale_fill_manual(values = cols)

class_special_grassland_top10 <- class_special_grassland[order(class_special_grassland$x, decreasing = T),][1:10,]
class_special_grassland_top10$type <- "Grassland"

class_special_grassland_top10$Group.1 <- factor(class_special_grassland_top10$Group.1, level = c("c__",rev(class_special_grassland_top10$Group.1[-2])) )
#class_special_grassland_top10$color <- cols[c(2,1,3:10)]

ggplot(class_special_grassland_top10, aes(fill=Group.1, y=x, x=type)) + 
  geom_bar( stat="identity") + 
  scale_fill_manual(values = cols)


## cumulative abundance
CUM_arable <- as.data.frame(OTU_RA[,Arable])
for (i in 1:length(colnames(CUM_arable)))
{
  CUM_arable[,i] <- CUM_arable[,i][order(CUM_arable[,i], decreasing = T)]
  
}

CUM_arable$mean <- rowMeans(CUM_arable)

## cumulative abundance
CUM_grassland <- as.data.frame(OTU_RA[,Grassland])
for (i in 1:length(colnames(CUM_grassland)))
{
  CUM_grassland[,i] <- CUM_grassland[,i][order(CUM_grassland[,i], decreasing = T)]
  
}
CUM_grassland$mean <- rowMeans(CUM_grassland)


CUM_arable <- CUM_arable[CUM_arable$mean!= 0,]
CUM_grassland <- CUM_grassland[CUM_grassland$mean!= 0,]

CUM_arable$cumulate <- CUM_arable$mean
for (i in 2:length(CUM_arable$cumulate))
  CUM_arable$cumulate[i] <- CUM_arable$mean[i]+CUM_arable$cumulate[i-1]


CUM_grassland$cumulate <- CUM_grassland$mean
for (i in 2:length(CUM_grassland$cumulate))
  CUM_grassland$cumulate[i] <- CUM_grassland$mean[i]+CUM_grassland$cumulate[i-1]


ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$cumulate)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$cumulate)), size = 1, color = col_vector[3])


CUM_trial <- CUM_grassland[1:200,]
CUM_trial$mean <- CUM_trial$mean * 100/ sum(CUM_trial$mean)
CUM_trial$cumulate <- CUM_trial$mean
for (i in 2:length(CUM_trial$cumulate))
  CUM_trial$cumulate[i] <- CUM_trial$mean[i]+CUM_trial$cumulate[i-1]

OTU_norm <- sqrt(OTU)
# select data for arable
CUM_arable <- as.data.frame(OTU_norm[,Arable])
# loop over all sites, for each site, order OTUs based on decreasing abundance
for (i in 1:length(colnames(CUM_arable)))
{
  CUM_arable[,i] <- CUM_arable[,i][order(CUM_arable[,i], decreasing = T)]
  
}

CUM_arable$mean <- rowMeans(CUM_arable)

CUM_grassland <- as.data.frame(OTU_norm[,Grassland])
for (i in 1:length(colnames(CUM_grassland)))
{
  CUM_grassland[,i] <- CUM_grassland[,i][order(CUM_grassland[,i], decreasing = T)]
  
}
CUM_grassland$mean <- rowMeans(CUM_grassland)

CUM_arable <- CUM_arable[CUM_arable$mean!= 0,]
CUM_grassland <- CUM_grassland[CUM_grassland$mean!= 0,]


CUM_arable$cumulate <- CUM_arable$mean
for (i in 2:length(CUM_arable$cumulate))
  CUM_arable$cumulate[i] <- CUM_arable$mean[i]+CUM_arable$cumulate[i-1]

CUM_grassland$cumulate <- CUM_grassland$mean
for (i in 2:length(CUM_grassland$cumulate))
  CUM_grassland$cumulate[i] <- CUM_grassland$mean[i]+CUM_grassland$cumulate[i-1]

ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$cumulate)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$cumulate)), size = 1, color = col_vector[3])


SP <- row.names(design[design$Country == "Spain",])

## cumulative abundance
CUM_arable <- as.data.frame(OTU_norm[,Arable[Arable %in% SP]])

for (i in 1:length(colnames(CUM_arable)))
{
  CUM_arable[,i] <- CUM_arable[,i][order(CUM_arable[,i], decreasing = T)]
  
}

CUM_arable$mean <- rowMeans(CUM_arable)

## cumulative abundance
CUM_grassland <- as.data.frame(OTU_norm[,Grassland[Grassland %in% SP]])
for (i in 1:length(colnames(CUM_grassland)))
{
  CUM_grassland[,i] <- CUM_grassland[,i][order(CUM_grassland[,i], decreasing = T)]
  
}
CUM_grassland$mean <- rowMeans(CUM_grassland)


CUM_arable <- CUM_arable[CUM_arable$mean!= 0,]
CUM_grassland <- CUM_grassland[CUM_grassland$mean!= 0,]

CUM_arable$cumulate <- CUM_arable$mean
for (i in 2:length(CUM_arable$cumulate))
  CUM_arable$cumulate[i] <- CUM_arable$mean[i]+CUM_arable$cumulate[i-1]


CUM_grassland$cumulate <- CUM_grassland$mean
for (i in 2:length(CUM_grassland$cumulate))
  CUM_grassland$cumulate[i] <- CUM_grassland$mean[i]+CUM_grassland$cumulate[i-1]




ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$cumulate)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$cumulate)), size = 1, color = col_vector[3])


ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$mean)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$mean)), size = 1, color = col_vector[3])


FR <- row.names(design[design$Country == "France",])

## cumulative abundance
CUM_arable <- as.data.frame(OTU_norm[,Arable[Arable %in% FR]])

for (i in 1:length(colnames(CUM_arable)))
{
  CUM_arable[,i] <- CUM_arable[,i][order(CUM_arable[,i], decreasing = T)]
  
}

CUM_arable$mean <- rowMeans(CUM_arable)

## cumulative abundance
CUM_grassland <- as.data.frame(OTU_norm[,Grassland[Grassland %in% FR]])
for (i in 1:length(colnames(CUM_grassland)))
{
  CUM_grassland[,i] <- CUM_grassland[,i][order(CUM_grassland[,i], decreasing = T)]
  
}
CUM_grassland$mean <- rowMeans(CUM_grassland)


CUM_arable <- CUM_arable[CUM_arable$mean!= 0,]
CUM_grassland <- CUM_grassland[CUM_grassland$mean!= 0,]

CUM_arable$cumulate <- CUM_arable$mean
for (i in 2:length(CUM_arable$cumulate))
  CUM_arable$cumulate[i] <- CUM_arable$mean[i]+CUM_arable$cumulate[i-1]


CUM_grassland$cumulate <- CUM_grassland$mean
for (i in 2:length(CUM_grassland$cumulate))
  CUM_grassland$cumulate[i] <- CUM_grassland$mean[i]+CUM_grassland$cumulate[i-1]


ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$cumulate)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$cumulate)), size = 1, color = col_vector[3])

plot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$mean)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$mean)), size = 1, color = col_vector[3])


CH <- row.names(design[design$Country == "Switzerland",])

## cumulative abundance
CUM_arable <- as.data.frame(OTU_norm[,Arable[Arable %in% CH]])

for (i in 1:length(colnames(CUM_arable)))
{
  CUM_arable[,i] <- CUM_arable[,i][order(CUM_arable[,i], decreasing = T)]
  
}

CUM_arable$mean <- rowMeans(CUM_arable)

## cumulative abundance
CUM_grassland <- as.data.frame(OTU_norm[,Grassland[Grassland %in% CH]])
for (i in 1:length(colnames(CUM_grassland)))
{
  CUM_grassland[,i] <- CUM_grassland[,i][order(CUM_grassland[,i], decreasing = T)]
  
}
CUM_grassland$mean <- rowMeans(CUM_grassland)


CUM_arable <- CUM_arable[CUM_arable$mean!= 0,]
CUM_grassland <- CUM_grassland[CUM_grassland$mean!= 0,]

CUM_arable$cumulate <- CUM_arable$mean
for (i in 2:length(CUM_arable$cumulate))
  CUM_arable$cumulate[i] <- CUM_arable$mean[i]+CUM_arable$cumulate[i-1]


CUM_grassland$cumulate <- CUM_grassland$mean
for (i in 2:length(CUM_grassland$cumulate))
  CUM_grassland$cumulate[i] <- CUM_grassland$mean[i]+CUM_grassland$cumulate[i-1]


ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$cumulate)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$cumulate)), size = 1, color = col_vector[3])

ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$mean)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$mean)), size = 1, color = col_vector[3])

```

#Germany
GE <- row.names(design[design$Country == "Germany",])

## cumulative abundance
CUM_arable <- as.data.frame(OTU_norm[,Arable[Arable %in% GE]])

for (i in 1:length(colnames(CUM_arable)))
{
  CUM_arable[,i] <- CUM_arable[,i][order(CUM_arable[,i], decreasing = T)]
  
}

CUM_arable$mean <- rowMeans(CUM_arable)

## cumulative abundance
CUM_grassland <- as.data.frame(OTU_norm[,Grassland[Grassland %in% GE]])
for (i in 1:length(colnames(CUM_grassland)))
{
  CUM_grassland[,i] <- CUM_grassland[,i][order(CUM_grassland[,i], decreasing = T)]
  
}
CUM_grassland$mean <- rowMeans(CUM_grassland)


CUM_arable <- CUM_arable[CUM_arable$mean!= 0,]
CUM_grassland <- CUM_grassland[CUM_grassland$mean!= 0,]

CUM_arable$cumulate <- CUM_arable$mean
for (i in 2:length(CUM_arable$cumulate))
  CUM_arable$cumulate[i] <- CUM_arable$mean[i]+CUM_arable$cumulate[i-1]


CUM_grassland$cumulate <- CUM_grassland$mean
for (i in 2:length(CUM_grassland$cumulate))
  CUM_grassland$cumulate[i] <- CUM_grassland$mean[i]+CUM_grassland$cumulate[i-1]




ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$cumulate)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$cumulate)), size = 1, color = col_vector[3])

ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$mean)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$mean)), size = 1, color = col_vector[3])


S_GE <- row.names(design[design$Country == "Germany" & design$Lat < 50,])

## cumulative abundance
CUM_arable <- as.data.frame(OTU_norm[,Arable[Arable %in% S_GE]])

for (i in 1:length(colnames(CUM_arable)))
{
  CUM_arable[,i] <- CUM_arable[,i][order(CUM_arable[,i], decreasing = T)]
  
}

CUM_arable$mean <- rowMeans(CUM_arable)

## cumulative abundance
CUM_grassland <- as.data.frame(OTU_norm[,Grassland[Grassland %in% S_GE]])
for (i in 1:length(colnames(CUM_grassland)))
{
  CUM_grassland[,i] <- CUM_grassland[,i][order(CUM_grassland[,i], decreasing = T)]
  
}
CUM_grassland$mean <- rowMeans(CUM_grassland)


CUM_arable <- CUM_arable[CUM_arable$mean!= 0,]
CUM_grassland <- CUM_grassland[CUM_grassland$mean!= 0,]

CUM_arable$cumulate <- CUM_arable$mean
for (i in 2:length(CUM_arable$cumulate))
  CUM_arable$cumulate[i] <- CUM_arable$mean[i]+CUM_arable$cumulate[i-1]


CUM_grassland$cumulate <- CUM_grassland$mean
for (i in 2:length(CUM_grassland$cumulate))
  CUM_grassland$cumulate[i] <- CUM_grassland$mean[i]+CUM_grassland$cumulate[i-1]

ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$cumulate)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$cumulate)), size = 1, color = col_vector[3])

ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$mean)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$mean)), size = 1, color = col_vector[3])


N_GE <- row.names(design[design$Country == "Germany" & design$Lat > 50,])

## cumulative abundance
CUM_arable <- as.data.frame(OTU_norm[,Arable[Arable %in% N_GE]])

for (i in 1:length(colnames(CUM_arable)))
{
  CUM_arable[,i] <- CUM_arable[,i][order(CUM_arable[,i], decreasing = T)]
  
}

CUM_arable$mean <- rowMeans(CUM_arable)

## cumulative abundance
CUM_grassland <- as.data.frame(OTU_norm[,Grassland[Grassland %in% N_GE]])
for (i in 1:length(colnames(CUM_grassland)))
{
  CUM_grassland[,i] <- CUM_grassland[,i][order(CUM_grassland[,i], decreasing = T)]
  
}
CUM_grassland$mean <- rowMeans(CUM_grassland)


CUM_arable <- CUM_arable[CUM_arable$mean!= 0,]
CUM_grassland <- CUM_grassland[CUM_grassland$mean!= 0,]

CUM_arable$cumulate <- CUM_arable$mean
for (i in 2:length(CUM_arable$cumulate))
  CUM_arable$cumulate[i] <- CUM_arable$mean[i]+CUM_arable$cumulate[i-1]


CUM_grassland$cumulate <- CUM_grassland$mean
for (i in 2:length(CUM_grassland$cumulate))
  CUM_grassland$cumulate[i] <- CUM_grassland$mean[i]+CUM_grassland$cumulate[i-1]


ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$cumulate)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$cumulate)), size = 1, color = col_vector[3])

ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$mean)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$mean)), size = 1, color = col_vector[3])

#Sweden
SW <- row.names(design[design$Country == "Sweden",])

## cumulative abundance
CUM_arable <- as.data.frame(OTU_norm[,Arable[Arable %in% SW]])

for (i in 1:length(colnames(CUM_arable)))
{
  CUM_arable[,i] <- CUM_arable[,i][order(CUM_arable[,i], decreasing = T)]
  
}

CUM_arable$mean <- rowMeans(CUM_arable)

## cumulative abundance
CUM_grassland <- as.data.frame(OTU_norm[,Grassland[Grassland %in% SW]])
for (i in 1:length(colnames(CUM_grassland)))
{
  CUM_grassland[,i] <- CUM_grassland[,i][order(CUM_grassland[,i], decreasing = T)]
  
}
CUM_grassland$mean <- rowMeans(CUM_grassland)


CUM_arable <- CUM_arable[CUM_arable$mean!= 0,]
CUM_grassland <- CUM_grassland[CUM_grassland$mean!= 0,]

CUM_arable$cumulate <- CUM_arable$mean
for (i in 2:length(CUM_arable$cumulate))
  CUM_arable$cumulate[i] <- CUM_arable$mean[i]+CUM_arable$cumulate[i-1]


CUM_grassland$cumulate <- CUM_grassland$mean
for (i in 2:length(CUM_grassland$cumulate))
  CUM_grassland$cumulate[i] <- CUM_grassland$mean[i]+CUM_grassland$cumulate[i-1]

ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$cumulate)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$cumulate)), size = 1, color = col_vector[3])

ggplot()+
  geom_point(aes(x=(c(1:length(CUM_grassland$cumulate))), y = (CUM_grassland$mean)), size = 1, color = col_vector[1])+
  geom_point(aes(x=(c(1:length(CUM_arable$cumulate))), y = (CUM_arable$mean)), size = 1, color = col_vector[3])


########## Alluvial Plots and Niche Breadth Analysis########################

library(tidyverse)
library(microbiome)
library(phyloseq)
library(metagMisc)
library(MiscMetabar)
library(EnhancedVolcano)
library(DESeq2)


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


# PCoA Endemic Calculations new analysis 08/7/2023 -------------

#### Niche Breadth Analysis
library(MicroNiche)
###### Niche Breadth
BacteriaNiche2<- as.data.frame(physeq.tab.2sites@otu_table) %>% 
  rownames_to_column(.,"OTU")

Bacteria_RareRanks2 <- as.data.frame(physeq.tab.2sites@tax_table) %>% 
  rownames_to_column(.,"OTU")

Prueba.NicheBreadth <- levins.Bn(BacteriaNiche2,2,
                                 design$Treatment)

specialist.2sites<- Prueba.NicheBreadth %>% filter(Bn<0.55 & P.adj < 0.05)

keep.specialist <- rownames(specialist.2sites)


specialist.2site.physeq<- prune_taxa(keep.specialist,physeq.tab.2sites)


##### Endemic Specialists for MI index


physeq2sites_plus_1.specialist <- transform_sample_counts(specialist.2site.physeq, 
                                                          function(x) x+1 )

TSS_Physeq.2.specialist <- metagMisc::phyloseq_standardize_otu_abundance(physeq2sites_plus_1.specialist, 
                                                                         method = "total")
#TSS_Physeq.2<- psadd::subset_samples_no_zero(TSS_Physeq.2,Treatment=='Arable')

bact.ord.2.specialist <- ordinate(TSS_Physeq.2.specialist ,method = "PCoA",distance="bray")

### Aitchison 

aitchison.dis.specialist<-dist_calc(specialist.2site.physeq,dist ="aitchison") %>% 
  ord_calc()

# unconstrained PCA ordination
unconstrained_aitchison_pca_rclr <- physeq2sites_plus_1.specialist %>%
  tax_transform(rank = "unique",trans = "rclr") %>%
  ord_calc()

unconstrained_aitchison_pca %>%
  ord_plot(color = "Treatment", shape = "Country.x", size = 2) +
  scale_colour_brewer(palette = "Dark2")

unconstrained_aitchison_pca@ord$Ybar
### Aitchison 
specialist.2site.physeq %>%
  tax_transform(rank = "unique", trans = "identity") %>% # don't transform!
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "Treatment", shape = "Country.x", size = 2) +
  scale_colour_brewer(palette = "Dark2")

specialist.2site.physeq %>%
  tax_transform("clr") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>%
  ord_plot(color = "Treatment", shape = "Country.x", size = 2) +
  scale_colour_brewer(palette = "Dark2")+
  theme_cowplot()


#### Testing aitchinson per Quagliariello et al. 2022

pp.hell.bray<- specialist.2site.physeq %>%
  tax_transform(rank = "unique", trans = "hellinger",zero_replace = 1) %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "PCoA"
  ) %>%  
  ord_plot(color = "Treatment", shape = "Country.x", size = 2) +
  scale_colour_brewer(palette = "Dark2")

pp.3<-pp.hell.bray$data
ggplot(pp12,aes(PC1,PC2))+geom_point(aes(colour=Treatment,size=5))+
  theme_minimal()


dist_calc(physeq2sites_plus_1.specialist,dist ="bray")-> bray.endemic.calc


as.data.frame(bray.endemic.calc@dist)-> prueba1


### 


bact.pcoa.DF.2.specialist <- plot_ordination(specialist.2sites,
                                             aitchison.dis.specialist,
                                             type = "sites", 
                                             color = "Treatment",
                                             shape = "Country",justDF = TRUE)

###### alpha endemic#####

physeq_table %>% microbiome::alpha(index = c("observed","shannon")) %>% 
  rownames_to_column("Sample_ID")->richness.endemic

left_join(richness.endemic,design %>% select(Sample_ID,Country.x,Treatment),
          by="Sample_ID")->combined.endemic.richness.samples


combined.endemic.richness.samples %>%filter(Treatment=="Arable") %>%  group_by(Treatment) %>%
  slice_sample(n = 61)-> Arable.61


combined.endemic.richness.samples %>%filter(Treatment=="Grassland")->Grassland.61


rbind(Arable.61,Grassland.61)->combined.61

combined.61$Country.x <-factor(combined.61$Country.x, 
                               levels =c("Spain","France","Switzerland","Germany","Sweden") )

combined.61  %>% 
  ggplot(aes(x=Country.x,y=observed,fill=Treatment))+
  geom_boxplot()+
  labs(y="Observed Richness",x="Country", title = "Overall Observed Richness")+
  scale_fill_manual(values = c("#FDC086","#7FC97F"))+
  ggpubr::stat_compare_means(label = "p.signif",hide.ns = TRUE,size=10)+
  theme_bw(16)+
  theme(legend.position = "bottom")

combined.61  %>% 
  ggplot(aes(x=Treatment,y=observed,fill=Treatment))+
  geom_boxplot()+
  labs(y="Observed Richness",x="Country")+
  scale_fill_manual(values = c("#FDC086","#7FC97F"))+
  theme_bw(16)+
  theme(legend.position = "none")+
  ggpubr::stat_compare_means(hide.ns = TRUE,size=5)


#### Cowplot test
library(cowplot)

main.plot <- 
  combined.61  %>% 
  ggplot(aes(x=Country.x,y=observed,fill=Treatment))+
  geom_boxplot()+
  labs(y="Observed Richness",x="Country", title = "Overall Observed Richness")+
  scale_fill_manual(values = c("#FDC086","#7FC97F"))+
  ggpubr::stat_compare_means(label = "p.signif",hide.ns = TRUE,size=10)+
  theme_cowplot(16)+
  theme(legend.position = "bottom")

inset.plot <- 
  combined.61  %>% 
  ggplot(aes(x=Treatment,y=observed,fill=Treatment))+
  geom_boxplot()+
  labs(y="",x="")+
  scale_fill_manual(values = c("#FDC086","#7FC97F"))+
  theme_cowplot()+
  theme(legend.position = "none")+
  ggpubr::stat_compare_means(label = "p.signif",hide.ns = TRUE,size=5)

plot.with.inset <-
  ggdraw() +
  draw_plot(main.plot) +
  draw_plot(inset.plot, x = 0.08, y = .63, width = .23, height = .33)






##### Proportion

psadd::subset_samples_no_zero(physeq_table,Country.x=="Spain") # 1524
psadd::subset_samples_no_zero(physeq_table,Country.x=="France") # 1717
psadd::subset_samples_no_zero(physeq_table,Country.x=="Switzerland") # 2546
psadd::subset_samples_no_zero(physeq_table,Country.x=="Germany") # 2438
psadd::subset_samples_no_zero(physeq_table,Country.x=="Sweden") # 1446


psadd::subset_samples_no_zero(specialist.2site.physeq, Country.x=="Spain") # 234
psadd::subset_samples_no_zero(specialist.2site.physeq, Country.x=="France") # 126
psadd::subset_samples_no_zero(specialist.2site.physeq, Country.x=="Switzerland") #362
psadd::subset_samples_no_zero(specialist.2site.physeq, Country.x=="Germany")# 279
psadd::subset_samples_no_zero(specialist.2site.physeq, Country.x=="Sweden") # 110



OTUbycountry.endemicprop<- data.frame(
  Country=c("Spain","France","Switzerland","Germany","Sweden"),
  Endemic.OTUs=c(234,126,362,279,110),
  Overall.OTUs=c(1524,1717,2546,2438,1446)) %>% mutate(Endemic=round((Endemic.OTUs/Overall.OTUs),3),
                                                       Overall=1-Endemic)

OTUbycountry.endemicprop$Country<-factor(OTUbycountry.endemicprop$Country,
                                         levels = c("Spain","France","Switzerland","Germany","Sweden")) 

OTUbycountry.endemicprop %>% 
  ggplot(aes(x=Country,y=Endemic))+geom_bar(fill="gray",stat = "identity",width = 0.7)+
  #scale_fill_manual(values = c("gray"))+
  labs(fill="Proportion",y="Proportion of Endemic Taxa")+
  ylim(0,.20)+
  geom_text(aes(label = Endemic), 
            position = position_dodge(width = 0.6),
            vjust=-0.5, color = "black",
            size = 5, 
            fontface = "bold")+
  theme_bw(16)


######## Scatter Plots ##########

### Index
mi <- readxl::read_excel("Subset_CropRichness_ManagementPractices2.xlsx")
mi2<- readxl::read_excel("MI Index for Sam.xlsx",sheet = "MI index")



#####
bact.pcoa.DF.2.specialist %>% filter(Treatment=='Arable') -> pcoa.bray.rareTaxa.2sites.Arable.specialist
pcoa.bray.rareTaxa.2sites.Arable.specialist %>% left_join(mi2)->axis.MI.endemic

PC.AXIS %>% select(Sample_ID,PC1) %>% left_join(mi2) %>% 
  drop_na() %>% 
  ggplot(aes(as.numeric(MI),PC1)) +  
  geom_point(size=5, 
             color="black", 
             alpha = .85,
             shape=21,
             fill="orange",
             position="dodge")+ 
  geom_smooth(method = lm,se=TRUE,colour="black", size=1.5, formula = y ~ x)+
  ggpubr::stat_cor(label.x = 0.65, 
                   label.y = 0.05, 
                   method = "spearman",label.sep = "\n", size = 10)

hell.bray.axis<-pp.3%>%  select(Sample_ID,MDS1) %>% left_join(mi2) %>% 
  drop_na() %>% 
  ggplot(aes(as.numeric(MI),MDS1)) +  
  geom_point(size=5, 
             color="black", 
             alpha = .85,
             shape=21,
             fill="gray",
             position="dodge")+ 
  geom_smooth(method = lm,se=TRUE,colour="black", size=1.5, formula = y ~ x)+
  ggpubr::stat_cor(label.x = 0.65, 
                   label.y = 0.05, 
                   method = "spearman",label.sep = "\n", size = 10)



aitchison.axis+labs(title = "Aitchison and PCoA") +clr.PCA.axis +labs(title = "CLR and PCA")+ hell.bray.axis
  
  

####
  axis.MI.endemic %>% 
  ggplot(aes(as.numeric(MI),Axis.1)) +  
  geom_point(size=5, 
             color="black", 
             alpha = .85,
             shape=21,
             fill="orange",
             position="dodge")+ 
  geom_smooth(method = lm,se=TRUE,colour="black", size=1.5, formula = y ~ x)+
  #geom_smooth(method=lm , se = TRUE, colour="black", size=1.5, formula = y ~  poly(x, 2))+
  #scale_color_manual(values = "#FDC086")+
  #facet_wrap(~name, scales = "free_y" , 
  #strip.position="left")+
  #ggpubr::stat_regline_equation(label.x =0, label.y = 0.40,)+
  ggpubr::stat_cor(label.x = 0.65, label.y = 0.25, method = "spearman",label.sep = "\n", size = 10) +
  theme_bw(24)+
  labs(x="Managemane Intensity",y="Endemic PCoA 1")+
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
  
  

specialist.2site.physeq %>%
    tax_transform(rank = "unique", trans = "hellinger",zero_replace = 1) %>%
    dist_calc(dist = "bray")->loco 

dist_get(loco) ->loco1
dist_df <- as.data.frame(as.matrix(loco1)

loco@dist

  
physeq2sites_plus_1.specialist %>% dist_calc(dist = "bray") %>% dist_get()
  
dist_df %>% rownames_to_column("Sample_ID")
  
dist_df %>% pivot_longer(2:218) %>% 
  left_join(mi2) %>% 
  left_join(mi %>% select(Sample_ID,CropRichness)) %>% 
  left_join(design %>% select(Sample_ID,Country.x)) %>% 
  drop_na() %>% 
  ggplot(aes(as.numeric(MI),value)) +  
  geom_point(size=2, 
             color="black", 
             alpha = .60,
             shape=21,
             fill="orange",
             position="dodge")+ 
  geom_smooth(method = lm,se=TRUE,colour="black", size=1.5, formula = y ~ x)+
  ggpubr::stat_cor(label.x = 0.65, 
                   label.y = 0.30, 
                   method = "spearman",label.sep = "\n", size = 10)+
  theme_bw()

dist_df %>% 
  rowwise() %>% 
  mutate(rowmean=mean(c(2:218))) %>% select(Sample_ID,rowmean) ->locotest


# 2023 Script BioDiversa version fungal traits included ----------------

## adding Verified column of guilds
#### Getting the verified funguilds

verified.guilds<-readxl::read_excel('~/Dropbox/Lennel/FunGuild.FunTraits.Genus.clean.verified.xlsx') %>% 
  column_to_rownames('OTUs') %>% select(Verified)

otu_taxa_new_guilds<-cbind(otu_taxa,verified.guilds)

## import design file
design <- read.csv("DATA/design.csv",row.names = 1, stringsAsFactors = F)

design<-left_join(design,soilProp,by="Sample_ID")
row.names(design) <- design$Sample_ID

## correct typos
design$Treatment[design$Treatment == "High "]<- "High"
design$Treatment[design$Treatment == "Low "]<- "Low"

design$Treatment[design$Treatment == "High" | design$Treatment == "Low" ] <- "Arable"

design[design$Treatment == "Organic"| design$Treatment == "Conventional","Treatment"] <- "Arable"
design_pre <- design
## create phyloseq object
physeq_table<- phyloseq(otu_table(OTU, taxa_are_rows=T),
                        tax_table(as.matrix(otu_taxa_new_guilds)),
                        sample_data(design))


### Filtering taxa 2 sites or less

physeq.tab.2sites <-filter_taxa(physeq = physeq_table,
                                function(x){sum(x > 0) <= 2}, prune = TRUE)

#### Niche Breadth Analysis
library(MicroNiche)
###### Niche Breadth
BacteriaNiche2<- as.data.frame(physeq.tab.2sites@otu_table) %>% 
  rownames_to_column(.,"OTU")

Bacteria_RareRanks2 <- as.data.frame(physeq.tab.2sites@tax_table) %>% 
  rownames_to_column(.,"OTU")




Prueba.NicheBreadth <- levins.Bn(BacteriaNiche2,2,
                    design$Treatment)


specialist.2sites<- Prueba.NicheBreadth %>% filter(Bn<0.55 & P.adj < 0.05)

keep.specialist <- rownames(specialist.2sites)


specialist.2site.physeq<- prune_taxa(keep.specialist,physeq.tab.2sites)


#### Composition of specialist 2 site or less. 

specialist.venn.2sites <- MicEco::ps_venn(specialist.2site.physeq,'Treatment',
                fill = c("#FDC086","#7FC97F"))


df.barplot.specialists <- data.frame(Treatment=c("Arable", "Grassland"),
                         Total_OTU=c(438, 692))

BarplotRareTaxa.specialists<-
  ggplot(df.barplot.specialists,
         aes(x=Treatment,y=Total_OTU,fill=Treatment))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#FDC086","#7FC97F"))+
  theme_bw(24)+labs(x="",y="Specialists total number of OTUs")+
  theme(legend.position = "none")



TotalOTusplotwVenn.specialists<-
  BarplotRareTaxa.specialists+
  inset_element(specialist.venn.2sites,
                left = 0.05,
                bottom = 0.74,
                right = 0.45,
                top = 0.99)

TotalOTusplotwVenn+TotalOTusplotwVenn.specialists

#### Composition
no_unidentified_order_depleted <-subset_taxa(genus.2sites.specialists.physeq,Genus != "g__")
ps1.com.fam_depleted <- microbiomeutilities::aggregate_top_taxa2(no_unidentified_order_depleted, 
                                                                 "Genus", top = 20)
ps2_depleted<- transform_sample_counts(ps1.com.fam_depleted, function(x) (x / sum(x))*100)


library(metagMisc)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

specialist.df<-phyloseq_to_df(ps1.com.fam_depleted,addtax = TRUE)


pivot_longer(specialist.df,cols = 4:220) %>% 
  left_join(design,by=c('name'='Sample_ID')) %>% 
  group_by(Treatment,Order) %>% 
  summarize(value=sum(value)) %>% 
  #group_by(Order) %>% 
  #arrange(desc(value)) %>% 
  ggplot(aes(Treatment,value,fill=Order))+
  #geom_boxplot()+
  #facet_wrap(~Order,scales = 'free')+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = getPalette(20))+
  theme_bw(16)+
  theme(legend.position = 'bottom')+
  guides(fill=guide_legend(ncol=2))+
  labs(y='Total reads',title = "Specialist (Levin's Bn < 0.55) in 2 sites or less",fill='Order')+
  coord_flip()-> composition.2sites.specialist.plot



composition.2sites.specialist.plot/composition.2sites.specialist.plot.guild -> combined.specialist.composition.new



TotalOTusplotwVenn+TotalOTusplotwVenn.specialists


specialist.venn.2sites/composition.2sites.specialist.plot


ggpubr::ggarrange(specialist.venn.2sites,composition.2sites.specialist.plot,
                  nrow = 2)


cleanTraits.Guilds$OTUs %in% keep.specialist

specialist.guild.traits<-subset(cleanTraits.Guilds, OTUs %in% keep.specialist)

##### Endemic Specialists for MI index


physeq2sites_plus_1.specialist <- transform_sample_counts(specialist.2site.physeq, 
                                               function(x) x+1 )

TSS_Physeq.2.specialist <- metagMisc::phyloseq_standardize_otu_abundance(physeq2sites_plus_1.specialist, 
                                                              method = "total")
#TSS_Physeq.2<- psadd::subset_samples_no_zero(TSS_Physeq.2,Treatment=='Arable')

bact.ord.2.specialist <- ordinate(TSS_Physeq.2.specialist ,method = "PCoA",distance="bray")

bact.pcoa.DF.2.specialist <- plot_ordination(TSS_Physeq.2.specialist,
                                  bact.ord.2.specialist,
                                  type = "sites", 
                                  color = "Treatment",
                                  shape = "Country",justDF = TRUE )




#### Scatter plots M.I index. 
bact.pcoa.DF.2.specialist %>% filter(Treatment=='Arable') -> pcoa.bray.rareTaxa.2sites.Arable.specialist


PCoA.1.MI.scatter.2sites.specialist<-alpha.endemic %>% left_join(mi) %>% 
  #drop_na(.) %>% 
  ggplot( aes(as.numeric(CropRichness),Shannon)) +  
  geom_point(size=5, 
             color="black", 
             alpha = .85,
             shape=21,
             fill="orange",
             position="dodge")+ 
  geom_smooth(method = lm,se=TRUE,colour="black", size=1.5, formula = y ~ x)+
  #geom_smooth(method=lm , se = TRUE, colour="black", size=1.5, formula = y ~  poly(x, 2))+
  #scale_color_manual(values = "#FDC086")+
  #facet_wrap(~name, scales = "free_y" , 
  #strip.position="left")+
  #ggpubr::stat_regline_equation(label.x =0, label.y = 0.40,)+
  ggpubr::stat_cor(label.x = 0.65, label.y = 2.5, method = "spearman",label.sep = "\n", size = 10) +
  theme_bw(24)+
  labs(x="Crop Richness",y="Endemic Shannon")+
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



genus.test<- phyloseq_to_df(ps2_depleted,addtax = TRUE)


combined.genus<-genus.test %>%
  pivot_longer(cols = 4:220) %>% 
  left_join(design,by=c("name"="Sample_ID")) 



combined.genus%>% 
  ggplot(aes(x=Organic_C_.,y=value))+geom_point()+facet_wrap(~Genus,scales = "free_y")


alpha.endemic <-estimate_richness(specialist.2site.physeq,measures = c("Observed","Shannon")) %>% rownames_to_column("Sample_ID")

#### funfun

funfunTable<-read_tsv("Results.tsv")


funfun.new <- funfunTable[,keep.specialist]

funfun.new %>% data_frame(funfunTable$Function,funfun.new)-> funfun.2


ko.hierchaies<-read_csv("~/Downloads/Table.KO.path.hierchachies.csv") %>% select(A:D) %>% 
  mutate(A=gsub("^A\\s+", "",A),
         B=gsub("^B\\s+", "",B),
         C=gsub("^C\\s+", "",C),
         D=gsub("^D\\s+", "",D)) %>% 
  select(A:C) %>% 
  na.omit() %>% 
  distinct(C,.keep_all = TRUE)




funfunTable %>% left_join(ko.hierchaies %>% select(A:C),by = c("Function"="C"),keep=FALSE) %>% 
  select(A,B,"Function",starts_with("OTU"))->test.new.overall
  








funfun.2$`funfunTable$Function`
install.packages("treemap")
library(treemap)


test.new.overall %>% pivot_longer(cols = 4:4176) %>%
  group_by(B) %>% 
  filter(!str_detect(`funfunTable$Function`,"PATH:")) %>% 
  dplyr::summarize(value = sum(value)) %>% 
  arrange(-value) %>% treemap(., 
                              index = "B", 
                              vSize="value",
                              type = "index", 
                              fontsize.labels = 10, 
                              fontcolor.labels = "black")


test.new.overall %>% pivot_longer(cols = 4:4176) %>%
  filter(name %in% get_unique.otus$Arable) %>% 
  filter(B == "09111 Xenobiotics biodegradation and metabolism") %>% 
  group_by(Function) %>%  
  dplyr::summarize(value = sum(value)) %>% 
  mutate(percentage = value / sum(value) * 100,
         csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos)) %>%
  ggplot(aes(x = "", y = value, fill = Function)) +
  geom_col(width = 1,color=1,alpha=0.7) +
  #geom_text(aes(label = paste0(round(percentage, 3), "%")), 
            #position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  ggrepel::geom_label_repel(aes(y = pos, label = paste0(round(percentage,2), "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE,alpha=0.7)+
  scale_fill_manual(values = c("#0000ee","#9933cc","#009999", "#ff0000","#ff9933",
                               "#ff6600","#3399ff","#ff6699","#00cc33","#cc3366","#ccaa99","gray"))+
  labs(title = "Arable",fill = "Functional category") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))->arable331.pie

arable331.pie+grassland535.pie+plot_layout(guides = "collect") & theme(legend.position = "bottom")

test.new %>% pivot_longer(cols = 4:976) %>%
  filter(name %in% get_unique.otus$Arable) %>% 
  filter(A == "09100 Metabolism") %>% 
  group_by(B) %>%  
  dplyr::summarize(value = sum(value)) %>% 
  mutate(percentage = value / sum(value) * 100) %>%
  ggplot(aes(x = "", y = value, fill = B)) +
  geom_col(width = 1) +
  geom_text(aes(label = paste0(round(percentage, 3), "%")), 
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  labs(title = "Grassland Endemic OTUs", fill = "Functional category") +
  theme_void() -> Grassland.pie




library(patchwork)
  
Arable.pie+Grassland.pie+plot_layout(guides = "collect")& theme(legend.position = "bottom")
  
  
  
  ggplot(aes(x=reorder(`funfunTable$Function`,-value),y=value))+
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





  ggplot(aes(x="", y=value, fill=`funfunTable$Function`)) + 
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0) + 
  geom_text(aes(label = paste0(round(value*100), "%")), 
            position = position_stack(vjust = 0.5))



summarise(count=sum(value)) %>% 
  mutate(prop=count/sum(value)) %>% 
  filter(prop > 0.02)->keep

temp <-
  otuL %>%
  mutate(newOrder=case_when(Order %in% keep$Order~Order,
                            TRUE~'other'))

temp %>% 
  group_by(SampleID) %>% 
  mutate(prop=count/sum(count))->rel.abun.ecm

rel.abun.ecm %>% left_join(.,meta,by="SampleID") %>% 
  ggplot(aes(x=newOrder,
             y=prop,
             fill=site_code))+
  geom_boxplot()

ggplot(temp,aes(x=site_code,y=count,fill=newOrder))+
  geom_bar(stat="identity",position=position_fill())+
  scale_fill_brewer(palette = 'Set3') -> plot



funfunTable %>% separate(Function, c('KO_number', 'Annotation'),extra = "merge") %>% 
  mutate(KO_number = ifelse(is.na(KO_number), "NA", paste0("ko", as.character(KO_number)))) ->test


test

##### Piechart dependent on sites


arable.specialist <-specialist.2site.physeq %>% psadd::subset_samples_no_zero(Treatment=="Arable")
grassland.specialist <- specialist.2site.physeq %>% psadd::subset_samples_no_zero(Treatment=="Grassland")

get_unique.otus<-MicEco::ps_venn(specialist.2site.physeq, group = "Treatment",plot = FALSE)


arable.onlyspecialistotus <-rownames(arable.specialist@otu_table)

grassland.onlyspecialistotus <-rownames(grassland.specialist@otu_table)

