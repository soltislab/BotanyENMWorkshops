# Point Based
# Ecological Analysis using Points
# This is based off scripts created by Hannah Owens and James Watling, and Anthony Melton.  
# Modified and created by ML Gaynor. 

# Load Packages
library(raster)
library(dplyr)
library(tidyr)
library(gtools)
library(factoextra)
library(FactoMineR)
library(multcompView)
library(ggsci)
library(gridExtra)
library(ggplot2)
library(ecospat)
library(dismo)

# Load function
source("functions/ggbiplot_copy.R")

# Load data file
alldf <- read.csv("data/cleaning_demo/maxent_ready/diapensiaceae_maxentready_20210625.csv")


# Raster layers
list <- list.files("data/climate_processing/PresentLayers/all", full.names = T, recursive = FALSE) 
list <- mixedsort(sort(list))
envtStack <- stack(list)


####################################################################################################################################################################################
# ENM - Realized Niche

## For each occurence record, extract the value for each bioclim variables using the package raster.
##### Extract value for each point
ptExtracted <- raster::extract(envtStack, alldf[3:2])

#### Convert to data frame
ptExtracteddf <- as.data.frame(ptExtracted)

#### Add species name
ptExtracteddf <- ptExtracteddf %>%
                 dplyr::mutate(name = as.character(alldf$name), x = alldf$long, y = alldf$lat)

#### Drop any NA
ptExtracteddf <- ptExtracteddf %>% 
                 tidyr::drop_na(bio_3, bio_7, bio_8, bio_9, bio_14, bio_15, bio_18, elev)

## PCA 
### Create two dataframes.
data.bioclim <- ptExtracteddf[, 1:8]
data.species <- ptExtracteddf[, 9]

### Using only the bioclim columns to run the principal components analysis.
data.pca <- prcomp(data.bioclim, scale. = TRUE) 

#### Understanding the PCA - Optional 
##### When you use the command prcomp your loading variables show up as rotational variables. 
##### Thanks to a really great answer on stack overflow: https://stackoverflow.com/questions/43407859/how-do-i-find-the-link-between-principal-components-and-raw-datas-variables
##### you can even convert the rotational 
##### variable to show the relative contribution.
loadings <- data.pca$rotation
summary(loadings)

##### There are two options to convert the loading to show the relative contribution, 
##### they both give the same answer so either can be used.
loadings_relative_A <- t(t(abs(loadings))/rowSums(t(abs(loadings))))*100
loadings_relative_A

loadings_relative_B <- sweep(x = abs(loadings), MARGIN = 2, STATS = colSums(abs(loadings)), FUN = "/")*100
loadings_relative_B

#### Plotting the PCA
##### First, I made a theme to change the background of the plot. Next, I changed the plot margins and the text size.
theme <- theme(panel.background = element_blank(),
               panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"), 
               axis.text = element_text(size = 12), 
               legend.text = element_text(size = 12), 
               legend.title = element_text(size = 12),
               text = element_text(size = 12))

##### Set colors
pal <- pal_locuszoom()(4)

##### Next, ggbiplot where obs.scale indicates the scale factor to apply to observation, 
##### var.scale indicates the scale factor to apply to variables, 
##### ellipse as TRUE draws a normal data ellipse for each group, 
##### and circle as TRUE draws a correlation circle.
g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, 
              groups = data.species, ellipse = TRUE, circle = TRUE) +
     scale_color_manual(name = '', values = pal) +
     theme(legend.direction = 'vertical', legend.position = 'bottom', 
           legend.text = element_text(size = 12, face = "italic")) +
     theme
 
g

####################################################################################################################################################################################
# ANOVA
## Simple function to run an ANOVA and a post-hoc Tukey-HSD test
stat.test <- function(data = ptExtracteddf, x = "name", y){
    bioaov <- aov(as.formula(paste0(y,"~",x)), data = data) 
    TH <- TukeyHSD(bioaov, "name")
    m <- multcompLetters(TH$name[,4])
    groups <- data.frame(groups = m$Letters, name = names(m$Letters))
    return(groups)
}

###### BIO3 ONLY ######
bio3aovplot <- ggplot(ptExtracteddf, aes(x = name, y = bio_3)) +
               geom_boxplot(aes(fill = name)) +
               scale_color_manual(name = '', values = pal) +
               geom_text(data = stat.test(y = "bio_3"), 
                        mapping = aes(x = name,
                                      y = max(ptExtracteddf["bio_3"]+1), 
                                      label = groups), 
                        size = 5, inherit.aes = FALSE) +
               theme(axis.text.x = element_text(angle = 90, size = 8, face = 'italic'))
bio3aovplot

## Loop through all variables
variablelist <- colnames(ptExtracteddf)[1:8]
plotlist <- c()
for(i in 1:8){
  bio <- variablelist[i]
  tempdf <- ptExtracteddf %>%
            dplyr::select(name, variablelist[i])
  plotlist[[i]] <- ggplot(tempdf, aes(x = name, y = tempdf[,2])) +
                    geom_boxplot(aes(fill = name)) +
                    scale_colour_manual(name = 'Species', values = pal) +
                    geom_text(data = stat.test(y = variablelist[i]), 
                              mapping = aes(x = name,
                                            y = max(tempdf[,2]+1), 
                                            label = groups), 
                              size = 5, inherit.aes = FALSE) +
                    scale_x_discrete(labels = c('G', 'Pba','Pbr', 'S')) +
                    ggtitle(label = paste0(variablelist[i])) +
                    ylab(paste0(variablelist[i])) +
                    theme(legend.position = "none")
}

gridExtra::grid.arrange(grobs = plotlist)

#################
# Ecospat Niche Overlap and Niche Equivalency
## Set up background points
bg1 <- randomPoints(mask = envtStack, n = 1000, p = alldf[,3:2])
bg1.env <- raster::extract(envtStack, bg1)
bg1.env <- data.frame(bg1.env)
allpt.bioclim <- rbind(bg1.env, data.bioclim)

## dudi.PCA to reduce variables
pca.env <- dudi.pca(allpt.bioclim,
                    center = TRUE, # Center by the mean
                    scannf = FALSE, # Don't plot
                    nf = 2) # Number of axis to keep 
## Pull out scores for each species
p1.score <- suprow(pca.env, dplyr::filter(ptExtracteddf, name == "Galax urceolata")[, 1:8])$li
p2.score <- suprow(pca.env, dplyr::filter(ptExtracteddf, name == "Pyxidanthera barbulata")[, 1:8])$li
p3.score <- suprow(pca.env, dplyr::filter(ptExtracteddf, name == "Pyxidanthera brevifolia")[, 1:8])$li
p4.score <- suprow(pca.env, dplyr::filter(ptExtracteddf, name == "Shortia galacifolia")[, 1:8])$li
scores.clim <- pca.env$li

## Visualize 
plot(scores.clim, pch = 16, asp = 1,
     col = adjustcolor(1, alpha.f = 0.2), cex = 2,
     xlab = "PC1", ylab = "PC2") 
points(p1.score, pch = 18, col = pal[1], cex = 2)
points(p2.score, pch = 18, col = pal[2], cex = 2)
points(p3.score, pch = 18, col = pal[3], cex = 2)
points(p4.score, pch = 18, col = pal[4], cex = 2)

## Kernel density estimates
### create occurrence density grids based on the ordination data
z1 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, p1.score, R = 100)
z2 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, p2.score, R = 100)
z3 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, p3.score, R = 100)
z4 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, p4.score, R = 100)
zlist  <- list(z1, z2, z3, z4)

## Niche Overlap 
### Schoener's D ranges from 0 to 1
### 0 represents no similarity between niche space
### 1 represents completely identical niche space
overlapD <- matrix(ncol = 2, nrow = 7)
n <- 1
for(i in 1:3){
  for(j in 2:4){
    if(i != j){
      overlapD[n, 1]<- paste0("z", i, "-", "z", j)
      overlapD[n, 2]<- ecospat.niche.overlap(zlist[[i]], zlist[[j]], cor = TRUE)$D
      n <- n + 1
    }
  }
}

overlapDdf <- data.frame(overlapD)
overlapDdf

## Niche Overlap Visualization 
par(mfrow=c(2,1))
ecospat.plot.niche.dyn(z1, z4, quant=0.25, interest = 1
                       , title= "Niche Overlap - Z1 top", name.axis1="PC1", name.axis2="PC2")
ecospat.plot.niche.dyn(z1, z4, quant=0.25, interest = 2
                       , title= "Niche Overlap - Z4 top", name.axis1="PC1", name.axis2="PC2")

## Niche Equivalency Test
### Based on Warren et al. 2008 - Are the two niche identical?
#### Hypothesis test for D, null based on randomization 
#### H1: the niche overlap is higher than expected by chance (or when randomized)
eq.test <- ecospat.niche.equivalency.test(z1, z4, rep = 10, alternative = "greater")
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")


## Niche Similarity Test
### Based on Warren et al. 2008 - Are the two niche similar?
### Can one species’ niche predicted the occurrences of a second species better than expected by chance.
sim.test <- ecospat.niche.similarity.test(z1, z4, rep = 10, alternative = "greater", rand.type=2)
ecospat.plot.overlap.test(sim.test, "D", "Similarity")

