rm(list=ls())
setwd("~/R/MRes winter project")

{library(rlang)
  library(ape)
  library(phylolm)
  library(phytools)
  library(stringr)
  library(dplyr)
  library(tibble)
  library(janitor)
  library(magrittr)
  library(caper)
  library(ggtree)
  library(ggtreeExtra)
  library(ggplot2)
  library(ggnewscale)
  library(readr)
  library(ggpubr)
  library(vioplot)
  library(GGally)
}

# Read the tree
all_birdtree <- read.tree("BirdzillaHackett1.tre")
# Read the trait data
passerine <- read_csv("Elev_Wing_data_passerines.csv")%>% clean_names()


################################ 1.DATA PREPARATION ############################

# Create new columns for overall max, min and mean elevation
passerine$all_max<- pmax(passerine$upper_limit, passerine$sp_max_elev_quintero, 
                         passerine$max_wb, passerine$new_max, na.rm=TRUE)
passerine$all_min<- pmin(passerine$lower_limit, passerine$sp_min_elev_quintero, 
                         passerine$min_wb, passerine$new_min, na.rm=TRUE)
passerine$all_mean<- (passerine$all_max + passerine$all_min)/2

# Group al index into binary values, low(0) and high(1)
passerine$al_binary <- ifelse(passerine$al_index == 1, 0, NA)
passerine$al_binary <- ifelse(passerine$al_index > 1, 1, passerine$al_binary)

# Group migration into binary values, non-migratory(0) and migratory(1)
passerine$migration_binary <- ifelse(passerine$migration < 3, 0, NA)
passerine$migration_binary <- ifelse(passerine$migration == 3, 1, passerine$migration_binary)

# Group habitat density into binary values, dense habitat(0) and open habitat(1)
passerine$habitat_binary <- ifelse(passerine$habitat_density == 1, 0, NA)
passerine$habitat_binary <- ifelse(passerine$habitat_density > 1, 1, passerine$habitat_binary)

# Make trophic level as numeric variables
passerine$trophic_level_num <- ifelse(passerine$trophic_level == 'Herbivore', 1, NA)
passerine$trophic_level_num <- ifelse(passerine$trophic_level == 'Omnivore', 2, passerine$trophic_level_num)
passerine$trophic_level_num <- ifelse(passerine$trophic_level == 'Carnivore', 3, passerine$trophic_level_num)

# Change jetz name
passerine$species <- gsub(" ", "_", passerine$species)

# Change row names to species name
row.names(passerine) <- passerine$species 

# Create high-certainty dataset
high <- passerine %>% filter(al_uncertainty < 3)


################################ 2.MODEL ANALYSIS ##############################

# Read 50 random phylogenetic trees
ran_tree <- sample(all_birdtree, size=50)

# Create empty dataframes to store analysis results
summary_hwi_max<-NA
summary_wl_max<-NA
summary_sl_max<-NA

# Loop 50 times
for (i in 1:50){
  # Select one tree each time
  one_tree <- ran_tree[[i]]
  tree <- drop.tip(one_tree, setdiff(one_tree$tip.label, passerine$species))
  # Create comparative data for PGLS
  comp<-comparative.data(data = passerine, phy = tree, names.col = 'species', na.omit=FALSE,warn.dropped = FALSE)
  
  # PGLS analysis for HWI
  ran_hwi_max <- phylolm(scale(sqrt(hand_wing_index)) ~ scale(sqrt(all_max)) + scale(log(mass)) + migration_binary + al_binary+ habitat_binary + scale(trophic_level_num), 
                     data = comp$data, phy = comp$phy, model='lambda')
  summary_hwi_max_temp<-as.data.frame(summary(ran_hwi_max)$coefficients)
  summary_hwi_max_temp$predictor <- rownames(summary_hwi_max_temp)
  summary_hwi_max<-rbind(summary_hwi_max,summary_hwi_max_temp)
  
  # PGLS analysis for wing length
  ran_wl_max <- phylolm(scale(log(wing_length)) ~ scale(sqrt(all_max)) + scale(log(mass)) + migration_binary + al_binary+ habitat_binary+ scale(trophic_level_num),
                    data = comp$data, phy = comp$phy, model='lambda')
  summary_wl_max_temp<-as.data.frame(summary(ran_wl_max)$coefficients)
  summary_wl_max_temp$predictor <- rownames(summary_wl_max_temp)
  summary_wl_max<-rbind(summary_wl_max,summary_wl_max_temp)
  
  # PGLS analysis for secondary wing length
  ran_sl_max <- phylolm(scale(log(secondary1)) ~ scale(sqrt(all_max)) + scale(log(mass)) + migration_binary + al_binary+ habitat_binary+ scale(trophic_level_num),
                    data = comp$data, phy = comp$phy, model='lambda')
  summary_sl_max_temp<-as.data.frame(summary(ran_sl_max)$coefficients)
  summary_sl_max_temp$predictor <- rownames(summary_sl_max_temp)
  summary_sl_max<-rbind(summary_sl_max,summary_sl_max_temp)
}

# Average out the results to get average coefficients
{averaged_hwi_max <- summary_hwi_max %>% group_by (predictor) %>%summarise_all(mean)
  averaged_wl_max <- summary_wl_max %>% group_by (predictor) %>%summarise_all(mean)
  averaged_sl_max <- summary_sl_max %>% group_by (predictor) %>%summarise_all(mean)
 
  averaged_hwi_max <- averaged_hwi_max %>% slice(-c(1,8)) # remove the first and last row
  averaged_wl_max <- averaged_wl_max %>% slice(-c(1,8)) # remove the first and last row
  averaged_sl_max <- averaged_sl_max %>% slice(-c(1,8)) # remove the first and last row
  
  averaged_hw_max$conf_low <- averaged_hwi_max$Estimate-1.96*averaged_hwi_max$StdErr
  averaged_hwi_max$conf_high <- averaged_hwi_max$Estimate+1.96*averaged_hwi_max$StdErr
  averaged_wl_max$conf_low <- averaged_wl_max$Estimate-1.96*averaged_wl_max$StdErr
  averaged_wl_max$conf_high <- averaged_wl_max$Estimate+1.96*averaged_wl_max$StdErr
  averaged_sl_max$conf_low <- averaged_sl_max$Estimate-1.96*averaged_sl_max$StdErr
  averaged_sl_max$conf_high <- averaged_sl_max$Estimate+1.96*averaged_sl_max$StdErr
  
  averaged_hwi_max$Significance <- NA
  averaged_hwi_max$Significance[averaged_hwi_max$p.value < 0.05] <- "Significant"
  averaged_hwi_max$Significance[averaged_hwi_max$p.value > 0.05] <- "Insignificant"
  
  averaged_wl_max$Significance <- NA
  averaged_wl_max$Significance[averaged_wl_max$p.value < 0.05] <- "Significant"
  averaged_wl_max$Significance[averaged_wl_max$p.value > 0.05] <- "Insignificant"
  
  averaged_sl_max$Significance <- NA
  averaged_sl_max$Significance[averaged_sl_max$p.value < 0.05] <- "Significant"
  averaged_sl_max$Significance[averaged_sl_max$p.value > 0.05] <- "Insignificant"}

# Re-run the analysis with high-certainty dataset and mean elevation

summary_hwi_mean<-NA
summary_wl_mean<-NA

for (i in 1:50){
  # Select one tree each time
  one_tree <- ran_tree[[i]]
  tree <- drop.tip(one_tree, setdiff(one_tree$tip.label, high$species))
  # Create comparative data for PGLS
  comp<-comparative.data(data = high, phy = tree, names.col = 'species', na.omit=FALSE,warn.dropped = FALSE)
  
  # PGLS analysis for HWI
  ran_hwi_mean <- phylolm(scale(sqrt(hand_wing_index)) ~ scale(sqrt(all_mean)) + scale(log(mass)) + migration_binary + al_binary+ habitat_binary + scale(trophic_level_num), 
                         data = comp$data, phy = comp$phy, model='lambda')
  summary_hwi_mean_temp<-as.data.frame(summary(ran_hwi_mean)$coefficients)
  summary_hwi_mean_temp$predictor <- rownames(summary_hwi_mean_temp)
  summary_hwi_mean<-rbind(summary_hwi_mean,summary_hwi_mean_temp)
  
  # PGLS analysis for wing length
  ran_wl_mean <- phylolm(scale(log(wing_length)) ~ scale(sqrt(all_mean)) + scale(log(mass)) + migration_binary + al_binary+ habitat_binary+ scale(trophic_level_num),
                        data = comp$data, phy = comp$phy, model='lambda')
  summary_wl_mean_temp<-as.data.frame(summary(ran_wl_mean)$coefficients)
  summary_wl_mean_temp$predictor <- rownames(summary_wl_mean_temp)
  summary_wl_mean<-rbind(summary_wl_mean,summary_wl_mean_temp)
}

# Average out the results to get average coefficients
{averaged_hwi_mean <- summary_hwi_mean %>% group_by (predictor) %>%summarise_all(mean)
  averaged_wl_mean <- summary_wl_mean %>% group_by (predictor) %>%summarise_all(mean)
 
 averaged_hwi_mean <- averaged_hwi_mean %>% slice(-c(1,8)) # remove the first and last row
  averaged_wl_mean <- averaged_wl_mean %>% slice(-c(1,8)) # remove the first and last row
 
  averaged_hw_meani$conf_low <- averaged_hwi_mean$Estimate-1.96*averaged_hwi_mean$StdErr
  averaged_hwi_mean$conf_high <- averaged_hwi_mean$Estimate+1.96*averaged_hwi_mean$StdErr
  averaged_wl_mean$conf_low <- averaged_wl_mean$Estimate-1.96*averaged_wl_mean$StdErr
  averaged_wl_mean$conf_high <- averaged_wl_mean$Estimate+1.96*averaged_wl_mean$StdErr

    averaged_hwi_mean$Significance <- NA
  averaged_hwi_mean$Significance[averaged_hwi_mean$p.value < 0.05] <- "Significant"
  averaged_hwi_mean$Significance[averaged_hwi_mean$p.value > 0.05] <- "Insignificant"
  
  averaged_wl_mean$Significance <- NA
  averaged_wl_mean$Significance[averaged_wl_mean$p.value < 0.05] <- "Significant"
  averaged_wl_mean$Significance[averaged_wl_mean$p.value > 0.05] <- "Insignificant"}


################################ 3.VISUALIZATION ###############################

##### (1) Genus-level phylogenetic tree (Figure 3) #####

# Copy a list of all the tips from the tree.
bird_tips <- one_tree$tip.label

# Split the labels into two strings where there's an underscore 
passerine_genus <- passerine$species %>% str_split(pattern = "_", simplify= TRUE)
colnames(passerine_genus) <- c("genus", "species")
passerine_genus <- as.data.frame(passerine_genus)

# Copy the trait data to the new data frame
passerine_genus$all_max <- c(passerine$all_max)
passerine_genus$hand_wing_index <- c(passerine$hand_wing_index)
passerine_genus$wing_length <- c(passerine$wing_length)
# Remove NAs
passerine_genus <- na.omit(passerine_genus) 

# Combine the columns and add back in the underscore so they match the labels in the tree
genera_tips <- paste(passerine_genus$genus, passerine_genus$species, sep="_")

# Remove all the species except one per genus.
genera_tree <- drop.tip(one_tree, setdiff(one_tree$tip.label, genera_tips))
genera_tree$tip.label <- bird_genera$genus

genus_scores <- passerine_genus %>% 
  group_by(genus) %>% 
  summarise(genus = first(genus), # Get the name for tree plotting.
            genus_size = length(all_mean), #Get the genus size
            mean_max = sum(as.numeric(all_max))/genus_size, # Get the mean elevation for each genus
            mean_hwi = sum(as.numeric(hand_wing_index))/genus_size, # Get the mean HWI score for each genus
            mean_wl = sum(as.numeric(wing_length))/genus_size) # Get the mean wing length score for each genus

# Drop non-passerine genus
genera_tree <- drop.tip(genera_tree, setdiff(genera_tree$tip.label, genus_scores$genus))
# Add a column that's the node number matching the tree
genus_scores$node <- nodeid(genera_tree, genus_scores$genus)
# Join our data and tree together
genus_plot_data <-  full_join(genera_tree, genus_scores, by = "node")

# Plot the tree
(genus_plot <- ggtree(genus_plot_data, layout="fan",  size = 0.4,aes(colour=mean_max)) + 
    xlim(0,240) +   # Branch colour represents max elevation
    scale_colour_gradient(low='lightskyblue', high='orangered') +
    theme(plot.margin=margin(-500,-500,-500,-500)) +
   geom_fruit(geom=geom_bar,  # Add bars to represent mean wing length 
               mapping=aes(y=node, x=mean_wl ),  
               pwidth=0.14,
               orientation="y", 
               offset=0.03,
               stat="identity", fill="navy", colour="navy", width=0.2) +
   geom_fruit(geom=geom_bar,  # Add bars to represent mean HWI
               mapping=aes(y=node, x=mean_hwi ),  
               pwidth=0.14,
               orientation="y", 
               offset=0,
               stat="identity", fill="lightgreen", colour="lightgreen", width=0.2))

# Get the legend for elevation
(genus_legend <- ggtree(genus_plot_data, layout="fan",  size = 0.4,aes(colour=mean_max)) + 
    xlim(0,240) +  
    scale_colour_gradient(low='lightskyblue', high='orangered')+
    theme(legend.title=element_blank(), text = element_text(size = 8),
          legend.position = c(0.56,0.49), legend.direction = "horizontal", legend.title.align = 1,
          legend.key.width = unit(1.8, "cm"), legend.key.height = unit(0.6, "cm"), 
          legend.text = element_text(size = 9.0), 
          legend.margin = NULL, plot.margin=margin(-500,-500,-500,-500)))
plot(get_legend(genus_legend))

##### (2) Forest plots (Figure 4 & S1) #####

# Plots for full dataset using maximum elevation
cust_hwi_max <- data.frame(
  term = c("Habitat density", "Max elevation",'Trophic level','Aerial lifestyle','Migration','Body mass'),
  estimate = c(averaged_hwi_max$Estimate[2],averaged_hwi_max$Estimate[5],averaged_hwi_max$Estimate[6],averaged_hwi_max$Estimate[1],averaged_hwi_max$Estimate[3],averaged_hwi_max$Estimate[4]),
  conf.low = c(averaged_hwi_max$conf_low[2],averaged_hwi_max$conf_low[5],averaged_hwi_max$conf_low[6],averaged_hwi_max$conf_low[1],averaged_hwi_max$conf_low[3],averaged_hwi_max$conf_low[4]),
  conf.high = c(averaged_hwi_max$conf_high[2],averaged_hwi_max$conf_high[5],averaged_hwi_max$conf_high[6],averaged_hwi_max$conf_high[1],averaged_hwi_max$conf_high[3],averaged_hwi_max$conf_high[4]),
  variable = c(averaged_hwi_max$Significance[2],averaged_hwi_max$Significance[5],averaged_hwi_max$Significance[6],averaged_hwi_max$Significance[1],averaged_hwi_max$Significance[3],averaged_hwi_max$Significance[4])
)
cust_hwi_max$term <- factor(cust_hwi_max$term, cust_hwi_max$term)

cust_wl_max <- data.frame(
  term = c("Habitat density", "Max elevation",'Trophic level','Aerial lifestyle','Migration','Body mass'),
  estimate = c(averaged_wl_max$Estimate[2],averaged_wl_max$Estimate[5],averaged_wl_max$Estimate[6],averaged_wl_max$Estimate[1],averaged_wl_max$Estimate[3],averaged_wl_max$Estimate[4]),
  conf.low = c(averaged_wl_max$conf_low[2],averaged_wl_max$conf_low[5],averaged_wl_max$conf_low[6],averaged_wl_max$conf_low[1],averaged_wl_max$conf_low[3],averaged_wl_max$conf_low[4]),
  conf.high = c(averaged_wl_max$conf_high[2],averaged_wl_max$conf_high[5],averaged_wl_max$conf_high[6],averaged_wl_max$conf_high[1],averaged_wl_max$conf_high[3],averaged_wl_max$conf_high[4]),
  variable = c(averaged_wl_max$Significance[2],averaged_wl_max$Significance[5],averaged_wl_max$Significance[6],averaged_wl_max$Significance[1],averaged_wl_max$Significance[3],averaged_wl_max$Significance[4])
)
cust_wl_max$term <- factor(cust_wl_max$term, cust_wl_max$term)

#Plot the forest plots
p1 <- ggcoef(cust_hwi_max, exponentiate = FALSE, mapping = aes(x = estimate, y = term, colour = variable),size=3)
p2 <- ggcoef(cust_wl_max, exponentiate = FALSE, mapping = aes(x = estimate, y = term, colour = variable),size=3)
p1 + theme_pubr()
p2 + theme_pubr()

# Plots for high-certianty dataset using mean elevation
cust_hwi_mean <- data.frame(
  term = c("Habitat density", "Mean elevation",'Trophic level','Aerial lifestyle','Migration','Body mass'),
  estimate = c(averaged_hwi_mean$Estimate[2],averaged_hwi_mean$Estimate[5],averaged_hwi_mean$Estimate[6],averaged_hwi_mean$Estimate[1],averaged_hwi_mean$Estimate[3],averaged_hwi_mean$Estimate[4]),
  conf.low = c(averaged_hwi_mean$conf_low[2],averaged_hwi_mean$conf_low[5],averaged_hwi_mean$conf_low[6],averaged_hwi_mean$conf_low[1],averaged_hwi_mean$conf_low[3],averaged_hwi_mean$conf_low[4]),
  conf.high = c(averaged_hwi_mean$conf_high[2],averaged_hwi_mean$conf_high[5],averaged_hwi_mean$conf_high[6],averaged_hwi_mean$conf_high[1],averaged_hwi_mean$conf_high[3],averaged_hwi_mean$conf_high[4]),
  variable = c(averaged_hwi_mean$Significance[2],averaged_hwi_mean$Significance[5],averaged_hwi_mean$Significance[6],averaged_hwi_mean$Significance[1],averaged_hwi_mean$Significance[3],averaged_hwi_mean$Significance[4])
)
cust_hwi_mean$term <- factor(cust_hwi_mean$term, cust_hwi_mean$term)

cust_wl_mean <- data.frame(
  term = c("Habitat density", "Mean elevation",'Trophic level','Aerial lifestyle','Migration','Body mass'),
  estimate = c(averaged_wl_mean$Estimate[2],averaged_wl_mean$Estimate[5],averaged_wl_mean$Estimate[6],averaged_wl_mean$Estimate[1],averaged_wl_mean$Estimate[3],averaged_wl_mean$Estimate[4]),
  conf.low = c(averaged_wl_mean$conf_low[2],averaged_wl_mean$conf_low[5],averaged_wl_mean$conf_low[6],averaged_wl_mean$conf_low[1],averaged_wl_mean$conf_low[3],averaged_wl_mean$conf_low[4]),
  conf.high = c(averaged_wl_mean$conf_high[2],averaged_wl_mean$conf_high[5],averaged_wl_mean$conf_high[6],averaged_wl_mean$conf_high[1],averaged_wl_mean$conf_high[3],averaged_wl_mean$conf_high[4]),
  variable = c(averaged_wl_mean$Significance[2],averaged_wl_mean$Significance[5],averaged_wl_mean$Significance[6],averaged_wl_mean$Significance[1],averaged_wl_mean$Significance[3],averaged_wl_mean$Significance[4])
)
cust_wl_mean$term <- factor(cust_wl_mean$term, cust_wl_mean$term)

#Plot the forest plots
p3 <- ggcoef(cust_hwi_mean, exponentiate = FALSE, mapping = aes(x = estimate, y = term, colour = variable),size=3)
p4 <- ggcoef(cust_wl_mean, exponentiate = FALSE, mapping = aes(x = estimate, y = term, colour = variable),size=3)
p3 + theme_pubr()
p4 + theme_pubr()

##### (3) Vioplots (Figure 5 & 6) #####

# Vioplot for HWI and aerial lifestyle index
vioplot(hand_wing_index~al_binary,data=passerine, col='lightgray',
        xlab='Aerial lifestyle index',
        ylab='HWI')

# Vioplot for HWI and migration
vioplot(hand_wing_index~migration_binary,data=passerine,col='lightgray',
        xlab='Migration',
        ylab='Hand wing index')

# Vioplot for HWI and habitat density
vioplot(hand_wing_index~habitat_binary,data=passerine,col='lightgray',
        xlab='Habitat density',
        ylab='Hand wing index')

# Vioplot for HWI and trophic level
vioplot(hand_wing_index~trophic_level_num,data=passerine,col='lightgray',
        xlab='Trophic level',
        ylab='Hand wing index')
