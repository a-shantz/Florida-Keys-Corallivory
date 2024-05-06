#Keys Corallivory Final Code
library(readxl)
library(readr)
library(tidyverse)
library(DHARMa)
library(fitdistrplus)
library(sjPlot)
library(effects)
library(lme4)
library(glmmTMB)
library(emmeans)
library(magrittr)
library(ggeffects)
library(car)
library(ggthemes)
library(vegan)
library(patchwork)


###########
#Coral Data
###########

#Import the data from your source directory
coral_dat <- read_csv("~/Keys_Coral_Comparison.csv")

#Aggregate total live coral area per transect
coral_dat$dummie <- 1

total_live_coral <- coral_dat %>%
  group_by(Site, Transect, Year) %>%
  summarise(Total_live_area = sum(Colony_area),
            Colony_count = sum(dummie)) %>%
  ungroup()

#aggregate live area for each group on the transect
total_group_area <- coral_dat %>%
  group_by(Site, Transect, Coral_group, Year) %>%
  summarise(Group_live_coral = sum(Colony_area)) %>%
  ungroup() %>%
  complete(Coral_group,
           nesting(Site, Transect, Year),
           fill = list(Group_live_coral = 0))

#Examine % cover for species/sites/etc
total_group_area %>%
  group_by(Coral_group, Year) %>%
  summarise(Average_Cover = mean(Group_live_coral/5000),
            SE = plotrix::std.error(Group_live_coral/5000)) %>%
  print(n = 36)

#Merge the two datasets
prop_dat <- merge(total_group_area, total_live_coral, by = c("Site", "Transect", "Year"))
prop_dat$Community_prop <- prop_dat$Group_live_coral/prop_dat$Total_live_area 
prop_dat$Site <- factor(prop_dat$Site)
levels(prop_dat$Site)

prop_dat <- prop_dat %>%
  mutate(Perc_c_unit = round(Group_live_coral/5000, 2))

#Drop one small Manicina
prop_dat <- subset(prop_dat, 
                   Coral_group != "Manicina areolata")

#Test for changes in total percent cover/live area by year IGNORING groups
coral_mod <- lm(log(Total_live_area+0.01) ~ Site * Year, data = total_live_coral) 

simulateResiduals(fittedModel = coral_mod, n = 250, plot = T)
car::Anova(coral_mod)
summary(coral_mod)


#Wilcoxon rank tests
#Excluded groups are: ACER (0 obs in 2010), CNAT (5 obs in 2010/1 in 2022),
#MMEA (3/4 obs), SBOU (0/2), MARE (1/0), MYSP (3/0)
APAL <- subset(prop_dat, Coral_group == "Acropora palmata") 
DSTO <- subset(prop_dat, Coral_group == "Dichocoenia stokesii")
EFAS <- subset(prop_dat, Coral_group == "Eusmilia fastigiata")
AGAR <- subset(prop_dat, Coral_group == "Agaricia spp.")
BPSP <- subset(prop_dat, Coral_group == "Branched Porites spp.")
DISP <- subset(prop_dat, Coral_group == "Diploria spp.")
MASP <- subset(prop_dat, Coral_group == "Madracis spp.")
MCAV <- subset(prop_dat, Coral_group == "Montastraea cavernosa")
ORSP <- subset(prop_dat, Coral_group == "Orbicella spp.")
PAST <- subset(prop_dat, Coral_group == "Porites astreoides")
SISP <- subset(prop_dat, Coral_group == "Siderastrea spp.")
SINT <- subset(prop_dat, Coral_group == "Stephanocoenia intersepta")

#one sided wilcox tests to see if changes were significant
wilcox.test(Perc_c_unit ~ Year, paired =F,
            alternative = "two.sided", data = APAL) 
wilcox.test(Perc_c_unit ~ Year, paired =F,
            alternative = "two.sided", data = DSTO) 
wilcox.test(Perc_c_unit ~ Year, paired =F,
            alternative = "two.sided", data = EFAS) 
wilcox.test(Perc_c_unit ~ Year, paired =F,
            alternative = "two.sided", data = AGAR) 
wilcox.test(Perc_c_unit ~ Year, paired =F,
            alternative = "two.sided", data = BPSP)
wilcox.test(Perc_c_unit ~ Year, paired =F,
            alternative = "two.sided", data = DISP) 
wilcox.test(Perc_c_unit ~ Year, paired =F,
            alternative = "two.sided", data = MASP) 
wilcox.test(Perc_c_unit ~ Year, paired =F,
            alternative = "two.sided", data = MCAV) 
wilcox.test(Perc_c_unit ~ Year, paired =F,
            alternative = "two.sided", data = ORSP) 
wilcox.test(Perc_c_unit ~ Year, paired =F,
            alternative = "two.sided", data = PAST) 
wilcox.test(Perc_c_unit ~ Year, paired =F,
            alternative = "two.sided", data = SISP) 
wilcox.test(Perc_c_unit ~ Year, paired =F,
            alternative = "two.sided", data = SINT) 

#House keeping
rm(APAL, AGAR, BPSP, EFAS, DISP, DSTO, MASP, MCAV, mod, ORSP, PAST, SINT, SISP)


##################
#NMDS - using proportion of data (relative contribution) to catch changes in the composition rather than absolute changes
prop_dat_wide <- prop_dat %>%
  dplyr::select(-c(Total_live_area, Perc_c_unit, Group_live_coral, Colony_count)) %>%
  pivot_wider(names_from = Coral_group, values_from = Community_prop) %>%
  arrange(Year, Site)


nmds_df <- prop_dat_wide[,c(4:20)]
Year <- factor(unlist(prop_dat_wide[,3]))
Site <- prop_dat_wide[,1]

coral_nmds <- metaMDS(nmds_df, distance = 'bray', k = 2) #Stress generally around 19
stressplot(coral_nmds)
coral_nmds
coral_sp_fit <- envfit(coral_nmds, nmds_df, permutations = 999)
head(coral_sp_fit)


plot(coral_nmds, type = "n")
ordihull(coral_nmds, groups = Year, draw = 'polygon', col = "grey90")
points(coral_nmds, display = "sites", pch = c(21, 16)[as.factor(Year)], col = "black")
plot(coral_sp_fit, p.max = 0.05, col = "black", cex = 0.7)


#Permanova
#Dispersion
dis <- betadisper(dist(nmds_df, method = "euclidean"), group = Year)
dis
TukeyHSD(dis, which = "group", ordered = FALSE,
         conf.level = 0.95)

#Permanova with Site as a nesting factor
perm_df <- prop_dat_wide[, c(1,3)]
perm <- how(nperm = 999)
setBlocks(perm) <- with(perm_df, Site)
adonis2(nmds_df ~ Year, data = perm_df, permutations = perm, method = "bray", 
        by = "margin")

#House keeping
rm(dis, perm, perm_df, nmds_df, prop_dat_wide, coral_nmds, coral_sp_fit, Site, Year)



############
#Parrotfish#
############

#load the dataset
pd_all <- read_csv("Final products/Keys_Parrotfish_Comparison.csv") #Update to your directory

#Drop non-corallivores using data from Burkepile et al 2019 (Coral reefs):
  #They show that coeruleus, chrysopterum, and aurofrenatum never fed on corals during follows...
  #but aurofrenatum were caught biting sids on video.
  #We will drop Coeruleus and chrysopterum

#Make species a factor
pd_all$Species <- factor(pd_all$Species) 
levels(pd_all$Species)

#Add zeros for transects with no observations
pd_all <- pd_all %>%
  complete(Site, Transect, Year, Species, fill =list(Count = 0, biomass_g = 0))

#Restrict to fish > 15 cm to match Burkepile 2012
corallivores <- pd_all %>%
  subset(Size > 15 & Species %in% c("viride", "vetula", "aurofrenatum","coelestinus", "guacamaia","iseri",
                                    "rubripinne", "taeniopterus"))

#Reset factor levels
corallivores$Species <- factor(corallivores$Species)

#Add zeros for transects with no observed fish
corallivores <- corallivores %>%
  complete(Site, Transect, Year, Species, fill =list(Count = 0, biomass_g = 0))

#Rename
corallivores$Species <- factor(corallivores$Species, levels = c("viride", "vetula", "taeniopterus","rubripinne", 
                                                                "iseri", "guacamaia", "coelestinus", "aurofrenatum"),
                               labels = c("Sp. viride", "Sc. vetula", "Sc. taeniopterus","Sp. rubripinne", 
                                          "Sc. iseri", "Sc. guacamaia", "Sc. coelestinus", "Sp. aurofrenatum"))

#Sum species level counts and biomass by site, transect, and year
pdat <- corallivores %>%
  group_by(Site, Transect, Year, Species) %>%
  summarise(Transect_total = sum(Count),
            Transect_bm = sum(biomass_g)) %>%
  ungroup()

#Total counts and biomass by transect/year
pdat2 <- pdat %>%
  group_by(Site, Transect, Year) %>%
  summarise(Parrot_count = sum(Transect_total),
            Parrot_biomass = round(sum(Transect_bm),0)) %>%
  ungroup()

#Has parrotfish abundance or biomass changed between years?
#Abundance/density
pd_mod <- lm(log(Parrot_count+1) ~ Year * Site, data = pdat2) 
simulateResiduals(fittedModel = pd_mod, n = 1000, plot = T)
car::Anova(pd_mod)

#Biomass
pb_mod <- lm(sqrt(Parrot_biomass) ~ Year * Site, data = pdat2) 
simulateResiduals(fittedModel = pb_mod, n = 250, plot = T)
car::Anova(pb_mod)

EMM <- emmeans(pb_mod, ~ "Year|Site") #
summary(pairs(EMM), by = NULL) 

#Token analysis of parrtofish changes with coral cover - with 8 reps and average values don't expect much!
avg_coral_area <- total_live_coral %>%
  group_by(Site, Year) %>%
  summarise(Average_live = mean(Total_live_area),
            live_se = plotrix::std.error(Total_live_area)) %>%
  ungroup()

coral_2010 <- subset(avg_coral_area, Year == 2010)
coral_2022 <- subset(avg_coral_area, Year == 2022)
coral_2010$change <- coral_2022$Average_live - coral_2010$Average_live

avg_parrots <- pdat2 %>%
  group_by(Site, Year) %>%
  summarize(Avg_biomass = mean(Parrot_biomass),
            Avg_count = mean(Parrot_count)) %>%
  ungroup()

fish_2010 <- subset(avg_parrots, Year == 2010)
fish_2022 <- subset(avg_parrots, Year == 2022)

coral_2010$parrot_bm_change <- fish_2022$Avg_biomass - fish_2010$Avg_biomass
coral_2010$parrot_count_change <- fish_2022$Avg_count - fish_2010$Avg_count

mod_bm <- lm(parrot_bm_change ~ change, weights = live_se, data = coral_2010)
plot(mod_bm) #lol well, we have 8 replicates...
Anova(mod_bm)

mod_den <- lm(parrot_count_change ~ change, weights = live_se, data = coral_2010)
plot(mod_den) #Same
Anova(mod_den)

#Clean up this mess
rm(EMM, pd_all, pdat, pb_mod, pd_mod, mod_bm, mod_den, fish_2010, fish_2022, coral_2010, coral_2022)

#############
#Corallivory#
#############

#Changes in corallivory patterns#
#Make a bites per meter squared live tissue column
cd_full <- coral_dat %>%
  mutate(BPM = (Bites/Colony_area)*10000)

#For proportional bite analyses
coral_transect_summary <- cd_full %>%
  group_by(Year, Site, Transect, Coral_group) %>%
  reframe(Total_group_area_cm2 = sum(Colony_area),
          Total_group_area_m2 = sum(Colony_area)/10000,
          Total_group_cover_perc = sum(Colony_area)/500000,
          Num_colon = sum(dummie),
          Num_bitten = sum(Bitten),
          Perc_bitten = sum(Bitten)/sum(dummie),
          Total_bites = sum(Bites)) %>%
  ungroup()

bites_explore <- left_join(cd_full, coral_transect_summary, by = c("Year", "Site", "Transect", "Coral_group"), 
                           relationship = "many-to-many") %>%
  mutate(Year = factor(Year))

pdat2$Year <- factor(pdat2$Year)

bitesNfish <- left_join(bites_explore, pdat2, by = c("Year", "Site", "Transect"), relationship = "many-to-many")

#Add a rounded BPM column
bitesNfish$Rounded_BPM <- round(bitesNfish$BPM,0)


#Does probability of being bitten differ between 2010 and 2022?
#First lets take a look at our counts
bitesNfish %>%
  group_by(Coral_group, Year) %>% #Site
  summarise(Num_bitten = sum(Bitten),
            Num_colonies = sum(dummie)) %>%
  print(n = 178)


#Excluding corals where no bites observed in either of the years or <= 8 total observations were made in either year.
bnf_limited <- bitesNfish %>%
  filter(Coral_group %in% c("Agaricia spp.", "Branched Porites spp.", "Madracis spp.", "Montastraea cavernosa",
                            "Orbicella spp.", "Porites astreoides", "Siderastrea spp.", "Stephanocoenia intersepta"))


##Probability of being bitten
bite_prob_mod <- glmmTMB(Bitten ~ Year * Coral_group + (1|Site:Transect),
                         family = binomial(), 
                         control = glmmTMBControl(optimizer = nlminb, optCtrl=list(iter.max=1e3,eval.max=1e3)),
                         data = bnf_limited)
simulateResiduals(bite_prob_mod, n = 250, plot = T)
summary(bite_prob_mod)
Anova(bite_prob_mod)

emmeans(bite_prob_mod, pairwise ~ "Year|Coral_group")

mod_emm <- emmeans(bite_prob_mod, ~Year * Coral_group, p.adjust = "mvt")
pairs(mod_emm, simple = "Year") 

#Changes in bite intensity
##ZINB model
bpm_zimod <- glmmTMB(Rounded_BPM ~ Year * Coral_group + (1|Site:Transect),
                     family = nbinom1(),
                     #control=glmmTMBControl(optimizer=optim,
                      #                      optArgs=list(method="BFGS")),
                     control = glmmTMBControl(optimizer = nlminb, 
                                              optCtrl=list(iter.max=1e3,eval.max=1e3)), 
                     ziformula = ~Year * Coral_group, 
                     data = bnf_limited)

#Check assumptions
plot(simulateResiduals(bpm_zimod, n = 250))
testDispersion(simulateResiduals(bpm_zimod, n = 250))
testZeroInflation(simulateResiduals(bpm_zimod, n = 250))
testOutliers(simulateResiduals(bpm_zimod, n = 250)) #ns

Anova(bpm_zimod)
summary(bpm_zimod)

#Who differed in main model?
emmeans(bpm_zimod, pairwise ~ "Year|Coral_group")


#convert to odds, calculate robust error and export table
#Variance-covariance matrix for robust standard errors
a <- vcov(bpm_zimod, full = FALSE) 

tab_model(bpm_zimod, 
          p.style = "numeric_stars", vcov.fun = a, show.se = T, transform = "exp")


####How does coral cover impact bites
total_live_coral$Year <- factor(total_live_coral$Year)

df <- left_join(bnf_limited, total_live_coral, by = c("Site", "Year", "Transect"),
                relationship = 'many-to-many')

df <- df %>%
  filter(Coral_group %in% c("Agaricia spp.", "Branched Porites spp.", "Orbicella spp.", "Porites astreoides",
                            "Siderastrea spp."))

#Rename stuff cause this is a mess
df$Total_percent_cover <- df$Total_live_area/50000
#Put percent cover on same scale for groups
df$Group_percent_cover <- df$Total_group_area_cm2/50000


#Model selection 
#Lets make a lot of dumb models!
SPGxG <- glmmTMB(Bitten ~ scale(Group_percent_cover) * Coral_group + 
                   Site +
                   Parrot_count +
                   (1|Year:Site:Transect), 
                 family = binomial(),
                 control = glmmTMBControl(optimizer = nlminb, 
                                          optCtrl=list(iter.max=1e3,eval.max=1e3)),
                 data = df) 

SGxG <- glmmTMB(Bitten ~ scale(Group_percent_cover) * Coral_group + 
                  Site +
                  (1|Year:Site:Transect), 
                family = binomial(),
                control = glmmTMBControl(optimizer = nlminb, 
                                         optCtrl=list(iter.max=1e3,eval.max=1e3)),
                data = df) 

PGxG <- glmmTMB(Bitten ~ scale(Group_percent_cover) * Coral_group + 
                  Parrot_count +
                  (1|Year:Site:Transect), 
                family = binomial(),
                control = glmmTMBControl(optimizer = nlminb, 
                                         optCtrl=list(iter.max=1e3,eval.max=1e3)),
                data = df) 
GxG <- glmmTMB(Bitten ~ scale(Group_percent_cover) * Coral_group + 
                 (1|Year:Site:Transect), 
               family = binomial(),
               control = glmmTMBControl(optimizer = nlminb, 
                                        optCtrl=list(iter.max=1e3,eval.max=1e3)),
               data = df) 

G <- glmmTMB(Bitten ~ Coral_group + 
               (1|Year:Site:Transect), 
             family = binomial(),
             control = glmmTMBControl(optimizer = nlminb, 
                                      optCtrl=list(iter.max=1e3,eval.max=1e3)),
             data = df) 

SPGxT <- glmmTMB(Bitten ~ scale(Total_percent_cover, center = TRUE) * Coral_group +
                   Site +
                   Parrot_count +
                   (1|Year:Site:Transect), 
                 family = binomial(),
                 control = glmmTMBControl(optimizer = nlminb, 
                                          optCtrl=list(iter.max=1e3,eval.max=1e3)),
                 data = df) 

SGxT <- glmmTMB(Bitten ~ scale(Total_percent_cover, center = TRUE) * Coral_group +
                  Site +
                  (1|Year:Site:Transect), 
                family = binomial(),
                control = glmmTMBControl(optimizer = nlminb, 
                                         optCtrl=list(iter.max=1e3,eval.max=1e3)),
                data = df) 

PGxT <- glmmTMB(Bitten ~ scale(Total_percent_cover, center = TRUE) * Coral_group +
                  Parrot_count +
                  (1|Year:Site:Transect), 
                family = binomial(),
                control = glmmTMBControl(optimizer = nlminb, 
                                         optCtrl=list(iter.max=1e3,eval.max=1e3)),
                data = df) 

GxT <- glmmTMB(Bitten ~ scale(Total_percent_cover, center = TRUE) * Coral_group +
                 (1|Year:Site:Transect), 
               family = binomial(),
               control = glmmTMBControl(optimizer = nlminb, 
                                        optCtrl=list(iter.max=1e3,eval.max=1e3)),
               data = df) 

SPG <- glmmTMB(Bitten ~ Coral_group + 
                 Site +
                 Parrot_count +
                 (1|Year:Site:Transect), 
               family = binomial(),
               control = glmmTMBControl(optimizer = nlminb, 
                                        optCtrl=list(iter.max=1e3,eval.max=1e3)),
               data = df) 

SG <- glmmTMB(Bitten ~ Coral_group + 
                Site +
                (1|Year:Site:Transect), 
              family = binomial(),
              control = glmmTMBControl(optimizer = nlminb, 
                                       optCtrl=list(iter.max=1e3,eval.max=1e3)),
              data = df) 

PG <- glmmTMB(Bitten ~ Coral_group + 
                Parrot_count +
                (1|Year:Site:Transect), 
              family = binomial(),
              control = glmmTMBControl(optimizer = nlminb, 
                                       optCtrl=list(iter.max=1e3,eval.max=1e3)),
              data = df) 

SP <- glmmTMB(Bitten ~
                Site +
                Parrot_count +
                (1|Year:Site:Transect), 
              family = binomial(),
              control = glmmTMBControl(optimizer = nlminb, 
                                       optCtrl=list(iter.max=1e3,eval.max=1e3)),
              data = df) 

S <- glmmTMB(Bitten ~ Site +
               (1|Year:Site:Transect), 
             family = binomial(),
             control = glmmTMBControl(optimizer = nlminb, 
                                      optCtrl=list(iter.max=1e3,eval.max=1e3)),
             data = df) 

P <- glmmTMB(Bitten ~
               Parrot_count +
               (1|Year:Site:Transect), 
             family = binomial(),
             control = glmmTMBControl(optimizer = nlminb, 
                                      optCtrl=list(iter.max=1e3,eval.max=1e3)),
             data = df) 


null <- glmmTMB(Bitten ~ 1 +
                  (1|Year:Site:Transect), 
                family = binomial(),
                control = glmmTMBControl(optimizer = nlminb, 
                                         optCtrl=list(iter.max=1e3,eval.max=1e3)),
                data = df) 



cand.models <- list(SPGxG, SGxG, PGxG, SPGxT, SGxT, PGxT, SPG, GxG, GxT, G, SP, SG, PG, P, S, null) 
names <- c('SPGxG', 'SGxG', 'PGxG', 'SPGxT', 'SGxT', 'PGxT', 'SPG', 'GxG', 'GxT', 'G', 'SP',
           'SG', 'PG', 'P', 'S', 'null') 
AICcmodavg::aictab(cand.set = cand.models, modnames = names, sort = T)

performance::r2_nakagawa(SPGxG)

#Check for multicollinearity in model - you need to remove interaction terms first because they will inflate VIF 
SPGxG_VIF_check <- glmmTMB(Bitten ~ scale(Group_percent_cover) + Coral_group + 
                             Site +
                             (1|Year:Site:Transect), 
                           family = binomial(),
                           control = glmmTMBControl(optimizer = nlminb, 
                                                    optCtrl=list(iter.max=1e3,eval.max=1e3)),
                           data = df) 
performance::check_collinearity(SPGxG_VIF_check, component = "count")  
summary(SPGxG)
a <- vcov(SPGxG, full = FALSE) 
tab_model(SPGxG, show.p = FALSE, vcov.fun = a, show.se = T,
          show.df = TRUE, transform = "exp", digits = 4) #"exp"


#Clean up
rm(SPGxG, SGxG, PGxG, SPGxT, SGxT, PGxT, SPG, GxG, GxT, G, SP, SG, PG, P, S, null, SPGxG_VIF_check) 

#Bite intensity models
#These take a lot longer. Get them started and go make yourself a drink...
SPGxG_i <- glmmTMB(Rounded_BPM ~ scale(Group_percent_cover) * Coral_group + 
                     Parrot_count +
                     Site+
                     (1|Year:Site:Transect), 
                   family = nbinom2(), ziformula = ~Colony_area * Coral_group,
                   control = glmmTMBControl(optimizer = nlminb, 
                                            optCtrl=list(iter.max=1e3,eval.max=1e3)),
                   data = df) 

SGxG_i <- glmmTMB(Rounded_BPM ~ scale(Group_percent_cover) * Coral_group + 
                    Site+
                    (1|Year:Site:Transect), 
                  family = nbinom2(), ziformula = ~Colony_area * Coral_group,
                  control = glmmTMBControl(optimizer = nlminb, 
                                           optCtrl=list(iter.max=1e3,eval.max=1e3)),
                  data = df) 

PGxG_i <- glmmTMB(Rounded_BPM ~ scale(Group_percent_cover, center = TRUE) * Coral_group +
                    Parrot_count +
                    (1|Year:Site:Transect), 
                  family = nbinom2(), ziformula = ~Colony_area * Coral_group,
                  control = glmmTMBControl(optimizer = nlminb, 
                                           optCtrl=list(iter.max=1e3,eval.max=1e3)),
                  data = df) 

GxG_i <- glmmTMB(Rounded_BPM ~ scale(Group_percent_cover, center = TRUE) * Coral_group +
                   (1|Year:Site:Transect), 
                 family = nbinom2(), 
                 ziformula = ~Colony_area * Coral_group,
                 control = glmmTMBControl(optimizer = nlminb, 
                                          optCtrl=list(iter.max=1e3,eval.max=1e3)),
                 data = df) 

G_i <- glmmTMB(Rounded_BPM ~ Coral_group +
                 (1|Year:Site:Transect), 
               family = nbinom2(), ziformula = ~Colony_area * Coral_group,
               control = glmmTMBControl(optimizer = nlminb, 
                                        optCtrl=list(iter.max=1e3,eval.max=1e3)),
               data = df) 


SPGxT_i <- glmmTMB(Rounded_BPM ~ scale(Total_percent_cover, center = TRUE) * Coral_group +
                     Parrot_count +
                     Site +
                     (1|Year:Site:Transect), 
                   family = nbinom2(), ziformula = ~Colony_area * Coral_group,
                   control = glmmTMBControl(optimizer = nlminb, 
                                            optCtrl=list(iter.max=1e3,eval.max=1e3)),
                   data = df)  

SGxT_i <- glmmTMB(Rounded_BPM ~ scale(Total_percent_cover, center = TRUE) * Coral_group +
                    Site +
                    (1|Year:Site:Transect), 
                  family = nbinom2(), ziformula = ~Colony_area * Coral_group,
                  control = glmmTMBControl(optimizer=optim,
                                           optArgs=list(method="BFGS")),
                  data = df)  

PGxT_i <-glmmTMB(Rounded_BPM ~ scale(Total_percent_cover, center = TRUE) * Coral_group +
                   Parrot_count +
                   (1|Year:Site:Transect), 
                 family = nbinom2(), ziformula = ~Colony_area * Coral_group,
                 control = glmmTMBControl(optimizer = nlminb, 
                                          optCtrl=list(iter.max=1e3,eval.max=1e3)),
                 data = df)  

GxT_i <- glmmTMB(Rounded_BPM ~ scale(Total_percent_cover, center = TRUE) * Coral_group +
                   (1|Year:Site:Transect), 
                 family = nbinom2(), 
                 ziformula = ~Colony_area * Coral_group,
                 control = glmmTMBControl(optimizer=optim,
                                         optArgs=list(method="BFGS")),
                 data = df) 

SPG_i <- glmmTMB(Rounded_BPM ~ Site + Parrot_count + Coral_group +
                   (1|Year:Site:Transect), 
                 family = nbinom2(), 
                 ziformula = ~Colony_area * Coral_group,
                 control = glmmTMBControl(optimizer = nlminb, 
                                          optCtrl=list(iter.max=1e3,eval.max=1e3)),
                 data = df) 


SG_i <- glmmTMB(Rounded_BPM ~ Coral_group +
                  Site +
                  (1|Year:Site:Transect), 
                family = nbinom2(), ziformula = ~Colony_area * Coral_group,
                control = glmmTMBControl(optimizer = nlminb, 
                                         optCtrl=list(iter.max=1e3,eval.max=1e3)),
                data = df) 


PG_i <- glmmTMB(Rounded_BPM ~ Coral_group +
                  Parrot_count +
                  (1|Year:Site:Transect), 
                family = nbinom2(), ziformula = ~Colony_area * Coral_group,
                control = glmmTMBControl(optimizer = nlminb, 
                                         optCtrl=list(iter.max=1e3,eval.max=1e3)),
                data = df) 


SP_i <- glmmTMB(Rounded_BPM ~ Parrot_count +
                  Site +
                  (1|Year:Site:Transect), 
                family = nbinom1(), ziformula = ~Colony_area * Coral_group,
                control = glmmTMBControl(optimizer = nlminb, 
                                         optCtrl=list(iter.max=1e3,eval.max=1e3)),
                data = df) 


S_i <- glmmTMB(Rounded_BPM ~ Site + 
                 (1|Year:Site:Transect), 
               family = nbinom2(), ziformula = ~Colony_area * Coral_group,
               control = glmmTMBControl(optimizer=optim,
                                        optArgs=list(method="BFGS")),
               
               data = df) 

P_i <- glmmTMB(Rounded_BPM ~ 
                 Parrot_count +
                 (1|Year:Site:Transect), 
               family = nbinom2(), ziformula = ~Colony_area * Coral_group,
               control = glmmTMBControl(optimizer = nlminb, 
                                      optCtrl=list(iter.max=1e3,eval.max=1e3)),
               data = df) 


null_i <- glmmTMB(Rounded_BPM ~ 1 +
                    (1|Year:Site:Transect), 
                  family = nbinom2(), 
                  ziformula = ~Colony_area * Coral_group,
                  control = glmmTMBControl(optimizer=optim,
                                           optArgs=list(method="BFGS")),
                  data = df) 

cand.models <- list(SPGxG_i, SGxG_i, PGxG_i, SPGxT_i, SGxT_i, PGxT_i, 
                    SPG_i, GxG_i, GxT_i, G_i, SP_i, SG_i, PG_i, P_i, S_i, null_i) 
names <- c('SPGxG', 'SGxG', 'PGxG', 'SPGxT', 'SGxT', 'PGxT', 'SPG', 'GxG', 'GxT', 'G', 'SP',
           'SG', 'PG', 'P', 'S', 'null') 
AICcmodavg::aictab(cand.set = cand.models, modnames = names, sort = T)

performance::r2_nakagawa(SGxG_i)

#Check for multicollinearity in model without interaction terms again
SGxG_i_VIF_check <- glmmTMB(Rounded_BPM ~ scale(Group_percent_cover) + Coral_group + 
                              Site +
                              (1|Year:Site:Transect), 
                            family = nbinom2(), ziformula = ~Colony_area + Coral_group,
                            control = glmmTMBControl(optimizer = nlminb, 
                                                     optCtrl=list(iter.max=1e3,eval.max=1e3)),
                            data = df) 
performance::check_collinearity(SGxG_i_VIF_check, component = "count") #Use component = "count" and "zi" to check both parts 

a <- vcov(SGxG_i, full = FALSE) 
tab_model(SGxG_i, show.p = FALSE, vcov.fun = a, show.se = T,
          show.df = TRUE, transform = "exp", digits = 4) #"exp"


#Clean up
rm(SPGxG_i, SGxG_i, PGxG_i, SPGxT_i, SGxT_i, PGxT_i, SGxG_i_VIF_check,
     SPG_i, GxG_i, GxT_i, G_i, SP_i, SG_i, PG_i, P_i, S_i, null_i) 











