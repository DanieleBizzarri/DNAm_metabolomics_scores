#Set your working directory
wd<-setwd("~/Desktop/LUMC/Metabocpg/git")

###############
## Libraries ##
library(foreach)
library(ggplot2)
library(haven)
library(survival)
library(survminer)

###################
## Load datasets ##

# Load the DNAm metabolomics scores
DNAm_metab_models <- readRDS("~/researchdrive/dbizzarri/metabo_cpgs/Full_NTR_data/unscaled_cpgs/All_models_coefficients/DNAm_metabolomics_models.rds")
#Calculated scores
DNAm_metab_scores<-readRDS(paste0(wd,"/Results/DNAm_metabolomics_features.rds"))

# Load the Horvath scores. These scores must be obtained by mailing to Lu et al: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6366976/
library(readr)
DNAm_clocks<- data.frame(read_csv("~/researchdrive/LLS_private/DNAm_scores_Horvarth/DNAm_predictors_in_LLS_PARTOFFS.csv"))
rownames(DNAm_clocks)<-DNAm_clocks$uuid

#Load covariates (at least age and sex)
cov<-readRDS(paste0(wd,"/data/covariates.rds"))

# Load the EpiScores
EpiScores <- read.table("~/researchdrive/dbizzarri/metabo_cpgs/EpiScores/bage_predictor/episcore_projections.tsv", sep="\t", 
                        row.names = 1, header = T, stringsAsFactors = F)

# Load bAge
bAge <- read.table("~/researchdrive/dbizzarri/metabo_cpgs/EpiScores/bage_predictor/bage_predictions.tsv", sep="\t", 
                   row.names = 1, header = T, stringsAsFactors = F)

#Multivariate mortality models trained in RS
# Model with DNAm MeetaboHealth and GrimAge
coxph_DNAm_MetaboHealth_GrimAge<-readRDS(paste0(wd,"/data/Stepwise_Cox_regression_DNAm_MetaboHealth_GrimAge.rds"))
# Model with DNAm met and GrimAge
coxph_DNAm_met_GrimAge<-readRDS(paste0(wd,"/data/Stepwise_Cox_regression_DNAm_metab_GrimAge.rds"))
# Model with DNAm met and GrimAge proteins
coxph_DNAm_met_prot<-readRDS(paste0(wd,"/data/Stepwise_Cox_regression_DNAm_prot_metab.rds"))
# Model with DNAm met, GrimAge proteins, and EpiScores
coxph_DNAm_met_prot_EpiScores<-readRDS(paste0(wd,"/data/Stepwise_Cox_regression_DNAm_prot_metab_EpiScores.rds"))

################################
## Obtaining the age residuals #
################################
Age_accel_DNAm_clocks<-DNAm_clocks[,c(2:4,20)]
colnames(Age_accel_DNAm_clocks)<-paste0("AgeAccel_", colnames(Age_accel_DNAm_clocks))
for(d in colnames(DNAm_clocks)[c(2:4,20)]){
  fitage<-lm(formula=paste0(d,"~","Age.x"), data=DNAm_clocks,na.action=NULL)
  Age_accel_DNAm_clocks[,paste0("AgeAccel_",d)]<-residuals(fitage)/sd(residuals(fitage))
}
DNAm_clocks<-cbind(DNAm_clocks, Age_accel_DNAm_clocks)
DNAm_clocks<-data.frame(DNAm_clocks[,c("uuid", "Horvath", "Hannum", "Levine", "DNAmGrimAge",
                                       "AgeAccel_Horvath", "AgeAccel_Hannum", "AgeAccel_Levine", "AgeAccel_DNAmGrimAge",
                                       "DNAmGDF_15","DNAmB2M", "DNAmCystatin_C", "DNAmTIMP_1","DNAmadm",        
                                       "DNAmpai_1", "DNAmleptin", "DNAmPACKYRS")])

#Age acceleration bAge
bAge$AgeAccel_bAge<-NA
fitage<-lm(formula="bAge~Age", data=bAge,na.action=NULL)
bAge$AgeAccel_bAge<-residuals(fitage)/sd(residuals(fitage))
bAge<-bAge[,c(1,2,8)]

#########################
## Zscaling the scores ##
#########################
DNAm_metab_scores<-scale(DNAm_metab_scores)
DNAm_clocks<-scale(DNAm_clocks[,-1])

EpiScores<-scale(EpiScores)
bAge<-scale(bAge)

#################################################
# Preparing data.frame with all DNAm surrogates #
#################################################
DNAm_scores<-merge(DNAm_clocks, DNAm_metab_scores, by=0)
rownames(DNAm_scores)<-DNAm_scores$Row.names
DNAm_scores<-merge(DNAm_scores[,-1], EpiScores, by=0)
rownames(DNAm_scores)<-DNAm_scores$Row.names
DNAm_scores<-merge(DNAm_scores[,-1], bAge, by=0)
rownames(DNAm_scores)<-DNAm_scores$Row.names
colnames(DNAm_scores)[1]<-"ID"

DNAm_mort<-merge(mort,DNAm_scores ,by="ID")
rownames(DNAm_mort)<- DNAm_mort$Id

#######################################
## Project the Biological Age scores ##
#######################################
DNAm_MetaboHealth_GrimAge <- data.frame(predict(coxph_DNAm_MetaboHealth_GrimAge, DNAm_mort, type="lp"))
rownames(DNAm_MetaboHealth_GrimAge)<-rownames(DNAm_mort)
colnames(DNAm_MetaboHealth_GrimAge)<-"DNAm_MetaboHealth_GrimAge"
DNAm_met_GrimAge <- data.frame(predict(coxph_DNAm_met_GrimAge, DNAm_mort, type="lp"))
rownames(DNAm_met_GrimAge)<-rownames(DNAm_mort)
colnames(DNAm_met_GrimAge)<-"DNAm_met_GrimAge"
DNAm_met_prot <- data.frame(predict(coxph_DNAm_met_prot, DNAm_mort, type="lp"))
rownames(DNAm_met_prot)<-rownames(DNAm_mort)
colnames(DNAm_met_prot)<-"DNAm_met_prot"
DNAm_met_prot_EpiScores <- data.frame(predict(coxph_DNAm_met_prot_EpiScores, DNAm_mort, type="lp"))
rownames(DNAm_met_prot_EpiScores)<-rownames(DNAm_mort)
colnames(DNAm_met_prot_EpiScores)<-"DNAm_met_prot_EpiScores"


