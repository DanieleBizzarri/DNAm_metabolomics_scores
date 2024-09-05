library(ggplot2)
library(foreach)
library(reshape2)
library(ggridges)

#Set your working directory
wd<-setwd("~/Desktop/LUMC/Metabocpg/git")

###############
## Load data ##
###############
#Load the RS Betas
load(paste0(wd,"betas.RData"))

#Load the DNAm metabolomics models
DNAm_metabolomics_models <- readRDS(paste0(wd, "/data/DNAm_metabolomics_models.rds"))

##################
## Missing cpgs ##
##################
#Missing cpgs not in the RS datasets
missing_cpgs<-setdiff(rownames(DNAm_metabolomics_models$models_betas)[-c(1,2)], rownames(beta_all))
#Find the amount of missingness per model
missing<-DNAm_metabolomics_models$models_betas[missing_cpgs,]
N_missing<-table(which(missing!=0,arr.ind=T)[,2])
#Amount of cpgs used per model
N_cpgs<-table(which(DNAm_metabolomics_models$models_betas!=0,arr.ind=T)[,2])
names(N_cpgs)<-c("DNAm_MetaboHealth",paste0("DNAm_",DNAm_metabolomics_models$metabolites_info$name))
N_cpgs<-data.frame(DNAm_scores=names(N_cpgs),N=N_cpgs)
N_cpgs$N_missing<-N_missing
N_cpgs<-N_cpgs[,c(1,3,4)]
#Combine all the numbers
N_cpgs2<-data.frame(DNAm_scores=c(N_cpgs$DNAm_scores,N_cpgs$DNAm_scores),
                    N=c(N_cpgs$N.Freq,N_cpgs$N_missing),
                    type=c(rep("N",65),rep("Nmiss",65)))
#Factorize
N_cpgs$DNAm_scores<-factor(N_cpgs$DNAm_scores,levels=N_cpgs$DNAm_scores[65:1])
N_cpgs2$DNAm_scores<-factor(N_cpgs2$DNAm_scores,levels=N_cpgs$DNAm_scores[65:1])
#Plot
p<-ggplot(data=N_cpgs2,aes(x=DNAm_scores,y=N, fill=type))+
  geom_bar(stat="identity")+scale_fill_manual(values=c("steelblue","red"))+
  geom_text(aes(label=N),hjust=-0.5,color="black",size=3.5)+theme_bw()+
  coord_flip()
ggsave("~/DNAm_models/Results/Images/Missing_cpgs.pdf", plot=p,
       width=10,height=8,units="in",device="pdf")

#Plot the coefficients of the missings
missing_melt<-melt(missing)
missing_melt<-missing_melt[-which(missing_melt$value==0),]
missing_melt$variable<-factor(missing_melt$variable,level=colnames(DNAm_metabolomics_models$models_betas))
p<-ggplot(missing_melt,aes(x=value,y=variable, fill=stat(x)))+
  geom_density_ridges_gradient(scale=3,rel_min_height=0.01)+scale_fill_viridis_c(option="C")+
  labs(title="Distributions of the coefficients of the missing cpgs in RS",
       x="models coefficients",y="models")
p

#Distributions of the median substituting the missings
Median_missing<-foreach(i=colnames(DNAm_metabolomics_models$models_betas),.combine="rbind")%do%{
  cpg<-rownames(missing)[which(missing[,i]!=0)]
  dat<-data.frame(CpGs=cpg,Median=DNAm_metabolomics_models$median_cpgs[cpg],model=i)
}
Median_missing$model<-factor(Median_missing$model,level=colnames(DNAm_metabolomics_models$models_betas))
p<-ggplot(Median_missing,aes(x=Median,y=model, fill=stat(x)))+
  geom_density_ridges_gradient(scale=3,rel_min_height=0.01)+scale_fill_viridis_c(option="C")+
  labs(title="Distributions of Median substituted to the missing cpgs in RS",
       x="Medians",y="models")

##############
## Imputing ##
##############
#Substitute the missing with zeros
#gold_missing<-data.frame(matrix(0, nrow=dim(beta_all)[2] , ncol= length(missing_cpgs)))
#Substitute the missing with median in BIOS
gold_missing <- DNAm_metabolomics_models$median_cpgs[missing_cpgs] #use all is better, but did not work as some were missing in gold
gold_missing <- as.numeric(gold_missing)
gold_missing<- t(replicate(ncol(beta_all), gold_missing))
rownames(gold_missing)<-colnames(beta_all)
colnames(gold_missing)<-missing_cpgs
gc()

#Create a unique data.frame
beta_all<-data.frame(t(beta_all),gold_missing,sentrix_position=cov$Position)
gc()
saveRDS(beta_all, file=paste0(wd,"/data/RS_betas_with_missing_to_median_in_BIOS.rds"))

#Selecting the betas useful for the models
i<-which(colnames(beta_all) %in% rownames(DNAm_metabolomics_models$models_betas))
b<-beta_all[,rownames(DNAm_metabolomics_models$models_betas)[-1]]

identical(colnames(b),rownames(DNAm_metabolomics_models$models_betas)[-1])
######################
## Apply the models ##
######################
#Apply the models
DNAm_metab_scores<-foreach(x=names(DNAm_metabolomics_models$models_betas), .combine="cbind") %do% {
  p<-apply_DNAm_metab(b, models=DNAm_metabolomics_models,selected_model = x)
}
colnames(DNAm_metab_scores)<-names(DNAm_metabolomics_models$models_betas)

identical(rownames(DNAm_metab_scores),rownames(cov))
rownames(DNAm_metab_scores)<-cov$ErgoID

#Save the models
saveRDS(DNAm_metab_scores, file=paste0(wd,"/Results/DNAm_metabolomics_features.rds"))










