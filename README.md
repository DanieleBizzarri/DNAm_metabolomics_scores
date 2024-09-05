# DNAm Metabolites Extraction

This repository contains R scripts to replicate the analysis described in the paper Bizzarri et al., 2024 at: [DNAm_metabolites](https://pubmed.ncbi.nlm.nih.gov/39154540/) .

- **DNAm Metabolites Extraction**: Obtaining DNAm metabolites from the dataset.
- **Multivariate Cox Regressions**: Performing multivariate Cox regression analyses.

## Usage
With the script apply_DNAm_models.R  it is possible to project the DNAm metabolites.
To be able to also project the cox regression the scores within GrimAge and the Episcores must be obtained.

To obtain the GrimAge scores and its components (DNAm proteins) please contact the authors of [GrimAge](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6366976/).

To obtain the Episcores you can use the code at [bAge](https://github.com/elenabernabeu/cage_bage/tree/main/bage_predictor)

## Repository Structure

## data (Folder containing input data files)
1) DNAm_metabolomics_models.rds    # File with the DNAm prediction models
2) Stepwise_Cox_regression_DNA_metab_GrimAge.rds   # File with the cox regression combining the DNAm metabolites and GrimAge
3) Stepwise_Cox_regression_DNAm_MetaboHealth_GrimAge.rds   # File with the cox regression combining the DNAm MetaboHealth and GrimAge
4) Stepwise_Cox_regression_DNAm_prot_metab_EpiScores.rds  # File with the cox regression combining the DNAm proteins from GrimAge, the DNAm metabolites, and the Episcores
5) Stepwise_Cox_regression_DNAm_prot_metab.rds  # File with the cox regression combining the DNAm metabolites and the DNAm proteins from GrimAge
### R (Folder containing R scripts)
1) apply_DNAm_models.R             # Script to obtain DNAm metabolites
2) apply_multivariate_cox_regressions.R   # Script for multivariate Cox regressions

   
## Contact
For any questions or issues, please contact Daniele Bizzarri at d.bizzarri@lumc.nl or Erik van den Akker at e.b.van_den_akker@lumc.nl.


