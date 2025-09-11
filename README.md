# MethylationBraak

This repo is designed to analyse methylation data alongside phenotypical data. It was first used to analyse the ROS/MAP cohort.

It was first used in my undergraduate project "non-linear methods reveal 267 novel CpG sites associated with Alzheimer's disease"

First, Regress.R is used to regress out covariates. Namely, age, sex and cell proportion

Second, generalised additive models are applied to each probe in GAMs.R is used to model each probes methylation values in response to the Braak stage of postmortem tissue samples

Third, MultipleTestingCorrections.R is used to control for the false discovery rate with a Benjamini-Hochberg correction, since there are 486,000+ probes in the ROS/MAP dataset

Fourth, an ensemble feature selection pipeline uses random forests, elastic net regression, Spearman's rank correlation, and gradient boosting machines to rank the relative influence of each probe's methylation on the postmortem tissues Braak stage in EFS.R

Other scripts in this Repo include a GAM cross validation error script and other general data visualisations tools using ggplot2 



