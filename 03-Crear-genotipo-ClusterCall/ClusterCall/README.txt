clusterCall-analysis-agrosavia.R

# It uses the ClusterCall package to generate genotype calls for autotetraploid samples from single
# nucleotide polymorphism (SNP) marker array data. From theta values for F1 populations (such as biparental
# families), genotype cluster positions are calibrated and with marker genotypes (training dataset)
# for two or more F1 populations, genotypes can be predicted for a population of samples of arbitrary 
# composition, such as a diversity panel.

# Train and prediction of genotypes from theta ratios
# It uses three families to train 
#a INPUTS: - Trainin datasets with six files with theta and radio files for the three families
#         - Input theta ratios for the population to predict the genotypes
# OUTPUT: - One file with predictions only for the population
# AUTHOR: Luis Garreta (lgarreta@gmail.com)
