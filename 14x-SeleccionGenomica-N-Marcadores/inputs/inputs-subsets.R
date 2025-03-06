#!/usr/bin/Rscript
# Create a geno and pheno for a N number of samples

N = 200

geno   = read.csv ('genotipo-ALL.csv', check.names=F)
phenos = read.csv ('fenotipos-ALL.csv', check.names=F)

# Select phenos
testPheno = phenos [1:100, c(1, 2, 3, 5, 6)]; 
write.csv (testPheno, sprintf ('fenotipos-%s.csv', N), row.names=F)

# Select geno
regs = testPheno [,1];
testGeno = geno [, colnames (geno) %in% c('Markers', regs)]; 
write.csv (testGeno, sprintf ('genotipo-%s.csv', N),  row.names=F)
