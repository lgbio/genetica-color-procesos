#!/usr/bin/Rscript
source ("lglib14.R")

USAGE="
Update genotipe file with new chromosomal positions from genome v6.1"

INPUTGENO  = "inputs/genotipo-AndigenaCCC-ClusterCall2020-ACGT.csv"
INPUTMAP   = "inputs/solcap_69k_SNPs_DM_v6_1_pos.txt"
OUTPUTFILE = "outputs/out-genotipo-AndigenaCCC-ClusterCall2023-ACGT.csv"

geno = read.csv (INPUTGENO, check.names=F); #view (geno)
map  = read.csv (INPUTMAP, sep="\t"); #view (map)
map$SNP_ID = gsub ("solcap_snp_", "", map$SNP_ID); #view (map)
map  = map [!duplicated (map$SNP_ID),]

message ("Merging geno and map files...")
genom = merge (geno, map, by.x="Markers", by.y="SNP_ID", all.y=F, sort=F, no.dups=F); #view (genom)
genom$Position = genom$SNP_POS; #view (genom)
columns = setdiff (names (genom), names (map))
genoUpd = genom [, columns]; #view (genoUpd)

createDir ("outputs")
write.csv (genoUpd, OUTPUTFILE, row.names=F)

