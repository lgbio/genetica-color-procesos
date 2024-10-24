Annotations of SNP markers resulting from GWAS or GS with SPUDDB 
Genome assembly (v6.1) annotations.
----------------------------------------------------------------

"Improved genome assembly and annotation (v6.1) for the doubled 
monoploid potato DM 1-3 516 R44"
at http://spuddb.uga.edu/dm_v6_1_download.shtml

Input files:
	1. File with markers (SNPs) to annotate: short SNP id instead large id  (e.g c1_7694 
	   instead solcap_snp_c1_7694)
	   Example file: "markers-scores.csv"

	2. File with markers map: Chrommosome and Position
	   Genome v6.1 file: "solcap_69k_SNPs_DM_v6_1_pos.txt"

	3. File with gene annotations: Gen id with starting and ending position
	   Genome v6.1 file: "DM_1-3_516_R44_potato.v6.1.hc_gene_models.gff3"
	
	4. File with GO term assignments: Gen id with GO terms
	   Genome v6.1 file: "DM_1-3_516_R44_potato.v6.1.working_models.iprscan_go_terms.txt"
	
	5. File with InterPro assignments: Gen id with annotations
	   Genome v6.1 file: "DM_1-3_516_R44_potato.v6.1.working_models.iprscan_output.tsv"

Output file:
	1. File with columns "SNP","REF","SNP_POS","SEQID","GOTERMS","SOURCE","ANNOT"

