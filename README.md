
# P_repens

Population genetics analysis of Protea repens

Code by Colin T. Kremer & Rachel Prunier

Written for:

"Isolation by distance and isolation by environment contribute to population differentiation in Protea repens (Proteaceae L.), a widespread South African species"

Submitted to AJB.


#####


metadata

##### Output of the TASSLE UNEAK pipeline ####


These files are needed to do quality control and come up with the final set of markers.  They are found only in dryad and not on github because they are too large for github.


file: HapMap.fas.txt
	description: fasta files for each SNP identified by the UNEAK SNP calling pipeline
	Structure: This file has one column, each SNP is represented by two fastas, one for each allele.  Each fasta has two rows, one with the name beginning with a ">", and one with the sequence of that allele at that marker. 
		For example: marker TP1 has two fasta files, one the "query" and one the "hit" (i.e. the two different alleles) and each sequence is 64 base long
			>TP1_query_64
			TGCAGAAAAAAAAAAAAAAAAACCCGCAGAAAAAAAAAGAGAAGGGCAAGAAGAAATGTTATCT
			>TP1_hit_64
			TGCAGAAAAAAAAAAAAAAAAACCCGCAGAAAAAAAAAGAGAAGGGCAAGAAGTAATGTTATCT

file: HapMap.hmn.txt
	description: This is a tab separated file which includes the allele calls for each individual for each polymorphic locus (SNP) identified by the UNEAK SNP calling pipeline
	structure: This file has 11 information columns and then one column for each of the 717 individuals genotyped. There is one row for each SNP, the data in each cell are 0 if the the individual was homozygous for the first allele listed in the allele column, 1 if the individual is heterozygous, and 2 if the individual is homozygous for the second allele listed in the allele column. Missing data points are indicated with a "."
		information columns: (none of these, other than the name, were used for our analyses)
			rs = the name of the SNP (e.g. TP1)
			alleles = nucleotides present at the SNP
			chrom = if there was a chromosome identified (none were), the chromosome on which the SNP falls
			pos	= a unique number for each SNP, reiterated in the name
			strand = which DNA strand the read is from 
			assembly = NA, the SNPs were discovered denovo
			center = pipeline	
			protLSID = NA	
			assayLSID = NA	
			panelLSID = NA	
			QCcode = whether the SNP passed internal QC
		data columns: one column for each individual. 
		rows: each row represents one SNP

file: HapMap.hmc.txt
	description: This is a tab separated file which includes the number of reads on which each allele call is based
	structure: This file has 1 information column and one column for each of the 717 individuals genotyped. There is one row for each SNP. Each cell has two numbers separated by a |.  The first number is the number of reads that had the sequence of the first allele for that SNP, the second number is the number of reads that had the sequence of the second allele of that SNP.
		information column:
			rs = the name of the SNP (e.g. TP1)
		data columns: one column for each individual
		rows: each row represents one SNP
		
#### processing SNPs to get final set of usable SNPS   ####		

file: clean_processing_hap_maps_v3.R 
	description: this R file contains the code that takes the output of the UNEAK pipeline resulting in 2066 markers for 663 individuals (loci20_v2.csv). It also takes that that set of markers and subsets out the outliers and non-outliers for separate analyses.  Non-outliers get used for phylogeographic analyses (MDS, TreeMix, Structure).  We assessed measures of diversity in both the high- outliers and non-outliers.
	structure: R code

####files required for phylogeographic analyses####		

file: loci_20_v2.csv
	description: This is a comma separated file which includes the allele calls for individuals and loci which passed quality control 
	application: imput file for outlier analysis, subsetted 
	structure: The columns are for each of the 2066 loci which passed QC. While row labels are included, there is no column title for that column - this must be taken into account when the data are read in to R
		rows: each row represents one SNO, row labels are the name of each of the 663 individuals that passed QC


file: clean_processing_hap_maps_v3.R 
	description: this R file contains the code that takes the output of the UNEAK pipeline resulting in 2066 markers for 663 individuals (loci20_v2.csv). It also takes that that set of markers and subsets out the outliers and non-outliers for separate analyses.  Non-outliers get used for phylogeographic analyses (MDS, TreeMix, Structure).  We assessed measures of diversity in both the high- outliers and non-outliers.
	structure: R code

file: outliers_columnnumbers.csv
	description:This comma separated file includes some of the output of the FST outlier analysis, markers that were identified  as outliers are included.
	structure: This file has six columns, each containing information on a locus that was identified as an outlier
		columns:
			marker - the name of the marker
			right column - the column number to which it belongs in the loci_20_v2 file
			FST - the fst for that locus estimated by the FST outlier analysis
			outlierness - the "p-value" of that outlier (see manuscript)
			level - whether the marker is a low or high FST outlier
		rows: one row for each marker that is identified as an outlier


File: MDS.R
	description: this R file contains the code that does the MDS analysis and produces figure 2
	structure: R code
	
#### file for BLAST analysis ####


File: pulling_out_fastas.R
	description: this R file contains the code that pulls the fasta files for the 2066 markers (high outliers, low outliers, and non-outliers, so they can be blasted against the P repens transcriptome 
	
	
#### file for Figure 4 analysis ####
	
IMput data sets:

	
File: P_repens_PopGen_analysis_code
	description: this R file contains code for the GLMM analyses of FST vs FST analysis figure 4 and supplemental figures X and Y

File: PCA code

File :allele frequencies


file: marker_fst_table.csv
	description: this comma separated file includes some of the output from the FST outlier analysis.  These results are the actual FST estimates for each marker.
	structure: This file has three columns, one information column, and two columns with information on each marker
		columns:
			marker: the name of each marker
			fst: the FST estimate for that marker
			level: whether or not the marker is an outlier (high or low) or nonoutlier.
		rows: one row for each marker that is identified as an outlier
	
file: repens_7_env_1997.csv
	description: This comma separated file includes environmental data for all protea repens populations in the CFR documented by the protea atlas project.  Environmental data was extracted from Schultze et al. 1997 and 2007.
	structure: one information column and 9 information columns
		columns:
			name: the name of the species identified (all PRREPE; Protea repens)
			LONDD: longitude in decimal degrees
			LATDD: latitude in decimal degrees
			pptcon: rainfall concentration (no units)
			MAP: mean annual precipitation (mms/year)
			tmax_jan_a: mean maxiumum temperature in January (summer) in celcius
			MAT: mean annual temperature (celcius)
			altitude: estimated altitude for the population (meters above sea level)
			tmin_july_: mean minimum temperature in July (winter) in celcius
			summer_rainfall: proportion of rainfall that falls in the summer (no units)
			
file:sample_env_goldblatt.csv
	description: This comma separated file includes environmental data for the protea repens populations sampled for this study.  Environmental data was extracted from Schultze et al. 1997 and 2007.
	structure: one information column and 9 information columns
		columns:
			name: the name of the species identified (all PRREPE; Protea repens)
			LONDD: longitude in decimal degrees
			LATDD: latitude in decimal degrees
			pptcon: rainfall concentration (no units)
			MAP: mean annual precipitation (mms/year)
			tmax_jan_a: mean maxiumum temperature in January (summer) in celcius
			MAT: mean annual temperature (celcius)
			altitude: estimated altitude for the population (meters above sea level)
			tmin_july_: mean minimum temperature in July (winter) in celcius
			summer_rainfall: proportion of rainfall that falls in the summer (no units)
			goldblatt: the name of the phytogeographic zone (as defined by Goldblatt) into which each population falls