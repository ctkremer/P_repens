
# Population genetics analysis of Protea repens

Code by Rachel Prunier, Colin T. Kremer, & Kent Holsinger

"Isolation by distance and isolation by environment contribute to population differentiation in Protea repens (Proteaceae L.), a widespread South African species." 2017. American Journal of Botany. doi:______

See also Dryad repository, doi:_______

Here we describe the various R scripts, input & output data files created for this project. Collectively, they support the results and figures presented in our AJB paper. To aid in understanding the workflow, we provide a .pdf flow chart. Details on each script and associated files follow, first as a short list of each R script and data file associated with this project. Subsequent sections provide a narrative of the flow of the analysis, touching on each R script, and providing detailed metadata on each data set as it is first mentioned. 


### R SCRIPTS
`clean_PCA_files.R` - Script running PCA analysis of environmental variables

`MDS.R` - MDS analysis

`P_repens_PopGen_analysis_code.R` - IBD vs. IBE regression analyses

`processing_hap_map_data.R` - Data processing and QC of original HapMap data

`pulling_out_fastas.R` - Creates fasta table for focal loci


### DATA FILES

More detailed meta-data provided after this list, as files are referenced.

`all_allele_frequencies_032417.csv` - Table of allele frequencies (population by marker)

`all_coefficient_estimates_combined.csv` - Table of all MCMCglmm coefficient estimates

`HapMap.fas.txt` - fasta files for each SNP identified by the UNEAK SNP calling pipeline

`HapMap.hmn.txt` - alleles for each individual for each polymorphic locus (SNP) identified by UNEAK

`HapMap.hmc.txt` - number of reads on which each allele call is based

`loci_20_v2.csv` - allele data for all individuals and loci (those which passed QC)

`marker_fst_table.csv` - Table of each loci/marker, indicating its Fst and outlier category

`noouts_data.csv` - allele data for all individuals and non-outlier loci

`outliers_columnnumbers.csv` - FST and column indices for focal markers

`pairwise_population_differences.csv` - Pairwise genetic and environmental distances among populations

`randomized_all_marker_deltas_for_analysis_v1.csv` - Randomized pairwise distances, take 1

`randomized_all_marker_deltas_for_analysis_v2.csv` - Randomized pairwise distances, take 2

`randomized_all_marker_deltas_for_analysis_v3.csv` - Randomized pairwise distances, take 3

`randomized_all_marker_deltas_for_analysis_v4.csv` - Randomized pairwise distances, take 4

`randomized_all_marker_deltas_for_analysis_v5.csv` - Randomized pairwise distances, take 5

`repens_7_env_1997.csv` - environmental variables from Protea atlas database

`sample_env_goldblatt.csv` - location, environmental variables, and goldblatt data for focal populations

`sample_PCA_and_covars.csv` - PCA coordinates (based on environ. variables) for each population



### Processing the output of the TASSLE UNEAK pipeline

Initial sequencing output produces data on a large number of loci, but many of these received low coverage, or were present in only a few individuals. The script `processing_hap_map_data.R`, combined with details in the text of the paper, describe how we performed quality control on these data. Specifically, this code takes the output of the UNEAK pipeline (`HapMap.fas.txt`,`HapMap.hmn.txt`,`HapMap.hmc.txt`) and conducts QC, resulting in 2066 markers for 663 individuals (`loci20_v2.csv`). These markers are further subdivided into outliers and non-outliers (based on the results of the outlier analysis, contained in `outliers_columnnumbers.csv`). Non-outliers get used for phylogeographic analyses (MDS, TreeMix, Structure), based on output file created in this script (`noouts_data.csv`). Finally, within this script we calculate allele frequencies for each marker at the level of each population, producing the output file `all_allele_frequencies_032417.csv`

This script requires as input the following files, which are archived on Dryad but not github, due to file size restrictions:

#### `HapMap.fas.txt`
	This file (not on github) has one column, each SNP is represented by two fastas, one for each allele.  Each fasta has two rows, one with the name beginning with a ">", and one with the sequence of that allele at that marker. 
		For example: marker TP1 has two fasta files, one the "query" and one the "hit" (i.e. the two different alleles) and each sequence is 64 base long
			>TP1_query_64
			TGCAGAAAAAAAAAAAAAAAAACCCGCAGAAAAAAAAAGAGAAGGGCAAGAAGAAATGTTATCT
			>TP1_hit_64
			TGCAGAAAAAAAAAAAAAAAAACCCGCAGAAAAAAAAAGAGAAGGGCAAGAAGTAATGTTATCT

#### `HapMap.hmn.txt`
	This file (not on github) has 11 information columns and then one column for each of the 717 individuals genotyped. There is one row for each SNP, the data in each cell are 0 if the the individual was homozygous for the first allele listed in the allele column, 1 if the individual is heterozygous, and 2 if the individual is homozygous for the second allele listed in the allele column. Missing data points are indicated with a "."
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

#### `HapMap.hmc.txt`
	This file (not on github) has 1 information column and one column for each of the 717 individuals genotyped. There is one row for each SNP. Each cell has two numbers separated by a |.  The first number is the number of reads that had the sequence of the first allele for that SNP, the second number is the number of reads that had the sequence of the second allele of that SNP.
		information column:
			rs = the name of the SNP (e.g. TP1)
		data columns: one column for each individual
		rows: each row represents one SNP

#### `outliers_columnnumbers.csv`
	This comma separated value file includes some of the output of the FST outlier analysis, markers that were identified  as outliers are included. This file has six columns, each containing information on a locus that was identified as an outlier
		columns:
			marker = the name of the marker
			right column = the column number to which it belongs in the loci_20_v2 file
			FST = the fst for that locus estimated by the FST outlier analysis
			outlierness = the "p-value" of that outlier (see manuscript)
			level = whether the marker is a low or high FST outlier
		rows: one row for each marker that is identified as an outlier

This script produces as output the following files:

#### `loci_20_v2.csv`
	This comma separated value file contains data on the genotypes of individuals for each of the 2066 loci which passed QC (See above). While row labels are included, there is no column title for that column - this must be taken into account when the data are read in to R. The data contained in each cell are: 0 if the the individual was homozygous for the first allele listed in the allele column, 1 if the individual is heterozygous, and 2 if the individual is homozygous for the second allele listed in the allele column. Missing data points are indicated with a "."
	  information columns: 
	    the first column contains names of individuals
	  data columns: each column contains data on one of 2066 loci.
		rows: each row represents one SNP, row labels are the name of each of the 663 individuals that passed QC
		
#### `all_allele_frequencies_032417.csv`
	This comma separated value file, produced by this R script, contains the population-level allele frequencies of each of the 2066 focal markers. 
    information column:
      pop = the name of each P. repens population considered
    data columns: one column for each marker.
    rows: populations

#### `noouts_data.csv`
	This comma separated value file is structured similarly to `loci_20_v2.csv`, but contains only the markers determined to be non-outliers within the outlier analysis part of our study. 
	  information column:
	    names = individual names
	  data columns: one column for each non-outlier marker
		rows: each row represents one SNP, row labels are the name of each of the 663 individuals that passed QC

	
	
### Processing environmental data

To analyze the effects of environment on population divergence, we needed to characterize key environmental variables and how they differ among populations. As many environmental variables co-vary, we used a PCA to generate independent environmental axes, as detailed in the script `clean_PCA_files.R`. This required two input data sets (`repens_7_env_1997.csv`, `sample_env_goldblatt.csv`) and produces one output data set (`sample_PCA_and_covars.csv`):

#### `repens_7_env_1997.csv`
	This comma separated file includes environmental data for all protea repens populations in the CFR documented by the protea atlas project.  Environmental data was extracted from Schultze et al. 1997 and 2007. This data set contains one information column and 9 data columns.
		information column:
			name = the name of the species identified (all PRREPE; Protea repens)
		data columns:
			LONDD = longitude in decimal degrees
			LATDD = latitude in decimal degrees
			pptcon = rainfall concentration (no units)
			MAP = mean annual precipitation (mms/year)
			tmax_jan_a = mean maxiumum temperature in January (summer) in celcius
			MAT = mean annual temperature (celcius)
			altitude = estimated altitude for the population (meters above sea level)
			tmin_july_ = mean minimum temperature in July (winter) in celcius
			summer_rainfall = proportion of rainfall that falls in the summer (no units)

			
#### `sample_env_goldblatt.csv`
	This comma separated file includes environmental data for the additional protea repens populations specifically sampled for this study.  Environmental data was extracted from Schultze et al. 1997 and 2007. This file contains one information column and 9 information columns.
		information column:
			name = the P. repens population code 
		data columns:
			LONDD = longitude in decimal degrees
			LATDD = latitude in decimal degrees
			pptcon = rainfall concentration (no units)
			MAP = mean annual precipitation (mms/year)
			tmax_jan_a = mean maxiumum temperature in January (summer) in celcius
			MAT = mean annual temperature (celcius)
			altitude = estimated altitude for the population (meters above sea level)
			tmin_july_ = mean minimum temperature in July (winter) in celcius
			summer_rainfall = proportion of rainfall that falls in the summer (no units)
			goldblatt = the name of the phytogeographic zone (as defined by Goldblatt) into which each population falls

#### `sample_PCA_and_covars.csv`
	This comma separated file includes PCA axes based on environmental data (see above), matched to each population, as well as its phytogeographic zone. This file contains one information column and 9 information columns.
    information column:
		  sample_names = the name of the species identified (all PRREPE; Protea repens)
	  data columns:
			LONDD = longitude in decimal degrees
			LATDD = latitude in decimal degrees
			Comp.1 = PCA axis 1 (dominant axis)
			Comp.2 = PCA axis 2 (second most dominant axis), etc....
			Comp.7 = 7th PCA axis
			goldblatt = the name of the phytogeographic zone (as defined by Goldblatt) into which each population falls
			




### Phylogeographic analyses

Explorations of the structure of P. repens populations are described in detail in the text, and include MDS, TreeMix, Structure, and arlequin analyses. Most of these are conducted using stand-alone software. The MDS analysis, however, is found in the R script `MDS.R`, and was used to produce figure 2. These analyses required input files including `noouts_data.csv`.

	
### BLAST analysis

Before blasting sequences against the P. repens transcriptome, we needed to extract fasta files containing only the focal 2066 markers (high, low, and non-outliers). This data processing is detailed in the script `pulling_out_fastas.R`. As input, it requires `HapMap.fas.txt` (see meta-data description above), and `marker_fst_table.csv`.

#### `marker_fst_table.csv`
	This comma separated value file lists the focal 2066 markers, along with their actual FST estimates and classification as high, low, or non-outliers. This file has three columns, one information column and two data columns on each marker.
		information column:
			marker = the name of each marker
		data column:
			fst = the fst for a given locus estimated during the FST outlier analysis
			level = whether the marker is a low or high FST outlier
		rows: one row for each marker of the 2066 markers



### IBD vs. IBE regression analyses (Fig. 4 & Supplemental Files)
	
For each of the 2066 markers, we fit generalized linear mixed models relating the pairwise genetic distances among populations to environmental and physical distances. We then related the effect sizes characterizing each of the environmental and physical covariates to the FST of each individual marker. (For more details on this novel analysis, see the main text). These analyses are contained in the R script: `P_repens_PopGen_analysis_code.R`. This analysis required several input data sets (`all_allele_frequencies_032417.csv`, `marker_fst_table.csv`, both described above). For convenience, we have also provided several of the intermediate data sets produced in this analysis, as several steps required significant computation time. Starting with these intermediate files, interested readers can jump through the analysis.

#### `pairwise_population_differences.csv`
	A data set of the pairwise genetic and environmental distances among populations, calculated as an intermediate data file before MCMCglmm regressions.
    information columns:
      pop1 = population 1 name
      pop2 = population 2 name (different from 1)
      physdist = physical distance between population pairs
      distpc1.a	= distance on PCA axis 1 between population pairs
      distpc2.a	= distance on PCA axis 2 between population pairs
      distpc3.a	= distance on PCA axis 3 between population pairs
      goldblattsame = are populations in the same (1) or different (0) phytogeographic zones?
    data columns:
      genetic distances between each of the 2066 focal markers.
    rows: one for each unique pair of populations (order independent)


#### `randomized_all_marker_deltas_for_analysis_v1.csv`
	There are a set of 5 of these files, which have the exact same structure as `pairwise_population_differences.csv` (see above). However, the corresonding environmental, physical, and genetic distances are the result of unique randomizations, breaking any possible empirical correlation between genetic and physical/environmental distances observed in the original data.


#### `all_coefficient_estimates_combined.csv`
	This output table contains the coefficient estimates produced by running MCMCglmm regressions on all 2066 markers, examining the effects of physical and environmental distance.
      information columns:
        run = indicates whether coefficients are for analysis of the original data set, or one of 5 randomized simulations.
        marker = name of individual marker
        cf = name of the coefficient (intercept, physical distance, etc.)
        fst = estimated FST for each marker
        level = category of each marker as a high, low, or non-outlier
      data columns:
        post.mean = mean of the posterior distribution of each coefficient from the MCMCglmm regression
        pMCMC = p-value estimate for difference of coefficient from 0
      rows: one for each individual coefficient estimated.
      