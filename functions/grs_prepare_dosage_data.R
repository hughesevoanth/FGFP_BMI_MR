#########################################
## A function to prepare dosage data
## to estimate a GRS
## by: David Hughes
## date: March 15th
#########################################
grs_prepare_dosage_data = function( 
	dosage_data , 
	snpstats_data , 
	grs_data , 
	grs_data_SNP_col_name , 
	grs_data_EffectAllele_col_name ,
	grs_data_AltAllele_col_name ,
	grs_data_EffectAllele_FREQ_col_name ,
	grs_data_BETA_col_name,
	name_of_grs_being_processed = NA ,
	EAF_delta = 0.1 ,
	Flipped_EAF_delta = 0.15 ,
	info_filter_value = 0.8 ,
	HWE_filter_pvalue =  2.5E-8  ){

	############################
	## SNP IDs that belong in the GRS
	############################
	snps = grs_data[, grs_data_SNP_col_name]
	all_grs_snps = snps ## a source list to annotate those that are excluded and why.
	unfiltered_grs_data = grs_data
	number_of_grs_snps = length(snps)

	############################
	## Use SNPSTATS to identify those SNPs not present in
	## the sample population and remove them
	############################
	w = which(!snps %in% snpstats_data$rsid)
	
	number_of_grs_snps_not_in_dosage = length(w)
	name_of_grs_snps_not_in_dosage = snps[w]

	if(length(w)>0){ 
		## If any are abscent remove from gwas snp list
		snps = snps[-w] 
		## If any are abscent remove them from the grs data frame
		grs_data = grs_data[-w, ]
	}
	
	############################
	## match to SNP IDs in dosage
	############################
	m = match(snps, dosage_data$rsid)
	## limit and order dosage data to match grs_data
	dosage_data = dosage_data[m, ]

	############################
	## match to rsid in snpstats
	############################
	m = match(snps, snpstats_data$rsid)
	## limit and order snpstats data to match grs_data
	snpstats_data = snpstats_data[m, ]

	###########################################
	##
	##  (I) Sort, QC|Filter SNPs
	##
	##	Filter and allele sort SNPs
	## 	- run an sapply to iterate over all SNPs
	##
	###########################################
	cat( paste0("I) Filtering SNPs for allele matching and effect allele MAF\n")  )
	###
	qc_dosage_data = t( sapply(1:nrow(grs_data), function(i){
		print(i)

		## the dosage data for this snp
		dos = dosage_data[i,] 
		## easy matches; NOTE: qctools (SNPTEST) the effect allele is alleleB
		allele_match = sum( grs_data[i, grs_data_EffectAllele_col_name ] == snpstats_data[i,"alleleB"] & grs_data[i, grs_data_AltAllele_col_name ] == snpstats_data[i,"alleleA"] )
		eaf_match = abs(grs_data[i, grs_data_EffectAllele_FREQ_col_name ] - snpstats_data[i, "alleleB_frequency"])
		if(allele_match == 1 & eaf_match <= EAF_delta ){
			# data out
			out = dos
		} else {
			### if the alleles are just swapped
			allele_flip = sum( grs_data[i, grs_data_EffectAllele_col_name ] == snpstats_data[i,"alleleA"] & grs_data[i, grs_data_AltAllele_col_name ] == snpstats_data[i,"alleleB"] )
			eaf_flip = abs( grs_data[i, grs_data_EffectAllele_FREQ_col_name ] - snpstats_data[i,"alleleA_frequency"] )
			if(allele_flip == 1 & eaf_flip <= Flipped_EAF_delta ){
				## edit dosage for flipped alleles
				dos[1, 7:length(dos)] = 2 - dos[ 1, 7:length(dos)]
				## flip the alleles
				a = dos[1, "alleleA"]
				b = dos[1, "alleleB"]
				dos[1, "alleleA"] = b
				dos[1, "alleleB"] = a
				## alter SNPID
				id = dos[1, "SNPID"]
				id = strsplit(id, split = "_")[[1]]
				id = paste0(id[1], "_",id[3], "_", id[2] )
				dos[1, "SNPID"] = id
				## data out
				out = dos

			} else {
				# cat(paste0("\tCan not make a match for ", dos[1, "rsid"] , " \n"))
				## data out
				out = rep(NA, length(dos))
				}
			}
			return( unlist( out ) )
		})
	)

	############################
	## Define as data frame
	############################
	qc_dosage_data = data.frame(qc_dosage_data)

	############################
	## remove any rows with NAs
	############################
	w = which( apply(qc_dosage_data, 1, function(x){ sum(is.na(x)) == length(x) }) == 1 )
	names_snp_allele_match_filtering = grs_data[w,"SNP"]
	
	if(length(w) > 0){ 
		qc_dosage_data = qc_dosage_data[-w, ] 
		grs_data = grs_data[-w, ]
		snpstats_data = snpstats_data[-w, ]
	}

	############################
	## how many SNP were removed during allele matching?
	############################
	snp_allele_match_filtering = length(names_snp_allele_match_filtering)

	#############################
	## quick allele match count
	#############################
	# m = match(qc_dosage_data$rsid, grs_data[, grs_data_SNP_col_name] )
	snpcount =  sum( grs_data[ , grs_data_EffectAllele_col_name ] == unlist( qc_dosage_data[,"alleleB"] ) & grs_data[ , grs_data_AltAllele_col_name ] == unlist( qc_dosage_data[,"alleleA"] ), na.rm = TRUE )
	number_of_allele_matched_snps_in_qc_dosage_data = snpcount


	############################################
	##
	## (II) Make All effects positive
	##
	## IF the BMI beta estimate in the grs_data  
	## file is negative flip the allele
	############################################
	cat( paste0("II) Making all effect allele trait|disease increasing alleles\n")  )

	qc_dosage_data_pos_beta = t( sapply(1:nrow(grs_data), function(i){
		dos = qc_dosage_data[i,]
		beta = grs_data[i, grs_data_BETA_col_name ]
		if(beta <  0){
			## edit dosage for flipped alleles
			dos[7:length(dos)] = 2 - as.numeric(dos[7:length(dos)])
			## flip the alleles
			a = dos["alleleA"]
			b = dos["alleleB"]
			dos["alleleA"] = b
			dos["alleleB"] = a
			## alter SNPID (from 1:1590521_G_A to 1:1590521_A_G)
			# id = dos["SNPID"]
			# id = strsplit(id, split = "_")[[1]]
			# id = paste0(id[1], "_",id[3], "_", id[2] )
			# dos["SNPID"] = id
			# data out
			out = dos
			} else {
				out = dos
			}
			return(unlist(dos))
		})
	)

	## define as data.frame
	qc_dosage_data_pos_beta = as.data.frame(qc_dosage_data_pos_beta)

	## set dosage to numeric
	for(i in 7:ncol(qc_dosage_data_pos_beta)){
		qc_dosage_data_pos_beta[,i] = as.numeric(as.character(qc_dosage_data_pos_beta[,i]))
	}

	######################################
	##
	## (III) Generate GRS
	##
	######################################
	cat( paste0("III) Estimating GRS\n")  )

	#############################
	## define the ABSOLUTE ( beta )
	#############################
		## match grs_data snps to dosage snps
	# m = match( qc_dosage_data_pos_beta$rsid, grs_data[,grs_data_SNP_col_name] )
	
	beta = abs( as.numeric( as.character( grs_data[ , grs_data_BETA_col_name] ) ) )

	#############################
	## estimate GRS
	#############################
	grs = t( apply(qc_dosage_data_pos_beta[, 7:ncol(qc_dosage_data_pos_beta) ], 2 , function(dos){
		un_w_grs = sum(dos, na.rm = TRUE)
		w_grs = sum( dos * beta , na.rm = TRUE)
		out = c(un_w_grs, w_grs)
		return(out) 
		})
	)
	###
	colnames(grs) = c("GRS","wGRS")

	########################################
	##
	## (IV) Generate FILTERED GRS
	##		 filter on imputation quality and HWE
	##
	########################################
	cat( paste0("IV) Estimating SNP QC'd GRS\n")  )

	#############################
	## Which SNPs, if any, should 
	## be filtered in the sample population?
	#############################
	## limit snpstats to those snps still in the dosage file
	# m = match(qc_dosage_data_pos_beta$rsid, snpstats_data$rsid)
	# snpstats_data = snpstats_data[m,]

	## which SNPs should be filtered?
	w = which(snpstats_data$impute_info < info_filter_value | snpstats_data$HW_exact_p_value < HWE_filter_pvalue)
	# snps_2_filter = snpstats_data$rsid[w]

	############################
	## how many SNP were removed during allele matching?
	############################
	snp_info_hwe_filtering = length(w)
	name_of_snp_info_hwe_filtering = snpstats_data$rsid[w]
	#############################
	## Identify snps to filter out in the GRS data frame
	#############################
	# w = which( qc_dosage_data_pos_beta$rsid %in% snps_2_filter)
	
	#############################
	## If any are present lets re-estiamte the GRS
	#############################
	if(length(w)> 0){

		#############################
		## re-define data
		#############################
		snpstats_data = snpstats_data[-w, ]
		qc_dosage_data_pos_beta = qc_dosage_data_pos_beta[-w, ]
		grs_data = grs_data[-w, ]

		#############################
		## define the ABSOLUTE ( beta )
		#############################
		## match grs_data snps to dosage snps
		# m = match( qc_dosage_data_pos_beta$rsid, grs_data[, grs_data_SNP_col_name ] )

		beta_filter = abs( as.numeric( as.character( grs_data[  , grs_data_BETA_col_name] ) ) )

		#############################
		## estiamte GRS
		#############################
		grs_filter = t( apply( qc_dosage_data_pos_beta[ , 7:ncol(qc_dosage_data_pos_beta) ] , 2 , function(dos){
			un_w_grs = sum(dos, na.rm = TRUE)
			w_grs = sum( dos * beta_filter , na.rm = TRUE)
			out = c(un_w_grs, w_grs)
			return(out) 
			})
		)
		###
		colnames(grs_filter) = c("GRS_info_hwe_filtered","wGRS_info_hwe_filtered")
		###
	}

	if(exists("grs_filter")){
		grs = cbind(grs, grs_filter)
	}

	if(!is.na(name_of_grs_being_processed)){
		colnames(grs) = paste0( name_of_grs_being_processed, "_", colnames(grs) )	
	}
	
	sumstats = t( data.frame(
		number_of_grs_snps = number_of_grs_snps,
		number_of_grs_snps_not_in_dosage = number_of_grs_snps_not_in_dosage,
		number_of_grs_snps_filtered_during_allele_matching = snp_allele_match_filtering,
		number_of_allele_matched_snps_in_grs = number_of_allele_matched_snps_in_qc_dosage_data,
		number_of_grs_snps_filtered_for_info_hwe = snp_info_hwe_filtering,
		number_of_allele_matched_snps_in_filtered_grs = number_of_allele_matched_snps_in_qc_dosage_data - snp_info_hwe_filtering ) )
	
	## What SNPs filtered and why ?
	filtered_snp_ids = data.frame(rsids = all_grs_snps, 
	                              snps_not_genotyped = 0,
	                              snps_w_allele_or_EAF_mismatches = 0,
	                              snps_filtered_for_info_or_hwe = 0
	                              )
	if( length(name_of_grs_snps_not_in_dosage) > 0 ){
	  w = which(all_grs_snps %in% name_of_grs_snps_not_in_dosage)
	  filtered_snp_ids$snps_not_genotyped[w] = 1
	}
	if( length(names_snp_allele_match_filtering) > 0 ){
	  w = which(all_grs_snps %in% names_snp_allele_match_filtering)
	  filtered_snp_ids$snps_w_allele_or_EAF_mismatches[w] = 1
	}
	if( length(name_of_snp_info_hwe_filtering) > 0 ){
	  w = which(all_grs_snps %in% name_of_snp_info_hwe_filtering)
	  filtered_snp_ids$snps_filtered_for_info_or_hwe[w] = 1
	}
	
	### Add filtered data to the unfiltered_grs_data table
	m = match(unfiltered_grs_data$SNP, filtered_snp_ids$rsids)
	unfiltered_grs_data = cbind(unfiltered_grs_data,filtered_snp_ids[m, -1] )
	

	return( list( grs = grs, sumstats = sumstats, grs_data_with_filtering_information = unfiltered_grs_data)  )

} ## END OF FUNCTION




