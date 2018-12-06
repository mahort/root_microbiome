# TODO: Add comment
# 
# Author: matt
###############################################################################

getBlups <- function(kingdom, analyticalMethod, speciesMatrix, speciesMatrixOffsets){

	require(lme4);

	cat("------------------------------\n");
	cat("Using: ", analyticalMethod, " to extract the blups from the ", kingdom, ".\n", sep="");
	
	ecotypeIds <- as.factor(determineSamples(speciesMatrix));
	cat("In the case of the", kingdom, ", there are: ", nlevels(ecotypeIds)," genotypes.\n", sep="");
	
	speciesMatrixWideFormat <- list();
	index <- 0;
	
	###############################################################################	
	## how do we line these up/sync these?
	## options:
	## a.) use glmer(family="poisson", with offsets, ...) and then work on those random effects 
	## b.) use lmer(log(x + 1), with offsets)... work on those random effects
	## nb: offsets make more sense for these non-normalized data... hellinger transformation will vary on community evenness and its accuracy will be heavily affected by differences in sampling effort.
	###############################################################################	
	for( j in 1:nrow(speciesMatrix)){
		otu_j <- speciesMatrix[j,];
		
		if( analyticalMethod == "lmer" ){
			selectedModel <- lmer(I(log(otu_j + 1)) ~ 1 + (1|ecotypeIds) + offset(log(speciesMatrixOffsets)));

		} else if( analyticalMethod == "glmer" ){
			selectedModel <- glmer(otu_j ~ 1 + (1|ecotypeIds) + offset(log(speciesMatrixOffsets)), family=poisson);

		} else if( analyticalMethod == "logit" ){
			selectedModel <- glmer(otu_j ~ 1 + (1|ecotypeIds), family=binomial);
		}

		varianceCheck <- as.data.frame(VarCorr(selectedModel));
		varianceCheck <- subset(varianceCheck, grp == "ecotypeIds");
		if( is.nan(varianceCheck[1, "vcov"]) | varianceCheck[1, "vcov"] < 1e-08 ){ ## why might it be zero?
			if( j %% 500 == 0 ){ cat("Skipping and on index:", j, "\n"); };
			next;
		}

		if( index %% 250 == 0 ){ cat("On element: ", j, ".\n", sep=""); }

		index <- index + 1;
		otu.rf <- ranef(selectedModel)$ecotypeIds;
		colnames(otu.rf) <- c(paste(ifelse(kingdom == "bacteria", "B", "F"), "_", rownames(speciesMatrix)[j], sep="")); ##"ecotype_id",
		speciesMatrixWideFormat[[index]] <- t(otu.rf);
	}

	speciesMatrixWideFormat <- t(do.call(rbind, speciesMatrixWideFormat)); ## cbind doesn't check names.
	speciesMatrixWideFormat <- data.frame(scale(speciesMatrixWideFormat));

	return(speciesMatrixWideFormat)
}

## using the default arguments, this function converts all of the values in a matrix to the positive set.
getDistanceMatrix <- function( dataMatrix, method, positivity=TRUE ){

	require(vegan); 

	if( positivity ){
		cat("Scaling distance matrix to positive numbers.\n");
		dataMatrix <- apply( dataMatrix, 2, function(x){ x + (abs(min(x)) + 1); });
	}

	if( method == "chord" ){
		distanceMatrix <- vegdist( decostand(dataMatrix, "norm"), method="euclidean");

	} else if( method == "squaredEuclidean" ){
		distanceMatrix <- vegdist( dataMatrix, method="euclidean")^2;

	} else if( method == "chi.square" ){
		distanceMatrix <- vegdist( decostand( dataMatrix, "chi.square"), method="euclidean"); 

	} else {
		distanceMatrix <- vegdist( dataMatrix, method=method)
	}

	return(distanceMatrix);
}

############################################################################################################################################
## glm based r2
############################################################################################################################################
pseudoR2 <- function(model0, model1){ ## need to update this for *mer models, this works with base/native glm models...
	n <- nrow(model0$model);
	chisq <- model0$deviance - model1$deviance;
	
	cox <- 1 - exp(-chisq/n);
	nagelkerke <- cox/(1 - exp(-model0$deviance/n));
	return(list("cox"=cox, "nagelkerke"=nagelkerke));
}

normalizeInSomeWay <- function( speciesMatrix, normalizationMethod, eliminateSingletons=FALSE ){

	require(vegan);

	################################################################################################################################################
	## various normalization protocols for the microbial species matrix.
	################################################################################################################################################
	if( normalizationMethod == "raw" ){
		cat("passing the raw data matrix through.\n");
		resultingMatrix <- speciesMatrix;

	} else if( normalizationMethod == "sqrt" ){
		cat("Square-rooting the raw values.\n");
		resultingMatrix <- sqrt(speciesMatrix);

	} else if( normalizationMethod == "pa" ){
		cat("Turning the species data into a presence/absence matrix.\n");
		resultingMatrix <- t(decostand(t(speciesMatrix), "pa"));

	} else if( normalizationMethod == "tlc" ){
		resultingMatrix <- t(decostand(t(speciesMatrix), "total"));

	} else if( normalizationMethod == "asin" ){
		cat("Performing arc-sin transformation.\n");
		resultingMatrix <- asin(sqrt(t(decostand(t(speciesMatrix), "total"))));

	} else if( normalizationMethod == "hlgr" ){
		cat("Performing hellinger transformation.\n");
		resultingMatrix <- t(decostand(t(speciesMatrix), "hellinger"));

	} else if( normalizationMethod == "logx+1" ){
		resultingMatrix <- log(speciesMatrix + 1);

	} else {
		stop("Not handled.\n");
	}

	if( eliminateSingletons ){
		dropouts <- which(apply(resultingMatrix, 1, function(x){ sum( x > 0 );}) <= 1);

		if( length(dropouts) > 0 ){
			cat("Dropping", length(dropouts), "singletons.\n");
			resultingMatrix <- resultingMatrix[-dropouts,];

		} else {
			cat("There are:", nrow(resultingMatrix), "species under consideration.\n");
		}
	}

	return(resultingMatrix);
}

applyThreshold <- function( otuMatrix, otuCountThreshold ){
	
	if( otuCountThreshold > 0 ){
		numberOfSamples_OTU_is_in <- apply(otuMatrix, 1, function(x){ sum( x > 0 );});
		indicesToOmit <- which(numberOfSamples_OTU_is_in < otuThreshold);

		if( length(indicesToOmit) > 0 ){
			otuMatrix <- otuMatrix[-indicesToOmit,];
			cat("Omitting:", length(indicesToOmit), "OTUs.\n");
		}
		
	} else {

		dropouts <- which(rowSums(otuMatrix) == 0 );
		if( length(dropouts) != 0 ){
			cat("Removing some dropouts holms.\n");
			otuMatrix <- otuMatrix[-dropouts,];
		}
	}

	cat("Eliminated OTUs with counts <=", otuThreshold, "\n");
	return(otuMatrix);
}

getLeaves <- function(matrix, speciesCountCheck=FALSE){
	
	cat("getLeaves assumes the matrix is organized species by site.\n");
	plates <- determinePlates(matrix);

	leafColumnNames <- unlist(lapply("^L", grep, plates));
	matrix <- matrix[,leafColumnNames]; ## retrieve the columns that start with the characters "L_"

	## optional: ensure that all of the taxa have counts in this habitat
	if( speciesCountCheck ){
		matrix <- matrix[which(rowSums(matrix) > 0),];
	}

	return( matrix );
}

getRoots <- function(matrix, speciesCountCheck=FALSE){
	
	cat("getRoots assumes the matrix is organized species by site.\n");
	plates <- determinePlates(matrix);

	rootColumnNames <- unlist(lapply("^R", grep, plates));
	matrix <- matrix[,rootColumnNames]; ## retrieve the columns that start with the characters "R_"

	## optional: ensure that all of the taxa have counts in this habitat
	if( speciesCountCheck ){
		matrix <- matrix[which(rowSums(matrix) > 0),];
	}
	
	return( matrix );
}

## need a few supplementary files to make this work.
determineBlocks <- function( matrix, silent=F){
	if( !silent ){
		cat("nb: This expects a table in species (rows) by samples (columns).\n");		
	}

	blocks <- as.matrix(read.table(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "data/blocks.txt"), header=T, sep="\t", as.is=T, stringsAsFactors=FALSE));
	
	plates <- determinePlates(matrix, silent=TRUE);
	
	if( "replicate_id" %in% rownames(matrix)){
		cat("Determine Blocks: We have to handle replicates.\n");
		plates <- cbind(colnames(matrix), matrix["ecotype_id",], matrix["replicate_id",], plates);
		colnames(plates) <- c("sample_id", "ecotype_id", "replicate_id", "plate_id");
		
	} else {
		ecotypeIds <- colnames(matrix);
		ecotypeIds <- strsplit(ecotypeIds, "_");
		ecotypeIds <- do.call(rbind, ecotypeIds)[,3]
		plates <- unique(cbind(colnames(matrix), ecotypeIds, plates));
		colnames(plates) <- c("sample_id", "ecotype_id", "plate_id");
	}

	blocks <- merge(plates, blocks, by="plate_id", all.x=T);
	blocks <- as.matrix(blocks);
	
	# if there are lines w/o block info, use a patch (problem: plate 1 had mixed blocks).
	indices <- which(is.na(blocks[,"block_id"]));
	if( length(indices) > 0 ){
		
		sampleNames <- blocks[indices, "sample_id"];
		lastlines <- read.table( paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "data/last.lines.covariates.txt"), header=T, sep="\t", as.is=T );
		
		for( sample_i in 1:length(sampleNames)){
			name <- sampleNames[sample_i];
			index <- which(lastlines[,1] == name);
			
			if( length(index) == 0 ){
				cat("Sample:", name, "wasn't found. Skipping.\n");
				next;
			}
			
			blocks[indices[sample_i], "block_id"] <- lastlines[index, "block_id"];
			blocks[indices[sample_i], "ptp_bacteria"] <- lastlines[index, "cycle454"];
			blocks[indices[sample_i], "ptp_fungi"] <- lastlines[index, "cycle454"]; # for these last 13, cycle454 applies to both bacteria and fungi.
		}
	}

	blocks[,"plate_id"] <- sapply(blocks[,"plate_id"], gsub, pattern="^L|^R", replacement="");
	indices <- lapply(colnames(matrix), grep, blocks[,"sample_id"], ignore.case=T);
	if( length(which(lapply(indices, length) == 0)) > 0 ){
		stop("Error in resorting the matrix!!!!\n");
	}

	blocks <- blocks[unlist(indices),];
	return( blocks );
}

getBlockData <- function(matrix, silent=F){
	blocks <- determineBlocks(matrix, silent);
	return(blocks[,"block_id"]);
}

getPtpIds <- function(matrix, kingdom, silent=F){
	targetPtpColumn <- paste("ptp_", kingdom, sep="");
	ptpIds <- determineBlocks(matrix, silent)[,targetPtpColumn];

	return(ptpIds);
}

getTaxonomy <- function(taxonomicData){
	categories <- c("kingdom", "phylum", "class", "order", "family", "genus", "species");
	
	taxa <- strsplit(taxonomicData[,"lineage"], ";");
	numberOfLevels <- max(unlist(lapply(taxa, length)));
	taxonomy <- c();
	for( i in 1:length(taxa)){
		taxon <- taxa[[i]];
		taxon <- c( taxonomicData[i,"id"], taxon, rep("", numberOfLevels - length(taxon)));
		taxonomy <- rbind(taxonomy, taxon)
	}
	
	colnames(taxonomy) <- c("otu_id", categories[1:numberOfLevels]);
	return(taxonomy);
}

## nb: I wrote this for various 'apply' function (e.g. over rows), which perhaps explains why I had chosen to do this at the names level instead of the colnames level. ;)
getLowestTaxonomicAssignment <- function(taxon){
	for( i in which(names(taxon) == "genus"):which(names(taxon) == "kingdom")){
		if( is.na(taxon[i])){
			return( "" );

		} else if( taxon[i] != "" ){
			return(taxon[i]);
		}	
	}
	
	return("unk");
}

determinePlates <- function(matrix, silent=FALSE){
	
	if( !silent ){
		cat("nb: determinePlates expects a table in species (rows) by samples (columns).\n");		
	}
	
	samples <- colnames(matrix);
	samples <- strsplit(samples,"_");
	
	tmp <- list();
	for( i in 1:length(samples)){
		sample <- samples[[i]];
		tmp[[i]] <- sample[1];
	}
	
	return( do.call(rbind, tmp));
}

determineSamples <- function(matrix, silent=FALSE){
	if( !silent ){
		cat("nb: determineSamples expects a table in species (rows) by samples (columns).\n");
	}
	
	samples <- colnames(matrix);
	samples <- strsplit(samples,"_");
	
	tmp <- list();
	for( i in 1:length(samples)){
		sample <- samples[[i]];
		tmp[[i]] <- sample[3];
	}
	
	return( do.call(rbind, tmp));
}

restrictSamples <- function(matrix, min_n=2, silent=FALSE){
	
	if( !silent ){
		cat("nb: restrictSamples expects a table in species (rows) by samples (columns).\n");
	}
	
	sampleSet <- determineSamples(matrix, silent);
	omits <- names(which(table(sampleSet) < min_n));
	
	if( length(omits) > 0 ){
		indices <- unlist(lapply(omits,grep,sampleSet));
		matrix <- matrix[,-indices];
	}
	
	## some of the OTUs may now be represented by NO samples, REMOVE THEM
	dropouts <- which(rowSums(matrix) == 0);
	if( length(dropouts) > 0 ){
		matrix <- matrix[-dropouts,];
	}
	
	return(matrix);
}

removeDropouts <- function(matrix, silent=F){ ## called by openData
	
	# drop out the zeros
	dropouts <- which(rowSums(matrix) == 0 );
	if(length(dropouts) > 0){
		if( !silent ){
			cat("Eliminating zeros now.\n");
			cat("Was this a subset of a larger matrix... to a particular set of accessions or tissue?\n");
			cat("Achtung! Is this table in the format that you expect (e.g. minimumSamplesPerOTU?)\n");
		}
		
		matrix <- matrix[-dropouts,];
	}
	
	return(matrix);
}

##
rarefyColumn <- function(sampleNames, sample, toSize=800){
	
	frame <- data.frame(names=sampleNames, count=sample, stringsAsFactors=FALSE);
	frame <- do.call(c, 
			apply(frame, 1, function(x){ rep(x[1], x[2]); }));
	
	indices <- sample(1:length(frame), size=toSize, repl=F);
	frame <- frame[indices];
	
	return(table(frame));
}

rarefyTable <- function(samples, sampleSize, nperm=1){
	
	cat("we expect a species (row) by site (col) table (not regular community matrix)\n");
	samples <- as.data.frame(samples);
	
	rarefiedMatrix <- matrix(nrow=nrow(samples), ncol=ncol(samples));
	otuNames <- rownames(rarefiedMatrix) <- rownames(samples);
	colnames(rarefiedMatrix) <- colnames(samples);
	
	for( i in 1:ncol(samples)){
		cat("Resampling from sample-number:", i, "\n");
		
		newSample <- list();
		for( j in 1:nperm ){
			newSample[[j]] <- rarefyColumn(sampleNames=otuNames, samples[,i], toSize=sampleSize);
		}
		
		targets <- sort(unique(names(unlist(newSample))));
		resultMatrix <- matrix(data=0, nrow=nperm, ncol=length(targets));
		colnames(resultMatrix) <- targets;
		
		for( k in 1:nperm ){
			result_k <- newSample[[k]];
			result_k <- result_k[sort(names(result_k))];
			resultMatrix[k, which(targets %in% names(result_k))]  <- result_k;
		}
		
		means <- apply(resultMatrix, 2, mean);

		## after correcting the otu_ids (OTU0 instead of previous 1-based version OTU0_1) this had to be updated to include the word boundary indicator '\\b' to ensure exact matching.
		indices <- unlist(lapply(paste0(names(means), "\\b"), grep, otuNames));
		rarefiedMatrix[indices, i] <- means;
	}
	
	return(rarefiedMatrix);
}

openLeafRootAlignmentFile <- function( marker ){
	
	if( marker == "ITS" ){
		combinedTissueFileName <- paste0( Sys.getenv('ROOT_MICROBIOME_DIR'), "data/ITS/fungi.combined.tissue.txt" );
		
	} else {
		combinedTissueFileName <- paste0( Sys.getenv('ROOT_MICROBIOME_DIR'), "data/16S/bacteria.combined.tissue.txt" );
	}
	
	alignmentFile <- read.table(combinedTissueFileName, header=T, sep="\t", as.is=T, stringsAsFactors=FALSE);
	
	# dropout singletons (one tissue or the other)
	noMatch <- which(is.na(alignmentFile[,"combined_tissue_id"]));
	if( length(noMatch) > 0 ){
		alignmentFile <- alignmentFile[-noMatch,];
	}

	alignmentFile <- alignmentFile[order(alignmentFile[,"combined_tissue_id"], alignmentFile[,"tissue"]),];
	return(alignmentFile);
}

alignLeavesWithRoots <- function( marker, leaves, roots, ignoreSingletons=FALSE){

	## get the alignment file map...
	map <- openLeafRootAlignmentFile( marker );

	## reduce the map to the samples provided by the user.
	names <- c(colnames(leaves), colnames(roots));
	map <- subset(map, sample_id %in% names);

	singletons <- which(table(map[,"combined_tissue_id"]) == 1);
	singletons <- which(map[,"combined_tissue_id"] %in% names(singletons));
	
	if( !ignoreSingletons & length(singletons) > 0 ){
		cat("Now excluding the", length(singletons), "samples that do not have corresponding matches in the other tissue.\n");
		map <- map[-singletons,];	
	}

	map <- map[order(map[,"combined_tissue_id"]),];
	if( is.numeric(map[,"ptp_id"])){
		## I anticipate this being treated as a factor...
		if( marker == "16S" ){
			ptpIdPrefix <- "B";

		} else {
			ptpIdPrefix <- "F";
		}

		map[,"ptp_id"] <- paste(ptpIdPrefix, map[,"ptp_id"], sep="");
	}

	# what do we need to return?
	# 1.) the map
	# 2.) the leaves and roots subsetted to the matrices of interest.
	results <- list();
	results[["map"]] <- map;
	
	## nb: we aren't subsetting to taxa where the rowcounts are greater than 0, in case we need the overlapping OTU-ids to calculate beta diversity. 
	results[["leaves"]] <- leaves[, subset(map, tissue == "L")[,"sample_id"]];
	results[["roots"]] <- roots[, subset(map, tissue == "R")[,"sample_id"]];

	return(results);
}

# tissue should be either 'leaf', 'root' or 'all'
openData <- function(qiimeOTU_Table,
		minimumReadsPerSite,
		minimumSamples=20,
		taxonomicStatus=F,
		tissue="leaf"){

	otus <- read.table(qiimeOTU_Table, header=T, sep="\t", as.is=T, stringsAsFactors=FALSE, skip=1, comment.char="", row.names=1);

	colnames(otus) <- gsub("\\.", "_", colnames(otus));

	## remove samples without ecotype-ids, if they exist...
	dropNAs <- grep("_NA\\b", colnames(otus));
	if( length(dropNAs) > 0 ){
		otus <- otus[,-dropNAs];		
	}

	## give the rows ids, using the characters OTU. :)
	tmp <- cbind(rep("OTU", nrow(otus)), rownames(otus));
	rownames(otus) <- apply(tmp, 1, paste, collapse="", sep="");

	## Omit the column storing the taxonomic data...
	if( taxonomicStatus ){
		fieldOfTaxonomicInformation <- "Consensus_Lineage";
		taxonomy <- data.frame(id=rownames(otus), lineage=otus[,fieldOfTaxonomicInformation], stringsAsFactors=FALSE);
		otus <- otus[,-which(colnames(otus) == fieldOfTaxonomicInformation)]; # omit the column that held the taxonomic data...
	}

	otus <- as.matrix(otus);

	match.arg(tissue, c("leaf", "root", "all"));
	if( tissue == "leaf" ){
		otus <- getLeaves(otus);

	} else if( tissue == "root" ){
		otus <- getRoots(otus);
	}

	# optionally: remove OTUs that aren't observed in a minimum # of samples.
	tables <- apply(otus, 1, function(x){ sum( x > 0 );});
	dropouts_belowMinimumRequiredSamples <- which(tables < minimumSamples);
	if( length(dropouts_belowMinimumRequiredSamples) > 0 ){
		otus <- otus[-dropouts_belowMinimumRequiredSamples,];		
	}

	# optionally: remove samples that have too few reads.
	counts <- colSums(otus);	
	dropouts_amplicon <- which(counts < minimumReadsPerSite);
	if( length(dropouts_amplicon) > 0 ){
		otus <- otus[,-dropouts_amplicon];	
	}

	cat("The minimum read count is now:", min(colSums(otus)),"\n");

	# if a sample is removed for failing to maintain a minimum # of reads,
	# it could unduly impact certain OTUs that are only found in these odd samples.
	otus <- removeDropouts(otus, silent=T);

	if( taxonomicStatus ){
		taxonomy <- taxonomy[which(taxonomy[,"id"] %in% rownames(otus)),];

	} else {
		taxonomy <- NULL;
	}

	return(list("data"=otus, "taxa"=taxonomy, "readcounts"=counts));                        
}

filterPossibleHostContaminants <- function( dataset, taxonomicInformation, minimumReadCount ){

	## first, handle the obvious plant material...
	indices <- which(taxonomicInformation[,"kingdom"] == "Viridiplantae");
	if( length(indices) > 0 ){
		cat("OTUs assigned to the Viridiplantae were detected. We are excluding these now...\n");
		dataset <- dataset[-indices,];
		taxonomicInformation <- taxonomicInformation[-indices,];
	}

	## next, handle the material that wasn't assigned at the kingdom level
	indices <- which(taxonomicInformation[,"kingdom"] == "Unassigned");
	if( length(indices) > 0 ){
		cat(length(indices), "OTUs lack taxonomic information; excluding these now.\n");
		dataset <- dataset[-indices,];
		taxonomicInformation <- taxonomicInformation[-indices,];
	}

	## Silva sometimes returns things that might be real or that might be mtDNA; let's exclude... 
	indices <- which(taxonomicInformation[,"family"] == "Mitochondria");
	if( length(indices) > 0 ){
		cat(length(indices), "OTUs assigned to Mitochondria were detected; excluding these now...\n");
		dataset <- dataset[-indices,];
		taxonomicInformation <- taxonomicInformation[-indices,];
	}

	## Silva sometimes returns things simultaneously labeled Chloroplast and Cyanobacteria; not sure whether it's real or not; let's exclude... 
	indices <- which(taxonomicInformation[,"class"] == "Chloroplast");
	if( length(indices) > 0 ){
		cat(length(indices), "OTUs were simultaneously labeled as Cyanobacteria/Chloroplast; excluding these now...\n");
		dataset <- dataset[-indices,];
		taxonomicInformation <- taxonomicInformation[-indices,];
	}

	# It's possible that some of the samples are now below the minimum specified read count. Exclude...
	indices <- which(colSums(dataset) < minimumReadCount );
	if( length(indices) > 0 ){
		cat("You are using a minimum read count threshold of:", minimumReadCount, "\n");
		cat("Excluding host-DNA resulted in:", length(indices), "sample(s) falling below the minimum read count threshold; these samples have been omitted.\n");					
		dataset <- dataset[,-indices];

		## were any OTUs (only) observed in the excluded samples?
		otuDropouts <- which(rowSums(dataset) == 0);
		if( length(otuDropouts) > 0 ){
			dataset <- dataset[-otuDropouts,];
			taxonomicInformation <- taxonomicInformation[-otuDropouts,];
			cat("And their rare OTUs have been excluded.\n");	
		}
	}

	return(list("dataset"=dataset, "taxa"=taxonomicInformation));
}

## this is the main function, and is responsible for opening the qiime OTU table, possibly rarefying the data, and reporting back the taxonomic information for each OTU.
## options include thinning on a minimum read count per sample, or a minimum # of samples per OTU, and excluding data that are suspect (that is, dna that looks like it came from the host's genome).
getData <- function(cutoff, marker="16S", combined=FALSE, tissue="root", minReadsPerSite=800, minSamplePerOTU=1, whetherToGetTaxa=TRUE, rarefy=TRUE, numberOfRarefactions=1, excludeHostDNA=TRUE, redo=FALSE){

	match.arg( marker, c("16S", "ITS"));
	match.arg( tissue, c("leaf", "root", "all"));

	if( combined ){
		otuTable <- paste0( Sys.getenv('ROOT_MICROBIOME_DIR'), "data/", marker, "/total/otus", cutoff, "/picked_otus", cutoff, "/all_tissue.otus.final.txt" );

	} else {
		otuTable <- paste0( Sys.getenv('ROOT_MICROBIOME_DIR'), "data/", marker, "/", tissue, "/otus", cutoff, "/picked_otus", cutoff, "/sep_tissue.otus.final.2.txt" );
	}

	cat("Core file:", otuTable, "\n");
	prefix <- paste( tissue, 
			marker, 
			minReadsPerSite, 
			paste0( minSamplePerOTU, "sampPerOTU"),
			sep="." );

	# swap the folder names around (e.g. picked_otus97 for rarefied).
	filenamePrefix <- paste0( gsub(paste0( "picked_otus", cutoff), "rarefied", dirname(otuTable)), "/", prefix );

	# if we ask for rarefaction, first check to see if we can open a pre-existing robject file...
	# if not, we will create one... (further below).
	if( rarefy & !redo & file.exists(paste(filenamePrefix, ".rare.r", numberOfRarefactions, ".robj", sep=""))){
		cat("####################################################\n");
		cat("Found a previously rarefied dataset. Opening it now!\n");
		cat("####################################################\n");

		prefix <- paste0(filenamePrefix, ".rare.r", numberOfRarefactions);
		load(paste0(prefix, ".robj"));

	} else {

		## open the data, optionally exclude host-dna (e.g. Viridplantae, Chloroplast, Mitochondria) and then - optionally - rarefy.
		openDataset <- openData(otuTable, 
				minimumReadsPerSite=minReadsPerSite, 
				minimumSamples=minSamplePerOTU, 
				taxonomicStatus=whetherToGetTaxa, 
				tissue=tissue);

		cat("Splitting the taxa into levels.\n");
		openDataset$taxa <- getTaxonomy(openDataset$taxa);

		if( excludeHostDNA ){ # this is rarefaction independent...
			resultSet <- filterPossibleHostContaminants( openDataset$data, openDataset$taxa, minReadsPerSite );
			openDataset$data <- resultSet$dataset;
			openDataset$taxa <- resultSet$taxa;
		}

		if( rarefy ){
			## this only makes sense with - at most - one resampling event.	However, we aim to please...		
			prefix <- paste0( filenamePrefix, ".rare.r", numberOfRarefactions );
			filename <- paste0( prefix, ".robj" );

			cat("####################################################\n");
			cat("Preparing to rarefy.\n");
			cat("####################################################\n");
			cat("Please note that rarefaction will be done using:", numberOfRarefactions, "permutation(s).\n");

			## if necessary, create the output directory.
			if( !dir.exists( dirname(filenamePrefix)) ){
				cat("Now creating the directory for the 'rarefied' robjects.\n");
				dir.create( dirname(filenamePrefix) ); 
			}

			# now, we rarefy...
			data <- rarefyTable(samples=openDataset$data, sampleSize=minReadsPerSite, nperm=numberOfRarefactions);
			data[which(is.na(data))] <- 0;

			# after rarefying, there will be zeros. exclude the problematic OTUs...
			otuDropouts <- which(rowSums(data) == 0);
			if( length(otuDropouts) > 0 ){
				cat( "Rarefaction led to the dropout of", length(otuDropouts), "OTUs.\n" );
				data <- data[-otuDropouts,];

				## these dropouts also affect the taxonomic information; update now...
				openDataset$taxa <- openDataset$taxa[-otuDropouts,];

			} else {
				cat("There were no dropouts after rarefaction.\nHow could this be?\n");
				stop("Please investigate.\n");
			}

			## the dataset has now been processed and is ready for the user.
			## persist the file so that we can open it faster next time.
			openDataset$data <- data;
			save(openDataset, file=filename);
			cat("After rarefaction, the table has:", nrow(openDataset$data), "OTUs.\n");

		} else {
			prefix <- paste0( filenamePrefix, ".raw" );
		}
	}

	openDataset$filename <- prefix;
	return(openDataset); # a list w/ $data and $taxa.
}

writePhenoAndBlockFile <- function( matrix, fileNamePrefix, transformation ){
	
	ids <- do.call(rbind, strsplit(rownames(matrix), "_"))[,3];
	replicate_ids <- numeric(length(ids));
	
	cat("Duplicate ids detected, splitting these into replicate ids now.\n");
	cat("\n");
	
	idsAsLevels <- as.numeric(unique(ids));
	
	for( i in idsAsLevels ){
		indices <- which(ids == i);
		replicate_ids[indices] <- 1:length(indices);
	}
	
	if( length(grep("glm", transformation)) > 0 ){
		offset <- data.frame(ecotype_id=ids, replicate_id=replicate_ids, offset=rowSums(matrix));
		
	} else {
		offset <- data.frame(ecotype_id=ids, replicate_id=replicate_ids);
	}
	
	matrix <- cbind(ids, replicate_ids, matrix);
	colnames(matrix)[1:2] <- c("ecotype_id", "replicate_id");
	offset <- offset[order(offset[,"ecotype_id"], offset[,"replicate_id"]),];
	matrix <- matrix[order(matrix[,"ecotype_id"], matrix[,"replicate_id"]),];
	
	ids <- matrix[,c(1,2)];
	
	if( ncol(matrix) > 3 ){
		phenotypes <- matrix[,-c(1,2)];
		
	} else if( ncol(matrix) == 3 ){ # there's just one phenotype.
		name <- colnames(matrix)[3];
		phenotypes <- as.numeric(matrix[,3]);
		dim(phenotypes) <- c(length(phenotypes),1);
		colnames(phenotypes)[1] <-  name # we lost the name of the phenotypes!
		
	} else {
		cat("That doesn't make sense.\n");
	}
	
	matrix.new <- data.frame(phenotype_id=c(),phenotype_name=c(),ecotype_id=c(),value=c(),replicate_id=c());
	for( i in 1:ncol(phenotypes)){
		if( i %% 100 == 0 ){ cat("Formatting phenotype file.\n"); }
		
		subset = phenotypes[,i];
		names = rep(colnames(phenotypes)[i], length(subset));
		
		matrix.new <- rbind(matrix.new,
				data.frame(
						phenotype_id=rep(i,length(subset)),
						phenotype_name=names,
						ecotype_id=ids[,"ecotype_id"],
						value=subset,
						replicate_id=ids[,"replicate_id"]));
	}
	
	fileName <- paste(fileNamePrefix,".",transformation,".txt",sep="");
	write.table( matrix.new, fileName, quote=F, sep=",", row.names=F );
	
	# get a list of blocks for this file.
	blocks <- determineBlocks( t(ids) );
	blocks <- merge(offset, blocks, all.x=T);
	blocks <- blocks[order(blocks[,"ecotype_id"], blocks[,"replicate_id"]),];
	
	blockFileName <- paste( fileNamePrefix, ".",transformation, ".blocks.txt", sep="");
	write.table( blocks, blockFileName, quote=F, sep="\t", row.names=F );
	
	# return the matrix for any downstream analyses/block creation.
	return(c(fileName, blockFileName));
}


