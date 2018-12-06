## This R-script 
## 
## Author: matt.horton
###############################################################################

rm(list=ls());

source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.microbiome_methods.R"));
source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.base_methods.R"));

require(vegan); require(lme4);

###############################################################################
## metagenomic variables:
###############################################################################
numberOfPCsMicrobes <- 5;
organ <- "root";
phylogeneticMarker <- "ITS";
otuCutoff <- 97;
minimumReads <- 1;
otuThreshold <- 2;

normalizationMethod <- "raw";
analyticalMethod <- "glmer";
ordinationMethod <- "pca";
###############################################################################
## end of metagenomic variables. 
###############################################################################

{
	stopifnot( normalizationMethod != "hlgr" ); ## double-negative
	possibleAnalyticalMethods <- c("lmer", "glmer", "logit" );
	analyticalMethod <- match.arg( analyticalMethod, choices=possibleAnalyticalMethods);
	ordinationMethod <- match.arg( ordinationMethod, c("pca", "cca", "dca"));

	if( normalizationMethod != "pa" & analyticalMethod == "logit" ){
		cat("Changing normalization method for ", phylogeneticMarker, ".\n");
		normalizationMethod <- "pa";
		stopifnot(minimumReads == 250);

	} else if( normalizationMethod != "raw" & analyticalMethod == "glmer" ){
		normalizationMethod <- "raw";
		cat("Changing normalization method for ", phylogeneticMarker, " to: ", normalizationMethod, "\n");
		stopifnot(minimumReads == 1);
	}

	###############################################################################
	## OPEN the dataset
	###############################################################################
	listOfData <- getData(
			cutoff=otuCutoff,
			combined=FALSE,
			marker=phylogeneticMarker,
			tissue=organ, 
			minReadsPerSite=minimumReads, 
			minSamplePerOTU=1, 
			whetherToGetTaxa=TRUE, 
			rarefy=(minimumReads > 1),
			numberOfRarefactions=1);

	speciesMatrix <- listOfData$data;
	filenamePrefix <- basename(listOfData$filename);

	###############################################################################
	## estimate the effort/offset
	###############################################################################
	data.offsets <- colSums( speciesMatrix );
	kingdom <- ifelse( phylogeneticMarker == "16S", "bacteria", "fungi" );
	cat("Analyzing:", nrow(speciesMatrix), "OTUs, surveyed in:", ncol(speciesMatrix), "samples\n");

	transformation <- ordinationMethod; # begin writing the suffix.

	###############################################################################
	## eliminate rare OTUs with a threshold...
	## z.B. doing this with the threshold == 2 eliminates OTUs that occur in only one sample (singletons).
	###############################################################################
	if( !is.na( normalizationMethod )){
		speciesMatrix <- normalizeInSomeWay( speciesMatrix, normalizationMethod = normalizationMethod, eliminateSingletons=(otuThreshold > 1));
	}

	filenamePrefix <- gsub("1sampPerOTU\\.", "", filenamePrefix);

	###############################################################################
	##
	###############################################################################
	numbersToTest <- sort(c(25, 50, 100 )); #ceiling(speciesCutoffs * nrow(speciesMatrix))));
	speciesMatrix <- speciesMatrix[order(rowSums(speciesMatrix), decreasing=T),];
	speciesMatrix <- speciesMatrix[1:max(numbersToTest),];

	###############################################################################	
	################################## GET BLUPS ##################################
	###############################################################################	
	{
		cat("------------------------------\n");
		speciesMatrix <- getBlups( kingdom, analyticalMethod, speciesMatrix, data.offsets );

		## if lots of tests fail, we have to reduce this variable to account for the shrinkage!!!
		if(ncol(speciesMatrix) < max(numbersToTest)){ 
			numbersToTest[which(numbersToTest > ncol(speciesMatrix))] <- ncol(speciesMatrix); ## in case some of the models didn't converge.
			numbersToTest <- unique(numbersToTest);
		}
	}

	###############################################################################	
	## Move to the output directory.
	###############################################################################
	setwd(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/gwas/otus", otuCutoff, "/", organ, "/phenotypes/")); ## 
	cat("Performing", ordinationMethod, "\n");

	formatEigenvectorsForGWA <- function( eigenvectorMatrix, numberOfEigenvectors=numberOfPCs){
		outputMatrix <- data.frame(phenotype_id=numeric(), phenotype_name=character(), ecotype_id=numeric(), value=numeric(), replicate_id=numeric());

		for( j in 1:numberOfEigenvectors ){
			phenotype_j <- mstack(eigenvectorMatrix[,j], newHeaders=c("ecotype_id", "value"), sorted=F);

			eig_j <- data.frame(
					phenotype_id=rep(j, nrow(phenotype_j)), 
					phenotype_name=rep(colnames(eigenvectorMatrix)[j], nrow(phenotype_j)), 
					ecotype_id=phenotype_j[,"ecotype_id"], value=phenotype_j[,"value"], replicate_id=rep(1, nrow(phenotype_j)), stringsAsFactors=FALSE);

			outputMatrix <- rbind(outputMatrix, eig_j);
		}

		return(outputMatrix);
	}

	###############################################################################	
	## now that things are in order: let us ordinate and write out the phenotype file (mixmogam file format)
	###############################################################################	
	for( j in 1:length(numbersToTest)){
		nspecies_j <- numbersToTest[j];
		cat("Conducting", ordinationMethod, "using the top:", nspecies_j, "'species' in the", kingdom, "\n");
		subset_j <- speciesMatrix[,1:nspecies_j];
		suffix <- paste0(analyticalMethod, ".", ordinationMethod, ".n", nspecies_j);

		if( ordinationMethod == "pca" ){
			pca.obj <- prcomp( subset_j, scale=TRUE, center=TRUE);
			u <- pca.obj$x[,1:numberOfPCsMicrobes];

		} else {
			stop(paste0( ordinationMethod, "is not supported yet.\n"));
		}

		colnames(u) <- paste0(filenamePrefix, ".tn", otuThreshold, ".", suffix, "_", colnames(u));

		outputFile <- formatEigenvectorsForGWA(u, numberOfEigenvectors=numberOfPCsMicrobes);
		outputFileName <- paste0(filenamePrefix, ".tn", otuThreshold, ".", suffix, ".txt");
		write.table(outputFile, file=outputFileName, quote=F, sep=",", row.names=F);
	}
}