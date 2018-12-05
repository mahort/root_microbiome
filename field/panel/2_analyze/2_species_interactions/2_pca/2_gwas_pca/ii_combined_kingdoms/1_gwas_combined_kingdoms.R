# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

rm(list=ls());

require(vegan); require(lme4);

source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.microbiome_methods.R"));
source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.base_methods.R"));

###############################################################################
## metagenomic variables:
###############################################################################
centerPCA <- scalePCA <- TRUE;
possibleAnalyticalMethods <- c("lmer", "glmer", "logit");

numberOfPCsMicrobes <- 5;
organ <- "root";

otuCutoff <- 97;
otuThreshold <- 2;
minReadsBacteria <- minReadsFungi <- 1;

ordinationMethod <- "pca";

normalizationMethod16S <- normalizationMethodITS <- "raw";
analyticalMethodBacteria <- "glmer";
analyticalMethodFungi <- "glmer"; ## logistic regression requires the argument `logit`, which is a binomial family instead of a poisson family model (glmer)

###############################################################################
## end of metagenomic variables:
###############################################################################

{
	analyticalMethodBacteria <- match.arg( analyticalMethodBacteria, choices=possibleAnalyticalMethods);
	if( normalizationMethod16S != "pa" & analyticalMethodBacteria == "logit" ){
		cat("Changing normalization method for bacteria.\n");
		normalizationMethod16S <- "pa";
		stopifnot(minReadsBacteria == 250);
		
	} else if( normalizationMethod16S != "raw" & analyticalMethodBacteria == "glmer" ){
		cat("Changing normalization method for bacteria.\n");
		normalizationMethod16S <- "raw";
		stopifnot(minReadsBacteria == 1);
	}
	
	analyticalMethodFungi <- match.arg( analyticalMethodFungi, choices=possibleAnalyticalMethods);	
	if( normalizationMethodITS != "pa" & analyticalMethodFungi == "logit" ){
		cat("Changing normalization method for fungi.\n");
		normalizationMethodITS <- "pa";
		minReadsFungi <- 250;

	} else if( normalizationMethodITS != "raw" & analyticalMethodFungi == "glmer" ){
		cat("Changing normalization method for fungi.\n");
		normalizationMethodITS <- "raw";
		stopifnot(minReadsBacteria == 1);
	}

	###############################################################################
	## BACTERIAL data
	###############################################################################
	bacteria <- getData(
			cutoff=otuCutoff,
			combined=FALSE,
			marker="16S",
			tissue=organ, 
			minReadsPerSite=minReadsBacteria, 
			minSamplePerOTU=1, 
			whetherToGetTaxa=TRUE, 
			rarefy=(minReadsBacteria > 1),
			numberOfRarefactions=1)$data;

	###############################################################################
	## FUNGAL data
	###############################################################################
	fungi <- getData(
			cutoff=otuCutoff,
			combined=FALSE,
			marker="ITS",
			tissue=organ, 
			minReadsPerSite=minReadsFungi, 
			minSamplePerOTU=1, 
			whetherToGetTaxa=TRUE, 
			rarefy=(minReadsFungi > 1),
			numberOfRarefactions=1)$data;

	## if we support minReads == 1, offsets need to be calculated before normalization
	bacterialOffsets <- colSums(bacteria);
	fungalOffsets <- colSums(fungi);

	###############################################################################
	## find the overlap of the two datasets!
	###############################################################################
	cat("There are:", nrow(bacteria), "bacterial OTUs, surveyed in:", ncol(bacteria), "samples\n");
	cat("There are:", nrow(fungi), "fungal OTUs, surveyed in:", ncol(fungi), "samples\n");
	cat("Now sorting by abundance.\n");

	names <- intersect(colnames(bacteria), colnames(fungi));
	bacteria <- bacteria[, names];
	bacteria <- bacteria[order(rowSums(bacteria), decreasing=T),];
	bacterialOffsets <- bacterialOffsets[names];

	fungi <- fungi[, names];
	fungi <- fungi[order(rowSums(fungi), decreasing=T),];
	fungalOffsets <- fungalOffsets[names];
	stopifnot(all.equal(colnames(bacteria), colnames(fungi)));

	###############################################################################
	## eliminate rare OTUs with a threshold...
	## z.B. doing this with the threshold == 2 eliminates OTUs that occur in only one sample (singletons).
	###############################################################################
	if( !is.na( normalizationMethod16S )){
		bacteria <- normalizeInSomeWay( bacteria, normalizationMethod = normalizationMethod16S, eliminateSingletons=(otuThreshold > 1));	
	}

	if( !is.na( normalizationMethodITS )){
		fungi <- normalizeInSomeWay( fungi, normalizationMethod = normalizationMethodITS, eliminateSingletons=(otuThreshold > 1));	
	}

	###############################################################################
	## get the blups for the bacteria and fungi
	###############################################################################
	numbersToTestBacteria <- 100; #ceiling(speciesCutoffs * nrow(bacteria))));
	bacteria <- bacteria[1:max(numbersToTestBacteria),];

	numbersToTestFungi <- c(25, 100); #ceiling(speciesCutoffs * nrow(fungi))));
	fungi <- fungi[1:max(numbersToTestFungi),];

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
	################################## BACTERIA ###################################
	###############################################################################	
	{
		cat("Getting the blups for the bacteria.\n");
		bacteria <- getBlups( "bacteria", analyticalMethodBacteria, bacteria, bacterialOffsets );

		## if tests fail, we may have to reduce this variable to account for the shrinkage!!!
		if(ncol(bacteria) < numbersToTestBacteria[length(numbersToTestBacteria)]){
			numbersToTestBacteria[which(numbersToTestBacteria > ncol(bacteria))] <- ncol(bacteria); ## in case some of the models didn't converge.
			numbersToTestBacteria <- unique(numbersToTestBacteria);
		}
	}

	###############################################################################	
	#################################### FUNGI ####################################
	###############################################################################	
	{
		cat("Getting the blups for the fungi.\n");
		fungi <- getBlups( "fungi", analyticalMethodFungi, fungi, fungalOffsets );

		## if tests fail, we may have to reduce this variable to account for the shrinkage!!!
		if(ncol(fungi) < numbersToTestFungi[length(numbersToTestFungi)]){
			numbersToTestFungi[which(numbersToTestFungi > ncol(fungi))] <- ncol(fungi); ## in case some of the models didn't converge.
			numbersToTestFungi <- unique(numbersToTestFungi);
		}
	}

	{
		################################################################################################################################################
		# specify the mapping directory.
		################################################################################################################################################
		setwd(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/gwas/otus", otuCutoff, "/", organ, "/phenotypes/")); ## 
		minimumReads <- paste( unique(c( minReadsBacteria, minReadsFungi)), sep="", collapse="_"); rm(minReadsBacteria); rm(minReadsFungi);

		###############################################################################
		## combine subsets of the bacterial and fungal communities
		## to conduct GWAS, write out the phenotype file.
		###############################################################################
		cat("------------------------------\n");
		cat("------------------------------\n");
		cat("Performing", ordinationMethod, "\n");
		cat("------------------------------\n");
		cat("------------------------------\n");

		for( j in 1:length( numbersToTestBacteria )){
			cat("------------------------------\n");
			nspecies_bacteria <- numbersToTestBacteria[j];
			bac.subset_j <- bacteria[,1:nspecies_bacteria];			
			cat("Working on bacterial sample size:", nspecies_bacteria, "\n");

			for( k in 1:length( numbersToTestFungi )){
				nspecies_fungi <- numbersToTestFungi[k];
				cat("Working on fungal sample size:", nspecies_fungi, "\n");

				fun.subset_k <- fungi[,1:nspecies_fungi];

				###############################################################################
				## now combine everything and THEN conduct PCA, and write out the phenotype file
				###############################################################################
				## cbind calls data.frame, which doesn't check names by default
				stopifnot(sum(rownames(bac.subset_j) == rownames(fun.subset_k)) == nrow(bac.subset_j));
				compositeCommunity <- cbind(bac.subset_j, fun.subset_k); 

				compositeCommunity.pca <- prcomp(compositeCommunity, center=centerPCA, scale=scalePCA);
				u <- compositeCommunity.pca$x[,1:numberOfPCsMicrobes];

				filenamePrefix <- paste0(organ, ".", minimumReads,  ".combined_", ordinationMethod, ".tn", otuThreshold, ".bac_n", nspecies_bacteria, "_", analyticalMethodBacteria, ".fun_n", nspecies_fungi, "_", analyticalMethodFungi, "");
				colnames(u) <- paste0(filenamePrefix, colnames(u));

				outputFile <- formatEigenvectorsForGWA(u, numberOfEigenvectors=numberOfPCsMicrobes);
				outputFileName <- paste0(filenamePrefix, ".txt");
				write.table(outputFile, file=outputFileName, quote=F, sep=",", row.names=F);
				cat("------------------------------\n");
			}
		}
	}
}
