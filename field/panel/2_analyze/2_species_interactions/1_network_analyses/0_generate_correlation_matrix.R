## This R-script creates a correlation matrix for the leaf OR root microbiome
## focusing on a subset (e.g. top-N species) of the bacterial and fungal communities.
## 
## Author: matt.horton
###############################################################################

rm(list=ls());

require(Hmisc);
source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.microbiome_methods.R"));

###############################################################################
## metagenomic variables:
###############################################################################
possibleAnalyticalMethods <- c("lmer", "glmer", "logit");
correlationMethod <- "pearson";

organ <- "root";
otuCutoff <- 97;
numbersToTest <- 100; 

normalizationMethod16S <- normalizationMethodITS <- "raw";
minReadsBacteria <- 400;
analyticalMethodBacteria <- "glmer";
minReadsFungi <- 400;
analyticalMethodFungi <- "glmer";

###############################################################################
## end of metagenomic variables:
###############################################################################

{
	analyticalMethodBacteria <- match.arg( analyticalMethodBacteria, choices=possibleAnalyticalMethods);

	if(( normalizationMethod16S != "pa" ) & ( analyticalMethodBacteria == "logit" )){
		cat("Changing normalization method for bacteria.\n");
		normalizationMethod16S <- "pa";
		stopifnot(minReadsBacteria != 1);

	} else if(( normalizationMethod16S != "raw" ) & ( analyticalMethodBacteria == "glmer" )){
		cat("Changing normalization method for bacteria.\n");
		normalizationMethod16S <- "raw";
		stopifnot(minReadsBacteria == 1);
	}

	analyticalMethodFungi <- match.arg( analyticalMethodFungi, choices=possibleAnalyticalMethods);	
	if(( normalizationMethodITS != "pa" ) & ( analyticalMethodFungi == "logit" )){
		cat("Changing normalization method for fungi.\n");
		normalizationMethodITS <- "pa";
		stopifnot(minReadsBacteria != 1);

	} else if(( normalizationMethodITS != "raw" ) & ( analyticalMethodFungi == "glmer" )){
		cat("Changing normalization method for fungi.\n");
		normalizationMethodITS <- "raw";
		stopifnot(minReadsBacteria == 1);
	}

	###############################################################################	
	## BACTERIAL data
	###############################################################################	
	bacterial.obj <- getData(
			cutoff=otuCutoff,
			combined=FALSE,
			marker="16S",
			tissue=organ, 
			minReadsPerSite=minReadsBacteria, 
			minSamplePerOTU=1, 
			whetherToGetTaxa=TRUE, 
			rarefy=(minReadsBacteria > 1),
			numberOfRarefactions=1);

	###############################################################################	
	## FUNGAL data
	###############################################################################	
	fungal.obj <- getData(
			cutoff=otuCutoff,
			combined=FALSE,
			marker="ITS",
			tissue=organ, 
			minReadsPerSite=minReadsFungi, 
			minSamplePerOTU=1, 
			whetherToGetTaxa=TRUE, 
			rarefy=(minReadsFungi > 1),
			numberOfRarefactions=1);

	###############################################################################
	## estimate the effort/offset
	###############################################################################
	bacteria <- bacterial.obj$data;
	fungi <- fungal.obj$data;
	bacterialOffsets <- colSums(bacteria);
	fungalOffsets <- colSums(fungi);

	###############################################################################	
	## find the overlap of the two datasets!
	###############################################################################	
	cat("There are:", nrow(bacteria), "bacterial OTUs, surveyed in:", ncol(bacteria), "samples\n");
	cat("There are:", nrow(fungi), "fungal OTUs, surveyed in:", ncol(fungi), "samples\n");
	cat("Now sorting by abundance.\n");

	names <- intersect(colnames(bacteria), colnames(fungi));
	bacteria <- bacteria[,names];
	bacteria <- bacteria[order(rowSums(bacteria), decreasing=T),];
	bacterialOffsets <- bacterialOffsets[names];

	fungi <- fungi[,names];
	fungi <- fungi[order(rowSums(fungi), decreasing=T),];
	fungalOffsets <- fungalOffsets[names];
	stopifnot(colnames(bacteria) == colnames(fungi));

	###############################################################################	
	## eliminate rare OTUs with a threshold...
	## z.B. doing this with the threshold == 2 eliminates OTUs that occur in only one sample (singletons).
	###############################################################################	
	if( !is.na( normalizationMethod16S )){
		bacteria <- normalizeInSomeWay( bacteria, normalizationMethod = normalizationMethod16S, eliminateSingletons=TRUE);	
	}

	if( !is.na( normalizationMethodITS )){
		fungi <- normalizeInSomeWay( fungi, normalizationMethod = normalizationMethodITS, eliminateSingletons=TRUE);	
	}

	###############################################################################	
	################################## BACTERIA ###################################
	###############################################################################	
	{
		bacteria <- bacteria[1:300,]; ## ask for more OTUs than necessary, in case models fail to converge.
		cat("Getting the blups for the bacteria.\n");
		bacterial.blups <- getBlups( "bacteria", analyticalMethodBacteria, bacteria, bacterialOffsets );
		bacteria <- t(bacteria); ## for the original, unaveraged samples...
		colnames(bacteria) <- paste0("B_", colnames(bacteria)); ##"ecotype_id",
	}

	###############################################################################	
	#################################### FUNGI ####################################
	###############################################################################	
	{
		fungi <- fungi[1:300,]; ## ask for more OTUs than necessary, in case models fail to converge.
		cat("Getting the blups for the fungi.\n");
		fungal.blups <- getBlups( "fungi", analyticalMethodFungi, fungi, fungalOffsets );
		fungi <- t(fungi); ## for the original, unaveraged samples...
		colnames(fungi) <- paste0("F_", colnames(fungi)); ##"ecotype_id",
	}

	###############################################################################
	## move to the (parent of the) output directory
	###############################################################################	
	outputDirectory <- paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "structure/", organ, "/otus", otuCutoff, "/");
	if( !dir.exists(outputDirectory)){
		result <- dir.create(outputDirectory, recursive=TRUE);
		stopifnot(result);
		cat('Created output directory:', outputDirectory, "\n");
	}

	setwd(outputDirectory);

	{
		###############################################################################	
		## combine subsets of the bacterial and fungal communities; create the correlation matrix 
		###############################################################################	
		for( j in 1:length( numbersToTest )){

			## in the current code, the number of species is the same for both bacteria and fungi...
			nspecies_j <- numbersToTest[j];

			localDirectory <- paste0("n_", nspecies_j);
			if( !dir.exists(localDirectory)){
				result <- dir.create(localDirectory);
				stopifnot(result);
				cat("Created local directory:", localDirectory, "\n");
			}

			setwd(localDirectory);
			cat("------------------------------\n");
			cat("Creating a correlation matrix for n =", nspecies_j, "OTUs from each kingdom.\n");
			cat("------------------------------\n");


			bac.subset_j <- bacterial.blups[,1:nspecies_j];
			fun.subset_j <- fungal.blups[,1:nspecies_j];
			stopifnot(all.equal(rownames(bac.subset_j), rownames(fun.subset_j)));

			###############################################################################
			## combine everything and THEN create the correlation matrix.
			###############################################################################
			compositeCommunity <- as.matrix(cbind(bac.subset_j, fun.subset_j)); ## cbind calls data.frame, which doesn't check names by default - treat it with caution.		

			## this code needs to:
			## 1.) calculate a correlation matrix for downstream analysis/plotting...
			rcorr.obj <- rcorr( compositeCommunity, type=correlationMethod );
			corrMatrix <- rcorr.obj$r;
			pvals <- rcorr.obj$P;

			minimumReads <- ifelse( minReadsBacteria == minReadsFungi, minReadsBacteria, paste0(minReadsBacteria, "_", minReadsFungi));
			outputFileName <- paste0( correlationMethod, ".", organ, ".", minimumReads, 
										".", analyticalMethodBacteria, 
										"_bac_n", nspecies_j, ".", analyticalMethodFungi, 
										"_fun_n", nspecies_j, ".cor.txt");

			write.table(corrMatrix, file=outputFileName, quote=F, sep="\t");
			outputFileName <- gsub(".\\bcor", ".pvals", outputFileName);
			write.table(pvals, file=outputFileName, quote=F, sep="\t");

			setwd("../");
		}
	}
}