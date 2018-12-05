## TODO: Add comment
## 
## Author: mhorto
###############################################################################

whittakersBeta <- function( pairOfSites ){ ## we assume this is in site (row) by species (columns) format.
	pairOfSites <- pairOfSites[,which( colSums(pairOfSites) > 0 )];
	gamma <- ncol(pairOfSites);
	alpha <- mean(apply(pairOfSites, 1, function(x){ sum(x > 0); })); ## omitting specnumber makes it unnec. to load the entire vegan library
	return((gamma/alpha) - 1);
}

## this estimates richness and offsets (in case there's no rarefaction) for each kingdom.
determineVanillaRichness <- function( bacterialSpeciesMatrix, fungalSpeciesMatrix, otuThreshold=1, mergeMatrices=FALSE){

	############################################################################################################################################
	## determine the offsets before estimating richness metrics.
	############################################################################################################################################
	bacterialOffsets <- apply(bacterialSpeciesMatrix, 2, sum);
	bacterialOffsets <- mstack(bacterialOffsets, newHeaders=c("sample_id", "bacterial_effort"), sorted=FALSE);
	
	fungalOffsets <- apply(fungalSpeciesMatrix, 2, sum);
	fungalOffsets <- mstack(fungalOffsets, newHeaders=c("sample_id", "fungal_effort"), sorted=FALSE);
	
	############################################################################################################################################
	## eliminate rare OTUs with a threshold... need to do this after determining the offsets.
	############################################################################################################################################
	bacterialSpeciesMatrix <- applyThreshold( bacterialSpeciesMatrix, otuThreshold );
	fungalSpeciesMatrix <- applyThreshold( fungalSpeciesMatrix, otuThreshold );
	
	############################################################################################################################################
	## determine richness in the bacterial community
	############################################################################################################################################
	bacterialRichness <- apply(bacterialSpeciesMatrix, 2, function(x){ sum( x > 0 )});
	bacterialRichness <- mstack(bacterialRichness, newHeaders=c("sample_id", "bacterial_richness"), sorted=FALSE);
	bacteria <- merge(bacterialRichness, bacterialOffsets, by="sample_id");
	bacteria[,"sample_id"] <- as.character(bacteria[,"sample_id"]);
	
	############################################################################################################################################
	## determine richness in the fungal community
	############################################################################################################################################
	fungalRichness <- apply(fungalSpeciesMatrix, 2, function(x){ sum( x > 0 )});
	fungalRichness <- mstack(fungalRichness, newHeaders=c("sample_id", "fungal_richness"), sorted=FALSE);
	fungi <- merge(fungalRichness, fungalOffsets, by="sample_id");
	fungi[,"sample_id"] <- as.character(fungi[,"sample_id"]);

	if( mergeMatrices ){ ## long format... usually a glm* in long format will be more effective, as it allows us to take into account subject specific covariates
		dataset <- merge(bacteria, fungi, by="sample_id");
		rownames(dataset) <- dataset[,"sample_id"];
		
		############################################################################################################################################
		## determine factors and combine the dataset.
		############################################################################################################################################
		ecotypeIds <- do.call(rbind, strsplit(rownames(dataset), "_"))[,3];
		blockIds <- getBlockData(t(dataset));
		plateIds <- determinePlates(t(dataset));
		dataset <- data.frame(ecotype_id=as.factor(ecotypeIds), block_id=as.factor(blockIds), plate_id=as.factor(plateIds), dataset, stringsAsFactors=F);
		return(dataset);

	} else {
		## return the separate matrices.
		rownames(bacteria) <- bacteria[,"sample_id"];
		covariates <- determineBlocks(t(bacteria));
		stopifnot(all.equal(covariates[,"sample_id"], bacteria[,"sample_id"]));

		bacteria <- data.frame(ecotype_id=as.factor(covariates[,"ecotype_id"]), 
							block_id=as.factor(covariates[,"block_id"]), 
							plate_id=as.factor(paste0("B", covariates[,"plate_id"])), 
							ptp_id=as.factor(paste0("B", covariates[,"ptp_bacteria"])), 
							bacteria, stringsAsFactors=F);

					
		rownames(fungi) <- fungi[,"sample_id"];
		covariates <- determineBlocks(t(fungi));
		stopifnot(all.equal(covariates[,"sample_id"], fungi[,"sample_id"]));
		
		fungi <- data.frame(ecotype_id=as.factor(covariates[,"ecotype_id"]), 
							block_id=as.factor(covariates[,"block_id"]), 
							plate_id=as.factor(paste0("F", covariates[,"plate_id"])), 
							ptp_id=as.factor(paste0("F", covariates[,"ptp_fungi"])), 
							fungi, stringsAsFactors=F);

		return(list("bacteria"=bacteria, "fungi"=fungi));
	}
}
