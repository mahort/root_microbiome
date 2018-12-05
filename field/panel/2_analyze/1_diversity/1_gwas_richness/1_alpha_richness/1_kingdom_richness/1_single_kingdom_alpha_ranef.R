## This R-script uses all data to estimate species richness across samples,
## using an offset to correct for sampling effort, and while fitting block and
## sequencing run (ptp_id) as technical covariates. The BLUPs are written out
## in mixmogam file format for GWAS.
## 
## Author: matt.horton
###############################################################################

rm(list=ls());

source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.microbiome_methods.R"));
source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.base_methods.R"));
source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.diversity_methods.R"));

require(lme4);

###############################################################################
## variables:
###############################################################################
otuCutoff <- 97;
organ <- "root";
minimumNumberOfReads <- 1; ## this determines whether or not the code rarefies (below); glmer allows you to specify, so we skip rarefaction for everything but simple visualization
otuThreshold <- 2;

{
	###############################################################################
	## BACTERIA: RICHNESS
	###############################################################################
	bacteria <- getData(
			cutoff=otuCutoff,
			combined=FALSE,
			marker="16S",
			tissue=organ, 
			minReadsPerSite=minimumNumberOfReads, 
			minSamplePerOTU=1, 
			whetherToGetTaxa=TRUE, 
			rarefy=(minimumNumberOfReads > 1),
			numberOfRarefactions=1)$data;

	###############################################################################
	## FUNGAL: RICHNESS
	###############################################################################
	fungi <- getData(
			cutoff=otuCutoff,
			combined=FALSE,
			marker="ITS",
			tissue=organ, 
			minReadsPerSite=minimumNumberOfReads, 
			minSamplePerOTU=1, 
			whetherToGetTaxa=TRUE, 
			rarefy=(minimumNumberOfReads > 1),
			numberOfRarefactions=1)$data;
	
	###############################################################################
	## estimate the sequencing effort, which affects the number of observed species... 
	###############################################################################
	setwd(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/gwas/otus", otuCutoff, "/", organ, "/"));
	if( !dir.exists("phenotypes")){
		dir.create("phenotypes");
		cat("Making the phenotypes directory.\n");
	}

	###############################################################################
	## determine richness/offsets per kingdom
	###############################################################################
	dataset <- determineVanillaRichness(bacteria, fungi, otuThreshold=otuThreshold);
	bacteria <- dataset$bacteria;
	fungi <- dataset$fungi;

	###############################################################################
	## using a mm/blup-approach
	###############################################################################
	bac.glm0 <- glmer( bacterial_richness ~ 1 + (1|ecotype_id) + offset(log(bacterial_effort)), data=bacteria, family=poisson);
	bac.glm1 <- glmer( bacterial_richness ~ 1 + block_id + (1|ecotype_id) + offset(log(bacterial_effort)), data=bacteria, family=poisson); 
	bac.glm2 <- glmer( bacterial_richness ~ 1 + block_id + ptp_id + (1|ecotype_id) + offset(log(bacterial_effort)), data=bacteria, family=poisson);
	anova(bac.glm0, bac.glm1, bac.glm2, test="Chisq");

	fun.glm0 <- glmer( fungal_richness ~ 1 + (1|ecotype_id) + offset(log(fungal_effort)), data=fungi, family=poisson);
	fun.glm1 <- glmer( fungal_richness ~ 1 + block_id + (1|ecotype_id) + offset(log(fungal_effort)), data=fungi, family=poisson);
	fun.glm2 <- glmer( fungal_richness ~ 1 + block_id + ptp_id + (1|ecotype_id) + offset(log(fungal_effort)), data=fungi, family=poisson);
	anova(fun.glm0, fun.glm1, fun.glm2, test="Chisq");

	{
		################################################################################
		### subject specific random effects (mit block + ptp)
		################################################################################
		bac.glm2.rf <- ranef(bac.glm2)$ecotype_id;
		bac.glm2.rf <- data.frame(ecotype_id=rownames(bac.glm2.rf), random_effects=bac.glm2.rf[,1], stringsAsFactors=F);
		outputFile <- cbind(phenotype_id=rep(1, nrow(bac.glm2.rf)), phenotype_name=rep(paste("bacterial_richness.", minimumNumberOfReads, ".tn", otuThreshold, ".ranef", sep=""), nrow(bac.glm2.rf)), 
				ecotype_id=bac.glm2.rf[,"ecotype_id"], value=bac.glm2.rf[,"random_effects"], replicate_id=rep(1, nrow(bac.glm2.rf)));
		outputFileName <- paste("phenotypes/bacterial_richness.", minimumNumberOfReads, ".tn", otuThreshold, ".ranef.txt", sep="");
		write.table(outputFile, outputFileName, row.names=F, sep=",", quote=F);

		## write out the residuals from glmer for fungi.
		fun.glm2.rf <- ranef(fun.glm2)$ecotype_id;
		fun.glm2.rf <- data.frame(ecotype_id=rownames(fun.glm2.rf), random_effects=fun.glm2.rf[,1], stringsAsFactors=F);
		outputFile <- cbind(phenotype_id=rep(1, nrow(fun.glm2.rf)), phenotype_name=rep(paste("fungal_richness.", minimumNumberOfReads, ".tn", otuThreshold, ".ranef", sep=""), nrow(fun.glm2.rf)), 
				ecotype_id=fun.glm2.rf[,"ecotype_id"], value=fun.glm2.rf[,"random_effects"], replicate_id=rep(1, nrow(fun.glm2.rf)));
		outputFileName <- paste("phenotypes/fungal_richness.", minimumNumberOfReads, ".tn", otuThreshold, ".ranef.txt", sep="");
		write.table(outputFile, outputFileName, row.names=F, sep=",", quote=F);
	}

} ## 
