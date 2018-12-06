# TODO: Add comment
# 
# Author: matt.horton
###############################################################################

rm(list=ls());

source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.microbiome_methods.R"));
source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.base_methods.R"));
source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.diversity_methods.R"));

require(lme4); 

{
	###############################################################################
	## variables:
	###############################################################################
	otuCutoff <- 97;
	organ <- "root";
	minimumNumberOfReads <- 1; ## this determines whether or not the code rarefies (below)
	otuThreshold <- 2;

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
	## restrict our focus to samples that have been scored for both bacteria and fungi...
	###############################################################################
	cat("Now reducing the samples to only those that overlap!\n");
	names <- intersect(colnames(bacteria), colnames(fungi));
	bacteria <- bacteria[,names];
	fungi <- fungi[,names];
	stopifnot(colnames(bacteria) == colnames(fungi));

	############################################################################################################################################
	## determine richness in the community, bringing in important covariates (sequencing effort, sequencing run, block, etc.) 
	############################################################################################################################################
	dataset <- determineVanillaRichness(bacteria, fungi, otuThreshold=otuThreshold);
	bacteria <- data.frame(dataset$bacteria, kingdom="bacteria", stringsAsFactors=FALSE);
	fungi <- data.frame(dataset$fungi, kingdom="fungi", stringsAsFactors=FALSE);

	colnames(bacteria) <- gsub("bacterial_", "", colnames(bacteria));
	colnames(fungi) <- gsub("fungal_", "", colnames(fungi));

	###############################################################################
	## merge the total community
	###############################################################################
	compositeCommunity <- rbind(bacteria, fungi);
	compositeCommunity[,"sample_id"] <- as.character(compositeCommunity[,"sample_id"]);	
	compositeCommunity <- transform(compositeCommunity, kingdom=as.factor(kingdom));

	###############################################################################
	## let's use a glmer/blup approach.
	###############################################################################
	glmer0 <- glmer(richness ~ 1 + block_id + (1|ecotype_id) + offset(log(effort)), family=poisson, data=compositeCommunity);
	glmer1 <- glmer(richness ~ 1 + block_id + kingdom + (1|ecotype_id) + offset(log(effort)), family=poisson, data=compositeCommunity); ##
	glmer2 <- glmer(richness ~ 1 + block_id + kingdom/ptp_id + (1|ecotype_id) + offset(log(effort)), family=poisson, data=compositeCommunity); ## this is robust against prepending the sequencing run ids with the name; in the current code (in which we do prepend), this is the same as a model without the 'kingdom/', as these sequencing runs are separated among bacteria and fungi - I leave it in here for future situations where it might not be true (our Fragaria work, z.B.). 
	anova(glmer0, glmer1, glmer2, test="Chisq");

	setwd(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/gwas/otus", otuCutoff, "/", organ, "/"));
	if( !file.exists("phenotypes")){
		command <- paste("mkdir phenotypes", sep="");
		system(command, wait=T);
		cat("Making the phenotypes directory.\n");
	}

	############################################################################################################################################
	## model average richness in the microbiome, accounting for differences in kingdom and sequencing runs (nested in kingdom)
	############################################################################################################################################
	{
		############################################################################################################################################
		## controlling for block, kingdom, and ptp_id as fixed effects, ecotype_id as a random effect...
		############################################################################################################################################	
		poissonRandomEffects <- ranef(glmer2)$ecotype_id;
		poissonRandomEffects <- data.frame(ecotype_id=rownames(poissonRandomEffects), value=poissonRandomEffects[,1], stringsAsFactors=FALSE);
		outputFile <- cbind(phenotype_id=rep(1, nrow(poissonRandomEffects)), 
				phenotype_name=rep(paste("total.", minimumNumberOfReads, ".tn", otuThreshold, ".ranef_glmer", sep=""), nrow(poissonRandomEffects)), 
				ecotype_id=poissonRandomEffects[,"ecotype_id"], value=poissonRandomEffects[,"value"], replicate_id=rep(1, nrow(poissonRandomEffects)));
		outputFileName <- paste("phenotypes/total.", minimumNumberOfReads, ".tn", otuThreshold, ".ranef_glmer.txt", sep="");
		write.table(outputFile, outputFileName, row.names=F, sep=",", quote=F);
	}
}