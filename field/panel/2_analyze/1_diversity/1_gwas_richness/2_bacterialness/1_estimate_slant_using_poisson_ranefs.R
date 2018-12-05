## TODO: Add comment
## 
## Author: matt.horton
###############################################################################
## the following code considers the community in a poisson model, to allow the individual richness estimates to be attached to the covariates more effectively \
## (to put that another way, a cbind'ed binomial obfuscates over the fact that the successes and failures may arise from difft. sampling distributions)
rm(list=ls());

source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.microbiome_methods.R"));
source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.base_methods.R"));
source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.diversity_methods.R"));

require(lme4); 

{
	#####################################################################################################################
	## variables:
	#####################################################################################################################
	otuCutoff <- 97;
	organ <- "root";
	minimumNumberOfReads <- 1;
	otuThreshold <- 2;

	############################################################################################################################################
	## BACTERIA: RICHNESS
	############################################################################################################################################
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

	############################################################################################################################################
	## FUNGAL: RICHNESS
	############################################################################################################################################
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

	############################################################################################################################################
	## restrict our focus to samples that have been scored for both bacteria and fungi...
	############################################################################################################################################
	names <- intersect(colnames(bacteria), colnames(fungi));
	bacteria <- bacteria[,names];
	fungi <- fungi[,names];
	stopifnot(all.equal(colnames(bacteria), colnames(fungi)));

	############################################################################################################################################
	## determine richness/offsets per kingdom
	############################################################################################################################################
	dataset <- determineVanillaRichness(bacteria, fungi, otuThreshold=otuThreshold);
	bacteria <- data.frame(dataset$bacteria, kingdom="bacteria", stringsAsFactors=FALSE);
	fungi <- data.frame(dataset$fungi, kingdom="fungi", stringsAsFactors=FALSE);

	colnames(bacteria) <- gsub("bacterial_", "", colnames(bacteria));
	colnames(fungi) <- gsub("fungal_", "", colnames(fungi));
	
	############################################################################################################################################
	## merge the total community
	############################################################################################################################################
	compositeCommunity <- rbind(bacteria, fungi);
	compositeCommunity[,"sample_id"] <- as.character(compositeCommunity[,"sample_id"]);	
	compositeCommunity <- transform(compositeCommunity, kingdom=as.factor(kingdom)); ## long format...

	############################################################################################################################################
	## let's use a blup approach with glmer.
	## Douglas Bates recommends the following formulation for this scenario (categorical variable kingdom):
	## https://stat.ethz.ch/pipermail/r-sig-mixed-models/2009q1/001736.html 
	############################################################################################################################################
	glmer0 <- glmer(richness ~ 1 + (0 + kingdom|ecotype_id) + offset(log(effort)), family=poisson, data=compositeCommunity);
	glmer1 <- glmer(richness ~ 1 + block_id + (0 + kingdom|ecotype_id) + offset(log(effort)), family=poisson, data=compositeCommunity);
	glmer2 <- glmer(richness ~ 1 + block_id + kingdom + (0 + kingdom|ecotype_id) + offset(log(effort)), family=poisson, data=compositeCommunity); ##
	glmer3 <- glmer(richness ~ 1 + block_id + kingdom/ptp_id + (0 + kingdom|ecotype_id) + offset(log(effort)), family=poisson, data=compositeCommunity); ## nest ptp in kingdom for clarity, but I now prepend the ptp ids with kingdom letters B and F making it unnecessary (but harmless)

	anova(glmer0, glmer1, glmer2, glmer3, test="Chisq"); ## 

	setwd(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/gwas/otus", otuCutoff, "/", organ, "/")); 
	if( !file.exists("phenotypes")){
		command <- paste("mkdir phenotypes", sep="");
		system(command, wait=T);
		cat("Making directory.\n");
	}

	############################################################################################################################################
	## write out the relative-ness of one kingdom or another
	############################################################################################################################################
	{
		############################################################################################################################################
		## random effects, block_ids and ptp_id as a fixed_effect
		############################################################################################################################################
		poissonRandomEffects <- ranef(glmer3)$ecotype_id;
		poissonRandomEffects <- data.frame(poissonRandomEffects, e1=exp(poissonRandomEffects[,1]), e2=apply(poissonRandomEffects, 1, function(x){ prod(exp(x)); }));
		poissonRandomEffects <- transform(poissonRandomEffects, slant=(e1-e2));
		poissonRandomEffects <- data.frame(ecotype_id=rownames(poissonRandomEffects), slant=poissonRandomEffects[,"slant"], stringsAsFactors=FALSE);
		outputFile <- cbind(phenotype_id=rep(1, nrow(poissonRandomEffects)), 
				phenotype_name=rep(paste("slant.", minimumNumberOfReads, ".tn", otuThreshold, ".ranef_diff_exps", sep=""), nrow(poissonRandomEffects)), 
				ecotype_id=poissonRandomEffects[,"ecotype_id"], value=poissonRandomEffects[,"slant"], replicate_id=rep(1, nrow(poissonRandomEffects)));
		outputFileName <- paste("phenotypes/slant.", minimumNumberOfReads, ".tn", otuThreshold, ".ranef_diff_exps.txt", sep=""); ## 
		write.table(outputFile, outputFileName, row.names=F, sep=",", quote=F);
	}
}

##
#> anova(glmer3, glm)
#Data: compositeCommunity
#Models:
#		glm: richness ~ 1 + block_id + kingdom/ptp_id + offset(log(effort))
#glmer3: richness ~ 1 + block_id + kingdom/ptp_id + (0 + kingdom | ecotype_id) + 
#		glmer3:     offset(log(effort))
#Df   AIC   BIC   logLik deviance Chisq Chi Df Pr(>Chisq)    
#glm     8 24893 24934 -12438.5    24877                            
#glmer3 11 19143 19199  -9560.5    19121  5756      3  < 2.2e-16 ***
#		---

