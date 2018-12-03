## This R-script uses a poisson glm to identify those taxa that differ
## in abundance across the leaves and roots of thaliana...
## 
## Author: matt.horton
###############################################################################

rm(list=ls());

source("/Volumes/projects/thaliana/code/hdr.microbiome_methods.R");
source("/Volumes/projects/code/hdr.base_methods.R");

#####################################################################################################################
## parameters
#####################################################################################################################
otuCutoff <- 97;
otuThreshold <- 2; ## eliminate singletons if you want... nb: this is really only valid for values of 1 or 2, given the way we treat it as a boolean below...
phylogeneticTarget <- "ITS"; ## one of either 16S or ITS
minReads <- 1;
minimumSumEitherOrgan <- 1000;

{
	############################################################################################################################################
	## open the leaf and root data...
	############################################################################################################################################
	leaf.robj <- getData(
			cutoff=otuCutoff,
			combined=FALSE,
			marker=phylogeneticTarget,
			tissue="leaf", 
			minReadsPerSite=minReads, 
			minSamplePerOTU=1, 
			whetherToGetTaxa=TRUE, 
			rarefy=(minReads > 1),
			numberOfRarefactions=1);

	root.robj <- getData(
			cutoff=otuCutoff,
			combined=FALSE,
			marker=phylogeneticTarget,
			tissue="root", 
			minReadsPerSite=minReads, 
			minSamplePerOTU=1, 
			whetherToGetTaxa=TRUE, 
			rarefy=(minReads > 1),
			numberOfRarefactions=1);

	################################################################################################################################################
	## the data and taxonomic info...
	################################################################################################################################################
	leaf.data <- leaf.robj$data;
	root.data <- root.robj$data;

	leaf.offsets <- data.frame( mstack(colSums(leaf.data), sorted=FALSE, newHeaders=c("sample_id", "offset")), organ="leaf", stringsAsFactors=FALSE );
	root.offsets <- data.frame( mstack(colSums(root.data), sorted=FALSE, newHeaders=c("sample_id", "offset")), organ="root", stringsAsFactors=FALSE );

	leaf.taxa <- leaf.robj$taxa;
	root.taxa <- root.robj$taxa;

	rownames(leaf.taxa) <- leaf.taxa[,1];
	rownames(root.taxa) <- root.taxa[,1];

	#####################################################################################################################
	## threshold OTUs/ASVs based on 
	#####################################################################################################################
	leaf.data <- applyThreshold( leaf.data, otuThreshold );
	root.data <- applyThreshold( root.data, otuThreshold );
	kingdom <- ifelse(phylogeneticTarget == "16S", "bacteria", "fungi");
	cat("There are:", nrow(leaf.data), " leaf ", kingdom, " using an otu-threshold of: ", otuThreshold, "\n", sep="");
	cat("There are:", nrow(root.data), " root ", kingdom, " using an otu-threshold of: ", otuThreshold, "\n", sep="");

	#####################################################################################################################
	## sum them up and calculate the differences...
	#####################################################################################################################
	leaf.data <- leaf.data[order(rowSums(leaf.data), decreasing=T),];
	root.data <- root.data[order(rowSums(root.data), decreasing=T),];
	leaf.taxa <- leaf.taxa[rownames(leaf.data),];
	root.taxa <- root.taxa[rownames(root.data),];

	formatTaxonomicLevel <- function( organ.indices, organ.dataset, organ.offset){
		## we need to handle 0 (empty), 1 (no columns), or more (well-behaved) indices...
		if( length(organ.indices) > 1 ){
			sums <- merge(mstack(colSums(organ.dataset[organ.indices,]), sorted=FALSE, newHeaders=c("sample_id", "sum")), organ.offset, by="sample_id");

		} else if( length(organ.indices) == 1 ){
			sums <- merge(mstack(organ.dataset[organ.indices,], sorted=FALSE, newHeaders=c("sample_id", "sum")), organ.offset, by="sample_id");

		} else if( length(organ.indices) == 0 ){
			sums <- data.frame(sum=0, organ.offset, stringsAsFactors=FALSE);
		}

		return(sums);
	}

	## from phylum to genus...
	taxonomicLevels <- c("phylum", "class", "order", "family", "genus");
	results <- data.frame(order=numeric(), kingdom=character(), level=character(), habitat=character(), name=character(), estimate_leaf=numeric(), estimate_root=numeric(), root_leaf_prop=numeric(), pvalue=numeric(), otus_in_leaf=numeric(), otus_in_root=numeric(), reads_leaf=numeric(), reads_root=numeric(), stringsAsFactors=FALSE);
	for( j in 1:length(taxonomicLevels)){

		level_j <- taxonomicLevels[j];
		columnIndex_leaf <- which(colnames(leaf.taxa) == level_j);
		columnIndex_root <- which(colnames(root.taxa) == level_j);
		stopifnot(columnIndex_leaf == columnIndex_root);

		uniqueTaxa <- unique(c(leaf.taxa[,columnIndex_leaf], root.taxa[,columnIndex_leaf]));
		cat("Investigating rank:", level_j, "\n");

		taxonomicResults <- results[0,]; ## data.frame(level=character(), name=character(), estimate_leaf=numeric(), estimate_root=numeric(), root_leaf_prop=numeric(), pvalue=numeric(), stringsAsFactors=FALSE);
		for( k in 1:length(uniqueTaxa)){
			taxon_k <- uniqueTaxa[k];
			leaf.indices <- which(leaf.taxa[,columnIndex_leaf] == taxon_k);
			leaf.sums <- formatTaxonomicLevel(organ.indices = leaf.indices, organ.dataset = leaf.data, organ.offset=leaf.offsets );
			root.indices <- which(root.taxa[,columnIndex_root] == taxon_k);
			root.sums <- formatTaxonomicLevel(organ.indices = root.indices, organ.dataset = root.data, organ.offset=root.offsets );
			if( sum(leaf.sums[,"sum"]) > minimumSumEitherOrgan | sum(root.sums[,"sum"]) > minimumSumEitherOrgan ){ 
				all.sums <- rbind(leaf.sums, root.sums);
				glm0 <- glm(sum ~ 1, data=all.sums, offset=log(offset), family=poisson);
				glm1 <- glm(sum ~ organ, data=all.sums, offset=log(offset), family=poisson);
				
				taxonomicResults[nrow(taxonomicResults) + 1, "kingdom"] <- kingdom;
				taxonomicResults[nrow(taxonomicResults), "level"] <- level_j;
				taxonomicResults[nrow(taxonomicResults), "name"] <- taxon_k;
				taxonomicResults[nrow(taxonomicResults), "estimate_leaf"] <- estimateLeaf <- exp(coefficients(glm1)[1]);
				taxonomicResults[nrow(taxonomicResults), "estimate_root"] <- estimateRoot <- prod(exp(coefficients(glm1)));
				taxonomicResults[nrow(taxonomicResults), "root_leaf_prop"] <- exp(coefficients(glm1)[2]);
				taxonomicResults[nrow(taxonomicResults), "pvalue"] <- anova(glm0, glm1, test="Chisq")[2, "Pr(>Chi)"];
				taxonomicResults[nrow(taxonomicResults), "otus_in_leaf"] <- length(leaf.indices);
				taxonomicResults[nrow(taxonomicResults), "otus_in_root"] <- length(root.indices);
				taxonomicResults[nrow(taxonomicResults), "reads_leaf"] <- sum(leaf.sums[,"sum"]);
				taxonomicResults[nrow(taxonomicResults), "reads_root"] <- sum(root.sums[,"sum"]);
				taxonomicResults[nrow(taxonomicResults), "order"] <- max(c(estimateLeaf, estimateRoot));
			};
		}

		taxonomicResults <- taxonomicResults[order(taxonomicResults[,"order"], decreasing=T),];
		taxonomicResults[,"habitat"] <- ifelse(taxonomicResults[,"root_leaf_prop"] > 1, "Root", "Leaf"); 
		taxonomicResults[,"order"] <- 1:nrow(taxonomicResults);
		results <- rbind( results, taxonomicResults);
	}

	setwd("/Users/mhorto/Google Drive/Horton_lab/thaliana/projects/root_metagenomics/tables/taxonomic_differences/");
	outputFileName <- paste("taxonomy_vs_organ.", minReads, ".", phylogeneticTarget, ".", minReads, ".tn", otuThreshold, ".min_one_organ", minimumSumEitherOrgan, ".txt", sep="");
	write.table(results, outputFileName, quote=F, sep="\t", row.names=F);
}
