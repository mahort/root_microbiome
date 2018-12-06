## This R-script creates a heatmap of the leaf OR root microbiome
## focusing on a subset (e.g. top-N species) of the bacterial and fungal communities.
## 
## Author: matt.horton
###############################################################################

rm(list=ls());

require(pheatmap); require(RColorBrewer);

source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.microbiome_methods.R"));
source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.base_methods.R"));

###############################################################################
## hard coded variables
###############################################################################
organ <- "root";
otuCutoff <- 97;

bacterialColor <- "tomato1";
fungalColor <- "cadetblue2";
showLabels <- TRUE;

minimumNumberOfReads <- 400;
sizeOfCommunity <- c(100);

analyticalMethod <- "glmer";
absoluteIt <- FALSE;
writePngFile <- TRUE;

minimumRho <- 0; ##
normalizeCentralityMeasure <- TRUE;

{
	###############################################################################\
	## format the taxonomic files.
	###############################################################################
	bacterialTaxonomyFileName <- paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "data/16S/", organ, "/otus", otuCutoff, "/", organ, ".16S.", minimumNumberOfReads, ".taxa.txt");
	bacterialTaxonomy <- read.table(bacterialTaxonomyFileName, header=T, sep="\t", as.is=T, stringsAsFactors=F);

	fungalTaxonomyFileName <- paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "data/ITS/", organ, "/otus", otuCutoff, "/", organ, ".ITS.", minimumNumberOfReads, ".taxa.txt");
	fungalTaxonomy <- read.table(fungalTaxonomyFileName, header=T, sep="\t", as.is=T, stringsAsFactors=F);
	taxa <- rbind(bacterialTaxonomy, fungalTaxonomy);

	###############################################################################\
	## move to the directory
	###############################################################################
	cat("Working with the community at size:", sizeOfCommunity, "\n");
	setwd(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "structure/", organ, "/otus", otuCutoff, "/n_", sizeOfCommunity));

	###############################################################################
	## pick the centrality file, the names of the pvalue and cor file just 
	## differ by a suffix, so we'll be able to automate their selection
	###############################################################################
	cat("Using rho threshold:", minimumRho, "\n");
	suffixMatch <- paste0("_minRho", minimumRho, ifelse( normalizeCentralityMeasure, "_norml", ""), "_centrality.txt$");
	fileName <- chooseFile(c(paste0(organ, ".", minimumNumberOfReads), analyticalMethod, suffixMatch)); ## 

	###############################################################################
	## open the centrality, cor, and pval files
	###############################################################################
	fileName <- fileName$filename;		
	centrality <- read.table(fileName, header=T, sep="\t", row.names=1); ## alt: sort the rhoThresholds, but this forces this to be clean & clear

	## the correlation file
	correlationMatrixFileName <- gsub(suffixMatch, "\\.cor.txt", fileName); ## at present, this doesn't handle non-blupped data (e.g. with the raw suffix).
	cor <- read.table(correlationMatrixFileName, header=T, sep="\t", row.names=1);
	diag(cor) <- NA;

	if( absoluteIt ){
		compositeCommunity.correlationMatrix <- as.matrix(abs(cor));

	} else {
		compositeCommunity.correlationMatrix <- as.matrix(cor);
	}

	###############################################################################
	## combine subsets of the bacterial and fungal communities; make the heatmap
	###############################################################################
	cat("------------------------------\n");
	cat("Creating a heatmap------------\n");
	cat("------------------------------\n");

	getLabelsForHeatmap <- function( speciesMatrix, fungalTaxonomicData, bacterialTaxonomicData, byColumn=TRUE ){

		bacterialTaxonomicData[,"otu_id"] <- paste0("B_", bacterialTaxonomicData[,"otu_id"]);
		fungalTaxonomicData[,"otu_id"] <- paste0("F_", fungalTaxonomicData[,"otu_id"]);
		combinedTaxonomicDataset <- rbind( bacterialTaxonomicData, fungalTaxonomicData );
		rownames(combinedTaxonomicDataset) <- combinedTaxonomicDataset[,"otu_id"];

		if( byColumn ){
			originalLabels <- colnames(speciesMatrix);

		} else {
			originalLabels <- rownames(speciesMatrix);
		}

		dropouts <- which(!originalLabels %in% combinedTaxonomicDataset[,"otu_id"]);
		stopifnot( length(dropouts) == 0 );

		combinedTaxonomicDataset <- combinedTaxonomicDataset[originalLabels,];
		return(apply(combinedTaxonomicDataset[,which(colnames(combinedTaxonomicDataset) == "genus"):which(colnames(combinedTaxonomicDataset) == "kingdom")], 1, getLowestTaxonomicAssignment));
	}

	## write out a heatmap for a quick overview
	groups <- data.frame(id=substr(rownames(compositeCommunity.correlationMatrix), 1, 1));
	rownames(groups) <- rownames(compositeCommunity.correlationMatrix); ## colnames and rownames are the same for the correlation matrix, der.
	group_colors <- list(id =  c("tomato1", "cadetblue2")) ## brewer.pal?
	names(group_colors$id) <- unique(groups[,"id"]);

	rowLabels <- getLabelsForHeatmap( compositeCommunity.correlationMatrix, fungalTaxonomy, bacterialTaxonomy );
	epsOutputFileName <- gsub("\\.txt$", paste0(".hmp", ifelse(showLabels, "", "_nolabs"), ".eps"), correlationMatrixFileName);
	setEPS(paper="special", horizontal=TRUE);
	postscript(epsOutputFileName, width=7, height=7, onefile=FALSE);

	pmap <- pheatmap(compositeCommunity.correlationMatrix, 
			col=brewer.pal(6, "RdBu"), 
			show_rownames=showLabels,
			labels_row=rowLabels,
			fontsize_row=5,
			show_colnames=showLabels,
			labels_col=rowLabels,
			fontsize_col=5,
			annotation_col=groups,
			annotation_row=groups,
			annotation_colors=group_colors,
			clustering_method="ward.D2"
	);

	dev.off();

	if( writePngFile ){
		png( gsub(".eps", ".png", epsOutputFileName), units="in", width=7, height=7, res=1200);
		grid::grid.newpage()
		grid::grid.draw(pmap$gtable)
		dev.off();
	}
}