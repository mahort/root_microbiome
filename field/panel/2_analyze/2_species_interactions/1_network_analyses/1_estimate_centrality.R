## This R-script estimates 4 measures of centrality:
## 1.) degree
## 2.) betweenness
## 3.) eigenvector centrality
## 4.) closeness
##
## and writes these measures to file. It also prepares a graph object for plotting (in a separate function).
## 
## Author: matt.horton
###############################################################################

rm(list=ls());

require(igraph);
source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.microbiome_methods.R"));
source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.base_methods.R"));

###############################################################################
## hard coded variables
###############################################################################
organ <- "root";
otuCutoff <- 97;

bacterialColor <- "tomato1";
fungalColor <- "cadetblue2";

minimumNumberOfReads <- 400;
sizeOfCommunity <- c(100);

minimumRho <- 0; ## earlier, rho was 0.05
normalizeCentralityMeasure <- TRUE;

{	
	## move to the parent of the output diretory
	setwd(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "structure/", organ, "/otus", otuCutoff, "/"));

	###############################################################################\
	## 1.a) begin - format the taxonomic files.
	###############################################################################
	bacterialTaxonomyFileName <- paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "data/16S/", organ, "/otus", otuCutoff, "/", organ, ".16S.", minimumNumberOfReads, ".taxa.txt");
	bacterialTaxonomy <- read.table(bacterialTaxonomyFileName, header=T, sep="\t", as.is=T, stringsAsFactors=F);
	bacterialTaxonomy[,"otu_id"] <- paste0("B_", do.call(rbind, strsplit(bacterialTaxonomy[,"otu_id"], "_"))[,1]);

	fungalTaxonomyFileName <- paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "data/ITS/", organ, "/otus", otuCutoff, "/", organ, ".ITS.", minimumNumberOfReads, ".taxa.txt");
	fungalTaxonomy <- read.table(fungalTaxonomyFileName, header=T, sep="\t", as.is=T, stringsAsFactors=F);
	fungalTaxonomy[,"otu_id"] <- paste0("F_", do.call(rbind, strsplit(fungalTaxonomy[,"otu_id"], "_"))[,1]);
	taxa <- rbind(bacterialTaxonomy, fungalTaxonomy);

	###############################################################################\
	## 1.b) determine the lowest possible assignment, ranging from genus to kingdom, over rows
	###############################################################################
	lowest <- apply(taxa[,which(colnames(taxa) == "genus"):which(colnames(taxa) == "kingdom")], 1, getLowestTaxonomicAssignment);
	taxa <- cbind(taxa, assignment=lowest, stringsAsFactors=F);
	rownames(taxa) <- taxa[,"otu_id"];
	
	## for the choosefile function:
	cat("Choose the file of pvalues.\n");
	for( j in 1:length(sizeOfCommunity)){

		sizeOfCommunity_j <- sizeOfCommunity[j];
		cat("Working with the community at size:", sizeOfCommunity_j, "\n");
		setwd(paste("n_", sizeOfCommunity_j, sep=""));
		fileName <- chooseFile(c("\\bpvals", "\\.txt$")); ## 

		fileName <- fileName$filename;		
		pvals <- read.table(fileName, header=T, sep="\t", row.names=1); ## alt: sort the rhoThresholds, but this forces this to be clean & clear
		diag(pvals) <- 1;

		## the names of the cor and pvals.txt file only differ in the suffix, so we can use the choice of the pvals file to find the corresponding correlation matrix file.
		correlationMatrixFileName <- gsub("\\bpvals", "cor", fileName);
		cor <- read.table(correlationMatrixFileName, header=T, sep="\t", row.names=1);
		pvals <- pvals[rownames(cor), colnames(cor)]; ## these matrices should always be in the same order - but to be sure...

		if( minimumRho > 0 ){
			cat("Using rho threshold:", minimumRho, "\n");			
			pvals[which((cor > (-1*minimumRho)) & (cor < minimumRho), arr.ind=T)] <- 1;
		}

		## we're using an adjacency matrix to create the graph, where 0s mean no relationship (edge) between nodes and higher values mean yes...
		## therefore, we have to flip these p-values into -log10(pvalues)
		minimumPvalue <- min(pvals[pvals != 0], na.rm=T);
		pvals[pvals == 0] <- minimumPvalue/2; ## 

		rowsToUse <- which(apply(pvals, 1, min) != 1);
		colsToUse <- which(apply(pvals, 1, min) != 1);
		if((length(rowsToUse) == 0 || length(colsToUse) == 0 ) || (rowsToUse != colsToUse) || ( rowsToUse == 0 )){ ## check the last arg, doubt it can be triggered
			stop("this isn't working.\n");
		}

		###############################################################################
		## make the graph off of the p-values, using the cor matrix to color the edges for interpretation
		## nb: we DO NOT make this graph directly off of the correlation matrix
		###############################################################################
		pvals <- pvals[rowsToUse, colsToUse];
		minusLog10Pvalues <- -1*log10(pvals); 

		minusLog10Pvalues[minusLog10Pvalues < 2] <- 0;
		minusLog10Pvalues <- as.matrix(minusLog10Pvalues);

		graph <- graph.adjacency(minusLog10Pvalues, mode="undirected", weighted=TRUE, diag=FALSE); ## 
		graph <- simplify(graph);
		autocurve.edges(graph, start=0.4);

		###############################################################################
		## determine the centrality of every node, using four different metrics:
		## 1.) degree
		## 2.) betweenness
		## 3.) eigenvector centrality
		## 4.) closeness
		###############################################################################
		cent.degree <- mstack(degree(graph, normalized=normalizeCentralityMeasure), newHeaders=c("otu_id", "degree"), sorted=FALSE); ## most popular/familiar
		cent.between <- mstack(betweenness(graph, normalized=normalizeCentralityMeasure), newHeaders=c("otu_id", "between"), sorted=FALSE);
		cent.evec <- mstack(evcent(graph)$vector, newHeaders=c("otu_id", "evec"), sorted=FALSE);
		cent.closeness <- mstack(closeness(graph, normalized=normalizeCentralityMeasure), newHeaders=c("otu_id", "closeness"), sorted=FALSE);
		centrality <- merge(merge(merge(cent.degree, cent.between, by="otu_id"), cent.evec, by="otu_id"), cent.closeness, by="otu_id");
		centrality <- merge(centrality, taxa, by="otu_id");
		rownames(centrality) <- centrality[,"otu_id"];
		print(tapply(centrality[,"degree"], centrality[,"kingdom"], mean));

		## get the taxa by OTU id...
		taxa_k <- taxa[V(graph)$name,];
		centrality <- centrality[rownames(taxa_k),];
		stopifnot(all.equal(rownames(taxa_k), as.character(centrality[,"otu_id"])));

		###############################################################################
		## record the centrality measures...
		###############################################################################
		centralityOutputFileName <- gsub(".\\bpvals", paste("_minRho", minimumRho, ifelse( normalizeCentralityMeasure, "_norml", ""), "_centrality", sep=""), fileName );
		write.table( centrality, centralityOutputFileName, quote=F, sep="\t", row.names=F );
		write.table( cor(centrality[,c("degree", "between", "closeness", "evec")]), file=gsub(".\\bpvals", "_centrality_correlation", fileName), quote=F, sep="\t" );

		V(graph)$kingdom <- centrality[,"kingdom"];
		V(graph)$color <- ifelse(V(graph)$kingdom == "Bacteria", bacterialColor, fungalColor);

		V(graph)$degree <- centrality[,"degree"];
		V(graph)$between <- centrality[,"between"];
		V(graph)$evec <- centrality[,"evec"];
		V(graph)$closeness <- centrality[,"closeness"];
		V(graph)$labels <- centrality[,"assignment"];
		print(list.vertex.attributes(graph));
		
		###############################################################################
		## determine if the edges/interactions are positive or negative
		## for coloring, but also for asking whether the measures of centrality are 
		## associated with negative/positive correlations
		## also investigate the individual metrics of centrality to see if they are predictive of heritability in the community
		## or if individual nodes are associated with species richness?!?
		###############################################################################
		edges <- get.edgelist(graph);
		for( k in 1:nrow(edges)){
			E(graph)$cor[k] <- cor[edges[k, 1], edges[k, 2]];
		}

		print(list.edge.attributes(graph));

		## plot(E(graph)$cor, col=E(graph)$color);
		write.graph( graph, file=gsub(".txt", ".gml", centralityOutputFileName), format="gml" );
		setwd("../");
	}
}