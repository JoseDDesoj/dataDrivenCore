require(R6)
ACTIONnet <- R6Class("ACTIONnet", lock_objects=FALSE, lock_class=FALSE,
    public = list(
      sce = NULL,
      network = NULL,
      metacell.network = NULL,
      backbone.network = NULL,
      print = function(...) {
        cat("<ACTIONnet> of ", dim(self$.__enclos_env__$private$items$counts)[1], " genes and ", dim(self$.__enclos_env__$private$items$counts)[2],  " cells \n", sep = "")
      invisible(self)
      },
      ################################################################
      run_ACTIONnet_analysis = function(..., IncolData=NULL, InrowData=NULL, geneFilter = 2, cellFilter = 200, PCAdim=30, Kmin=2, Kmax = 20, NThreads=4, Knn = 30, Ntop=10, Ngenes=NULL, AllAttrib=TRUE, AttribVec="", iters2D = 100, repForcesStrength2D = 5, Iters3Dapp=0, Negweigths3Dapp=1) {
        self$build_sce(IncolData = IncolData, InrowData = InrowData, geneFilter = geneFilter, cellFilter = cellFilter)
        self$normalize_dge()
        self$reduce_dge(PCAdim = PCAdim)
        self$decompose_sce(Kmin = Kmin, Kmax = Kmax, NThreads = NThreads)
        self$build_cellnetwork(Kmax = Kmax, Knn = Knn, NThreads = NThreads)
        self$extract_archetypal_patterns()
        self$compute_gene_metacell_contribution()
        self$refine_metacell_space()
        self$extract_top_contributing_genes(Ntop = Ntop)
        self$extract_metacell_gene_ranked_profiles(Ngenes = Ngenes)
        self$annotate_cell_network(AllAttrib=AllAttrib, AttribVec=AttribVec)
        self$annotate_2D_coordinates_to_network(...,iters = iters2D, repForcesStrength = repForcesStrength2D)
        #self$annotate_3D_coordinates_to_network_aprox(Iters3Dapp = Iters3Dapp, Negweigths3Dapp=Negweigths3Dapp)
        self$build_metacell_network()
        self$build_network_backbone()
      },
      ################################################################
      #initialize = function(...) private$items <- list(...),
      initialize = function(...) private$items <- list(...),
      ################################################################
      build_sce = function(..., IncolData=NULL, InrowData=NULL, geneFilter=NULL, cellFilter=NULL) {
        require(SingleCellExperiment)
          self$sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=as(self$.__enclos_env__$private$items$counts, "sparseMatrix")))
          #if(is.null(IncolData)) SingleCellExperiment::colData(self$sce) <- S4Vectors::DataFrame(data.frame(Xindex=1:ncol(self$sce)))
          if(is.null(IncolData)) colData(self$sce) <- S4Vectors::DataFrame(data.frame(Xindex=1:ncol(self$sce)))
          if(!is.null(IncolData)) colData(self$sce) <- S4Vectors::DataFrame(data.frame(IncolData))
          if(!is.null(InrowData)) colData(self$sce) <- S4Vectors::DataFrame(data.frame(InrowData))
          if(!is.null(geneFilter)) self$sce <- self$sce[Matrix::rowSums(self$sce@assays[["counts"]]>0)>=geneFilter,]
          if(!is.null(cellFilter)) self$sce <- self$sce[,Matrix::colSums(self$sce@assays[["counts"]]>0)>=cellFilter]
          self$.__enclos_env__$private$items$counts <- self$sce@assays[["counts"]]
      },
      ################################################################
      normalize_dge = function(...) {
          X <- self$sce@assays[["counts"]]
          lib_size = Matrix::colSums(X)
          norm.out <- t(t(X)/lib_size * median(lib_size))
          norm.out <- log(norm.out+1)
          self$sce@assays[["logcounts"]] <- as(norm.out, "sparseMatrix")
      },
      ################################################################
      addQCMetrics = function(...) {
        #scater::calculateQCMetrics(self$sce, feature_controls = list(MIT=rownames(sce)[grep("MT-", rownames(sce))][-1]))
        scater::calculateQCMetrics(self$sce, ...)
      },
      ################################################################
      reduce_dge = function(..., PCAdim=30) {
        x <- ACTION::reduce_GeneExpression(self$sce@assays[["logcounts"]], PCA_dim=PCAdim)
        self$.__enclos_env__$private$items$metagene.cell <- x$cell_signatures
        self$.__enclos_env__$private$items$gene.metagene <- x$V
      },
      ################################################################
      decompose_sce = function(..., Kmin=2, Kmax = 20, NThreads=4) {
        x <- ACTION::run_ACTION(self$.__enclos_env__$private$items$metagene.cell, k_min=Kmin, k_max = Kmax, numThreads = NThreads)
        self$.__enclos_env__$private$items$ACTION.output <- x
        self$.__enclos_env__$private$items$metacell.cell.trace <- rlist::list.remove(x$CellProfiles, 1)
        self$.__enclos_env__$private$items$metagene.metacell.trace <- lapply(rlist::list.remove(x$SmoothingMats, 1), function(i) self$.__enclos_env__$private$items$metagene.cell%*%i)
      },
      ################################################################
      build_cellnetwork = function(..., Kmax = 20, Knn = 30, NThreads=4) {
        G <- ACTION::build_ACTIONnet(self$.__enclos_env__$private$items$ACTION.output$CellProfiles[1:Kmax], Knn, NThreads)
        self$network <- igraph::graph.adjacency(G, mode = "undirected", weighted = T)
      },
      ################################################################
      extract_archetypal_patterns = function(...) {
        names(self$.__enclos_env__$private$items$metacell.cell.trace) <-  paste0("D", 1:length(self$.__enclos_env__$private$items$metacell.cell.trace))
        for(i in 1:length(self$.__enclos_env__$private$items$metacell.cell.trace)) rownames(self$.__enclos_env__$private$items$metacell.cell.trace[[i]]) <- paste0("D", i, "_", 1:nrow(self$.__enclos_env__$private$items$metacell.cell.trace[[i]]))
        self$.__enclos_env__$private$items$metacell.cell.assignments.trace <- lapply(self$.__enclos_env__$private$items$metacell.cell.trace, function(i) apply(i, 2, which.max))
        self$.__enclos_env__$private$items$metacell.cell.assignments.trace <- lapply(self$.__enclos_env__$private$items$metacell.cell.trace, function(i) apply(i, 2, which.max))
        self$.__enclos_env__$private$items$multilevel.metacell.space <- do.call(rbind, self$.__enclos_env__$private$items$metacell.cell.trace)
        self$.__enclos_env__$private$items$landmarks.idx <- apply(self$.__enclos_env__$private$items$multilevel.metacell.space, 1, which.max)
        self$.__enclos_env__$private$items$landmarks <- as.data.frame(self$.__enclos_env__$private$items$landmarks.idx)
        colnames(self$.__enclos_env__$private$items$landmarks) <- "cell.index"
        self$.__enclos_env__$private$items$cell.multilevel.metacell.confidence <- apply(self$.__enclos_env__$private$items$multilevel.metacell.space, 2, max)
      },
      ################################################################
      compute_gene_metacell_contribution = function(...) {
        out <- self$.__enclos_env__$private$items$gene.metagene%*%self$.__enclos_env__$private$items$metagene.cell%*%t(self$.__enclos_env__$private$items$multilevel.metacell.space)
        rownames(out) <- rownames(self$sce)
        self$.__enclos_env__$private$items$multilevel.gene.metacell.contribution <- as.data.frame(out)
      },
      ################################################################
      refine_metacell_space = function(...) {
          temp <- ACTION::refine_solution(self$.__enclos_env__$private$items$metagene.cell, self$.__enclos_env__$private$items$metagene.cell%*%do.call(cbind, args = self$.__enclos_env__$private$items$ACTION.output$SmoothingMats))
          self$.__enclos_env__$private$items$refined.decomposition <- temp
          self$.__enclos_env__$private$items$multilevel.metacell.space.refined <- temp$Refined_H
      },
      ################################################################
      extract_top_contributing_genes = function(..., Ntop=10) {
        out <- lapply(1:ncol(self$.__enclos_env__$private$items$multilevel.gene.metacell.contribution), function(i)  X = rownames(self$.__enclos_env__$private$items$multilevel.gene.metacell.contribution)[order(self$.__enclos_env__$private$items$multilevel.gene.metacell.contribution[,i], decreasing = T)][1:Ntop])
        names(out) <- colnames(self$.__enclos_env__$private$items$multilevel.gene.metacell.contribution)
        self$.__enclos_env__$private$items$metacell.top.contributing.genes <- out
      },
      ################################################################
      extract_metacell_gene_ranked_profiles = function(..., Ngenes=NULL) {
         metacell.ranked.genes.list <- lapply(names(self$.__enclos_env__$private$items$multilevel.gene.metacell.contribution), function(i) rownames(self$.__enclos_env__$private$items$multilevel.gene.metacell.contribution)[order(self$.__enclos_env__$private$items$multilevel.gene.metacell.contribution[,i], decreasing = T)])
        names(metacell.ranked.genes.list) <- names(self$.__enclos_env__$private$items$multilevel.gene.metacell.contribution)
        if(!is.null(Ngenes)) metacell.ranked.genes.list <- lapply(metacell.ranked.genes.list, function(i) i[1:Ngenes])
        self$.__enclos_env__$private$items$metacell.ranked.genes.list <- metacell.ranked.genes.list
      },
      ################################################################
      annotate_cell_network = function(..., AllAttrib=TRUE, AttribVec="") {
        tempNet <- self$network
        if(AllAttrib) AttribVec <- colnames(self$sce@colData)
        Attrs <- DataFrame(self$sce@colData[,AttribVec])
        colnames(Attrs) <- colnames(self$sce@colData)
        for(i in AttribVec) if(is.factor(Attrs[,i])) Attrs[,i] <- as.character(Attrs[,i])
        landmark <- rep("other", igraph::vcount(tempNet))
        landmark[self$.__enclos_env__$private$items$landmarks.idx] <- "landmark"
        igraph::V(tempNet)$landmark <- landmark
        for(i in AttribVec) tempNet <-  igraph::set_vertex_attr(tempNet, name = i, value=as.character(self$sce@colData[,i]))
        self$network <- tempNet
      },
      ################################################################
      annotate_2D_coordinates_to_network = function(..., iters = 100, repForcesStrength = 5) {
        tempNet <- self$network
        coor <-  ACTION::layout_ACTIONet(igraph::get.adjacency(tempNet, attr = "weight"), iters, repForcesStrength = repForcesStrength,  ...)
        igraph::V(tempNet)$x <- coor[,1]
        igraph::V(tempNet)$y <- coor[,2]
        self$network <- tempNet
      },
      ################################################################
      annotate_3D_coordinates_to_network = function(..., iters = 100) {
        tempNet <- self$network
        coor <- igraph::layout_with_fr(tempNet, dim = 3, niter = iters)
        igraph::V(tempNet)$X <- coor[,1]
        igraph::V(tempNet)$Y <- coor[,2]
        igraph::V(tempNet)$Z <- coor[,3]
        self$network <- tempNet
      },
      ################################################################
      annotate_3D_coordinates_to_network_aprox = function(..., Iters3Dapp = 0, Negweigths3Dapp=1) {
        tempNet <- self$network
        coor <- ACTION::layout_ACTIONet_3D(igraph::get.adjacency(tempNet, attr = "weight"), iters = Iters3Dapp, negate_weights = Negweigths3Dapp)
        igraph::V(tempNet)$X <- coor[,1]
        igraph::V(tempNet)$Y <- coor[,2]
        igraph::V(tempNet)$Z <- coor[,3]
        self$network <- tempNet
      },
      ################################################################
      build_metacell_network = function(...) {
        CC <- cor(t(self$.__enclos_env__$private$items$multilevel.metacell.space))
        diag(CC) <- 0
        CC[CC<0] <- 0
        #CC <- glasso::glasso(CC, rho=0.2)
        self$metacell.network <- igraph::graph.adjacency(CC, weighted = T, mode = "undirected")
      },
      ################################################################
      build_network_backbone = function(...) {
        self$backbone.network <- igraph::graph.adjacency(ACTION::construct_backbone(self$.__enclos_env__$private$items$multilevel.metacell.space), weighted = T, mode = "undirected")
      },
      ################################################################
      plot_ACTIONnet = function(..., AttrName=NULL, CPall="npg", Vsize=3, PlotLegend=TRUE, Gene=NULL, ColGroups=FALSE, Basecolor="tomato", LegendPos="bottomleft", TransScore=NULL) {
        #require(ggpubr)
        InNet <- self$network
        Insce <- self$sce
        if(!is.null(AttrName)) ColGroups <- TRUE
        if(!ColGroups) PlotLegend <- FALSE
        #if(!is.null(Gene)) if(is.null(Insce)) return(print("sce object missing"))
        if(ColGroups) Ncols <- factor(ggpubr::get_palette(CPall, length(unique(igraph::get.vertex.attribute(InNet, AttrName)))))[factor(igraph::get.vertex.attribute(InNet, AttrName))] else Ncols <- rep(Basecolor, igraph::vcount(InNet))
        if(!is.null(TransScore)) Ncols <- ggplot2::alpha(Ncols, TransScore)
        if(is.null(Gene)) {
        plot(InNet, vertex.color=as.character(Ncols), edge.width=scale(igraph::E(InNet)$weight), vertex.label="", vertex.size=Vsize, layout=cbind(as.numeric(igraph::get.vertex.attribute(InNet, "x")),as.numeric(igraph::get.vertex.attribute(InNet, "y"))), edge.color="grey", ...)
        } else {
        temp.scores <- igraph::page_rank(InNet, directed = F, damping = 0.85, personalized = as.numeric(Insce@assays[["logcounts"]][Gene,]), weights = igraph::E(InNet)$weight)$vector
        temp.scores <- temp.scores/max(temp.scores)
        plot(InNet, vertex.color=as.character(Ncols), edge.width=scale(igraph::E(InNet)$weight), vertex.label="", vertex.size=as.numeric(temp.scores)*10, layout=cbind(as.numeric(igraph::get.vertex.attribute(InNet, "x")),as.numeric(igraph::get.vertex.attribute(InNet, "y"))), edge.color="grey", ...)
        }
        if(PlotLegend==TRUE) legend(LegendPos, legend=unique(igraph::get.vertex.attribute(InNet, AttrName)), col = as.character(unique(factor(ggpubr::get_palette(CPall, length(unique(igraph::get.vertex.attribute(InNet, AttrName)))))[factor(igraph::get.vertex.attribute(InNet, AttrName))])), bty = "n", pch=20 , pt.cex = 3, cex = 1, horiz = FALSE, inset = c(0.1, 0.1))
      },
      ################################################################
      plot_ACTIONnet_3D_points = function(..., AttrName=NULL, CPall="npg", Vsize=0.5, Gene=NULL, ColGroups=FALSE, Basecolor="tomato", PropSizeWeight=0.5, Bg="black") {
        InNet <- self$network
        Insce <- self$sce
        if(!is.null(AttrName)) ColGroups <- TRUE
        if(!is.null(Gene)) {
          temp.scores <- igraph::page_rank(InNet, directed = F, damping = 0.85, personalized = as.numeric(Insce@assays[["logcounts"]][Gene,]), weights = igraph::E(InNet)$weight)$vector
          temp.scores <- temp.scores/max(temp.scores)
          Vsize <- as.numeric(temp.scores)*PropSizeWeight
        }
        if(ColGroups) {
          if(ColGroups) Ncols <- factor(ggpubr::get_palette(CPall, length(unique(igraph::get.vertex.attribute(InNet, AttrName)))))[factor(igraph::get.vertex.attribute(InNet, AttrName))] else Ncols <- rep(Basecolor, igraph::vcount(InNet))
            threejs::scatterplot3js(x = igraph::V(InNet)$X, y = igraph::V(InNet)$Y, z = igraph::V(InNet)$Z, axis.scales = FALSE, size = Vsize, axis = F, grid = F, color = Ncols, bg=Bg, ...)
          } else {
            threejs::scatterplot3js(x = igraph::V(InNet)$X, y = igraph::V(InNet)$Y, z = igraph::V(InNet)$Z, axis.scales = FALSE, size = Vsize, axis = F, grid = F, color = Basecolor, bg=Bg, ...)
          }
      },
      ################################################################
      plot_ACTIONnet_3D_net = function(..., AttrName=NULL, CPall="npg", Vsize=0.05, Gene=NULL, ColGroups=FALSE, Basecolor="tomato", PropSizeWeight=0.5, Sparsify=FALSE) {
        require(ggpubr)
        require(threejs)
        InNet <- self$network
        Insce <- self$sce
        if(Sparsify) InNet <- self$.__enclos_env__$private$axFunctions$sparsify_net_edges(InNet)
        if(!is.null(AttrName)) ColGroups <- TRUE
        if(!is.null(Gene)) {
          temp.scores <- igraph::page_rank(InNet, directed = F, damping = 0.85, personalized = as.numeric(Insce@assays[["logcounts"]][Gene,]), weights = igraph::E(InNet)$weight)$vector
          temp.scores <- temp.scores/max(temp.scores)
          Vsize <- as.numeric(temp.scores)*PropSizeWeight
        }
        if(ColGroups) {
            if(ColGroups) Ncols <- factor(ggpubr::get_palette(CPall, length(unique(igraph::get.vertex.attribute(InNet, AttrName)))))[factor(igraph::get.vertex.attribute(InNet, AttrName))] else Ncols <- rep(Basecolor, igraph::vcount(InNet))
            threejs::graphjs(InNet, layout = cbind(igraph::V(InNet)$X, igraph::V(InNet)$Y, igraph::V(InNet)$Z), vertex.size = Vsize, vertex.color = Ncols, edge.alpha = 0.5, ...)
        } else {
            threejs::graphjs(InNet, layout = cbind(igraph::V(InNet)$X, igraph::V(InNet)$Y, igraph::V(InNet)$Z), vertex.size = Vsize, vertex.color = Basecolor, edge.alpha = 0.5, ...)
        }
      },
      ################################################################
      pick_landmark_cell = function(..., cex=1) {
        #X11()
        InNet <- self$network
        quartz()
        Cex <- rep(cex, igraph::vcount(InNet))
        Cex[igraph::V(InNet)$landmark=="landmark"] <- cex*2
        Alphas <- ifelse(igraph::V(InNet)$landmark=="landmark", 1, 0.25)
        Cols <- ifelse(igraph::V(InNet)$landmark=="landmark", "red", 1)
        plot(igraph::V(InNet)$x, igraph::V(InNet)$y, col=ggplot2::alpha(Cols, Alphas), pch=20, cex=Cex, axes = FALSE, xlab="", ylab="", ...)
        temp <- identify(igraph::V(InNet)$x, igraph::V(InNet)$y)
        return(temp)
      },
      ################################################################
      plot_cell_score_across_network_2D = function(..., InScore, Cpall="Blues") {
        InNetwork <- self$network
        score.color <- scales::col_numeric(Cpall, domain = NULL)(as.numeric(InScore))
        plot(InNetwork, vertex.color=score.color, edge.width=scale(igraph::E(InNetwork)$weight), vertex.label="", layout=cbind(igraph::V(InNetwork)$x, igraph::V(InNetwork)$y), edge.color=ggplot2::alpha("grey", 0.5), vertex.size=InScore/max(InScore)*10, vertex.frame.color=ggplot2::alpha("black", 0.8), ...)
      },
      ################################################################
      compute_gene_propagation = function(..., InGenes) {
        Insce <- self$sce
        InNetwork <- self$network
        InGenes <- InGenes[InGenes%in%rownames(Insce)]
        sce.genes <- self$.__enclos_env__$private$axFunctions$split_sce_rows(Insce, InRowNames=InGenes)
        #sce.genes <- sce.genes[rowSums(exprs(sce.genes))>0,]
        Genes.propagation <- do.call(cbind, lapply(rownames(sce.genes), function(i) igraph::page_rank(InNetwork, directed = F, damping = 0.85, personalized = sce.genes[i,]@assays[["logcounts"]], weights = igraph::E(InNetwork)$weight)$vector))
        colnames(Genes.propagation) <- rownames(sce.genes)
        Genes.propagation <- apply(Genes.propagation, 2, function(i) i/max(i))
        return(Genes.propagation)
      },
      ################################################################
      annotate_cells_based_on_geneSets = function(..., InMarkerList, AttrName="annot.action") {
        InNet <- self$network
        Insce <- self$sce
        if(AttrName%in%names(igraph::vertex_attr(InNet))) print("WARNING: AttrName already there...")
        InMarkerList <- lapply(InMarkerList, function(i) i[i%in%rownames(Insce)])
        #out <- Marker.based.cell.annotation(Insce, InNet, InMarkerList)
        #MarkerScoreMatrix <- self$compute_gene_propagation(InGenes = unlist(InMarkerList))
        #GroupScoreMatrix <- do.call(cbind, lapply(names(InMarkerList), function(i) apply(MarkerScoreMatrix[,InMarkerList[[i]]], 1, self$.__enclos_env__$private$axFunctions$geom_mean)))

        temp.scores.list <- lapply(InMarkerList, function(i) self$compute_gene_propagation(InGenes = i))
        GroupScoreMatrix <- do.call(cbind, lapply(temp.scores.list, function(i) apply(as.matrix(i), 1, self$.__enclos_env__$private$axFunctions$geom_mean)))
        colnames(GroupScoreMatrix) <- names(InMarkerList)
        Label <- colnames(GroupScoreMatrix)[apply(GroupScoreMatrix, 1, which.max)]
        Marker.confidence <- apply(GroupScoreMatrix, 1, max)
        InNet <- igraph::set_vertex_attr(graph = InNet, name = AttrName, value = Label)
        InNet <- igraph::set_vertex_attr(graph = InNet, name = paste0(AttrName, "_conf"), value = Marker.confidence)
        self$network <- InNet
      }
      ################################################################
    ),
  private = list(
    items = list(),
    axFunctions = list(
        sparsify_net_edges = function(InNet) {
          return(igraph::delete_edges(InNet, E(InNet)[scale(E(InNet)$weight)<0]))
        },
        ###
        geom_mean = function (x, na.rm = TRUE) {
           if (is.null(nrow(x))) {
            exp(mean(log(x), na.rm = TRUE))
          }
          else {
            exp(apply(log(x), 2, mean, na.rm = na.rm))
          }
        },
        ###
        split_sce_cols = function(InSCE, InClass) {
          temp.sce.split <- lapply(as.character(sort(unique(colData(InSCE)[,InClass]))), function(i) InSCE[,colData(InSCE)[,InClass]==i])
          names(temp.sce.split) <- as.character(sort(unique(colData(InSCE)[,InClass])))
          return(temp.sce.split)
        },
        ###
        split_sce_rows = function(InSCE, InRowNames) return(InSCE[match(InRowNames, rownames(InSCE)),]),
        ###
        compute_cell_group_propagation = function(InNetwork, ConditionName, WhichCondition) {
          Score.vec <- igraph::page_rank(InNetwork, directed = F, damping = 0.85, personalized = ifelse(igraph::get.vertex.attribute(InNetwork, ConditionName)==WhichCondition, 1, 0), weights = igraph::E(InNetwork)$weight)$vector
          Score.vec <- Score.vec/max(Score.vec)
          return(Score.vec)
        },
        ###
        propagate_cell_score = function(InNetwork, InScore) {
          Score.vec <- igraph::page_rank(InNetwork, directed = F, damping = 0.85, personalized = InScore, weights = igraph::E(InNetwork)$weight)$vector
          Score.vec <- Score.vec/max(Score.vec)
          return(Score.vec)
        }
        ###
    )
  )
)

Plot.filled.after.propagation <- function(Insce, InNet, Ingene, Cex=0.25, detectionCut=1) {
    temp <- split(as.numeric(Propagate.gene(Insce, InNet, InGenes = Ingene)),as.numeric(Insce@assays[["counts"]][Ingene,]))
    Insce@colData$x <- V(InNet)$x
    Insce@colData$y <- V(InNet)$y
    Insce@colData$temp.annot <- as.character(Insce@assays[["counts"]][Ingene,])
    sce.splitted <- Split.sce.cols(Insce, "temp.annot")
    sce.splitted$`0`@colData$temp.annot <- Binarize(temp$`0`)

    par(mfrow=c(2,2), mar=c(0,0,2,0))
      plot(sce.splitted$`0`@colData$x, sce.splitted$`0`@colData$y, pch=20, cex=Cex, axes=F, xlab="", ylab="", main=paste(Ingene, "undetected"))

      plot(sce.splitted$`0`@colData$x, sce.splitted$`0`@colData$y, pch=20, cex=Cex, col=ifelse(sce.splitted$`0`@colData$temp.annot==1,"red", "black"), axes=F, xlab="", ylab="", main=paste(Ingene, "undetected (after propagation - red)"))

      plot(Insce@colData$x[as.numeric(Insce@colData$temp.annot)>detectionCut], Insce@colData$y[as.numeric(Insce@colData$temp.annot)>detectionCut], pch=20, cex=Cex, axes=F, xlab="", ylab="", main=paste(Ingene, "detected"))

      plot(sce.splitted[[length(sce.splitted)]]@colData$x, sce.splitted[[length(sce.splitted)]]@colData$y, pch=20, cex=Cex, axes=F, xlab="", ylab="", main=paste(Ingene, "high counts"))
    par(mfrow=c(1,1))
}

Plot.filled.after.propagation.boxplot <- function(Insce, InNet, Ingene, Cex=0.2) {
  coords <- Extract.2D.layout(InNet)
   Insce@colData$x <- coords[,1]
   Insce@colData$y <- coords[,2]
  x.counts <- Insce@assays[["counts"]][Ingene,]
  PropScore <- Propagate.gene(Insce, InNet, Ingene)
  PropScore.splitted <- split(PropScore,x.counts)
  Insce@colData$temp.annot <- as.character(x.counts)
  sce.splitted <- Split.sce.cols(Insce, "temp.annot")
  sce.splitted$`0`@colData$temp.annot <- Binarize(PropScore.splitted$`0`)

  #pdf("temp.pdf")
  #layout(matrix(c(1,1,2,2)), heights = c(5,1))
  layout(matrix(c(1,2,1,2,3,3), 2, 3), heights = c(5,1), widths = c(1,1,3))
  par(mar=c(4,4,1,0))
  boxplot(PropScore.splitted, axes=F, ylab="propagation score", main=Ingene, pars=list(outcol=c("red", rep("grey", length(PropScore.splitted)-1)), outpch=20, outcex=c(1.5, rep(1, length(PropScore.splitted)-1))))
  #box()
  #axis(side = 1, labels = NA, tck=0)
  axis(side = 2)
  barplot(sapply(PropScore.splitted, length), ylab="cell count", xlab="read counts", col="black")
  plot(sce.splitted$`0`@colData$x, sce.splitted$`0`@colData$y, pch=20, cex=Cex, col=ifelse(sce.splitted$`0`@colData$temp.annot==1,"red", "black"), axes=F, xlab="", ylab="", main="propagated scores")
  #dev.off()
  #Openfile("temp.pdf")
}
