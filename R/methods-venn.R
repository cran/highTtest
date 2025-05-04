if(!isGeneric("vennD")){
  setGeneric(name="vennD", 
             def=function(x, gamma, ...){standardGeneric("vennD")})
}

setMethod("vennD",
  signature = c(x="highTtest", gamma="numeric"),
  definition = function(x, gamma, ...){

    tst <- x@gammas - gamma
    igamma <- which(tst > -1e-8 & tst < 1e-8)

    if(length(igamma) == 0) {
      cat("Requested gamma value not included in object provided.\n")
      return(NULL)
    }

    if(is.null(x@BH) && is.null(x@ST)){
      cat("A Venn Diagram method is not available for 1 set.\n")
      return(NULL)
    } else if(!is.null(x@BH) && !is.null(x@ST)){
      venn3(x, igamma, ...)
    } else {
      venn2(x, igamma, ...)
    }
  }
)

venn2_orig <- function(x, igamma, ...){

  if(!is.null(x@BH)){
    area <- cbind(x@CK[,igamma], x@BH[,igamma])
    category <- c("CK","BH")
  } else if(!is.null(x@ST)){
    area <- cbind(x@CK[,igamma], x@ST[,igamma])
    category <- c("CK","ST")
  }

  same <- sum(rowSums(area)==2)
  area <- colSums(area)

  cgy <- c(paste(category[1], " (", area[1], ")", sep=""),
           paste(category[2], " (", area[2], ")", sep=""))

  area[1] <- area[1] - same
  area[2] <- area[2] - same

  vec <- c(area, same)
  names(vec) <- c("10","01","11")

  args <- list(...)

  if(is.null(args$Title)){
    args$Title <- paste(paste(category,collapse=" & "),
                        "at level", x@gammas[igamma], sep=" ")
  }

  if(is.null(args$Colors)){
    args$Colors <- c("red","yellow")
  }

  args$labels <- cgy
  args$reverseLabelOrdering <- FALSE
  args$x <- vec

  plot.new()

#  do.call(colorfulVennPlot::plotVenn2d,args)
}

venn2 <- function(x, igamma, ...) {
  if (!is.null(x@BH)) {
    area_mat <- cbind(x@CK[, igamma], x@BH[, igamma])
    category <- c("CK", "BH")
  } else if (!is.null(x@ST)) {
    area_mat <- cbind(x@CK[, igamma], x@ST[, igamma])
    category <- c("CK", "ST")
  }
  
  A_only <- sum(area_mat[, 1] == 1 & area_mat[, 2] == 0)
  B_only <- sum(area_mat[, 2] == 1 & area_mat[, 1] == 0)
  both   <- sum(rowSums(area_mat) == 2)
  
  A_total <- A_only + both
  B_total <- B_only + both
  
  cat_labels <- paste(category, "(", c(A_total, B_total), ")", sep = "")
  
  args <- list(...)
  if (is.null(args$Title)) {
    args$Title <- paste(paste(category, collapse = " & "),
                        "at level", x@gammas[igamma], sep = " ")
  }
  if (is.null(args$Colors)) {
    args$Colors <- c("red", "yellow")
  }
  
  # Final sanity check
  if (both > min(A_total, B_total)) {
    warning("inconsistent Venn counts; skipping plot.")
    return(invisible(NULL))
  }
  
  grid.newpage()
  venn_plot <- VennDiagram::draw.pairwise.venn(
    area1 = A_total,
    area2 = B_total,
    cross.area = both,
    category = cat_labels,
    fill = args$Colors,
    scaled = TRUE,
    lty = "blank",
    cex = 1.2,
    cat.cex = 1.2,
    cat.pos = c(-20, 20),
    cat.dist = c(0.05, 0.05),
    main = args$Title
  )
  
  grid.draw(venn_plot)
}


venn3_orig <- function(x, igamma, ...){

  area <- cbind(x@CK[,igamma], x@BH[,igamma], x@ST[,igamma])

  same12 <- sum(rowSums(area[,c(1,2)])==2 & area[,3]==0)
  same13 <- sum(rowSums(area[,c(1,3)])==2 & area[,2]==0)
  same23 <- sum(rowSums(area[,c(2,3)])==2 & area[,1]==0)

  same123 <- sum(rowSums(area) == 3)

  area <- colSums(area)

  category <- c(paste("CK(", area[1], ")", sep=""),
                paste("BH(", area[2], ")", sep=""),
                paste("ST(", area[3], ")", sep=""))

  area[1] <- area[1] - same12 - same13 - same123
  area[2] <- area[2] - same12 - same23 - same123
  area[3] <- area[3] - same13 - same23 - same123

  vec <- c(area, same12, same23, same13, same123)
  names(vec) <- c("100","010","001","110","011","101","111")

  args <- list(...)
  if(is.null(args$Title)){
    args$Title <- paste("CK, BH, & ST at level", x@gammas[igamma], sep=" ")
  }
  if(is.null(args$Colors)){
    args$Colors <- c("red","yellow","orange")
  }

  args$labels <- category
  args$x <- vec

  plot.new()

#  do.call(colorfulVennPlot::plotVenn3d,args)


}

venn3 <- function(x, igamma, ...) {
  
  mat <- cbind(x@CK[, igamma], x@BH[, igamma], x@ST[, igamma])
  colnames(mat) <- c("CK", "BH", "ST")
  
  # Each region
  only1   <- sum(mat[,1] == 1 & mat[,2] == 0 & mat[,3] == 0)  # 100
  only2   <- sum(mat[,1] == 0 & mat[,2] == 1 & mat[,3] == 0)  # 010
  only3   <- sum(mat[,1] == 0 & mat[,2] == 0 & mat[,3] == 1)  # 001
  both12  <- sum(mat[,1] == 1 & mat[,2] == 1 & mat[,3] == 0)  # 110
  both23  <- sum(mat[,1] == 0 & mat[,2] == 1 & mat[,3] == 1)  # 011
  both13  <- sum(mat[,1] == 1 & mat[,2] == 0 & mat[,3] == 1)  # 101
  all123  <- sum(rowSums(mat) == 3)                           # 111
  
  # Total for each set
  A_total <- only1 + both12 + both13 + all123
  B_total <- only2 + both12 + both23 + all123
  C_total <- only3 + both13 + both23 + all123
  
  category <- c("CK", "BH", "ST")
  labels <- paste(category, "(", c(A_total, B_total, C_total), ")", sep = "")
  
  args <- list(...)
  if (is.null(args$Title)) {
    args$Title <- paste("CK, BH, & ST at level", x@gammas[igamma], sep = " ")
  }
  if (is.null(args$Colors)) {
    args$Colors <- c("red", "yellow", "orange")
  }
  
  grid.newpage()
  venn_plot <- VennDiagram::draw.triple.venn(
    area1 = A_total,
    area2 = B_total,
    area3 = C_total,
    n12 = both12 + all123,
    n23 = both23 + all123,
    n13 = both13 + all123,
    n123 = all123,
    category = labels,
    fill = args$Colors,
    scaled = TRUE,
    lty = "blank",
    cex = 1.2,
    cat.cex = 1.2,
    cat.pos = c(-20, 20, 0),
    cat.dist = c(0.05, 0.05, 0.05),
    main = args$Title
  )
  
  grid.draw(venn_plot)
}


