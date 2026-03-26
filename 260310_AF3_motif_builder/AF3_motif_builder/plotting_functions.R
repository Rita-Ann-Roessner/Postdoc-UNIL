### Modified MixTCRviz plotting functions ###
TCRgene2aes <- MixTCRviz:::TCRgene2aes

# Compute the weighted average of motifs for different V/J usage
# This is important since the choice of V/J has huge influence on the CDR3 motif.
# By considering V/J usage in baseline, it helps identifying what is specific to the input TCR
#' @export
weighted_countCDR3 <- function(countCDR3.VJL.baseline, countVJ.L.es){
  
  countCDR3 <- list()
  L <- names(countVJ.L.es)
  for(lg.c in L){
    
    lg <- as.numeric(unlist(strsplit(lg.c, split="_"))[2])
    
    countCDR3[[lg.c]] <- matrix(0, nrow=N.aa, ncol=lg)
    rownames(countCDR3[[lg.c]]) <- aa.list
    if(sum(countVJ.L.es[[lg.c]]) > 0){
      w <- countVJ.L.es[[lg.c]] / sum(countVJ.L.es[[lg.c]])
      for(v in rownames(countVJ.L.es[[lg.c]])){
        for(j in colnames(countVJ.L.es[[lg.c]])){
          if(w[v,j] > 0){
            s <- paste(v, j, sep="_")
            if(lg.c %in% names(countCDR3.VJL.baseline[[s]])){
              if(length(countCDR3.VJL.baseline[[s]][[lg.c]]) > 0){
                countCDR3[[lg.c]] <- countCDR3[[lg.c]] +
                  scale(
                    countCDR3.VJL.baseline[[s]][[lg.c]],
                    center = FALSE,
                    scale = colSums(countCDR3.VJL.baseline[[s]][[lg.c]])
                  ) * w[v,j]
              }
            }
          }
        }
      }
    }
    if(colSums(countCDR3[[lg.c]])[1] != 0){
      countCDR3[[lg.c]] <- scale(
        countCDR3[[lg.c]],
        center = FALSE,
        scale = colSums(countCDR3[[lg.c]])
      )
    } else {
      countCDR3[[lg.c]][1:N.aa, 1:lg] <- 0.05
    }
  }
  return(countCDR3)
}

# Compute the weighted average of CDR3 length for the observed V/J usage
# This is important since CDR3 length is primarily determined by the length of the V and J segments
#' @export
weighted_countL <- function(countL.VJ.baseline, countVJ.es){
  
  countL <- rep(0, times = Lmax - Lmin + 1)
  names(countL) <- paste("L", Lmin:Lmax, sep = "_")
  
  countVJ.es <- countVJ.es / sum(countVJ.es)
  for(v in rownames(countVJ.es)){
    for(j in colnames(countVJ.es)){
      if(countVJ.es[v, j] > 0){
        s <- paste(v, j, sep = "_")
        if(s %in% names(countL.VJ.baseline)){
          if(length(countL.VJ.baseline[[s]]) > 0){
            countL.VJ.baseline[[s]] <-
              countL.VJ.baseline[[s]] / sum(countL.VJ.baseline[[s]]) * countVJ.es[v, j]
            for(lc in names(countL.VJ.baseline[[s]])){
              countL[lc] <- countL[lc] + countL.VJ.baseline[[s]][lc]
            }
          }
        }
      }
    }
  }
  return(countL)
}

# Infer the short name for the MHC from the long name (HLA-A*02:01 -> A0201)
find_mhc <- function(m){
  mhc <- unlist(lapply(1:dim(m)[1], function(x){
    if(m[x, "MHC"] == "HLA-A:01") {
      h <- "A0101"
    } else {
      if(m[x, "species"] == "HomoSapiens"){
        a <- unlist(strsplit(m[x, "MHC"], split = "-", fixed = TRUE))
        h <- gsub("*", "", a[2], fixed = TRUE)
        h <- gsub(":", "", h, fixed = TRUE)
        h <- gsub("/", "_", h, fixed = TRUE)
      } else {
        h <- m[x, "MHC"]
      }
    }
    return(h)
  }))
}

#' Plot the comparison between V or J usage.
#'
#' @param pType Switch telling the type of plot to use. See plot.VJ.switch
#' @param species Tells the species, which is used for gene colors in the plots.
#' @param ret.resList when TRUE indicates to return a list of the results instead of the plots directly
#' @param combined.resList When this isn't NULL, use these combined results instead of building from scratch
#' @param label.neg If TRUE, also show labels of genes depleted in input1
#' @param label.diag Print labels on the diagonal above a certain value
#' @param label.min.fr Region of the lower-left corner with no labels
#' @export
plotVJ <- function(count.es, count.rep, plddt.es = NULL, sd.es = NULL, sd.rep = NULL,
                   distr.es = NULL, distr.rep = NULL, info = NULL, comp.baseline = TRUE,
                   pType = 1, species = "HomoSapiens", ret.resList = FALSE,
                   combined.resList = NULL, label.neg = FALSE,
                   ZscoreVJ.thresh = 0, FoldChangeVJ.thresh = 1.25,
                   label.diag = 0.3, label.min.fr = c(0.05, 0.05),
                   print.size = TRUE, plot.sd = TRUE, verbose = 1, show.plddt.legend = FALSE){
  
  if (is.null(combined.resList)) {
    if (length(count.es) == 0) {
      if (ret.resList) {
        return(NULL)
      } else {
        return(ggplot())
      }
    }
    
    segment <- info["gene"]
 
    n.es <- sum(count.es)
    n.rep <- sum(count.rep)
    
    nm <- unique(c(names(count.rep), names(count.es)))
    cn <- c(info["input1.name"], info["baseline.name"])
    
    count <- matrix(0, nrow = length(nm), ncol = 2)
    rownames(count) <- nm
    colnames(count) <- cn
    
    for(i in seq_along(count.es)){
      count[names(count.es[i]), info["input1.name"]] <- count.es[i]
    }
    for(i in seq_along(count.rep)){
      count[names(count.rep[i]), info["baseline.name"]] <- count.rep[i]
    }
    
    count <- scale(count, center = FALSE, scale = colSums(count))
    count.df <- data.frame(count)
    colnames(count.df) <- c("Y", "X")
    count.df$name <- rownames(count.df)
    count.df$gene <- gsub("\\*.*$", "", count.df$name)
    
    count.df$plddt <- NA_real_
    if (!is.null(plddt.es)) {
      idx_name <- match(count.df$name, names(plddt.es))
      matched_name <- !is.na(idx_name)
      if (any(matched_name)) {
        count.df$plddt[matched_name] <- as.numeric(plddt.es[idx_name[matched_name]])
      }
      
      idx_gene <- match(count.df$gene, names(plddt.es))
      matched_gene <- is.na(count.df$plddt) & !is.na(idx_gene)
      if (any(matched_gene)) {
        count.df$plddt[matched_gene] <- as.numeric(plddt.es[idx_gene[matched_gene]])
      }
    }
    
    if (!is.null(sd.es)) {
      if (n.es > 1.1) sd.es <- sd.es / n.es
      count.df$SD_es <- 0
      for(n in nm){
        if(!is.na(sd.es[n])) count.df[n, "SD_es"] <- sd.es[n]
      }
    }
    
    if (!is.null(sd.rep)) {
      if (n.rep > 1.1) sd.rep <- sd.rep / n.rep
      count.df$SD_rep <- 0
      for(n in nm){
        if(!is.na(sd.rep[n])) count.df[n, "SD_rep"] <- sd.rep[n]
      }
    }
    
    if (is.null(sd.es) || !plot.sd) {
      lim.y <- max(count.df[, "Y"]) * 1.4
    } else {
      lim.y <- max(max(count.df[, "Y"]) * 1.4, count.df[, "Y"] + count.df[, "SD_es"])
    }
    
    if (is.null(sd.rep) || !plot.sd) {
      lim.x <- max(count.df[, "X"]) * 1.3
    } else {
      lim.x <- max(max(count.df[, "X"]) * 1.3, count.df[, "X"] + count.df[, "SD_rep"])
    }
    
    ratio <- (count.df[, "Y"] + 0.0001) / (count.df[, "X"] + 0.0001)
    logFC <- log2(ratio)
    names(logFC) <- nm
    count.df$logFC <- as.character(signif(logFC, 2))
    
    if (!is.null(sd.rep)) {
      Zscore <- sapply(nm, function(n){
        if(count.df[n, "SD_rep"] > 0){
          (count.df[n, "Y"] - count.df[n, "X"]) / count.df[n, "SD_rep"]
        } else {
          if(abs(count.df[n, "Y"] - count.df[n, "X"]) > 0.0001) 5 else 0
        }
      })
    } else {
      Zscore <- setNames(rep(NA, length(nm)), nm)
    }
    count.df$Zscore <- as.character(signif(Zscore, 2))
    
    if (!is.null(distr.rep)) {
      Pval <- sapply(nm, function(n){
        if(n %in% rownames(distr.rep)){
          Ppos <- length(which(distr.rep[n, ] >= count.df[n, "Y"])) / length(distr.rep[n, ])
          Pneg <- length(which(distr.rep[n, ] <= count.df[n, "Y"])) / length(distr.rep[n, ])
          max(1 / dim(distr.rep)[2], min(Ppos, Pneg))
        } else {
          if(count.df[n, "Y"] > 0) 1 / dim(distr.rep)[2] else 1
        }
      })
    }
    
    label <- nm
    names(label) <- nm
    label.all <- label
    
    if (!label.neg) {
      label[logFC < 0] <- NA
    }
    
    label[count.df[, "Y"] < label.min.fr[1] & count.df[, "X"] < label.min.fr[2]] <- NA
    
    if (!is.null(sd.rep) > 0) {
      label[abs(Zscore) < ZscoreVJ.thresh] <- NA
    }
    
    n_lab_max <- ifelse(floor(pType) == 2, 12, 8)
    
    min.logFC <- log2(c(FoldChangeVJ.thresh, 1.5, 2, 3))
    min.fr <- c(label.diag, 0.2, 0.1, 0.1)
    
    label[which((abs(logFC) < min.logFC[1]) &
                  count.df[, "X"] < min.fr[1] &
                  count.df[, "Y"] < min.fr[1])] <- NA
    
    t <- 2
    while(length(which(!is.na(label))) > n_lab_max && t <= length(min.fr)){
      label[which((abs(logFC) < min.logFC[t]) &
                    count.df[, "X"] < min.fr[t] &
                    count.df[, "Y"] < min.fr[t])] <- NA
      t <- t + 1
    }
    
    if(length(which(!is.na(label))) > n_lab_max){
      if(verbose >= 1){
        print(paste("Too many labels for", segment, ", selection will be based primarily on logFC:"))
        print("  Consider augmenting values in label.min.fr if visualisation is not good")
      }
      logFC.sort <- sort(abs(logFC))
      t <- 1
      while(length(which(!is.na(label))) > n_lab_max){
        label[names(logFC.sort[t])] <- NA
        t <- t + 1
      }
    }
    
    if(label.diag){
      ind <- which(count.df[, "Y"] > label.diag &
                     count.df[, "X"] > label.diag &
                     abs(logFC) < 1.5)
      if(length(ind) < 10){
        label[ind] <- label.all[ind]
      } else {
        print("Warning: the number of labels to be shown along the diagonal is very large. They will not be shown. Consider increasing label.diag")
      }
    }
    
    if(print.size){
      ylab <- paste(info["input1.name"], " (", n.es, ")", sep = "")
    } else {
      ylab <- info["input1.name"]
    }
    
    if(comp.baseline || !print.size){
      xlab <- info["baseline.name"]
    } else {
      xlab <- paste(info["baseline.name"], " (", n.rep, ")", sep = "")
    }
    
    if (floor(pType) == 1) {
      count.df$label <- label
    } else {
      count.df$label <- count.df$name
    }
    
    namesToKeep <- setdiff(label, NA)
    count.df$baselineName <- xlab
    count.df$log2FC <- log2(ratio + 1e-5)
    count.df$model <- paste0(info["model"], " (", n.es, ")")
    
    if (ret.resList) {
      return(list(count.df = count.df, namesToKeep = namesToKeep, segment = segment))
    }
    
  } else {
    count.df <- combined.resList$count.df
    count.df$model <- factor(count.df$model, levels = rev(unique(count.df$model)))
    namesToKeep <- combined.resList$namesToKeep
    segment <- combined.resList$segment
  }
  
  if (pType == 1) {
    colorScale <- TCRgene2aes[[species]][[segment]]$color1
  } else {
    colorScale <- TCRgene2aes[[species]][[segment]]$color2
  }
  
  if (floor(pType) == 1) {
    if (!is.null(combined.resList)) {
      stop("Didn't implement the use of combined.ResList when floor(pType) == 1")
    }
    
    if (info["gene"] == "TRAV" || info["gene"] == "TRBV") {
      CDR1_seq <- gsub("g", "-", cdr123[[species]][[substr(info["gene"], 1, 3)]][count.df$name, "CDR1"])
      CDR2_seq <- gsub("g", "-", cdr123[[species]][[substr(info["gene"], 1, 3)]][count.df$name, "CDR2"])
      CDR3_seq <- cdr123[[species]][[substr(info["gene"], 1, 3)]][count.df$name, "CDR3"]
      
      count.df$CDR1_seq <- CDR1_seq
      count.df$CDR2_seq <- CDR2_seq
      count.df$CDR3_seq <- CDR3_seq
      
      count.plot <- ggplot(
        count.df,
        aes(
          x = X, y = Y,
          label = label, label2 = name, label3 = logFC, label4 = Zscore,
          label5 = CDR1_seq, label6 = CDR2_seq, label7 = CDR3_seq
        )
      ) +
        geom_abline(col = "orange", linetype = "dashed", linewidth = 1)
      
    } else if (info["gene"] == "TRAJ" || info["gene"] == "TRBJ") {
      CDR3_seq <- Jseq[[species]][[substr(info["gene"], 1, 3)]][count.df$name, "CDR3"]
      count.df$CDR3_seq <- CDR3_seq
      
      count.plot <- ggplot(
        count.df,
        aes(
          x = X, y = Y,
          label = label, label2 = name, label3 = logFC, label4 = Zscore,
          label5 = CDR3_seq
        )
      ) +
        geom_abline(col = "orange", linetype = "dashed", linewidth = 1)
    }
    
    if (!is.null(sd.es) && plot.sd) {
      count.plot <- count.plot +
        geom_errorbar(
          aes(ymax = Y + SD_es, ymin = sapply(Y - SD_es, function(x) max(0.001, x))),
          width = 0.015 * lim.x, linewidth = 0.4, color = "grey50"
        )
    }
    
    if (!is.null(sd.rep) && plot.sd) {
      count.plot <- count.plot +
        geom_errorbarh(
          aes(xmax = X + SD_rep, xmin = sapply(X - SD_rep, function(x) max(0.001, x))),
          height = 0.015 * lim.y, linewidth = 0.4, color = "grey50"
        )
    }
    
    if (pType == 1.3) {
      count.plot <- count.plot + geom_point()
      
    } else if (pType == 1.2) {
      outerWidth <- 0.5
      
      if (!is.null(plddt.es)) {
        count.plot <- count.plot +
          geom_point(
            aes(fill = plddt),
            shape = 21, color = "gray20",
            stroke = outerWidth, size = 2.5
          ) +
          scale_fill_gradient2(
            low = "red",
            mid = "white",
            high = "blue",
            midpoint = 0.75,
            limits = c(0.5, 1),
            oob = scales::squish,
            na.value = "grey85",
            name = "pLDDT",
            guide = if (show.plddt.legend) "colourbar" else "none"
          )
      } else {
        count.plot <- count.plot +
          geom_point(
            aes(fill = gene),
            shape = 21, color = "gray20",
            stroke = outerWidth, size = 2.5
          ) +
          scale_fill_manual(values = colorScale, guide = "none")
      }
      
    } else {
      outerWidth <- 0.5
      
      if (!is.null(plddt.es)) {
        count.plot <- count.plot +
          geom_point(
            aes(fill = plddt, shape = gene, color = gene),
            stroke = outerWidth, size = 2.5
          ) +
          scale_fill_gradient2(
            low = "red",
            mid = "white",
            high = "blue",
            midpoint = 0.75,
            limits = c(0.5, 1),
            oob = scales::squish,
            na.value = "grey85",
            name = "pLDDT",
            guide = if (show.plddt.legend) "colourbar" else "none"
          )
      } else {
        count.plot <- count.plot +
          geom_point(
            aes(fill = gene, shape = gene, color = gene),
            stroke = outerWidth, size = 2.5
          ) +
          scale_fill_manual(values = colorScale, guide = "none")
      }
      
      count.plot <- count.plot +
        scale_color_manual(
          values = TCRgene2aes[[species]][[segment]]$outerColor1,
          guide = "none"
        ) +
        scale_shape_manual(
          values = TCRgene2aes[[species]][[segment]]$shape1,
          guide = "none"
        )
    }
    
    if (!is.null(plddt.es)) {
      count.plot <- count.plot +
        ggrepel::geom_label_repel(
          data = count.df[!is.na(count.df$label), ],
          aes(x = X, y = Y, label = label, fill = plddt),
          size = 3,
          nudge_y = 0.02,
          box.padding = 0.25,
          show.legend = FALSE,
          na.rm = TRUE,
          color = "black"
        )
    } else {
      count.plot <- count.plot +
        ggrepel::geom_label_repel(
          data = count.df[!is.na(count.df$label), ],
          aes(x = X, y = Y, label = label),
          size = 3,
          nudge_y = 0.02,
          box.padding = 0.25,
          show.legend = FALSE,
          na.rm = TRUE,
          fill = "white",
          color = "black"
        )
    }
    
    count.plot <- count.plot +
      ggtitle(segment) +
      theme_bw() +
      scale_x_continuous(limits = c(0, lim.x), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, lim.y), expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      xlab(xlab) +
      ylab(ylab) +
      theme(
        plot.title = element_text(size = 14, hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_blank()
      )
    
    if (show.plddt.legend) {
      count.plot <- count.plot +
        theme(
          legend.position = c(0.05, 0.95),   # (x, y) inside panel
          legend.justification = c(0, 1),    # anchor top-left
          legend.background = element_rect(fill = "white", color = "white"),
          legend.key.height = unit(0.4, "cm"),
          legend.key.width  = unit(0.2, "cm")
        )
    }
    
  } else {
    count_other <- count.df[!count.df$name %in% namesToKeep, , drop = FALSE] %>%
      dplyr::summarize(
        Y = sum(Y), X = 0, name = "Other", gene = "Other",
        label = "Other", log2FC = NA, baselineName = count.df$baselineName[1],
        .by = model
      )
    
    count.df <- count.df[count.df$name %in% namesToKeep, , drop = FALSE]
    count.df <- count.df[order(count.df$log2FC, decreasing = TRUE), , drop = FALSE]
    
    if (pType == 2.1) {
      count.df <- dplyr::bind_rows(count.df, count_other)
    }
    
    count.df$label <- factor(count.df$label, levels = rev(unique(count.df$label)))
    
    labYlen <- max(nchar(levels(count.df$label)))
    YaxSize <- ifelse(labYlen < 8, 14,
                      ifelse(labYlen < 10, 13,
                             ifelse(labYlen < 12, 12,
                                    ifelse(labYlen < 14, 11, 10))))
    
    legLen <- max(nchar(levels(count.df$model)))
    legSize <- ifelse(legLen < 15, 12,
                      ifelse(legLen < 20, 10,
                             ifelse(legLen < 25, 8, 6)))
    
    if (is.null(combined.resList)) {
      count.plot <- ggplot(count.df, aes(x = Y, y = label, fill = gene)) +
        geom_col(color = "gray10", linewidth = 1)
      figTitle <- paste0(segment, " (", n.es, ")")
    } else {
      count.plot <- ggplot(count.df, aes(x = Y, y = label, fill = gene, color = model)) +
        geom_col(position = "dodge", linewidth = 1, width = 0.8) +
        scale_color_manual(
          values = set_model_colPals(rev(levels(count.df$model))),
          guide = guide_legend(ncol = 1, order = 1, reverse = TRUE)
        )
      figTitle <- segment
    }
    
    count.plot <- count.plot +
      scale_fill_manual(values = colorScale, guide = "none") +
      ggtitle(figTitle) +
      xlab("Frequency") +
      ylab(NULL) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.05)), n.breaks = 3) +
      geom_col(
        aes(x = X, alpha = baselineName),
        position = "dodge",
        width = ifelse(is.null(combined.resList), 0.9, 0.8),
        linewidth = 0.5, color = "gray60", fill = "gray80"
      ) +
      scale_alpha_manual(values = 0.7, guide = guide_legend(order = 2)) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 14),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = YaxSize, face = "bold", color = "black"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = legSize)
      )
  }
  
  return(count.plot)
}

# ret.resList when TRUE indicates to return a list of the results instead of the plot directly
# ret.resList when TRUE indicates to return a list of the results instead of the plot directly
plotLD <- function(countL.es, countL.rep, plddtL.es = NULL, info = NULL,
                   sd.es = NULL, sd.rep = NULL, plot.oneline = 0,
                   ret.resList = FALSE, combined.resList = NULL,
                   comp.baseline = TRUE, print.size = TRUE, plot.sd = TRUE,
                   pseudo = 1e-4, show_baseline_ref = TRUE) {
  
  if (is.null(combined.resList)) {
    
    L.all <- Lmin:Lmax
    cn <- paste("L", L.all, sep = "_")
    
    ld.es <- rep(0, length(L.all)); names(ld.es) <- cn
    ld.rep <- rep(0, length(L.all)); names(ld.rep) <- cn
    n.es <- sum(countL.es)
    n.rep <- sum(countL.rep)
    
    for (lc in cn) {
      if (!is.na(countL.es[lc])) ld.es[lc] <- countL.es[lc]
      if (!is.na(countL.rep[lc])) ld.rep[lc] <- countL.rep[lc]
    }
    
    if (sum(ld.es) > 0) ld.es <- ld.es / sum(ld.es)
    if (sum(ld.rep) > 0) ld.rep <- ld.rep / sum(ld.rep)
    
    lds.es <- rep(0, length(L.all)); names(lds.es) <- cn
    lds.rep <- rep(0, length(L.all)); names(lds.rep) <- cn
    
    if (!is.null(sd.es)) {
      if (n.es > 1.1) sd.es <- sd.es / n.es
      for (lc in cn) {
        if (!is.na(sd.es[lc])) lds.es[lc] <- sd.es[lc]
      }
    }
    
    if (!is.null(sd.rep)) {
      if (n.rep > 1.1) sd.rep <- sd.rep / n.rep
      for (lc in cn) {
        if (!is.na(sd.rep[lc])) lds.rep[lc] <- sd.rep[lc]
      }
    }
    
    # normalize input by baseline
    ld.norm <- (ld.es + pseudo) / (ld.rep + pseudo)
    
    # optional approximate propagated SD for the ratio
    lds.norm <- rep(NA_real_, length(L.all))
    names(lds.norm) <- cn
    if (!is.null(sd.es) || !is.null(sd.rep)) {
      for (lc in cn) {
        A <- ld.es[lc]
        B <- ld.rep[lc]
        sA <- lds.es[lc]
        sB <- lds.rep[lc]
        if ((B + pseudo) > 0) {
          lds.norm[lc] <- sqrt(
            (sA / (B + pseudo))^2 +
              (((A + pseudo) * sB) / (B + pseudo)^2)^2
          )
        }
      }
    }
    
    ld.df <- data.frame(
      v1 = L.all,
      v2 = as.numeric(ld.norm),
      SD = as.numeric(lds.norm),
      group = "Input/Baseline"
    )
    
    # add pLDDT for normalized input points
    ld.df$plddt <- NA_real_
    if (!is.null(plddtL.es)) {
      len_names <- paste0("L_", ld.df$v1)
      match_idx <- match(len_names, names(plddtL.es))
      ld.df$plddt <- as.numeric(plddtL.es[match_idx])
    }
    
    if (ret.resList) {
      if (length(info) < 4) {
        stop("The 'info' vector given in plotLD input should have 4 elements when ret.resList is TRUE.")
      }
      ld.df$model <- paste0(info["model"], gsub(".*( \\(\\d+\\))", "\\1", info["input1.name"]))
      return(list(ld.df = ld.df, info = info["chain"]))
    }
    
    legend.size <- 12
    if (plot.oneline != 0) {
      if (nchar(info["input1.name"]) > 23) legend.size <- 11
      if (nchar(info["input1.name"]) > 25) legend.size <- 10
    }
    
    ld.plot <- ggplot(ld.df, aes(x = v1, y = v2))
    
    if (show_baseline_ref) {
      ld.plot <- ld.plot +
        geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.5)
    }
    
    if (plot.sd && any(!is.na(ld.df$SD))) {
      ld.plot <- ld.plot +
        geom_errorbar(
          aes(
            ymax = v2 + SD,
            ymin = pmax(v2 - SD, 0)
          ),
          width = 0.3,
          linewidth = 0.4,
          color = "grey60"
        )
    }
    
  } else {
    ld.df <- combined.resList$ld.df
    ld.df$model <- factor(ld.df$model, levels = unique(ld.df$model))
    info <- combined.resList$info
    legend.size <- 10
    
    ld.plot <- ggplot(ld.df, aes(x = v1, y = v2, color = model)) +
      facet_grid(rows = vars(model), scales = "free_y") +
      scale_color_manual(values = set_model_colPals(levels(ld.df$model))) +
      guides(color = "none")
    
    if (show_baseline_ref) {
      ld.plot <- ld.plot +
        geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.5)
    }
  }
  
  size <- ifelse(plot.oneline == 0, 2.5, 2)
  
  if (is.null(combined.resList)) {
    
    # normalized curve in grey
    ld.plot <- ld.plot +
      geom_line(color = "grey70", linewidth = 0.7)
    
    # points colored by pLDDT if available
    if (!is.null(plddtL.es)) {
      ld.plot <- ld.plot +
        geom_point(
          aes(fill = plddt),
          shape = 21,
          color = "black",
          size = size,
          stroke = 0.6
        ) +
        scale_fill_gradient2(
          low = "red",
          mid = "white",
          high = "blue",
          midpoint = 0.75,
          na.value = "grey85",
          name = "pLDDT",
          guide = "none"
        )
    } else {
      ld.plot <- ld.plot +
        geom_point(size = size, shape = 16, color = "black")
    }
    
  } else {
    ld.plot <- ld.plot +
      geom_line(linewidth = 0.7) +
      geom_point(size = size)
  }
  
  ld.plot <- ld.plot +
    theme_bw() +
    theme(
      legend.key.size = unit(0.65, "cm"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = legend.size)
    ) +
    xlab(paste("Length_CDR3", info["chain"], sep = "")) +
    ylab("Input / Baseline") +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 15, hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  return(ld.plot)
}

make_plddt_logo <- function(pwm_mat, plddt_mat, title = "", logo.type = "bits",
                            axisTextSizeX = 10, axisTextSizeY = 8,
                            plddt_limits = c(0.5, 1.0)) {
  
  p0 <- ggseqlogoMOD::ggseqlogoMOD(
    data = pwm_mat,
    additionaAA = additionalAA,
    axisTextSizeX = axisTextSizeX,
    axisTextSizeY = axisTextSizeY,
    methods = logo.type
  )
  
  dat <- p0$layers[[1]]$data
  
  letter_col <- intersect(c("letter", "aa", "AA", "label"), names(dat))[1]
  dat$letter <- toupper(trimws(as.character(dat[[letter_col]])))
  rownames(plddt_mat) <- toupper(trimws(rownames(plddt_mat)))
  
  dat$.rowid <- seq_len(nrow(dat))
  dat$pos <- if ("position" %in% names(dat)) {
    as.integer(dat$position)
  } else if ("pos" %in% names(dat)) {
    as.integer(dat$pos)
  } else {
    round(dat$x)
  }
  
  dat$poly_id <- interaction(dat$pos, dat$letter, drop = TRUE)
  
  poly_meta <- dat[!duplicated(dat$poly_id), c("poly_id", "letter", "pos")]
  poly_meta$plddt <- mapply(
    function(letter, pos) {
      if (!is.na(letter) &&
          !is.na(pos) &&
          letter %in% rownames(plddt_mat) &&
          pos >= 1 &&
          pos <= ncol(plddt_mat)) {
        as.numeric(plddt_mat[letter, pos])
      } else {
        NA_real_
      }
    },
    poly_meta$letter,
    poly_meta$pos
  )
  
  dat$plddt <- as.numeric(poly_meta$plddt[match(dat$poly_id, poly_meta$poly_id)])
  dat <- dat[order(dat$poly_id, dat$.rowid), ]
  
  gb <- ggplot_build(p0)
  xr <- gb$layout$panel_params[[1]]$x.range
  yr <- gb$layout$panel_params[[1]]$y.range
  
  n_pos <- ncol(pwm_mat)
  
  ggplot() +
    geom_polygon(
      data = dat,
      aes(x = x, y = y, group = poly_id, fill = plddt),
      color = 'grey85',
      inherit.aes = FALSE
    ) +
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0.75,
      limits = plddt_limits,
      oob = scales::squish,
      na.value = "grey85",
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = 1:n_pos,
      labels = 1:n_pos,
      expand = c(0, 0)
    ) +
    coord_cartesian(xlim = xr, ylim = yr, expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      #axis.line.x = element_line(linewidth = 0.8, colour = "black"),
    #axis.line.y = element_line(linewidth = 0.8, colour = "black"),
      #axis.ticks.x = element_line(linewidth = 0.8, colour = "black"),
      #axis.ticks.y = element_line(linewidth = 0.8, colour = "black"),
      #axis.text.x = element_text(size = axisTextSizeX, colour = "black"),
      #axis.text.y = element_text(size = axisTextSizeY, colour = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(l=25, b=20, t=5, r=5)
    )
}

plotCDR3 <- function(countL.es, countL.rep, countCDR3.es, countCDR3.rep, plddtCDR3.L.es = NULL,info = NULL,
                     comp.baseline = TRUE, plot.oneline = 0, plot.all.length = FALSE,
                     logo.type = "bits", plot.cdr3.subtract.baseline = 0,
                     set.cdr3.length = NA, print.size = TRUE){
  
  L.es <- as.numeric(lapply(names(countL.es), function(x) unlist(strsplit(x, split = "_"))[2]))
  L.rep <- as.numeric(lapply(names(countL.rep), function(x) unlist(strsplit(x, split = "_"))[2]))
  
  L.TR <- intersect(L.es, L.rep)
  
  pwm.rep <- list()
  pwm.es <- list()
  
  logo.CDR3.L.es <- list()
  logo.CDR3.L.rep <- list()
  
  if(length(L.TR) > 0){
    
    tl <- countL.es[paste("L", L.TR, sep = "_")]
    
    if(is.na(set.cdr3.length)){
      lmax <- as.numeric(unlist(strsplit(names(tl[which.max(tl)]), split = "_"))[2])
    } else {
      if(set.cdr3.length %in% L.TR){
        lmax <- set.cdr3.length
      } else {
        lmax <- as.numeric(unlist(strsplit(names(tl[which.max(tl)]), split = "_"))[2])
        print(paste("set.cdr3", info["chain"], ".length=", set.cdr3.length,
                    " is incompatible with the input data. Default value of ", lmax,
                    " will be used.", sep = ""))
      }
    }
    
    if(plot.oneline != 0){
      axis.size.max <- ifelse(lmax < 15, 8, 7)
      title.size <- 11
    } else {
      axis.size.max <- 10
      title.size <- 12
    }
    
    ylab <- ""
    
    if(!plot.all.length){
      L.TR <- c(lmax)
    }
    
    for(l in L.TR){
      
      lc <- paste("L", l, sep = "_")
      
      pwm.es[[lc]] <- scale(countCDR3.es[[lc]], center = FALSE, scale = colSums(countCDR3.es[[lc]]))
      pwm.rep[[lc]] <- scale(countCDR3.rep[[lc]], center = FALSE, scale = colSums(countCDR3.rep[[lc]]))
      
      if(plot.cdr3.subtract.baseline == 1){
        size.es <- apply(pwm.es[[lc]], 2, function(x){
          ind <- which(x != 0)
          IC <- log(N.aa)/log(2) + sum(x[ind] * log(x[ind]) / log(2))
          IC * x
        })
        size.rep <- apply(pwm.rep[[lc]], 2, function(x){
          ind <- which(x != 0)
          IC <- log(N.aa)/log(2) + sum(x[ind] * log(x[ind]) / log(2))
          IC * x
        })
        x.norm <- size.es - size.rep
        y.inc <- 1
      }
      
      if(plot.cdr3.subtract.baseline == 2){
        rc <- max(5, 0.1 * countL.es[[lc]])
        pseudo <- 1 / N.aa * rc / countL.es[[lc]]
        x.baseline <- pwm.rep[[lc]] + pseudo
        x.es <- pwm.es[[lc]] + pseudo
        x.norm <- x.es / x.baseline
        x.norm <- scale(x.norm, center = FALSE, scale = colSums(x.norm))
        y.inc <- 4
      }
      
      title <- info["input1.name"]
      if(print.size) title <- paste(title, " (", countL.es[[lc]], ")", sep = "")
      title <- paste(title, ", CDR3", info["chain"], "_", l, sep = "")
      
      if (!is.null(plddtCDR3.L.es) && lc %in% names(plddtCDR3.L.es)) {
        logo.CDR3.L.es[[lc]] <- make_plddt_logo(
          pwm_mat = pwm.es[[lc]],
          plddt_mat = plddtCDR3.L.es[[lc]],
          title = title,
          logo.type = logo.type,
          axisTextSizeX = 12,
          axisTextSizeY = 8,
          plddt_limits = c(0.5, 1.0)
        )
      } else {
        logo.CDR3.L.es[[lc]] <-
          ggseqlogoMOD::ggseqlogoMOD(
            data = pwm.es[[lc]], additionaAA = additionalAA,
            axisTextSizeX = 12, axisTextSizeY = 8, methods = logo.type
          ) +
          labs(title = title) +
          ylab(ylab) +
          theme(plot.title = element_text(size = 15, hjust = 0.5))
      }
      
      title.baseline <- info["baseline.name"]
      if(!comp.baseline && print.size){
        title.baseline <- paste(title.baseline, " (", countL.rep[[lc]], ")", sep = "")
      }
      
      if(plot.cdr3.subtract.baseline == 0){
        title.baseline <- paste(title.baseline, ", CDR3", info["chain"], "_", l, sep = "")
      } else if(plot.cdr3.subtract.baseline == 1){
        title.baseline <- paste(title.baseline, " subtract, CDR3", info["chain"], "_", l, sep = "")
      } else if(plot.cdr3.subtract.baseline == 2){
        title.baseline <- paste(title.baseline, " renorm, CDR3", info["chain"], "_", l, sep = "")
      }
      
      if(plot.cdr3.subtract.baseline == 0){
        logo.CDR3.L.rep[[lc]] <-
          ggseqlogoMOD::ggseqlogoMOD(
            data = pwm.rep[[lc]], additionaAA = additionalAA,
            axisTextSizeX = 12, axisTextSizeY = 8, methods = logo.type
          ) +
          labs(title = title.baseline) +
          ylab(ylab) +
          theme(plot.title = element_text(size = 15, hjust = 0.5))
        
      } else if(plot.cdr3.subtract.baseline == 1){
        y.min <- min(apply(x.norm, 2, function(x) sum(x[x < 0])))
        y.max <- max(apply(x.norm, 2, function(x) sum(x[x > 0])))
        y.min <- max(-log(N.aa)/log(2), y.inc * y.min)
        y.max <- log(N.aa)/log(2)
        
        logo.CDR3.L.rep[[lc]] <-
          ggseqlogoMOD::ggseqlogo(data = x.norm, method = "custom") +
          labs(title = title.baseline) +
          ylim(y.min, y.max) +
          ylab(ylab) +
          theme(plot.title = element_text(size = title.size, hjust = 0.5)) +
          theme(legend.position = "none")
        
      } else if(plot.cdr3.subtract.baseline == 2){
        IC.max <- max(unlist(apply(x.norm, 2, function(x){
          ind <- which(x != 0)
          log(N.aa)/log(2) + sum(x[ind] * log(x[ind]) / log(2))
        })))
        y.max <- min(IC.max * y.inc, log(N.aa)/log(2))
        
        logo.CDR3.L.rep[[lc]] <-
          ggseqlogoMOD::ggseqlogoMOD(
            data = x.norm, additionaAA = additionalAA,
            axisTextSizeX = 12, axisTextSizeY = 8,
            ylim = c(0, y.max), methods = logo.type
          ) +
          labs(title = title.baseline) +
          ylab(ylab) +
          theme(plot.title = element_text(size = 15, hjust = 0.5))
      }
      
      if(l == lmax){
        title <- info["input1.name"]
        if(print.size) title <- paste(title, " (", countL.es[[lc]], ")", sep = "")
        title <- paste(title, ", CDR3", info["chain"], "_", l, sep = "")
        
        title.baseline <- info["baseline.name"]
        if(!comp.baseline && print.size){
          title.baseline <- paste(title.baseline, " (", countL.rep[[lc]], ")", sep = "")
        }
        
        if(plot.cdr3.subtract.baseline == 0){
          title.baseline <- paste(title.baseline, ", CDR3", info["chain"], "_", l, sep = "")
        } else if(plot.cdr3.subtract.baseline == 1){
          title.baseline <- paste(title.baseline, " subtract, CDR3", info["chain"], "_", l, sep = "")
        } else if(plot.cdr3.subtract.baseline == 2){
          title.baseline <- paste(title.baseline, " renorm, CDR3", info["chain"], "_", l, sep = "")
        }
        
        if(plot.oneline != 0 && (nchar(title) > 26 || nchar(title.baseline) > 26)){
          title <- paste("CDR3", info["chain"], "_", l, " ", info["input1.name"], sep = "")
          if(print.size) title <- paste(title, "\n(", countL.es[[lc]], ")", sep = "")
          
          if(plot.cdr3.subtract.baseline == 0){
            title.baseline <- paste("CDR3", info["chain"], "_", l, " ", info["baseline.name"], "\n", sep = "")
          } else if(plot.cdr3.subtract.baseline == 1){
            title.baseline <- paste("CDR3", info["chain"], "_", l, " subtract \n", info["baseline.name"], sep = "")
          } else if(plot.cdr3.subtract.baseline == 2){
            title.baseline <- paste("CDR3", info["chain"], "_", l, " renorm \n", info["baseline.name"], sep = "")
          }
          
          if(!comp.baseline && print.size){
            title.baseline <- paste(title.baseline, "(", countL.rep[[lc]], ")", sep = "")
          }
        }
        
        if (!is.null(plddtCDR3.L.es) && lc %in% names(plddtCDR3.L.es)) {
          logo.CDR3.L.es.max <- make_plddt_logo(
            pwm_mat = pwm.es[[lc]],
            plddt_mat = plddtCDR3.L.es[[lc]],
            title = title,
            logo.type = logo.type,
            axisTextSizeX = axis.size.max,
            axisTextSizeY = 8,
            plddt_limits = c(0.5, 1.0)
          )
        } else {
          logo.CDR3.L.es.max <-
            ggseqlogoMOD::ggseqlogoMOD(
              data = pwm.es[[lc]], additionaAA = additionalAA,
              axisTextSizeX = axis.size.max, axisTextSizeY = 8, methods = logo.type
            ) +
            labs(title = title) +
            ylab(ylab) +
            theme(plot.title = element_text(size = title.size, hjust = 0.5))
        }
        
        if(plot.cdr3.subtract.baseline == 0){
          logo.CDR3.L.rep.max <-
            ggseqlogoMOD::ggseqlogoMOD(
              data = pwm.rep[[lc]], additionaAA = additionalAA,
              axisTextSizeX = axis.size.max, axisTextSizeY = 8, methods = logo.type
            ) +
            labs(title = title.baseline) +
            ylab(ylab) +
            theme(plot.title = element_text(size = title.size, hjust = 0.5))
          
        } else if(plot.cdr3.subtract.baseline == 1){
          y.min <- min(apply(x.norm, 2, function(x) sum(x[x < 0])))
          y.max <- max(apply(x.norm, 2, function(x) sum(x[x > 0])))
          y.min <- max(-log(N.aa)/log(2), 2 * y.min)
          y.max <- log(N.aa)/log(2)
          
          logo.CDR3.L.rep.max <-
            ggseqlogoMOD::ggseqlogo(data = x.norm, method = "custom") +
            ggtitle(title.baseline) +
            ylim(y.min, y.max) +
            ylab(ylab) +
            theme(plot.title = element_text(size = title.size, hjust = 0.5)) +
            theme(legend.position = "none")
          
        } else if(plot.cdr3.subtract.baseline == 2){
          IC.max <- max(unlist(apply(x.norm, 2, function(x){
            ind <- which(x != 0)
            log(N.aa)/log(2) + sum(x[ind] * log(x[ind]) / log(2))
          })))
          y.max <- min(IC.max * y.inc, log(N.aa)/log(2))
          
          logo.CDR3.L.rep.max <-
            ggseqlogoMOD::ggseqlogoMOD(
              data = x.norm, additionaAA = additionalAA,
              axisTextSizeX = axis.size.max, axisTextSizeY = 8,
              ylim = c(0, y.max), methods = logo.type
            ) +
            labs(title = title.baseline) +
            ylab(ylab) +
            theme(plot.title = element_text(size = title.size, hjust = 0.5))
        }
      }
    }
    
    ls <- list(
      logo.CDR3.L.es, logo.CDR3.L.rep, L.TR, lmax,
      logo.CDR3.L.es.max, logo.CDR3.L.rep.max
    )
    names(ls) <- c("ES", "Baseline", "length", "lmax", "ES_max", "Baseline_max")
    
  } else {
    ls <- list(list(), list(), c(), 0, ggplot(), ggplot())
    names(ls) <- c("ES", "Baseline", "length", "lmax", "ES_max", "Baseline_max")
  }
  
  return(ls)
}

#' Defines the color palette for plots with models shown together.
set_model_colPals <- function(models){
  colPal_modelsCombined <- c(
    palette.colors(palette = "Okabe-Ito"),
    palette.colors(palette = "Tableau 10"),
    palette.colors(palette = "Polychrome 36")
  )
  if (length(models) > length(colPal_modelsCombined)){
    colPal_modelsCombined <- c(
      colPal_modelsCombined,
      rainbow(n = length(models) - length(colPal_modelsCombined), alpha = NULL)
    )
  }
  colPal_modelsCombined <- colPal_modelsCombined[1:length(models)]
  names(colPal_modelsCombined) <- models
  if (any(grepl("Trash", models))){
    colPal_modelsCombined[grep("Trash", models)] <- "gray70"
  }
  return(colPal_modelsCombined)
}

create_interactive_plots <- function(countV.plot, countJ.plot, ld.plot, CDR3, plot.oneline){
  countV.plot_not_title <- countV.plot + labs(title = NULL)
  
  p1 <- plotly::ggplotly(
    countV.plot_not_title,
    tooltip = c("name", "logFC", "Zscore", "CDR1_seq", "CDR2_seq", "CDR3_seq")
  )
  
  p1$x$data <- lapply(p1$x$data, function(trace) {
    if (!is.null(trace$marker$size)) trace$marker$size <- 11
    if (!is.null(trace$error_x)){
      trace$error_x$width <- 3
      trace$error_x$color <- rgb(102/255, 102/255, 102/255, alpha = 0.4)
    }
    trace
  })
  
  p1 <- p1 %>%
    plotly::layout(
      title = NULL,
      margin = list(t = 80),
      annotations = list(
        list(
          text = unname(countV.plot$labels$title),
          x = 0.5, y = 1.04,
          xref = "paper", yref = "paper",
          xanchor = "center", yanchor = "bottom",
          font = list(size = 20),
          showarrow = FALSE
        )
      ),
      showlegend = FALSE
    )
  
  countJ.plot_not_title <- countJ.plot + labs(title = NULL)
  
  p2 <- plotly::ggplotly(
    countJ.plot_not_title,
    tooltip = c("name", "logFC", "Zscore", "CDR3_seq")
  )
  
  p2$x$data <- lapply(p2$x$data, function(trace) {
    if (!is.null(trace$marker$size)) trace$marker$size <- 11
    if (!is.null(trace$error_x)){
      trace$error_x$width <- 3
      trace$error_x$color <- rgb(102/255, 102/255, 102/255, alpha = 0.4)
    }
    trace
  })
  
  p2 <- p2 %>%
    plotly::layout(
      title = NULL,
      margin = list(t = 80),
      annotations = list(
        list(
          text = unname(countJ.plot$labels$title),
          x = 0.5, y = 1.04,
          xref = "paper", yref = "paper",
          xanchor = "center", yanchor = "bottom",
          font = list(size = 20),
          showarrow = FALSE
        )
      ),
      showlegend = FALSE
    )
  
  p3 <- plotly::ggplotly(ld.plot, tooltip = "none") %>%
    plotly::layout(
      legend = list(
        orientation = "h",
        x = 0.5,
        xanchor = "center",
        y = 1.125,
        yanchor = "bottom",
        title = list(text = NULL)
      ),
      margin = list(t = 50)
    )
  
  if(plot.oneline != 2){
    p4 <- plotly::ggplotly(CDR3$ES_max, tooltip = "none")
    p5 <- plotly::ggplotly(CDR3$Baseline_max, tooltip = "none")
    bottom_subplot <- manipulateWidget::combineWidgets(p4, p5, ncol = 1, title = NULL)
  }
  
  if(plot.oneline == 0){
    combined_plots <- manipulateWidget::combineWidgets(
      p1, p2, p3, bottom_subplot,
      ncol = 2, title = NULL
    )
  }
  
  if(plot.oneline == 1){
    combined_plots <- manipulateWidget::combineWidgets(
      p1, p2, p3, bottom_subplot,
      ncol = 4, nrow = 1, title = NULL
    )
  }
  
  if(plot.oneline == 2){
    combined_plots <- manipulateWidget::combineWidgets(
      p1, p2, p3,
      ncol = 3, nrow = 1, title = NULL
    )
  }
  
  return(combined_plots)
}