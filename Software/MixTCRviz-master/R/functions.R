
#####
# Define some function
#####

#' @export
build_stat <- function(input, chain="AB", species="HomoSapiens", comp.VJL=0){

  # comp.VJL>=1 means we are computing length distributions and motifs knowing P(VJ)
  # It takes some time, but still reasonable.

  chain.list <- paste("TR",strsplit(chain,"")[[1]], sep="")

  es <- list()
  es$species <- species
  if(comp.VJL==0){
    stat.list <- c("L", "countL", "countV", "countJ", "countV.L", "countJ.L", "countCDR1", "countCDR2", "countCDR3.L", "countVJ", "countVJ.L")
  } else {
    stat.list <- c("L", "countL", "countV", "countJ", "countV.L", "countJ.L", "countL.VJ", "countCDR1", "countCDR2", "countCDR3.L", "countCDR3.VL", "countCDR3.JL", "countCDR3.VJL", "countVJ", "countVJ.L")
  }
  for(s in stat.list){
    es[[s]] <- list()
  }

  input[input==""] <- NA # Useful if using build_stat outside of MixTCRviz

  for(ch in chain.list){

    Vn <- paste(ch,"V", sep="")
    Jn <- paste(ch,"J", sep="")
    cdr3 <- paste("cdr3_",ch,sep="")

    es$countV[[ch]] <- table(input[,Vn])
    es$countJ[[ch]] <- table(input[,Jn])

    es$countL[[ch]] <- table(nchar(input[,cdr3]))
    es$L[[ch]] <- as.numeric(names(es$countL[[ch]]))
    if(length(es$countL[[ch]])>0){
      names(es$countL[[ch]]) <- paste("L",names(es$countL[[ch]]),sep="_")
    }
    es$countVJ[[ch]] <- table(input[,Vn], input[,Jn])

    es$countV.L[[ch]] <- list()
    es$countJ.L[[ch]] <- list()
    es$countVJ.L[[ch]] <- list()
    for(lg in es$L[[ch]]){
      ind <- which(nchar(input[,cdr3])==lg)
      lg.c <- paste("L",lg,sep="_")
      es$countV.L[[ch]][[lg.c]] <- table(input[ind,Vn])
      es$countJ.L[[ch]][[lg.c]] <- table(input[ind,Jn])
      es$countVJ.L[[ch]][[lg.c]] <- table(input[ind,Vn],input[ind,Jn])
    }

    if(comp.VJL>=1){

      for(V in names(es$countV[[ch]])){
        indv <- which(input[,Vn]==V)
        es$countCDR3.VL[[ch]][[V]] <- count_aa(input[indv,cdr3], keep.gap=0)
      }
      for(J in names(es$countJ[[ch]])){
        indj <- which(input[,Jn]==J)
        es$countCDR3.JL[[ch]][[J]] <- count_aa(input[indj,cdr3], keep.gap=0)
      }
      for(V in names(es$countV[[ch]])){
        indv <- which(input[,Vn]==V)

        for(J in names(es$countJ[[ch]])){
          indj <- which(input[indv,Jn]==J)
          ind <- indv[indj]
          s <- paste(V,J, sep="_")
          if(length(ind)>0){
            es$countL.VJ[[ch]][[s]] <- table(nchar(input[ind,cdr3]))
            if(length(es$countL.VJ[[ch]][[s]])>0){
              names(es$countL.VJ[[ch]][[s]]) <- paste("L", names(es$countL.VJ[[ch]][[s]]),sep="_")
            }
            es$countCDR3.VJL[[ch]][[s]] <- count_aa(input[ind,cdr3], keep.gap=0)
          } else {
            es$countL.VJ[[ch]][[s]] <- table(NA)
            es$countCDR3.VJL[[ch]][[s]] <- table(NA)
          }
        }
      }
    }
    if(length(es$countV[[ch]])>0){
      es$countCDR1[[ch]] <- count_aa(cdr123[[species]][[ch]][input[,Vn],"CDR1"], keep.gap=1)
      es$countCDR2[[ch]] <- count_aa(cdr123[[species]][[ch]][input[,Vn],"CDR2"], keep.gap=1)
    }
    es$countCDR3.L[[ch]] <- count_aa(input[,cdr3], keep.gap=0)


  }

  return(es)

}

#Compute the counts of each aa at each position.
#Would be better to have the option of treating gaps in CDR1/2 as separate amino acids, while missing data can be treated as 'unspecific'
count_aa <- function(cdr.seq, keep.gap=0){

  if(keep.gap == 0){    #All gaps, including "g" are treated as unspecific data. This can be useful for visualisation.
    tgap <- c(gap,"g")
    taa.list <- aa.list
  } else {   # "x" are treated as additonal aa, other 'gaps' are discarded. This is more correct for modelling CDR1/CDR2 loops.
    tgap <- c()
    taa.list <- c(aa.list,"g")
  }

  #First get the list of length
  l.seq <- nchar(cdr.seq)
  L <- sort(unique(l.seq))
  L <- L[L>0]  # This ensures that the cases of empty sequences are never counted.
  m.list <- list()
  for(lg in L){
    ind <- which(l.seq==lg)
    tcdr.seq <- cdr.seq[ind]
    m <- matrix(0, nrow=length(taa.list), ncol=lg)
    rownames(m) <- taa.list
    for(p in 1:lg){
      s <- substr(tcdr.seq,p,p)
      tb <- table(s)
      for(a in names(tb)){
        if(keep.gap==0){
          if(a %in% tgap){
            m[,p] <- m[,p] + tb[a]/20
          } else {
            m[a,p] <- m[a,p]+tb[a]
          }
        } else {  #Here "x" are treated as a separate amino acid and missing data (e.g., "*", "X",etc.) are not included
          if(a %in% gap==F){
            m[a,p] <- m[a,p]+tb[a]
          }
        }
      }
    }
    lc <- paste("L",lg,sep="_")
    m.list[[lc]] <- m
  }
  return(m.list)
}

#This assumes the matrix includes gaps (21 rows)
build_cdr12_motif <- function(cdr.seq, keep.gap=0){

  if(keep.gap==0){
    g <- cdr.seq["g",]
    for(a in aa.list){
      cdr.seq[a,] <- cdr.seq[a,]+g/length(aa.list)
    }
    cdr.seq <- cdr.seq[aa.list,]
  }
  cdr.seq <- scale(cdr.seq, center=F, scale=colSums(cdr.seq))
  return(cdr.seq)
}


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
    if(sum(countVJ.L.es[[lg.c]])>0){
      w <- countVJ.L.es[[lg.c]]/sum(countVJ.L.es[[lg.c]])
      for(v in rownames(countVJ.L.es[[lg.c]])){
        for(j in colnames(countVJ.L.es[[lg.c]])){
          if(w[v,j]>0){
            s <- paste(v,j,sep="_")
            if(lg.c %in% names(countCDR3.VJL.baseline[[s]])){
              if(length(countCDR3.VJL.baseline[[s]][[lg.c]]) > 0){
                countCDR3[[lg.c]] <- countCDR3[[lg.c]]+scale(countCDR3.VJL.baseline[[s]][[lg.c]], center=F, scale=colSums(countCDR3.VJL.baseline[[s]][[lg.c]]))*w[v,j]
              }
            }
          }
        }
      }
    }
    if(colSums(countCDR3[[lg.c]])[1]!=0){
      countCDR3[[lg.c]] <- scale(countCDR3[[lg.c]],center=F, scale=colSums(countCDR3[[lg.c]]))
    } else {
      countCDR3[[lg.c]][1:N.aa,1:lg] <- 0.05 #This is the case if the were no data for all V-J for this length in the repertoire. Ideally, we should indicate this in the title of the plot
    }
  }
  return(countCDR3)
}

# Compute the weighted average of CDR3 length for the observed V/J usage
# This is important since CDR3 length is primarily determined by the length of the V and J segments
#' @export
weighted_countL <- function(countL.VJ.baseline, countVJ.es){


  countL <- rep(0,times=Lmax-Lmin+1)
  names(countL) <- paste("L", Lmin:Lmax, sep="_")

  countVJ.es <- countVJ.es/sum(countVJ.es)
  for(v in rownames(countVJ.es)){
    for(j in colnames(countVJ.es)){
      if(countVJ.es[v,j]>0){
        s <- paste(v,j,sep="_")
        if(s %in% names(countL.VJ.baseline)){
          if(length(countL.VJ.baseline[[s]])>0){
            countL.VJ.baseline[[s]] <- countL.VJ.baseline[[s]]/sum(countL.VJ.baseline[[s]])*countVJ.es[v,j]
            for(lc in names(countL.VJ.baseline[[s]])){
              countL[lc] <- countL[lc]+countL.VJ.baseline[[s]][lc]
            }
          }
        }
      }
    }
  }
  return(countL)
}



#Infer the short name for the mhc from the long name (HLA-A*02:01 -> A0201)
#The input should be a matrix/dataframe with the MHC columns and a specices column, since it only works for human
#This is only needed for the MixTCRpred1.0 input.
find_mhc <- function(m){
  mhc <- unlist(lapply(1:dim(m)[1], function(x){
    if(m[x,"MHC"]=="HLA-A:01") {h <- "A0101"}
    else{
      if(m[x,"species"] == "HomoSapiens"){
        a <- unlist(strsplit(m[x,"MHC"], split="-", fixed=T))
        h <- gsub('*', '', a[2], fixed=T)
        h <- gsub(':', '', h, fixed=T)
        h <- gsub("/", "_", h, fixed=T)
      } else {
        h <- m[x,"MHC"]
      }
    }
    return(h)
  }))
}

#' Plot the comparison between V or J usage.
#'
#' We should have different criteria for selecting the labels for the ES versus Pred.
#'
#' @param pType Switch telling the type of plot to use. See plot.VJ.switch
#'     option from MixTCRviz.
#' @param species Tells the species, which is used for gene colors in the plots.
#' @param ret.resList when T indicates to return a list of the results instead of
#'    the plots directly (in this case, the "info" should have a 4th element,
#'   indicating the model name).
#' @param combined.resList When this isn't NULL, we'll use the results from
#'    this list to plot the results (from multiple models combined together).
#' @param label.neg (default=F): If T, show also the labels of the genes most depleted in input1
#' @param label.diag (default=0.3): Print label on the diagonal above a certain value for both x and y axis
#' @param label.min.fr (default=c(0.05, 0.05)): Region (rectangle) of the left corner of V/J plots with no gene label
#' @export
plotVJ <- function(count.es, count.rep, sd.es=NULL, sd.rep=NULL, distr.es=NULL, distr.rep=NULL, info=NULL, comp.baseline=T, pType=1,
                   species="HomoSapiens", ret.resList=F, combined.resList=NULL, label.neg=F,
                   ZscoreVJ.thresh=0, FoldChangeVJ.thresh=1.25,
                   label.diag=0.3, label.min.fr=c(0.05, 0.05), print.size=T, plot.sd=T, verbose=1){

  if (is.null(combined.resList)){
    if (length(count.es) == 0){
      # Can directly return an empty plot when this doesn't contain any data.
      if (ret.resList){
        return(NULL)
      } else {
        return(ggplot())
      }
    }
    #count.es is on the Y axis, count.rep on the X
    segment <- info["gene"]  # e.g., TRAV
    n.es <- sum(count.es)

    n.rep <- sum(count.rep)
    v <- c(names(count.rep), names(count.es))
    nm <- unique(v)
    cn <- c(info["input1.name"], info["baseline.name"])
    count <- matrix(0, nrow=length(nm), ncol=2)
    rownames(count) <- nm
    colnames(count) <- cn
    for(i in 1:length(count.es)){ count[names(count.es[i]),info["input1.name"]] <- count.es[i] }
    for(i in 1:length(count.rep)){ count[names(count.rep[i]),info["baseline.name"]] <- count.rep[i] }
    count <- scale(count, center=F, scale=colSums(count))
    count.df <- data.frame(count)
    colnames(count.df) <- c("Y","X")
    count.df$name <- rownames(count.df)
    count.df$gene <- gsub("\\*.*$", "", count.df$name)
    # gsub is used to possible remove the allele information from the 'name'.

    #Now create the sd to show as error bars
    if(!is.null(sd.es)){

      if(n.es>1.1){  #This means that sd.rep should have the same normalisation as count.rep
        sd.es <- sd.es/n.es
      }
      cn <- colnames(count.df)
      count.df <- cbind(count.df,0)
      colnames(count.df) <- c(cn,"SD_es")
      for(n in nm){
        if(!is.na(sd.es[n])){
          count.df[n,"SD_es"] <- sd.es[n]
        }
      }

    }
    #Now create the sd to show as error bars
    if(!is.null(sd.rep)){

      if(n.rep>1.1){  #This means that sd.rep should have the same normalisation as count.rep
        sd.rep <- sd.rep/n.rep
      }
      cn <- colnames(count.df)
      count.df <- cbind(count.df,0)
      colnames(count.df) <- c(cn,"SD_rep")
      for(n in nm){
        if(!is.na(sd.rep[n])){
          count.df[n,"SD_rep"] <- sd.rep[n]
        }
      }
    }


    if(is.null(sd.es) | !plot.sd){
      lim.y <- max(count.df[,"Y"] )*1.4
    } else {
      lim.y <- max(max(count.df[,"Y"] )*1.4, count.df[,"Y"] + count.df[,"SD_es"])
    }
    if(is.null(sd.rep) | !plot.sd){
      lim.x <- max(count.df[,"X"] )*1.3
    } else {
      lim.x <- max(max(count.df[,"X"] )*1.3, count.df[,"X"] + count.df[,"SD_rep"])
    }


    ratio <- (count.df[,"Y"]+0.0001)/(count.df[,"X"]+0.0001)
    logFC <- log2(ratio);
    names(logFC) <- nm
    count.df$logFC <- as.character(signif(logFC,2))


    if(!is.null(sd.rep)){
      Zscore <- sapply(nm, function(n){
        if(count.df[n,"SD_rep"]>0){
          Z <- (count.df[n,"Y"] - count.df[n,"X"])/count.df[n,"SD_rep"]
        } else {
          #If there is no sd in baseline
          if(abs(count.df[n,"Y"] - count.df[n,"X"])>0.0001){
            Z <- 5 #This a max value fixed by default
          } else {
            #This is when there is no sd, but es and rep have the same value (typically 0 everywhere)
            Z <- 0
          }
        }
        return(Z)
      })
    } else {
      Zscore <- setNames(rep(NA, length(nm)), nm)
    }
    count.df$Zscore <- as.character(signif(Zscore,2))

    #print(count.df)

    #Compute an experimental P-value (i.e., %rank by comparing with repertoire)
    #This is currently not used, since it would require many more samples to estimated baseline variability
    if(!is.null(distr.rep)){
      Pval <- sapply(nm, function(n){
        if(n %in% rownames(distr.rep)){
          Ppos <- length(which(distr.rep[n,] >= count.df[n,"Y"]))/length(distr.rep[n,])
          Pneg <- length(which(distr.rep[n,] <= count.df[n,"Y"]))/length(distr.rep[n,])
          P <- max(1/dim(distr.rep)[2], min(Ppos, Pneg))
        } else {
          if(count.df[n,"Y"]>0){
            P <- 1/dim(distr.rep)[2]   # Min value, corresponding to cases where the gene was never seen in baseline
          } else {
            P <- 1
          }
        }
        return(P)
      })
    }

    print.example <- F
    if(print.example){
      gn <- "TRBJ2-7"
      if(gn %in% nm){
        print(count.df[gn,])
        if(!is.null(sd.rep)){
          print(Zscore[gn])
        }
        if(!is.null(distr.rep)){
          print(Pval[gn])
        }
      }

    }
    #Now decide which label to show.

    label <- nm; names(label) <- nm
    label.all <- label

    #If label.neg==F, put labels of all genes with negative logFC to NA
    if(!label.neg){
      label[logFC < 0] <- NA
    }

    #Hide labels for points in the lower left corner, based on label.min.fr criteria
    label[ count.df[,"Y"] < label.min.fr[1] & count.df[,"X"] < label.min.fr[2] ] <- NA

    #If Z-scores can be computed, hide labels for cases with Z-scores smaller than ZscoreVJ.thresh (default=0, nothing hidden)
    if(!is.null(sd.rep)>0){
      label[abs(Zscore) < ZscoreVJ.thresh] <- NA
    }

    n_lab_max <- ifelse(floor(pType)==2, 12, 8)

    #Hide labels for points with low fold change and not very high frequencies
    #Do this iteratively
    min.logFC <- log2(c(FoldChangeVJ.thresh,1.5,2,3))
    min.fr <- c(label.diag, 0.2, 0.1, 0.1) #This means that genes with higher frequency than min.fr in input1 will always be labelled, irrespective of their logFC

    #Test different stringency on the logFC thresholds (bar plot can
    #accomodate more labels before it gets too confusing visually).

    label[which( (abs(logFC) < min.logFC[1]) & count.df[,"X"] < min.fr[1] & count.df[,"Y"] < min.fr[1])] <- NA
    t <- 2
    while(length(which(!is.na(label))) > n_lab_max & t <= length(min.fr)){
      label[which( (abs(logFC) < min.logFC[t]) & count.df[,"X"] < min.fr[t] & count.df[,"Y"] < min.fr[t])] <- NA
      t <- t+1
    }


    #If we still have too many labels, take those with the top logFC


    if(length(which(!is.na(label))) > n_lab_max){
      if(verbose>=1){
        print(paste("Too many labels for ",segment,", selection will be based primarily on logFC:",sep=""))
        print("  Consider augmenting values in label.min.fr if visualisation is not good")
      }
      logFC.sort <- sort(abs(logFC))
      t <- 1

      while(length(which(!is.na(label))) > n_lab_max){
        label[names(logFC.sort[t])] <- NA
        t <- t+1
        #print(t)
      }
    }
    if(label.diag){
      #Force to show points along the diagonal (this can be useful when comparing two similar repertoires, with some enriched V or J)
      ind <- which(count.df[,"Y"] > label.diag & count.df[,"X"] > label.diag & abs(logFC) < 1.5)
      if(length(ind)<10){
        label[ind] <- label.all[ind]
      } else {
        print("Warning: the number of labels to be shown along the diagonal is very large. They will not be shown. Consider increasing label.diag")
      }
    }
    if(print.size){
      ylab <- paste(info["input1.name"]," (",n.es,")", sep="")
    } else {
      ylab <- info["input1.name"]
    }
    if(comp.baseline | !print.size){
      xlab <- info["baseline.name"]
    } else {
      xlab <- paste(info["baseline.name"]," (",n.rep,")", sep="")
    }
    if (floor(pType) == 1){
      count.df$label <- label
    } else {
      count.df$label <- count.df$name
    }
    # When using scatter plots, the labels not to show have NAs, while for the
    # bar plots, we keep all labels here but will subset count.df below to only
    # keep genes of interest (we keep all labels as needed for the
    # combined.resList case).

    # # If we don't want to show the TRAV/TRAJ/... info in the label, we can
    # # uncomment below line:
    # count.df$label <- gsub(segment, "", count.df$label)
    # # The 'segment' is just TRAV, TRAJ, ..., while count.df$name and
    # # count.df$gene are TRAV5-4 for example (with allele name possibly present
    # # in $name). We could thus remove the TRAV, ... info from the label to
    # # only keep the digits/code following these names.

    # And add/define some information needed for the bar plot.
    namesToKeep <- setdiff(label, NA)
    count.df$baselineName <- xlab
    # 'baselineName' will be used to show the baseline value from the genes in the
    # bar plot.
    count.df$log2FC <- log2(ratio+1e-5)
    count.df$model <- paste0(info["model"], " (", n.es, ")")

    if (ret.resList){
      return(list(count.df=count.df, namesToKeep=namesToKeep, segment=segment))
    }
  } else {
    count.df <- combined.resList$count.df
    count.df$model <- factor(count.df$model, levels=rev(unique(count.df$model)))
    # Make it as a factor so that the order in which the models were obtained
    # is kept.
    namesToKeep <- combined.resList$namesToKeep
    segment <- combined.resList$segment
  }

  if (pType == 1){
    colorScale <- TCRgene2aes[[species]][[segment]]$color1
  } else {
    colorScale <- TCRgene2aes[[species]][[segment]]$color2
  }

  #Plot the comparison between input and repertoires
  if (floor(pType) == 1){
    if (!is.null(combined.resList)){
      stop("Didn't implement the use of combined.ResList when floor(pType) == 1")
    }

    if(info["gene"] == "TRAV" | info["gene"] == "TRBV" ){
      CDR1_seq <- gsub("g","-",cdr123[[species]][[substr(info["gene"],1,3)]][count.df$name,"CDR1"])
      CDR2_seq <- gsub("g","-",cdr123[[species]][[substr(info["gene"],1,3)]][count.df$name,"CDR2"])
      CDR3_seq <- cdr123[[species]][[substr(info["gene"],1,3)]][count.df$name,"CDR3"]
      count.df$CDR1_seq <- CDR1_seq
      count.df$CDR2_seq <- CDR2_seq
      count.df$CDR3_seq <- CDR3_seq
      count.plot <- ggplot(count.df, aes(x=X, y=Y, label=label, label2=name, label3=logFC, label4=Zscore, label5=CDR1_seq, label6=CDR2_seq, label7=CDR3_seq)) +
        geom_abline(col="orange",linetype="dashed",linewidth=1)

    } else if(info["gene"] == "TRAJ" | info["gene"] == "TRBJ" ){
      CDR3_seq <- Jseq[[species]][[substr(info["gene"],1,3)]][count.df$name,"CDR3"]
      count.df$CDR3_seq <- CDR3_seq
      count.plot <- ggplot(count.df, aes(x=X, y=Y, label=label, label2=name, label3=logFC, label4=Zscore, label5=CDR3_seq)) +
        geom_abline(col="orange",linetype="dashed",linewidth=1)
    }





    if(!is.null(sd.es) & plot.sd){
      count.plot <- count.plot + geom_errorbar(aes(ymax=Y+SD_es, ymin=sapply(Y-SD_es, function(x){max(0.001,x)})), width=0.015*lim.x, linewidth=0.4, color="grey50")
    }
    if(!is.null(sd.rep) & plot.sd){
      count.plot <- count.plot + geom_errorbarh(aes(xmax=X+SD_rep, xmin=sapply(X-SD_rep, function(x){max(0.001,x)})), height=0.015*lim.y, linewidth=0.4, color="grey50")
      #count.plot <- count.plot + geom_errorbar(aes(xmax=X+SD_rep, xmin=sapply(X-SD_rep, function(x){max(0.001,x)})), orientation = "y", width=0.015*lim.y, linewidth=0.4, color="grey50")
    }

    if (pType == 1.3){
      count.plot <- count.plot + geom_point()
    } else {
      # Define parameters to have a thin line around the points in scatter plot.
      if (pType == 1.2){
        shape_color <- "common"
        shapeScale <- 21
        outerColorScale <- "gray20"
      } else if (pType == 1){
        shape_color <- count.df$gene
        shapeScale <- TCRgene2aes[[species]][[segment]]$shape1
        outerColorScale <- TCRgene2aes[[species]][[segment]]$outerColor1
      }
      outerWidth <- 0.5

      count.plot <- count.plot + geom_point(aes(fill=gene, shape=shape_color,
                                                color=shape_color), stroke=outerWidth, size=2.5) +
        scale_color_manual(values=outerColorScale, guide="none") +
        scale_shape_manual(values=shapeScale, guide="none")
    }
    count.plot <- count.plot +
      ggtitle(segment) +
      xlim(0, lim.x) + ylim(0,lim.y) + theme_bw() +
      theme(plot.title = element_text(size = 14, hjust=0.5),
            axis.text=element_text(size=10), axis.title=element_text(size=14)) +
      geom_label_repel(aes(fill=gene), size = 3, nudge_y=0.02, box.padding = 0.25,
                       show.legend=F, na.rm=T) +
      scale_fill_manual(values=colorScale, guide="none") +
      xlab(xlab) + ylab(ylab) + theme(panel.grid.minor = element_blank())

  } else {
    # Show results as bar plots. Will only keep most significant genes, possibly
    # summing together all the other and rework a bit the data.
    count_other <- count.df[!count.df$name %in% namesToKeep,,drop=F] %>%
      dplyr::summarize(Y=sum(Y), X=0, name="Other", gene="Other",
        label="Other", log2FC=NA, baselineName=count.df$baselineName[1],
        .by=model)
    # All baselineName should be the same value (and we put a value of 0
    # for the 'baseline' of these other, so that no baseline value is showed
    # there).

    count.df <- count.df[count.df$name %in% namesToKeep,,drop=F]
    count.df <- count.df[order(count.df$log2FC, decreasing = T),,drop=F]
    # Order genes based on the log2FC between input and baseline to show
    # most important ones on top.
    if (pType == 2.1){
      # Adding the summed frequency of all remaining genes if pType was 2.1 only.
      count.df <- dplyr::bind_rows(count.df, count_other)
    }

    count.df$label <- factor(count.df$label, levels=rev(unique(count.df$label)))
    # Use a factor with given levels to keep the order constructed above
    # (otherwise it'd order those labels alphabetically) - we use rev to make
    # the most important genes on top of the plot.
    labYlen <- max(nchar(levels(count.df$label)))
    YaxSize <- ifelse(labYlen < 8, 14, ifelse(labYlen < 10, 13,
      ifelse(labYlen < 12, 12, ifelse(labYlen < 14, 11, 10))))
    # Choose different sizes for the y-axis labels, a trade-off between reading
    # these labels well and seeing well the values (because otherwise the
    # long labels take up all the space from the plot and we cannot see the
    # bars).
    legLen <- max(nchar(levels(count.df$model)))
    legSize <- ifelse(legLen < 15, 12, ifelse(legLen < 20, 10,
      ifelse(legLen < 25, 8, 6)))
    # And also for the legend indicating the model names.

    if (is.null(combined.resList)){
      count.plot <- ggplot(count.df, aes(x=Y, y=label, fill=gene)) +
        geom_col(color="gray10", linewidth=1)

      figTitle <- paste0(segment, " (", n.es, ")")
    } else {
      count.plot <- ggplot(count.df, aes(x=Y, y=label, fill=gene, color=model)) +
        geom_col(position="dodge", linewidth=1, width=0.8) +
        scale_color_manual(values=set_model_colPals(rev(levels(count.df$model))),
                           guide=guide_legend(ncol=1, order=1, reverse=T))
      figTitle <- segment
    }

    count.plot <- count.plot +
      scale_fill_manual(values=colorScale, guide="none") +
      ggtitle(figTitle) + xlab("Frequency") + ylab(NULL) +
      scale_x_continuous(expand=expansion(mult=c(0, 0.05)), n.breaks=3) +
      # Make the x-axis isn't expanded on the left and is expanded as usual on
      # the right.
      geom_col(aes(x=X, alpha=baselineName),
               position="dodge", width=ifelse(is.null(combined.resList), 0.9, 0.8),
               linewidth=0.5, color="gray60", fill="gray80") +
      scale_alpha_manual(values=0.7, guide=guide_legend(order=2)) +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, hjust=0.5),
        axis.title=element_text(size=14),
        panel.grid.major.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=YaxSize, face="bold", color="black"),
        legend.position="top", legend.title=element_blank(),
        legend.text=element_text(size=legSize))
  }

  return(count.plot)
}


# ret.resList when T indicates to return a list of the results instead of
# the plots directly (in this case, the "info" should have a 4th element,
# indicating the model name).
# And when combined.resList isn't NULL, we'll use the results from this list
# to plot the results (from multiple models combined together).
plotLD <- function(countL.es, countL.rep, info=NULL, sd.es=NULL, sd.rep=NULL, plot.oneline=0, ret.resList=F,
                   combined.resList=NULL, comp.baseline=T, print.size=T, plot.sd=T){

  if (is.null(combined.resList)){

    L.all <- Lmin:Lmax
    cn <- paste("L",L.all,sep="_")

    ld.es <- rep(0,length(L.all)); names(ld.es) <- cn
    ld.rep <- rep(0,length(L.all)); names(ld.rep) <- cn
    n.es <- sum(countL.es)
    n.rep <- sum(countL.rep)

    for(lc in cn){
      if(!is.na(countL.es[lc])){ ld.es[lc] <- countL.es[lc]}
      if(!is.na(countL.rep[lc])){ld.rep[lc] <- countL.rep[lc]}
    }
    if(sum(ld.es)>0){
      ld.es <- ld.es/sum(ld.es)
    }
    if(sum(ld.rep)>0){
      ld.rep <- ld.rep/sum(ld.rep)
    }

    lds.es <- rep(0,length(L.all)); names(lds.es) <- cn
    lds.rep <- rep(0,length(L.all)); names(lds.rep) <- cn
    #Now create the sd to show as error bars
    if(!is.null(sd.es)){
      if(n.es>1.1){  #This implies that sd.rep should have the same normalisation as count
        sd.es <- sd.es/n.es
      }
      for(lc in cn){
        if(!is.na(sd.es[lc])){ lds.es[lc] <- sd.es[lc]}
      }
    }
    #Now create the sd to show as error bars
    if(!is.null(sd.rep)){
      if(n.rep>1.1){  #This means that sd.rep should have the same normalisation as count.rep
        sd.rep <- sd.rep/n.rep
      }
      for(lc in cn){
        if(!is.na(sd.rep[lc])){ lds.rep[lc] <- sd.rep[lc]}
      }
    }

    #Plot the comparison for length distribution

    v1 <- c(L.all,L.all);
    v2 <- c(ld.es, ld.rep);
    v3 <- c( rep(info["input1.name"], length(L.all)), rep(info["baseline.name"], length(L.all))) ;
    v4 <- c(lds.es,lds.rep)
    ld.df <- data.frame(v1,v2,v3,SD=v4)
    ld.df$v3 <- factor(ld.df$v3, levels=c(info["input1.name"], info["baseline.name"]))

    if (ret.resList){
      if (length(info) < 4){
        stop("The 'info' vector given in plotLD input should have 4 elements ",
             "when ret.resList is T.")
      }
      ld.df$model <- paste0(info["model"], gsub(".*( \\(\\d+\\))", "\\1", info["input1.name"]))
      levels(ld.df$v3) <- gsub("(.*) \\(\\d+\\)", "\\1", levels(ld.df$v3))
      # The number of sequences is inputted in info["input1.name"] (saved in v3 column of
      # lf.df). We instead indicate this information in the model column when
      # combining results from multiple models.
      return(list(ld.df=ld.df, info=info["chain"]))
      # Return the data.frame, as well as info (but only 1st element as the other
      # are already encode in ld.df).
    }

    legend.size <- 12
    if(plot.oneline!=0){
      if(nchar(info["input1.name"])>23){legend.size=11}
      if(nchar(info["input1.name"])>25){legend.size=10}
    }

    ld.plot <-  ggplot(ld.df, aes(x=v1, y=v2, color=v3, shape=v3)) +
      guides(color = guide_legend(ncol = 1, order=1), shape=guide_legend(ncol = 1, order=1))

    if(plot.sd & (!is.null(sd.es) | !is.null(sd.rep))){
      ld.plot <- ld.plot +
        geom_errorbar(aes(ymax=v2+SD, ymin=sapply(v2-SD, function(x){max(0.001,x)})), width=0.3, linewidth=0.4)
    }

    # The rest of the plot is the same if combined.resList was NULL or if showing
    # the results from multiple models combined, so we'll draw it below.
  } else {
    ld.df <- combined.resList$ld.df
    ld.df$model <- factor(ld.df$model, levels=unique(ld.df$model))
    info <- combined.resList$info
    legend.size <- 10
    ld.plot <-  ggplot(ld.df, aes(x=v1, y=v2, color=model, linetype=v3)) +
      facet_grid(rows=vars(model), scales="free_y") +
      scale_color_manual(values=set_model_colPals(levels(ld.df$model))) +
      guides(linetype=guide_legend(nrow=1, order=2), color="none")
  }

  if(plot.oneline==0){
    size <- 2.5
  } else {
    size <- 2
  }

  ld.plot <- ld.plot + geom_point(size=size) + geom_line() + theme_bw() +
    theme(legend.key.size = unit(0.65, 'cm'), legend.position="top",
          legend.title=element_blank(), legend.text=element_text(size=legend.size)) +
    xlab(paste("Length_CDR3",info["chain"],sep="")) + ylab("") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14),
          plot.title = element_text(size=15,hjust = 0.5)) + theme(panel.grid.minor = element_blank())



  return(ld.plot)
}


plotCDR3 <- function(countL.es, countL.rep, countCDR3.es, countCDR3.rep, info=NULL,
                     comp.baseline=T, plot.oneline=0, plot.all.length=F, logo.type = "bits",
                     plot.cdr3.subtract.baseline=0, set.cdr3.length=NA, print.size=T){

  L.es <- as.numeric(lapply(names(countL.es), function(x){unlist(strsplit(x,split="_"))[2]}))
  L.rep <- as.numeric(lapply(names(countL.rep), function(x){unlist(strsplit(x,split="_"))[2]}))

  L.TR <- intersect(L.es,L.rep)

  pwm.rep <- list()
  pwm.es <- list()

  logo.CDR3.L.es <- list()
  logo.CDR3.L.rep <- list()

  if(length(L.TR)>0){

    tl <- countL.es[paste("L",L.TR,sep="_")]

    if(is.na(set.cdr3.length)){
      lmax <- as.numeric(unlist(strsplit(names(tl[which.max(tl)]), split="_"))[2])
    } else {
      if(set.cdr3.length %in% L.TR){
        lmax <- set.cdr3.length
      } else {
        lmax <- as.numeric(unlist(strsplit(names(tl[which.max(tl)]), split="_"))[2])
        print(paste("set.cdr3",info["chain"],".length=",set.cdr3.length," is incompatible with the input data. Default value of ",lmax," will be used.", sep=""))
      }
    }

    if(plot.oneline!=0){
      if(lmax<15){
        axis.size.max <- 8
      } else {
        axis.size.max <- 7
      }
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

      lc <- paste("L",l,sep="_")

      pwm.es[[lc]] <- scale(countCDR3.es[[lc]], center=F, scale=colSums(countCDR3.es[[lc]]))
      pwm.rep[[lc]] <- scale(countCDR3.rep[[lc]], center=F, scale=colSums(countCDR3.rep[[lc]]))

      if(plot.cdr3.subtract.baseline==1){

        #Compute the logo representing the difference in frequencies renormalized by information content
        #This is currently not supported by ggseqlogoMOD...
        #Compute the matrix representing the size of each letter
        size.es <- apply(pwm.es[[lc]], 2, function(x){ ind <- which(x!=0); IC <- log(N.aa)/log(2)+sum(x[ind]*log(x[ind])/log(2)); return(IC*x) })
        size.rep <- apply(pwm.rep[[lc]], 2, function(x){ ind <- which(x!=0); IC <- log(N.aa)/log(2)+sum(x[ind]*log(x[ind])/log(2)); return(IC*x) })
        x.norm <- size.es-size.rep
        y.inc <- 1
      }
      if(plot.cdr3.subtract.baseline==2){

        #Compute the logo based on normalised fold-change
        rc <- max(5,0.1*countL.es[[lc]])
        #rc <- 5
        pseudo <- 1/N.aa*rc/countL.es[[lc]]
        x.baseline <- pwm.rep[[lc]]+pseudo
        x.es <- pwm.es[[lc]]+pseudo
        x.norm <- x.es/x.baseline
        x.norm <- scale(x.norm,center = F, scale=colSums(x.norm))
        y.inc <- 4
      }
      title <- info["input1.name"]
      if(print.size){ title <- paste(title, " (",countL.es[[lc]],")", sep="")  }
      title <- paste(title,", CDR3", info["chain"],"_",l, sep="")

      logo.CDR3.L.es[[lc]] <- ggseqlogoMOD::ggseqlogoMOD(data=pwm.es[[lc]], additionaAA=additionalAA,  axisTextSizeX = 12, axisTextSizeY = 8, methods = logo.type) +
        labs(title=title) + ylab(ylab) + theme(plot.title=element_text(size=15, hjust=0.5))

      title.baseline <- info["baseline.name"]
      if(!comp.baseline & print.size){title.baseline <- paste(title.baseline, " (",countL.rep[[lc]],")", sep="")}

      if(plot.cdr3.subtract.baseline==0){
        title.baseline <- paste(title.baseline,", CDR3", info["chain"],"_",l, sep="")
      } else if(plot.cdr3.subtract.baseline==1){
        title.baseline <- paste(title.baseline," subtract, CDR3", info["chain"],"_",l, sep="")
      } else if(plot.cdr3.subtract.baseline==2){
        title.baseline <- paste(title.baseline," renorm, CDR3", info["chain"],"_",l, sep="")
      }

      if(plot.cdr3.subtract.baseline==0){
        logo.CDR3.L.rep[[lc]] <- ggseqlogoMOD::ggseqlogoMOD(data=pwm.rep[[lc]], additionaAA=additionalAA,  axisTextSizeX = 12, axisTextSizeY = 8, methods = logo.type) +
          labs(title=title.baseline) + ylab(ylab) + theme(plot.title=element_text(size=15, hjust=0.5))
      } else if(plot.cdr3.subtract.baseline==1){
        y.min <- min(apply(x.norm, 2, function(x){ sum(x[x<0]) }))
        y.max <- max(apply(x.norm, 2, function(x){ sum(x[x>0]) }))
        y.min <- max(-log(N.aa)/log(2), y.inc*y.min)
        y.max <- log(N.aa)/log(2) # min(log(N.aa)/log(2), y.inc*y.max)
        logo.CDR3.L.rep[[lc]] <- ggseqlogoMOD::ggseqlogo(data=x.norm, method='custom') +
          labs(title=title.baseline) + ylim(y.min,y.max) + ylab(ylab) +
          theme(plot.title=element_text(size=title.size, hjust=0.5)) + theme(legend.position = 'none')
      } else if(plot.cdr3.subtract.baseline==2){
        IC.max <- max(unlist(apply(x.norm, 2, function(x){ ind <- which(x!=0); IC <- log(N.aa)/log(2)+sum(x[ind]*log(x[ind])/log(2)); return(IC) })))
        y.max <- min(IC.max*y.inc, log(N.aa)/log(2))
        logo.CDR3.L.rep[[lc]] <- ggseqlogoMOD::ggseqlogoMOD(data=x.norm, additionaAA=additionalAA,  axisTextSizeX = 12, axisTextSizeY = 8, ylim=c(0, y.max), methods = logo.type) +
          labs(title=title.baseline) + ylab(ylab) + theme(plot.title=element_text(size=15, hjust=0.5))
      }
      #For the special case where l==lmax, build the logo with different graphical parameters, depending on the plot.oneline
      #So far, we redo everything, since the graphical outline has to be a little bit different,
      # but this is not optimal since any change has to be performed multiple times
      if(l==lmax){
        title <- info["input1.name"]
        if(print.size){ title <- paste(title, " (",countL.es[[lc]],")", sep="")  }
        title <- paste(title,", CDR3", info["chain"],"_",l, sep="")

        title.baseline <- info["baseline.name"]
        if(!comp.baseline & print.size){title.baseline <- paste(title.baseline, " (",countL.rep[[lc]],")", sep="")}

        if(plot.cdr3.subtract.baseline==0){
          title.baseline <- paste(title.baseline,", CDR3", info["chain"],"_",l, sep="")
        } else if(plot.cdr3.subtract.baseline==1){
          title.baseline <- paste(title.baseline," subtract, CDR3", info["chain"],"_",l, sep="")
        } else if(plot.cdr3.subtract.baseline==2){
          title.baseline <- paste(title.baseline," renorm, CDR3", info["chain"],"_",l, sep="")
        }

        if(plot.oneline!=0 & (nchar(title)>26 | nchar(title.baseline)>26)){
          title <- paste("CDR3", info["chain"],"_",l," ",info["input1.name"], sep="")
          if(print.size){ title <- paste(title,"\n(",countL.es[[lc]],")", sep="")}
          if(plot.cdr3.subtract.baseline==0){
            title.baseline <- paste("CDR3", info["chain"],"_",l," ",info["baseline.name"],"\n", sep="")
          } else if(plot.cdr3.subtract.baseline==1){
            title.baseline <- paste("CDR3", info["chain"],"_",l," subtract \n",info["baseline.name"], sep="")
          } else if(plot.cdr3.subtract.baseline==2){
            title.baseline <- paste("CDR3", info["chain"],"_",l," renorm \n",info["baseline.name"], sep="")
          }
          if(!comp.baseline & print.size){
            title.baseline <- paste(title.baseline, "(",countL.rep[[lc]],")", sep="")
          }
        }

        logo.CDR3.L.es.max <- ggseqlogoMOD::ggseqlogoMOD(data=pwm.es[[lc]], additionaAA=additionalAA,  axisTextSizeX = axis.size.max, axisTextSizeY = 8, methods = logo.type) +
          labs(title=title) + ylab(ylab) + theme(plot.title=element_text(size=title.size, hjust=0.5))

        if(plot.cdr3.subtract.baseline==0){
          logo.CDR3.L.rep.max <- ggseqlogoMOD::ggseqlogoMOD(data=pwm.rep[[lc]], additionaAA=additionalAA,  axisTextSizeX = axis.size.max, axisTextSizeY = 8, methods = logo.type) +
            labs(title=title.baseline) + ylab(ylab) + theme(plot.title=element_text(size=title.size, hjust=0.5))
        } else if(plot.cdr3.subtract.baseline==1){
          y.min <- min(apply(x.norm, 2, function(x){ sum(x[x<0]) }))
          y.max <- max(apply(x.norm, 2, function(x){ sum(x[x>0]) }))
          y.min <- max(-log(N.aa)/log(2), 2*y.min)
          y.max <- log(N.aa)/log(2) # min(log(N.aa)/log(2), 2*y.max)
          logo.CDR3.L.rep.max <- ggseqlogoMOD::ggseqlogo(data=x.norm, method='custom') +
            ggtitle(title.baseline) + ylim(y.min,y.max) + ylab(ylab) +
            theme(plot.title=element_text(size=title.size, hjust=0.5)) + theme(legend.position = 'none')
        } else if (plot.cdr3.subtract.baseline==2){
          IC.max <- max(unlist(apply(x.norm, 2, function(x){ ind <- which(x!=0); IC <- log(N.aa)/log(2)+sum(x[ind]*log(x[ind])/log(2)); return(IC) })))
          y.max <- min(IC.max*y.inc, log(N.aa)/log(2))
          logo.CDR3.L.rep.max <- ggseqlogoMOD::ggseqlogoMOD(data=x.norm, additionaAA=additionalAA,  axisTextSizeX = axis.size.max, axisTextSizeY = 8, ylim=c(0, y.max), methods = logo.type) +
            labs(title=title.baseline) + ylab(ylab) + theme(plot.title=element_text(size=title.size, hjust=0.5))
        }
      }
    }
    ls <- list(logo.CDR3.L.es, logo.CDR3.L.rep, L.TR, lmax, logo.CDR3.L.es.max, logo.CDR3.L.rep.max)
    names(ls) <- c("ES", "Baseline", "length", "lmax", "ES_max", "Baseline_max")

  }  else {
    ls <- list(list(), list(), c(), 0, ggplot(), ggplot())
    names(ls) <- c("ES", "Baseline", "length", "lmax", "ES_max", "Baseline_max")

  }

  return(ls)
}


#' @export
check_input <- function(input, chain="AB", name="input1", species.default="HomoSapiens", infer.VJ=F,infer.CDR3=F,
                        model.default="Model_default", input.list=F, build.clones=F, verbose=1){
  #Check if some columns are missing, and add them with default values

  chain.list <- paste("TR",strsplit(chain,split="")[[1]], sep="")

  if (is.data.frame(input)){
    input <- as.data.frame(input)
    # Some code in this function doesn't handle tibbles that are a type of
    # data.frame (as.data.frame drops these tibble class).
  }

  if(is.data.frame(input) & !input.list){

    format <- determine.format(input)

    #Do some corrections specific for VDJdb
    if(format=="VDJdb"){
      if("Species" %in% colnames(input)){
        input$species <- input$Species
      }
      #If not model, use the MHC_Epitope, or only the Epitope
      if(length(intersect(colnames(input),"model"))==0){
        if("Epitope" %in% colnames(input)){
          print("WARNING: Epitopes will be used as models - this is risky if the same epitope is restricted to multiple MHCs.")
          print("  Consider adding a \'model\' column to your data.")
          input$model <- input$Epitope
        }
      }
    }

    #Handle the format with alpha and beta chains on different rows (e.g., with clone.id)
    if(format %in% names(clone.format.col)){
      if(!build.clones){
        input <- stack_clones(input, format)
      } else {
        input <- merge_clones(input, format)
      }
    }


    col <- as.character(sapply(chain.list, function(x){c(paste(x,"V",sep=""), paste(x,"J",sep=""),paste("cdr3_",x,sep=""))}))

    #Deal with cases where column name do not follow the default of MixTCRviz
    #Do different corrections for single chain (e.g., allow V or CDR3) and paired data (allow only Va or CDR3a).
    for(i in 1:length(colnames(input))){
      cl <- colnames(input)[i]
      #Check if the col.names is in the mapping
      if(cl %in% names(mapping.colnames[[chain]])){
        mp <- mapping.colnames[[chain]][cl]
        #Make sure another entry does not already have the corrected name and add the column
        if( !mp %in% colnames(input)){
          input[[mp]] <- input[[cl]]
        }
      }
    }

    #Check missing input
    for(cl in col){
      if(cl %in% colnames(input) == F){
        cn <- colnames(input)
        input <- cbind(input,"")
        colnames(input) <- c(cn, cl)
        if(cl %in% c("cdr3_TRA", "cdr3_TRB") & !infer.CDR3){
          print(paste("Missing",cl,"information in",name))
        }
        if(cl %in% c("TRAV", "TRAJ" ,"TRBV", "TRBJ") & !infer.VJ){
          print(paste("Missing",cl,"information in",name))
        }
      } else {
        input[,cl] <- as.character(input[,cl])
      }
    }
    if(infer.VJ | infer.CDR3){
      #This is the special case where the V and J and/or CDR3 are inferred from TCRa and TCRb columns
      col <- as.character(sapply(chain.list, function(x){paste0("TCR",tolower(substr(x,3,3)))}))
      for(cl in col){
        if(! cl %in% colnames(input)){
          stop(paste("Missing",cl,"information - this is incompatible with infer.VJ=TRUE or infer.CDR3=TRUE"))
        }
      }
    }
    #If the "species" column is not provided, we add a column with species.default
    #This is a bit suboptimal, but ok for now
    if("species" %in% colnames(input) == F){
      cn <- colnames(input)
      input <- cbind(input,species.default)
      colnames(input) <- c(cn, "species")
      if(verbose>0){
        print(paste("Using",species.default,"as species for all entries"))
      }
    }
    if("model" %in% colnames(input) == F){
      cn <- colnames(input)
      input <- cbind(input,model.default)
      colnames(input) <- c(cn, "model")
      if(verbose>0){
        print(paste("Using",model.default,"as model for all entries"))
      }
    }

    #Force all columns describing the TCR to be characters


  } else if(input.list) {
    nm.list <- c("L","countL","countV","countJ", "countVJ.L", "countCDR3.L", "countVJ",  "countV.L", "countJ.L") # the last 4 are actually not very useful
    for(nm in nm.list){
      if(is.null(input[[nm]])){
        stop(paste("Missing feature in input object:",nm,sep=" "))
      } else {
        for(ch in chain.list){
          if(is.null(input[[nm]][[ch]])){
            stop(paste("Missing feature in input object:",nm,ch,sep=" "))
          }
        }
      }
    }
    #Check the 'model' field.
    if(is.null(input$model)){
      input$model <- model.default
      print(paste("Using",model.default,"as model"))
    }
    if(is.null(input$species)){
      input$species <- species.default
      print(paste("Using",species.default,"as species"))
    }

  } else {
    stop("Issues with input1 format")
  }
  lst <- list()
  lst[["data"]] <- input
  return(lst)

}

#' @export
clean_input <- function(input, use.allele=F, correct.gene.names=T, use.mouse.strain=F,
                        chain="AB", species.default="HomoSapiens", check.cdr3.mode=1,
                        keep.incomplete.chain=T, start.lg=1, end.lg=2, seq.protocol="Default",
                        merge.ambiguous=T, verbose=1){


  # Clean the input by removing CDR3 with weird characters, longer than Lmax or shorter than Lmin
  # Correct VJ genes based on our dictionary
  # species.default is only used if input does not contain the "species" column
  # merge.ambiguous should be set to FALSE ONLY is clean_input is used outside of MixTCRviz
  # this will prevent mapping TRBV6-2 to TRBV6-2/6-3

  #print("Start")

  if(is.data.frame(input)){
    input <- as.data.frame(input)
    # Some code used here don't handle tibbles that are a type of data.frame
    # (as.data.frame drops these tibble class).
  }

  if("species" %in% colnames(input)){
    species.list <- unique(input[,"species"])
    use.species.default <- F
  } else {
    species.list <- c(species.default)
    use.species.default <- T
  }

  if(chain=="AB"){
    col <- c("TRAV","TRAJ","cdr3_TRA","TRBV","TRBJ","cdr3_TRB")
    segment.list <- c("TRAV","TRAJ","TRBV","TRBJ")
    cdr3.list <- c("cdr3_TRA","cdr3_TRB")
  }
  if(chain=="A"){
    col <- c("TRAV","TRAJ","cdr3_TRA")
    segment.list <- c("TRAV","TRAJ")
    cdr3.list <- c("cdr3_TRA")
  }
  if(chain=="B"){
    col <- c("TRBV","TRBJ","cdr3_TRB")
    segment.list <- c("TRBV","TRBJ")
    cdr3.list <- c("cdr3_TRB")
  }



  #Replace empty values by NA
  for(i in col){
    input[which(input[,i] == '' | input[,i] == ""),i] <- NA
  }

  #Set to NA CDR3 sequences with incompatible lengths or weird characters

  #print("Checking non-aa characters")
  for(cdr3 in cdr3.list){
    nc <- nchar(input[,cdr3])
    ind <- which( nc < Lmin | nc > Lmax | grepl('[^ACDEFGHIKLMNPQRSTVWY]', input[,cdr3]) == T)
    input[ind,cdr3] <- NA
  }

  #Remove anything that comes after parenthesis (this is the case for MiXCR data for instance)
  for(s in segment.list){
    ind <- grep("\\(",input[,s])
    input[ind,s] <- sapply(ind,function(i){strsplit(input[i,s], split="\\(")[[1]][1]})
  }

  #If multiple V or J genes (separated by "," or "\" or ";" or " or "), keep only the first one
  for(s in segment.list){
    ind <- grep("\\\\|,|;| or ",input[,s])
    if(length(ind)>0){
      cor <- sapply(ind,function(i){strsplit(input[i,s], split="\\\\|,|;| or ")[[1]][1]})
      if(verbose>=1) {
        entry <- ifelse(length(ind)==1,"entry","entries")
        print(paste(length(ind)," ",entry," with multiple ",s," segments (only the first segment will be kept):", sep=""))

        if(verbose==1 | verbose==2){
          print("Use verbose=3 to see them all")
        } else {
          mt <- cbind(input[ind,s], cor)
          colnames(mt) <- c("Original", "Corrected")
          print(mt)
        }
      }
      input[ind,s] <- cor
    }
  }

  #Remove spaces
  for(s in segment.list){
    input[,s] <- gsub(" ","",input[,s])
  }

  #Remove or add alleles
  if(!use.allele){
    #print("Removing alleles")
    for(s in segment.list){
      ind <- grep("*",input[,s], fixed=T)
      input[ind,s] <- gsub("\\*0[0-9]/0[0-9]", "", input[ind,s]) #This are entries with ambiguous allele assignment
      input[ind,s] <- gsub("\\*0[0-9]", "", input[ind,s])
    }
  } else{
    #print("Adding alleles")
    for(s in segment.list){
      ind <- which(!grepl("*",input[,s], fixed=T) & !is.na(input[,s]) )
      if(use.species.default){
        al <- allele.default[[seq.protocol]][[species.default]][input[ind,s]]
      } else {
        al <- sapply(ind, function(i){allele.default[[seq.protocol]][[input[i,"species"]]][input[i,s]]})
      }
      al[is.na(al)] <- "01" #This happens in case of wrong gene names, since gene names were not yet corrected
      input[ind,s] <- paste(input[ind,s], al, sep="*")
    }
  }

  if(merge.ambiguous){
    #Do a manual correction for TRBV6-2 and TRBV6-3 -> TRBV6-2/6-3
    if(chain=="B" | chain=="AB"){
      if(seq.protocol=="Default" | seq.protocol=="SEQTR"){
        ind <- which((grepl("TRBV6-2",input[,"TRBV"], fixed=T) | grepl("TRBV6-02",input[,"TRBV"], fixed=T) | grepl("TRBV6-3",input[,"TRBV"], fixed=T) | grepl("TRBV6-03",input[,"TRBV"], fixed=T))
                     & input[,"TRBV"] != "TRBV6-2/6-3" & input[,"TRBV"] != "TRBV6-2/6-3*01")
        if(length(ind)>0){
          if(verbose>0){
            print("Mapping all TRBV6-2 and TRBV6-3 to TRBV6-2/6-3, since they cannot be distinguished at the sequencing level")
          }
          if(use.allele){
            input[ind,"TRBV"] <- "TRBV6-2/6-3*01"
          } else {
            input[ind,"TRBV"] <- "TRBV6-2/6-3"
          }
        }
      }
    }

    #Do a manual correction for TRBV12-3 or TRBV12-4 -> TRBV12-3/12-4

    if(chain=="B" | chain=="AB"){
      if(seq.protocol=="SEQTR"){
        ind <- which((grepl("TRBV12-3",input[,"TRBV"], fixed=T) | grepl("TRBV12-03",input[,"TRBV"], fixed=T) | grepl("TRBV12-4",input[,"TRBV"], fixed=T) | grepl("TRBV12-04",input[,"TRBV"], fixed=T)) &
                       input[,"TRBV"] != "TRBV12-3/12-4" & input[,"TRBV"] != "TRBV12-3/12-4*01")

        if(length(ind)>0){
          if(verbose>0){
            print("Mapping all TRBV12-3 or TRBV12-4 to TRBV12-3/12-4, since they cannot be distinguished with SEQTR protocol")
          }
          if(use.allele){
            input[ind,"TRBV"] <- "TRBV12-3/12-4*01"
          } else {
            input[ind,"TRBV"] <- "TRBV12-3/12-4"
          }
        }
      } else {
        ind <- grep("TRBV12-3/12-4",input[,"TRBV"])
        if(length(ind)>0){
          if(verbose>0){
            print("*** WARNING: TCRs contain TRBV12-3/12-4 entries.")
            print("    If you are using data generated with SEQTR protocol, make sure to specify it with seq.protocol=\"SEQTR\"")
          }
        }
      }
    }
  }

  # Correct gene names
  # If alleles, it will correct the gene name, and keep the allele. If the allele cannot be found, it will remove it
  # If genes, it will correct the gene name
  # If gene name cannot be corrected, it gives NA

  if(correct.gene.names){
    #print("Check V/J names")
    input <- correct.VJnames(input=input, segment.list=segment.list, species.default=species.default,
                             use.allele=use.allele, seq.protocol=seq.protocol,verbose=verbose)
  } else {
    for(species in species.list){

      if(!use.species.default){
        ind.species <- which(input[,"species"]==species)
      } else {
        ind.species <- 1:dim(input)[1]
      }
      for(s in segment.list){
        if(use.allele){
          name.list <- gene.allele.list[[seq.protocol]][[species]][substr(gene.allele.list[[seq.protocol]][[species]],1,4)==s]
        } else {
          name.list <- gene.list[[seq.protocol]][[species]][substr(gene.list[[seq.protocol]][[species]],1,4)==s]
        }
        ind <- which(input[ind.species,s] %in% name.list == F & !is.na(input[ind.species,s]) )
        if(length(ind)>=1 & verbose != 0){
          missing <- sort(unique(input[ind.species[ind],s]))
          nm <- ifelse(length(missing)==1,"name","names")
          cat("\n")
          print(paste("*** ",length(missing)," ",s," ",nm," in ",length(ind)," entries absent from IMGT in ",species," ***",sep=""))
          print(missing)
          cat("\n")
        }
        input[ind.species[ind],s] <- NA
      }
    }
  }

  if(check.cdr3.mode > 0){
    #print("Checking CDR3")
    input <- check_cdr3(input=input, chain=chain, species.default=species.default, check.cdr3.mode=check.cdr3.mode, start.lg=start.lg, end.lg=end.lg, verbose=verbose)
  }


  # Do an extra correction for mouse entries, where only gene level analyses are allowed
  # and TRAV genes can be merged


  if(!use.species.default){
    ind <- which(input[,"species"]=="MusMusculus")
  } else {
    if(species.default=="MusMusculus"){
      ind <- 1:dim(input)[1]
    } else {
      ind <- c()
    }
  }
  if(length(ind)>0){

    if(use.allele){
      #Remove the alleles (if(use.allele==F), this was done before)
      for(s in segment.list){
        input[ind,s] <- unlist(lapply(input[ind,s], function(x){unlist(strsplit(x,split="*", fixed=T))[1]}))
      }
    }

    if(!use.mouse.strain & chain != "B"){
      input[ind,] <- merge_mouse_TRAV(input[ind,])  #WARNING: This only works if alleles have been removed (so far always the case in mouse)
    }
  }

  if(!keep.incomplete.chain){
    chain.list <- paste("TR",strsplit(chain,split="")[[1]], sep="")
    for(ch in chain.list){
      cl <- c(paste(ch,"V",sep=""),paste(ch,"J",sep=""),paste("cdr3_",ch,sep=""))
      ind <- apply(input[,cl],1,function(x){any(is.na(x))})
      input[ind,cl] <- NA
    }
  }
  #Remove empty lines (No longer since it's convenient to write them in processed_data)
  #ind <- apply(es.all,1,function(x){ s <- length(which(is.na(x[col])==F)); return(s)})
  #es.all <- es.all[which(ind>0),]

  return(input)

}

#' @export
check_cdr3 <- function(input, chain="AB", species.default="HomoSapiens", check.cdr3.mode=1, start.lg=1, end.lg=2, verbose=1){

  # Clean the CDR3 based on the V and J usage.
  # This should be applied after correcting the gene names, and adding the species if needed
  # species.default is only used if es.all does not contain the "species" column
  # If the allele is given in the gene name, the allele will be used.

  use.species.default <- F
  if("species" %in% colnames(input)){
    species.list <- unique(input[,"species"])
  } else {
    species.list <- c(species.default)
    use.species.default <- T
  }

  chain.list <- paste("TR",strsplit(chain,split="")[[1]], sep="")


  for(ch in chain.list){

    V <- paste(ch,"V",sep="")
    J <- paste(ch,"J",sep="")
    cdr3 <- paste("cdr3_",ch,sep="")

    for(species in species.list){

      if(!use.species.default){
        ind.species <- which(input[,"species"]==species)
      } else{
        ind.species <- 1:dim(input)[1]
      }

      if(check.cdr3.mode==0){
        ind.first <- c()
        ind.last <- c()
      }

      ind.traj38 <- c()

      if(check.cdr3.mode==1){

        #Correct TRAJ38 in human (e.g., issue with some 10X data)
        if(species=="HomoSapiens" & ch=="TRA"){

          ind.traj38 <- which((input[ind.species,J]=="TRAJ38" | input[ind.species,J]=="TRAJ38*01") & str_sub(input[ind.species,cdr3],start=-2)=="LI" & nchar(input[ind.species,cdr3]) < Lmax)
          input[ind.species[ind.traj38],cdr3] <- paste(input[ind.species[ind.traj38],cdr3], "W", sep="")

        }

        #Extract the first (start.lg) and last (end.lg) amino acids
        first <- substr(input[ind.species,cdr3], 1, start.lg)
        last <- str_sub(input[ind.species,cdr3], start=-end.lg)

        #print(first)
        #print(last)

        #Find cases incompatible with the reference (allowing matching to any allele)
        nm <- input[ind.species,V]
        rf <- sapply(ref.cdr3.first[[species]][[ch]], function(x){unique(substr(x,1,start.lg))}) #This includes all non-redundant allelic variants

        diff.first <- sapply(1:length(first), function(i){
          diff <- F
          if(!is.na(nm[i]) & !is.na(first[i])){
            diff <- T
            #Check if 'first' matches one of the possible allele of the V gene
            for(st in rf[[nm[i]]]){
              if(nchar(first[i])==nchar(st)){
                if(first[i] == st ){
                  diff <- F
                }
              } else if (nchar(first[i]) > nchar(st)) {
                #This is the special case were 'first' is actually longer than some rf[[nm[i]]]
                #Typically, this is the case when using large values for start.lg
                if(substr(first[i],1,nchar(st)) == st){
                  diff <- F
                }
              }
            }
          }
          return(diff)
        })
        ind.first <- (1:length(first))[diff.first]

        #Find cases incompatible with the reference (allowing matching to any allele)
        nm <- input[ind.species,J]
        rf <- sapply(ref.cdr3.last[[species]][[ch]], function(x){unique(str_sub(x,start=-end.lg))})

        diff.last <- sapply(1:length(last), function(i){
          diff <- F
          if(!is.na(nm[i]) & !is.na(last[i])){
            diff <- T
            #Check if 'last' matches one of the possible allele of the V gene
            for(st in rf[[nm[i]]]){
              if(nchar(last[i])==nchar(st)){
                if(last[i] == st ){
                  diff <- F
                }
              } else if (nchar(last[i]) > nchar(st)) {
                #This is the special case were 'last' is actually longer than some rf[[nm[i]]]
                #Typically, this is the case when using large values for start.lg
                if(str_sub(last[i],start=-nchar(st)) == st){
                  diff <- F
                }
              }
            }
          }
          return(diff)
        })
        ind.last <- (1:length(last))[diff.last]
      }

      if(verbose>0){

        if(length(ind.traj38)>0 & ch=="TRA"){
          print("Adding W on CDR3a for TRAJ38 entries:")
          if(verbose<=2){
            n <- min(10,length(ind.traj38))
            if(n<length(ind.traj38)){
              print("Examples  (use verbose > 2 to see them all):")
            }
          }
          if(verbose > 2){
            n <- length(ind.traj38)
          }
          ti <- ind.species[ind.traj38[1:n]]
          print(input[ti,c(cdr3,J)])
        }

        if(length(ind.first)>0){

          nt <- length(which(!is.na(input[ind.species,V]) & !is.na(input[ind.species,cdr3])))

          entry <- ifelse(length(ind.first)==1,"entry","entries")
          print(paste("*** Likely inconsistencies between ",ch,"V gene and CDR3",chain.small[ch]," in ",length(ind.first)," ",entry," (out of ",nt,") in ",species," ***",sep=""))
          if(verbose==1){
            n <- min(10,length(ind.first))
            if(length(ind.first)>n){
              print("Examples  (use verbose >= 2 to see them all):")
            }
          }
          if(verbose > 1){
            n <- length(ind.first)
          }
          if(verbose>=1){
            ti <- ind.species[ind.first[1:n]]
            sg <- input[ti,V]
            m.prob <- cbind(input[ti,c(V,cdr3)], sapply(sg,function(x){paste(ref.cdr3.first[[species]][[ch]][[x]], collapse=" / ")}))
            colnames(m.prob) <- c(V,cdr3,"Ref_CDR3_start")
            print(m.prob)
          }
          cat("\n")
        }
        if(length(ind.last)>0){

          nt <- length(which(!is.na(input[ind.species,J]) & !is.na(input[ind.species,cdr3])))

          entry <- ifelse(length(ind.first)==1,"entry","entries")
          print(paste("*** Likely inconsistencies between ",ch,"J gene and CDR3",chain.small[ch]," in ",length(ind.last)," ",entry," (out of ",nt,") in ",species," ***",sep=""))
          if(verbose==1){
            n <- min(10,length(ind.last))
            if(length(ind.last)>10){
              print("Examples (use verbose>=2 to see them all):")
            }
          }
          if(verbose>1){
            n <- length(ind.last)
          }
          if(verbose>=1){
            ti <- ind.species[ind.last[1:n]]
            sg <- input[ti,J]
            m.prob <- cbind(input[ti,c(J,cdr3)], sapply(sg,function(x){paste(ref.cdr3.last[[species]][[ch]][[x]], collapse=" / ")}))
            colnames(m.prob) <- c(J,cdr3,"Ref_CDR3_end")
            print(m.prob)
          }
          cat("\n")
        }
      }

      input[ind.species[ind.first],c(V,cdr3)] <- NA
      input[ind.species[ind.last],c(J,cdr3)] <- NA

    }
  }

  return(input)
}


correct.VJnames <- function(input, segment.list=c("TRAV","TRAJ","TRBV","TRBJ"), species.default="HomoSapiens",
                            use.allele=F, seq.protocol="Default", verbose=1){

  if("species" %in% colnames(input)){
    species.list <- unique(input[,"species"])
  } else {
    species.list <- c(species.default)
  }


  for(species in species.list){
    for(s in segment.list){
      if(use.allele){
        name.list <- gene.allele.list[[seq.protocol]][[species]][substr(gene.allele.list[[seq.protocol]][[species]],1,4)==s]
      } else {
        name.list <- gene.list[[seq.protocol]][[species]][substr(gene.list[[seq.protocol]][[species]],1,4)==s]
      }
      if("species" %in% colnames(input)){
        ind <- which(!(input[,s] %in% name.list) & input[,"species"]==species & !is.na(input[,s]))
      } else {
        ind <- which(!(input[,s] %in% name.list) & !is.na(input[,s]))
      }

      if(length(ind)>0){

        nm <- strsplit(input[ind,s], split="*", fixed=T)
        gene <- unlist(lapply(nm, function(x){x[1]}))
        allele <- unlist(lapply(nm, function(x){x[2]}))

        ga <- unlist(lapply(1:length(gene), function(x){ clean.name.allele(gene=gene[x], allele=allele[x], species=species,
                                                                           use.allele=use.allele, seq.protocol=seq.protocol)}))

        #Set to NA cases where the gene names comes from another segment
        #This is because the clean.name.allele does not check if the gene is in the right column (e.g., TRAV12-2 in TRAJ column)
        ga[substr(ga,1,4) != s] <- NA

        if(verbose>0){
          i <- which(input[ind,s] != ga & is.na(ga)==F)
          if(length(i)>0){

            m.cor <- data.frame(original.name = input[ind[i],s], corrected.name = ga[i],row.names = NULL)
            m.cor <- m.cor[!duplicated(m.cor),]

            entry <- ifelse(length(i)==1,"entry","entries")
            nm <- ifelse(dim(m.cor)[1]==1,"name","names")
            verb <- ifelse(dim(m.cor)[1]==1,"was","were")
            print(paste("*** ",dim(m.cor)[1]," ",s," ",nm," in ",length(i)," ",entry," ",verb," corrected in ",species, "***",sep=""))
            if(verbose==1 | verbose==2){
              print("Use verbose=3 to see them")
            }
            if(verbose==3){
              print(m.cor)
            }

            cat("\n")
          }

          #Check the cases where the segment was not NA, but was put to NA (i.e., mapping of gene name failed)
          i <- which(!is.na(input[ind,s]) & input[ind,s] != "" & is.na(ga)==T)
          if(length(i)>0){
            v <- unique(input[ind[i],s])
            v <- v[!is.na(v)]
            print(paste("*** ",length(v), " ", s, " names in ",length(i)," entries could not be corrected in ",species," ***", sep=""))
            if(verbose==1){
              n <- min(10,length(v))
              if(length(v)>n){
                print("Examples  (use verbose >= 2 to see them all):")
              }
              print(v[1:n])
            }
            if(verbose>1){
              print(v)
            }
            cat("\n")
          }
        }

        input[ind,s] <- ga
      }
    }
  }
  return(input)
}

#' @export
merge_mouse_TRAV <- function(input){

  # This has to be run after alleles have been removed and genes have been corrected
  # If the "species" field is present, it takes only "MusMusculus" entries
  # If not, it assumes all entries are "MusMusculus"
  # It also assumes that alleles have been removed

  if("TRAV" %in% colnames(input)){

    if("species" %in% colnames(input)){
      ind <- which(input[,"species"]=="MusMusculus")
    } else {
      ind <- 1:dim(input)[1]
    }

    v.cor <- as.character(unlist(lapply(input[["TRAV"]][ind], function(y){  # WARNING: I don't understand why not using es[ind,"TRAV"]
      if (y %in% names(merge.mouse.TRAV)){
        y <- merge.mouse.TRAV[y]
      }
      return(y)
    })))

    input[ind,"TRAV"] <- v.cor
  }
  return(input)

}

#' Take the gene + allele.
#' If allele is empty, try to correct the gene if needed, and return only the gene.
#' If allele is not empty, try to correct the gene\*allele.
#' If the gene can be corrected, but the gene\*allele does not exist, return gene\*default.allele
#' If the gene cannot be corrected, return NA.
clean.name.allele <- function(gene=gene, allele=allele, species="HomoSapiens", seq.protocol="Default", use.allele=F){

  if(species != "HomoSapiens" & species != "MusMusculus"){
    print("Undefined species: ",species)
  }

  if(!use.allele){
    allele <- ""
  }

  if(is.na(gene) | gene==""){  #This is not needed in MixTCRviz, but can be useful in other cases
    ga <- NA
  } else {

    #Check if the gene needs to be corrected
    if(gene %in% gene.list[[seq.protocol]][[species]] == F){
      #Do a few automatic corrections
      gene <- gsub("TCR","TR",gene)
      gene <- gsub("–","-", gene)
      gene <- gsub("-0","-", gene)
      if(species=="HomoSapiens"){ gene <- gsub("hTR","TR",gene) }
      if(species=="MusMusculus"){ gene <- gsub("mTR","TR",gene) }
      st0 <- substr(gene,1,2)
      st <- substr(gene,3,4)
      if(st0=="TR" & st %in% c("AV","AJ","BV","BJ")){
        gene <- gsub(paste("TR",st,"0",sep=""),paste("TR",st,sep=""),gene)
      }
      if(substr(gene,1,4)=="TRAJ" & grepl("-",gene)){
        gene <- unlist(strsplit(gene,"-"))[1]
      }

      if(gene %in% gene.list[[seq.protocol]][[species]] == F ){
        #The gene is still not ok, try using our manual dictionary
        if(!is.na(map[[species]][gene])){
          gene <- map[[species]][gene] #The gene was wrong, but can be corrected
        } else {
          gene <- NA #The gene was wrong and could not be corrected
        }
      }
    }

    if(use.allele){
      if(!is.na(gene)){
        #The gene was ok or could be corrected
        v <- paste(gene,allele,sep="*")
        if(v %in% gene.allele.list[[seq.protocol]][[species]] == T){
          ga <- v #the gene*allele is ok
        } else {
          ga <- paste(gene,allele.default[[seq.protocol]][[species]][gene],sep="*") #If not, Use the gene*default.allele
        }
      } else {
        #The gene could not be corrected
        ga <- NA
      }
    } else {
      ga <- gene
    }
  }
  return(ga)
}

#No longer used...
add_alleles <- function(TCR, segment.list=c("TRAV", "TRAJ", "TRBV", "TRBJ"), species.default="HomoSapiens"){

  # If allele is missing, add the default one (or "01" is default is not known, which can happen if people use non-standard V/J names)
  # Important: this function does not attempt to correct V/J names

  if("species" %in% names(TCR)){
    species <- as.character(TCR["species"])
  } else {
    species <- species.default
  }

  for(s in segment.list){
    if(s %in% names(TCR)){
      a <- unlist(strsplit(as.character(TCR[s]),split="*", fixed=T))

      if(length(a)==1){
        if(a[1] %in% names(allele.default[[seq.protocol]][[species]])){
          TCR[s]=paste(a[1],allele.default[[seq.protocol]][[species]][a[1]], sep="*")
        } else if(!is.na(a[1])) {
          TCR[s]=paste(a[1],"01",sep="*")
        }
      }
    }
  }
  return(TCR)
}


#' @export
inferVJ <- function(input, species.default="HomoSapiens", chain="AB", verbose=1){

  # This function uses all residues N-terminal (i.e., before C) of CDR3 for each V as unique barcodes
  # This function uses all residues C-terminal (i.e., including the last F/W) of CDR3 for each J as unique barcodes

  # For a few J segments, the seond-to-last amino acid in J segment is also used for disambiguation
  # For a few J segments, two additional amino acids in the CDR3 are used  for disambiguation

  #We could also provide a summary of the unsuccessful cases

  chain.list <- paste0("TR",strsplit(chain,split="")[[1]])
  chain.small <- tolower(strsplit(chain,split="")[[1]])
  names(chain.small) <- chain.list

  if("species" %in% colnames(input)){
    #This is always the case in MixTCRviz, but not if the function is used independently
    species.list <- unique(input$species)
  } else {
    species.list <- species.default
  }

  # Transform to simple data.frame in case a tibble due to code
  # incompatibilities below
  if (is_tibble(input)){
    input <- as.data.frame(input)
    if (verbose > 0){
      print(paste0("WARNING: input was a tibble, but has been transformed to ",
        "data.frame for compatibility with inferVJ function."))
    }
  }

  #Add the required columns
  for(ch in chain.list){
    if(paste0(ch,"V") %in% colnames(input) & verbose>0){
      if(!all(is.na(input[,paste0(ch,"V")])) && any( input[,paste0(ch,"V")]!="" ) ){
        print(paste0("WARNING: using infer.VJ=T will erase all data in ",ch,"V column"))
      }
    }
    input[[paste0(ch,"V")]] <- NA
    if(paste0(ch,"J") %in% colnames(input) & verbose>0){
      if(!all(is.na(input[,paste0(ch,"J")])) && any(input[,paste0(ch,"J")]!="" )){
        print(paste0("WARNING: using infer.VJ=T will erase all data in ",ch,"J column"))
      }
    }
    input[[paste0(ch,"J")]] <- NA
  }


  for(sp in species.list){

    if("species" %in% colnames(input)){
      ind.sp <- which(input$species==sp)
    } else {
      ind.sp <- 1:dim(input)[1]
    }

    #This part may be slow for large datasets... Things could be done MUCH faster with more efficient string matching algorithms.
    for(ch in chain.list){

      tc <- paste0("TCR",chain.small[ch])

      #Get the V gene
      nm <- names(inferV[[sp]][[ch]])
      pos.V <- lapply(input[ind.sp,tc], function(x){
        which(sapply(nm, function(n){grepl(n,x)}))})
      V <- sapply(pos.V, function(x){v <- unique(inferV[[sp]][[ch]][x]);if(length(v)!=1){return(NA)} else{ return(v)}})
      input[ind.sp,paste0(ch,"V")] <- V

      #Get the J gene
      nm <- names(inferJ[[sp]][[ch]])
      pos.J <- lapply(input[ind.sp,tc], function(x){
        which(sapply(nm, function(n){grepl(n,x)}))})
      J <- sapply(pos.J, function(x){j <- unique(inferJ[[sp]][[ch]][x]);if(length(j)!=1){return(NA)} else{ return(j)}})
      input[ind.sp,paste0(ch,"J")] <- J

    }
  }
  return(input)

}


#' @export
inferCDR3 <- function(input, species.default="HomoSapiens", chain="AB", verbose=1){

  #The starting and end position of the CDR3 are inferred based on the 10 residues before the CDR3 and all the residues after the CDR3
  #This function does not require V/J knowledge
  #This means the inferred CDR3 may not be compatible with V/J data given in input.

  chain.list <- paste0("TR",strsplit(chain,split="")[[1]])
  chain.small <- tolower(strsplit(chain,split="")[[1]])
  names(chain.small) <- chain.list

  if("species" %in% colnames(input)){
    #This is always the case in MixTCRviz, but not if the function is used independently
    species.list <- unique(input$species)
  } else {
    species.list <- species.default
  }

  #Add the required columns
  for(ch in chain.list){
    if(paste0("cdr3_",ch) %in% colnames(input) & verbose>0){
      if(any(input[,paste0("cdr3_",ch)]!="" )){
        print(paste0("WARNING: using infer.CDR3=T will erase all data in cdr3_",ch," column"))
      }
    }
    input[[paste0("cdr3_",ch)]] <- NA
  }


  for(sp in species.list){

    if("species" %in% colnames(input)){
      ind.sp <- which(input$species==sp)
    } else {
      ind.sp <- 1:dim(input)[1]
    }

    #This part may be slow for large datasets... Things could be done MUCH faster with more efficient string matching algorithms.
    for(ch in chain.list){

      tc <- paste0("TCR",chain.small[ch])

      #Take all possible last 10 amino acids before CDR3
      tind <- which(cdr123[[sp]][[ch]][,"CDR3"]!="" & cdr123[[sp]][[ch]][,"full"]!="") #Exclude weird genes/alleles with nothing in the CDR3 of in the full V
      V.out <- unique(apply(cdr123[[sp]][[ch]][tind,],1,function(x){
        substr(as.character(x["full"]), nchar(as.character(x["full"]))-nchar(as.character(x["CDR3"]))-9, nchar(as.character(x["full"]))-nchar(as.character(x["CDR3"])))
      }))
      #Take all possible amino acids after CDR3
      J.out <- unique(apply(Jseq[[sp]][[ch]],1,function(x){
        substr(as.character(x["full"]), nchar(as.character(x["CDR3"]))+1, nchar(as.character(x["full"])))
      }))

      #Get the CDR3
      CDR3 <- unlist(lapply(input[ind.sp,tc], function(x){

        #Check if some V.out can be found
        V.out.grep <- unique(unlist(lapply(V.out, function(s){p <- regexpr(s,as.character(x)); if(length(p)==1){return(p)} else{return(NA)}})))
        V.out.grep <- V.out.grep[V.out.grep>0]
        J.out.grep <- unique(unlist(lapply(J.out, function(s){p <- regexpr(s,as.character(x)); if(length(p)==1){return(p)} else{return(NA)}})))
        J.out.grep <- J.out.grep[J.out.grep>0]
        if(length(V.out.grep)==1 & length(J.out.grep)==1){
          cdr3 <- substr(x,V.out.grep[1]+10,J.out.grep[1]-1)
        } else {
          cdr3 <- NA
        }
        return(cdr3)
      }))

      input[ind.sp,paste0("cdr3_",ch)] <- CDR3
    }
  }
  return(input)

}

#' Defines the color palette for plots with models showed together.
#'
#' Little function that defines the color palette to use for the different
#' models that are combined together when plot.modelsCombined is T.
#' These are the colors used around the bars and for the length distribution plots.
#  In case a 'Trash' model is found, its color is set as light gray.
#'
#' @param models A vector of the models that are currently used in the figure.
#' @return A named vector of the colors to use for each model.
set_model_colPals <- function(models){
  colPal_modelsCombined <- c(palette.colors(palette="Okabe-Ito"),
                             palette.colors(palette="Tableau 10"), palette.colors(palette="Polychrome 36"))
  # Make it large, but will likely never use all its colors
  if (length(models) > length(colPal_modelsCombined)){
    colPal_modelsCombined <- c(colPal_modelsCombined,
                               rainbow(n=length(models)-length(colPal_modelsCombined), alpha=NULL))
  }
  colPal_modelsCombined <- colPal_modelsCombined[1:length(models)]
  names(colPal_modelsCombined) <- models
  if (any(grepl("Trash", models))){
    colPal_modelsCombined[grep("Trash", models)] <- "gray70"
  }
  return(colPal_modelsCombined)
}

create_interactive_plots <- function(countV.plot,countJ.plot,ld.plot,CDR3,plot.oneline){
  # Turn off legends in the first two plots
  countV.plot_not_title <- countV.plot + labs(title = NULL)

  p1 <- plotly::ggplotly(countV.plot_not_title, tooltip = c("name", "logFC", "Zscore", "CDR1_seq", "CDR2_seq", "CDR3_seq"))

  p1$x$data <- lapply(p1$x$data, function(trace) {
    # Set marker size and mode for points

    if (!is.null(trace$marker$size)){trace$marker$size <- 11}
    if (!is.null(trace$error_x)){
      trace$error_x$width <- 3
      trace$error_x$color <- rgb(102 / 255, 102 / 255, 102 / 255, alpha = 0.4)
    }

    return(trace)
  })

  # p1 <- p1 %>%
  #   plotly::layout(
  #     title = list(
  #       text    = countV.plot$labels$title,
  #       x       = 0.5,         # Center the title horizontally
  #       y       = 1.02,         # Position the title at the top edge of the plotting area
  #       xanchor = "center",
  #       yanchor = "bottom",    # With yanchor = "bottom", the bottom of the title text aligns with y = 1.0
  #       pad     = list(b = 0)  # Remove extra bottom padding for the title
  #     ),
  #     margin = list(t = 60),   # Reduce top margin so the gap between title and plot is smaller
  #     showlegend = FALSE
  #   )
  #
  p1 <- p1 %>%
    plotly::layout(
      title = NULL,           # no default title
      margin = list(t = 80),  # minimal top margin
      annotations = list(
        list(
          text      = unname(countV.plot$labels$title),
          x         = 0.5,     # center horizontally
          y         = 1.04,    # still in "paper" coords, so 1.05 is just above the plot
          xref      = "paper",
          yref      = "paper",
          xanchor   = "center",
          yanchor   = "bottom",
          font = list(size = 20),
          showarrow = FALSE
        )
      ),
      showlegend = FALSE
    )


  countJ.plot_not_title <- countJ.plot + labs(title = NULL)

  p2 <- plotly::ggplotly(countJ.plot_not_title, tooltip = c("name", "logFC" ,"Zscore", "CDR3_seq"))

  p2$x$data <- lapply(p2$x$data, function(trace) {
    # Set marker size and mode for points

    if (!is.null(trace$marker$size)){trace$marker$size <- 11}
    if (!is.null(trace$error_x)){
      trace$error_x$width <- 3
      trace$error_x$color <- rgb(102 / 255, 102 / 255, 102 / 255, alpha = 0.4)
    }

    return(trace)
  })


  # p2 <- p2 %>%
  #   plotly::layout(
  #     title = list(
  #       text = countJ.plot$labels$title
  #       #y = 0.9    # Move title a bit from the left, adjust as desired
  #     ),
  #     margin = list(t = 80) , # Increase top margin to prevent overlap
  #     showlegend = FALSE
  #   )


  p2 <- p2 %>%
    plotly::layout(
      title = NULL,           # no default title
      margin = list(t = 80),  # minimal top margin
      annotations = list(
        list(
          text      = unname(countJ.plot$labels$title),
          x         = 0.5,     # center horizontally
          y         = 1.04,    # still in "paper" coords, so 1.05 is just above the plot
          xref      = "paper",
          yref      = "paper",
          xanchor   = "center",
          yanchor   = "bottom",
          font = list(size = 20),
          showarrow = FALSE
        )
      ),
      showlegend = FALSE
    )

  # Adjust legend in the third plot and remove legend title
  p3 <- plotly::ggplotly(ld.plot, tooltip = "none") %>%
    plotly::layout(
      legend = list(
        orientation = "h",
        x = 0.5,
        xanchor = "center",
        y = 1.125,
        yanchor = "bottom",
        title = list(text = NULL)  # Remove legend title
      ),
      margin = list(t = 50)
    )



  # Prepare the last two plots
  if(plot.oneline!=2){
    p4 <- plotly::ggplotly(CDR3$ES_max, tooltip = "none")
    p5 <- plotly::ggplotly(CDR3$Baseline_max, tooltip = "none")
    bottom_subplot <- manipulateWidget::combineWidgets(p4, p5, ncol = 1, title = NULL)
  }


  if(plot.oneline==0){
    combined_plots <- manipulateWidget::combineWidgets(
      p1, p2,
      p3, bottom_subplot,
      ncol = 2,
      title = NULL
    )
  }
  if(plot.oneline==1){
    combined_plots <- manipulateWidget::combineWidgets(
      p1, p2,
      p3, bottom_subplot,
      ncol = 4,
      nrow = 1,
      title = NULL
    )
  }

  if(plot.oneline==2){
    combined_plots <- manipulateWidget::combineWidgets(
      p1, p2,p3,
      ncol = 3,
      nrow = 1,
      title = NULL
    )
  }

  return(combined_plots)
}


verify.chain <- function(input, chain){

  # Check cases where people leave the chain="AB",
  # but actually provide single chain data, including with ambiguous colnames like V,J,CDR3_seq

  format <- determine.format(input, verbose=0)
  #Only check for the format not based on clone.id (cases with formats based on clone.id will always be treated as alpha+beta, evn if one chain is empty)
  if(!format %in% names(clone.format.col)){
    if(chain=="AB"){

      map.A <- names(mapping.colnames[["AB"]][which(mapping.colnames[["AB"]] %in% c("TRAV", "TRAJ", "cdr3_TRA"))])
      map.B <- names(mapping.colnames[["AB"]][which(mapping.colnames[["AB"]] %in% c("TRBV", "TRBJ", "cdr3_TRB"))])
      inter.A <- intersect(colnames(input), c("TRAV", "TRAJ", "cdr3_TRA", map.A))
      inter.B <- intersect(colnames(input), c("TRBV", "TRBJ", "cdr3_TRB", map.B))

      if(length(inter.A)==0 | length(inter.B)==0){
        if(length(inter.A)==0 & length(inter.B)>0){
          print("Missing columns for alpha chain, only beta chain will be considered")
          chain <- "B"
        } else if(length(inter.A)>0 & length(inter.B)==0){
          print("Missing columns for beta chain, only alpha chain will be considered")
          chain <- "A"
        } else if(length(inter.A)==0 & length(inter.B)==0){
          # Check if the a single chain can be inferred.
          # This is the case for instance when providing data in AIRR format without specifying the chain
          map.unk <- names(mapping.colnames[["A"]][which(mapping.colnames[["A"]] == "TRAV")])
          if(length(intersect(map.unk, colnames(input)))==1){
            p <- intersect(map.unk, colnames(input))[1]
            if("TRA" %in% unique(substr(input[,p],1,3)) & ! "TRB" %in% unique(substr(input[,p],1,3))){
              print("Missing data for beta chain, only alpha chain will be considered")
              chain <- "A"
            } else if("TRB" %in% unique(substr(input[,p],1,3)) & ! "TRA" %in% unique(substr(input[,p],1,3))){
              print("Missing data for alpha chain, only beta chain will be considered")
              chain <- "B"
            } else {
              chain <- ""
            }
          }
        }
      }
    } else if(chain=="A"){

      map.A <- names(mapping.colnames[["A"]][which(mapping.colnames[["A"]] %in% c("TRAV", "TRAJ", "cdr3_TRA"))])
      inter.A <- intersect(colnames(input), c("TRAV", "TRAJ", "cdr3_TRA", map.A))
      if(length(inter.A)==0){
        chain <- ""
      }

    } else if(chain=="B"){

      map.B <- names(mapping.colnames[["B"]][which(mapping.colnames[["B"]] %in% c("TRBV", "TRBJ", "cdr3_TRB"))])
      inter.B <- intersect(colnames(input), c("TRBV", "TRBJ", "cdr3_TRB", map.B))
      if(length(inter.B)==0){
        chain <- ""
      }

    }
  }
  return(chain)
}

#Determine the format
determine.format <- function(input, verbose=1){

  col <- c("TRAV","TRAJ","cdr3_TRA","TRBV","TRBJ","cdr3_TRB")
  col.A <- col[1:3]
  col.B <- col[4:6]

  format.list <- names(clone.format.col)
  format <- "custom"
  if(length(intersect(colnames(input), col.A))==3 | length(intersect(colnames(input), col.B))==3){
    format <- "MixTCRviz"
  } else {
    for(f in format.list){
      #Check that the three field for specific formats are present, and none of the standard colnames in MixTCRviz
      if(length(intersect(colnames(input), clone.format.col[[f]]))==3 & length(intersect(colnames(input),col))==0){
        format <- f
        if(verbose==1){
          print(paste("Inferred format:",f))
        }
        break
      }
    }
  }
  return(format)
}
#Build input files with stacked alpha and beta chains on different rows if the data are provided in a format based on clone IDs.
#In this case, the actual clones are not reconstructed, and data are treated as unpaired
stack_clones <- function(input, format){

  if(format %in% names(clone.format.col)){

    col <- clone.format.col[[format]]
    other.col <- setdiff(colnames(input), col)

    input.f <- apply(input,1,function(x){
      if(substr(x[col[1]],1,3)=="TRA"){
        v <- c(x[col],NA,NA,NA,x[other.col])
      } else if(substr(x[col[1]],1,3)=="TRB"){
        v <- c(NA,NA,NA,x[col],x[other.col])
      } else if(format=="VDJdb" & "Gene" %in% names(x)){
        if(x["Gene"]=="TRA"){
          v <- c(x[col],NA,NA,NA,x[other.col])
        } else if(x["Gene"]=="TRB"){
          v <- c(NA,NA,NA,x[col],x[other.col])
        }
      } else if("chain" %in% names(x)){
        if(x["chain"]=="TRA"){
          v <- c(x[col],NA,NA,NA,x[other.col])
        } else if(x["chain"]=="TRB"){
          v <- c(NA,NA,NA,x[col],x[other.col])
        }
      } else {
        #In this case, inferring the chain failed
        v <- c(NA,NA,NA,NA,NA,NA,x[other.col])
      }
      return(v)
    })
    input.f <- t(input.f)
    colnames(input.f) <- c("TRAV","TRAJ","cdr3_TRA","TRBV","TRBJ","cdr3_TRB",other.col)
    input.f <- as.data.frame(input.f)
    if(format=="VDJdb" & "complex.id" %in% colnames(input.f)){
      input.f[,"complex.id"] <- as.numeric(input.f[,"complex.id"])
    }
  } else {
    input.f <- input
  }
  return(input.f)
}

#Build actual clones if the data are provided in a format based on clone IDs.
merge_clones <- function(input, format){

  col <- c("TRAV","TRAJ","cdr3_TRA","TRBV","TRBJ","cdr3_TRB")

  if(length(intersect(colnames(input), clone.id))==1){

    TCR.col <- clone.format.col[[format]]
    Vn <- TCR.col[1]
    barcode.label <- intersect(colnames(input), clone.id)

    #other.col <- setdiff(colnames(input), TCR.col) #keep all columns (can be complex if different values are used with the same barcode)
    other.col <- intersect(colnames(input),c(barcode.label,"model","species")) #keep only the barcode, the model and the species (if present)
    input.f <- as.data.frame(matrix(nrow=0, ncol=6+length(other.col)))

    if("model" %in% colnames(input)){
      barcode.list <- paste(input[,barcode.label], input[,"model"]) #Include the possibility of having redundant barcodes for different models.
    } else {
      barcode.list <- input[,barcode.label]
    }

    input.t <- input[order(barcode.list),]
    barcode.list <- sort(barcode.list)

    #Go through all clone IDs
    i <- 1
    while(i<=dim(input.t)[1]){
      seq.A <- list()
      seq.B <- list()
      ct.A <- 0
      ct.B <- 0
      #This is the special case of the '0' complex.id for VDJdb which contains all single chain entries
      if(input.t[i,barcode.label]==0 & format=="VDJdb"){
        ind.A <- which(input.t[,barcode.label]==0 & substr(input.t$V,1,3)=="TRA")
        ind.B <- which(input.t[,barcode.label]==0 & substr(input.t$V,1,3)=="TRB")
        if(length(ind.A)>0){
          input.f <- rbind(input.f,cbind(input.t[ind.A,TCR.col],NA,NA,NA,input.t[ind.A,other.col]))
        }
        if(length(ind.B)>0){
          input.f <- rbind(input.f,cbind(NA,NA,NA,input.t[ind.B,TCR.col],input.t[ind.B,other.col]))
        }
        i <- i+length(ind.A)+length(ind.B)

      } else {
        for(j in 0:(dim(input.t)[1]-i)){
          if(barcode.list[i+j]==barcode.list[i]){
            if(substr(input.t[i+j,Vn],1,3)=="TRA"){
              ct.A <- ct.A+1
              seq.A[[ct.A]] <- as.character(input.t[i+j,TCR.col])
            } else if(substr(input.t[i+j,Vn],1,3)=="TRB"){
              ct.B <- ct.B+1
              seq.B[[ct.B]] <- as.character(input.t[i+j,TCR.col])
            }
          } else {
            break
          }
        }


        if(length(seq.A)>0 & length(seq.B)>0){
          for(sA in seq.A){
            for(sB in seq.B){
              seq.all <- c(sA,sB,unlist(input.t[i,other.col]))
              input.f <- rbind(input.f, seq.all)
            }
          }
        } else if(length(seq.A)>0 & length(seq.B)==0){
          for(sA in seq.A){
            seq.all <- c(sA,c(NA,NA,NA),unlist(input.t[i,other.col]))
            input.f <- rbind(input.f, seq.all)
          }
        } else if(length(seq.A)==0 & length(seq.B)>0){
          for(sB in seq.B){
            seq.all <- c(c(NA,NA,NA),sB,unlist(input.t[i,other.col]))
            input.f <- rbind(input.f, seq.all)
          }
        }
        i <- i+ct.A+ct.B
      }
    }
    colnames(input.f) <- c(col, other.col)
    if(format=="VDJdb" & "complex.id" %in% colnames(input.f)){
      input.f[,"complex.id"] <- as.numeric(input.f[,"complex.id"])
      print(input.f[,"complex.id"])
    }
    return(input.f)

  } else {
    print(paste("WARNING: Unable to reconstruct clones. Clone_id should be exactly one element in",paste(clone.id, collapse=" / ")))
    print("  Try using build.clones = F")
    return(input)
  }


}



