
library(MixTCRviz)
library(stringr)
library(umap)
library(ggplot2)

library(bio3d)
Sys.setenv(PATH = paste("/Users/davidgfeller/miniconda3/bin/", Sys.getenv("PATH"), sep = .Platform$path.sep))
Sys.setenv(PATH = paste("/usr/local/bin/", Sys.getenv("PATH"), sep = .Platform$path.sep))


pA <- which(!grepl("*",rownames(MixTCRviz::cdr123$HomoSapiens$TRA), fixed=T))
pB <- which(!grepl("*",rownames(MixTCRviz::cdr123$HomoSapiens$TRB), fixed=T) & 
              rownames(MixTCRviz::cdr123$HomoSapiens$TRB)!="TRBV6-3" &
              rownames(MixTCRviz::cdr123$HomoSapiens$TRB)!="TRBV6-2/6-3" &
              rownames(MixTCRviz::cdr123$HomoSapiens$TRB)!="TRBV12-3/12-4")

cdr1 <- list()
cdr1[["TRA"]] <- setNames(gsub("g","",MixTCRviz::cdr123$HomoSapiens$TRA$CDR1[pA]), rownames(MixTCRviz::cdr123$HomoSapiens$TRA)[pA])
cdr1[["TRB"]] <- setNames(gsub("g","",MixTCRviz::cdr123$HomoSapiens$TRB$CDR1[pB]), rownames(MixTCRviz::cdr123$HomoSapiens$TRB)[pB])

cdr2 <- list()
cdr2[["TRA"]] <- setNames(gsub("g","",MixTCRviz::cdr123$HomoSapiens$TRA$CDR2[pA]), rownames(MixTCRviz::cdr123$HomoSapiens$TRA)[pA])
cdr2[["TRB"]] <- setNames(gsub("g","",MixTCRviz::cdr123$HomoSapiens$TRB$CDR2[pB]), rownames(MixTCRviz::cdr123$HomoSapiens$TRB)[pB])

cdr3_pre <- list()
cdr3_pre[["TRA"]] <- setNames(unlist(apply(MixTCRviz::cdr123$HomoSapiens$TRA[pA,], 1, 
                                           function(x){lg <- nchar(x["full"]); lg.cdr3 <- nchar(x["CDR3"]);
                                           substr(x["full"],lg-lg.cdr3-15, lg-lg.cdr3)})),
                              rownames(MixTCRviz::cdr123$HomoSapiens$TRA)[pA])
cdr3_pre[["TRB"]] <- setNames(unlist(apply(MixTCRviz::cdr123$HomoSapiens$TRB[pB,], 1, 
                                           function(x){lg <- nchar(x["full"]); lg.cdr3 <- nchar(x["CDR3"]);
                                           substr(x["full"],lg-lg.cdr3-15, lg-lg.cdr3)})),
                              rownames(MixTCRviz::cdr123$HomoSapiens$TRB)[pB])


fA <- which(!grepl("*",rownames(MixTCRviz::Jseq$HomoSapiens$TRA), fixed=T))
fB <- which(!grepl("*",rownames(MixTCRviz::Jseq$HomoSapiens$TRB), fixed=T))

cdr3_post <- list()
cdr3_post[["TRA"]] <- setNames(unlist(apply(MixTCRviz::Jseq$HomoSapiens$TRA[fA,], 1, 
                                            function(x){lg <- nchar(x["full"]); lg.cdr3 <- nchar(x["CDR3"]);
                                            substr(x["full"],lg.cdr3+1, lg)})),
                               rownames(MixTCRviz::Jseq$HomoSapiens$TRA)[fA])

cdr3_post[["TRB"]] <- setNames(unlist(apply(MixTCRviz::Jseq$HomoSapiens$TRB[fB,], 1, 
                                            function(x){lg <- nchar(x["full"]); lg.cdr3 <- nchar(x["CDR3"]);
                                            substr(x["full"],lg.cdr3+1, lg)})),
                               rownames(MixTCRviz::Jseq$HomoSapiens$TRB)[fB])



find.V <- function(s, chain){
  
  V.1 <- names(cdr1[[chain]])[sapply(cdr1[[chain]], function(x){grepl(x,s)})]
  V.2 <- names(cdr2[[chain]])[sapply(cdr2[[chain]], function(x){grepl(x,s)})]
  
  V <- NA
  if(length(V.1)>=1 & length(V.2)>=1 & length(intersect(V.1, V.2))==1){
    V <- intersect(V.1, V.2)
  } else {
    print("Issues")
  }
  return(V)
}



find.cdr <- function(s,chain,V,J){
  
  p1 <- str_locate_all(pattern = cdr1[[chain]][V],s)[[1]]
  p2 <- str_locate_all(pattern = cdr2[[chain]][V],s)[[1]]
  
  p3_pre <- str_locate_all(pattern = cdr3_pre[[chain]][[V]],s)[[1]]
  p3_post <- str_locate_all(pattern=cdr3_post[[chain]][[J]],s)[[1]]
  
  if(length(p3_post)==0){
    t <- str_locate_all(pattern="EQYV",s)[[1]]
    p3_post <- c(t[2]+1,t[2]+1)
  }
  
  p3 <- c(p3_pre[2]+1, p3_post[1]-1)
  
  return(list(p1[1,],p2[1,], p3))
  
}

path <- "LAU5013/YF_LAU5013_sc_WT/"
#path <- "Public_Data/YF_public_pairedData_20251010/"

m.data <- read.table(paste0("data/",path, strsplit(path,split="/",fixed=T)[[1]][2],"_AF3.tsv" ), header=T, row.names = 1)
ind <- which(m.data$TRBV=="TRBV6-3")
if(length(ind)>0){
  m.data[ind,"TRBV"] <- "TRBV6-2"
}
tcr.list <- rownames(m.data)[which(!is.na(m.data$TEMPOrank) & !is.na(m.data$AF3_iptm_pair_mean) & 
                                     nchar(m.data$cdr3_TRA)>=10 & nchar(m.data$cdr3_TRB)>=10)]

align <- F
if(align){
  
  pdb1 <- read.pdb(paste0("data/",path,"/model_pdb/",tcr.list[1],"_model.cif.pdb"))
  system(paste0("cp data/",path,"/model_pdb/",tcr.list[1],"_model.cif.pdb data/",path,"model_pdb_align/",tcr.list[1],".pdb"))
  
  for(tcr in tcr.list[2:length(tcr.list)]){
    print(tcr)
    pdb2 <- read.pdb(paste0("data/",path,"model_pdb/",tcr,"_model.cif.pdb"))
    
    system("rm align_MHC/*.pdb")
    
    aln_structure <- struct.aln(pdb1, pdb2, write.pdbs = T, 
                                fixed.inds=atom.select(pdb1, chain="A"), 
                                outpath="align_MHC")
    
    system(paste0("cp align_MHC/mobile_",length(aln_structure$rmsd)-1,".pdb data/",path,"model_pdb_align/",tcr,".pdb"))
  }
  
} 

#For each aligned pdb, extract the CDR coordinates


chain.list <- c("TRA", "TRB")

cdr.coord <- list()
for(tcr in tcr.list){
  
  print(tcr)
  
  pdb <- read.pdb(paste0("data/",path,"model_pdb_align/",tcr,".pdb"))
  
  ind.C <- which(pdb$atom$elety=="C")
  pdb.seq <- paste(pdbseq(pdb), collapse="")
  pdb.chain <- pdb$atom$chain[ind.C]
  
  cdr.coord[[tcr]] <- list()
  
  for(chain in chain.list){
    
    V <- m.data[tcr,paste0(chain,"V")]
    J <- m.data[tcr,paste0(chain,"J")]
    
    cdr.pos <- find.cdr(pdb.seq, chain=chain, V=V, J=J)
    
    cdr.coord[[tcr]][[chain]] <- list()
    for(i in 1:3){
      coord <- pdb$atom[pdb$atom$elety=="C",]
      coord <- coord[cdr.pos[[i]][1]:cdr.pos[[i]][2],]
      cdr.coord[[tcr]][[chain]][[paste0("cdr",i)]] <- coord[,c("resid","x","y","z")]
    } 
  }
}

#Find the minimal length

cdr.list <- c("cdr1", "cdr2", "cdr3")

min.cdr <- list()
for(chain in chain.list){
  min.cdr[[chain]] <- list()
  
  min.cdr[[chain]][["cdr1"]] <- min(nchar(cdr1[[chain]][m.data[tcr.list,paste0(chain,"V")]]))
  min.cdr[[chain]][["cdr2"]] <- min(nchar(cdr2[[chain]][m.data[tcr.list,paste0(chain,"V")]]))
  min.cdr[[chain]][["cdr3"]] <- min(nchar(m.data[tcr.list,paste0("cdr3_",chain)]))
}

coord.pca <- data.frame(matrix(nrow=0, ncol=3*(min.cdr$TRA$cdr1+min.cdr$TRA$cdr2+min.cdr$TRA$cdr3+
                                                 min.cdr$TRB$cdr1+min.cdr$TRB$cdr2+min.cdr$TRB$cdr3)))
xyz <- c("x","y","z")

for(tcr in tcr.list){
  p <- c()
  for(chain in chain.list){
    for(cdr in cdr.list[1:2]){
      #For CDR1/2, take the first two residues, and then the last N (= 2 or 3 depending on the minimal CDR1/2 length)
      d <- dim(cdr.coord[[tcr]][[chain]][[cdr]])[1]
      v <- c(1:2, (d-(min.cdr[[chain]][[cdr]]-2)+1):d)
      p <- c(p,as.numeric(unlist(t(cdr.coord[[tcr]][[chain]][[cdr]][v,xyz]))))
    }
    #For CDR1/2, take the first five residues, and then the last N (= 5 or 6 depending on the minimal CDR3 length)
        d <- dim(cdr.coord[[tcr]][[chain]]$cdr3)[1]
    v <- c(1:5, (d-(min.cdr[[chain]]$cdr3-5)+1):d)
    p <- c(p,as.numeric(unlist(t(cdr.coord[[tcr]][[chain]]$cdr3[v,xyz]))))
    
  }
  coord.pca <- rbind(coord.pca,p)
}

cn <- c()
for(chain in chain.list){
  for(cdr in cdr.list){
    for(p in 1:min.cdr[[chain]][[cdr]]){
      cn <- c(cn,c(paste0(toupper(cdr), substr(chain,3,3),"_",p,"_X"),
                   paste0(toupper(cdr), substr(chain,3,3),"_",p,"_Y"),
                   paste0(toupper(cdr), substr(chain,3,3),"_",p,"_Z")))
    } 
  }
}
colnames(coord.pca) <- cn

col <- c("TRAV", "TRAJ", "cdr3_TRA", "TRBV", "TRBJ", "cdr3_TRB")


pca <- prcomp(coord.pca)
df <- data.frame(PC1=pca$x[,"PC1"], PC2=pca$x[,"PC2"], AF3=m.data[tcr.list,"AF3_iptm_pair_mean"], TEMPO=-m.data[tcr.list,"TEMPOrank"])
g <- ggplot(df, aes(x=PC1, y=PC2, color=AF3)) + geom_point(size=2) + theme_bw()
ggsave(g, file=paste0("plots/",path,"/PCA_all.pdf"), height = 3, width = 4)
g <- ggplot(df, aes(x=PC1, y=PC2, color=TEMPO)) + geom_point(size=2) + theme_bw()
ggsave(g, file=paste0("plots/",path,"/PCA_all_TEMPO.pdf"), height = 3, width = 4)

if(path=="LAU5013/YF_LAU5013_sc_WT/"){
  
  g2 <- ggplot(df, aes(x=PC1, y=PC2, color=AF3)) + geom_point(size=2) + theme_bw() + ylim(c(-40,25)) + xlim(c(-30,40))
  ggsave(g2, file=paste0("plots/",path,"/PCA_zoom.pdf"), height = 3, width = 4)
  g2 <- ggplot(df, aes(x=PC1, y=PC2, color=TEMPO)) + geom_point(size=2) + theme_bw() + ylim(c(-40,25)) + xlim(c(-30,40))
  ggsave(g2, file=paste0("plots/",path,"PCA_zoom_TEMPO.pdf"), height = 3, width = 4)

  #Check the variance per position
  ind <- which(pca$x[,"PC1"] > -30 & pca$x[,"PC1"]<40)
  tvar <- apply(coord.pca[ind,],2,var)
  var <- c()
  for(i in 1:(dim(coord.pca)[2]/3)){
    var <- setNames(c(var,sum(tvar[(i*3-2):(i*3)])), c(names(var), substr(names(tvar[i*3-2]), 1, nchar(names(tvar[i*3-2]))-2)))
  }
  df <- data.frame(val=var, pos=factor(names(var), levels = names(var)))
  g <- ggplot(df, aes(x=pos,y=val)) + geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

}

t.noise <- tcr.list[pca$x[,1] < -30 | pca$x[,1] > 60]
motif <- MixTCRviz(m.data[t.noise,col], output.path = paste0("plots/",path,"Motifs"), model.default = "Noise_PCA_TRA", plot.oneline = 1, chain="A")
motif <- MixTCRviz(m.data[t.noise,col], output.path = paste0("plots/",path,"Motifs"), model.default = "Noise_PCA")

t.good <- tcr.list[pca$x[,1] > -30 & pca$x[,1] < 60]
motif <- MixTCRviz(m.data[t.good,], output.path = "plots/",path,"Motifs", model.default = "Good_PCA_TRA", plot.oneline = 1, chain="A")
motif <- MixTCRviz(m.data[t.good,], output.path = "plots/",path,"Motifs", model.default = "Good_PCA")

set.seed(101)
umap <- umap(coord.pca)
df <- data.frame(UMAP1=umap$layout[,1], UMAP2=umap$layout[,2], AF3=m.data[tcr.list,"AF3_iptm_pair_mean"], TEMPO=m.data[tcr.list,"TEMPOrank"])
g <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=AF3)) + geom_point() + theme_bw() 
ggsave(g, file=paste0("plots/",path,"UMAP_all.pdf"), height = 3, width = 4)
g <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=TEMPO)) + geom_point() + theme_bw() + 
  scale_colour_gradient(name = waiver(),high = "#132B43",low = "#56B1F7")
ggsave(g, file=paste0("plots/",path,"UMAP_all_TEMPO.pdf"), height = 3, width = 4)

t.noise <- tcr.list[umap$layout[,2] < -5]
motif <- MixTCRviz(m.data[t.noise,col], output.path = paste0("plots/",path,"Motifs"), model.default = "Noise_UMAP")

t.good <- tcr.list[umap$layout[,2] > -5]
motif <- MixTCRviz(m.data[t.good,col], output.path = paste0("plots/",path,"Motifs"), model.default = "Good_UMAP")


#Other plots

ind <- which(m.data[tcr.list,"TEMPOrank"]>0.1)
MixTCRviz(m.data[tcr.list[ind],col], output.path = "Motifs", model.default = "Low_TEMPO")
ind <- which(m.data[tcr.list,"AF3_iptm_pair_mean"]<0.7)
MixTCRviz(m.data[tcr.list[ind],col], output.path = "Motifs", model.default = "Low_AF3")




