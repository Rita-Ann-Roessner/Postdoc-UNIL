
devtools::load_all("/Users/davidgfeller/Research/LICR/Github/MixTCRviz/")

#Make sure all correct V/J names are part of the list

for(s in c("TRAV", "TRAJ", "TRBV", "TRBJ")){

  print(s)
  mp <- read.csv(paste("mapping_",s,".csv",sep=""))
    
  for(species in c("HomoSapiens", "MusMusculus")){
    print(species)
    tmp <- mp[which(mp$species==species & mp$Done=="Yes"),]
    ind <- which(!tmp$Vcorrect %in% gene.list[[species]])
    if(length(ind)>0){
      print(tmp[ind,])
    }
    
  }
  
}
