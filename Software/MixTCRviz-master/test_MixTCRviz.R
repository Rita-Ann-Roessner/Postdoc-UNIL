
#Run it from the MixTCRviz/ folder

#Do this if you do not want to install the package
devtools::load_all(".")

#Do this if you have installed the package
#library(MixTCRviz)

#Compare input TCRs (specific for A0201_LLWNGPMAV) with baseline repertoire
MixTCRviz(input1="test/test.csv", output.path="test/out/", output.stat = T)

new.data <- F

if(new.data){
  
  system("cp test/out/A0201_LLWNGPMAV* test/out_compare/")
  system("cp test/out/stats/A0201_LLWNGPMAV.rds test/out_compare/stats/A0201_LLWNGPMAV.rds")
  
} else {
  
  list_1 <- readRDS("test/out/stats/A0201_LLWNGPMAV.rds")
  list_2 <- readRDS("test/out_compare/stats/A0201_LLWNGPMAV.rds")
  
  comp <- identical(list_1, list_2)
  if(comp){
    print("No problem detected")
  } else {
    stop("There were some issues...")
  }
  
}
