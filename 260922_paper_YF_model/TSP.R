# Computing motifs for YF-binding modes
# Rita Ann Roessner (ritaann.roessner@unil.ch)
# 2026-02-17

library(MixTCRviz)

# pca clusters
df = read.csv('df_pca_clusters.csv')

for (c in 0:1){
  tmp = df[df$cluster == c,c('TRAV','TRAJ','cdr3_TRA','TRBV','TRBJ','cdr3_TRB','model')]
  out_path = paste0("TSPs")
  print(out_path)
  MixTCRviz(input1=tmp, output.path = out_path)
}