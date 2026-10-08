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

# additional TSP: major binding mode (cluster 0) with the CDR3 motif shown
# after SUBTRACTING the baseline repertoire (plot.cdr3.norm = 1)
maj = df[df$cluster == 0, c('TRAV','TRAJ','cdr3_TRA','TRBV','TRBJ','cdr3_TRB','model')]
MixTCRviz(input1 = maj,
          output.path    = "TSPs_major_cdr3subtract",
          plot.cdr3.norm = 1,            # 1 = subtract baseline from the CDR3 motif
          logo.type      = "probability")  # 'probability' = columns not scaled by information content (vs 'bits')

