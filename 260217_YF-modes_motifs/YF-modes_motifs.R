# Computing motifs for YF-binding modes
# Rita Ann Roessner (ritaann.roessner@unil.ch)
# 2026-02-17

library(MixTCRviz)

# pca clusters
df = read.csv('pca_CDR1A_cluster.csv')

for (c in 0:1){
  tmp = df[df$cluster == c,c('TRAV','TRAJ','cdr3_TRA','TRBV','TRBJ','cdr3_TRB')]
  out_path = paste0("out/", c)
  print(out_path)
  MixTCRviz(input1=tmp, output.path = out_path)
}

# outliers by AF3 confidence
df = read.csv('low_confidence_models.csv')
out_path = "out/low_confidence"

for (i in c(60,70,80)){
  tmp = df[df$lowest_conf < i,]
  file_name_out = paste0(out_path, 'plddt_smaller_than_', i)
  MixTCRviz(input1=tmp, output.path = out_path, filename.output=file_name_out)
}

# outliers by PCA coordinates
df = read.csv('outlier_models.csv')
out_path = "out/outlier"

tmp = df[df$outlier == 'yes',]
head(tmp)

MixTCRviz(input1=tmp, output.path = out_path)
