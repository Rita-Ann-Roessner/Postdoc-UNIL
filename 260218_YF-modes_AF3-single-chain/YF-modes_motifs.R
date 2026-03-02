# Computing motifs for YF-binding modes
# Rita Ann Roessner (ritaann.roessner@unil.ch)
# 2026-02-17

library(MixTCRviz)

df = read.csv('chainA_rmsd_low.csv')
MixTCRviz(input1=df, output.path = './', filename.output='chainA_rmsd_low.pdf')

df = read.csv('chainA_rmsd_high.csv')
MixTCRviz(input1=df, output.path = './', filename.output='chainA_rmsd_high.pdf')
