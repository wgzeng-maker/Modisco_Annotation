# Modisco_ClusterAnalysis
Downstream analysis pipeline on Jacob Schreiber's Modisco-lite
Modisco-Lite generates seqlets clusters (called Patterns) and match each cluster to known mofit from the public motif dataset JASPAR. It also generate subclusters (called Subpattern) within each cluster. This repository aims to analyze and visualize the clusters and subclusters heterogeneity.
Umap to visulaize the Patterns:
<img width="1950" height="1650" alt="all_patterns_seqlet_umap" src="https://github.com/user-attachments/assets/9f1141e5-e951-42fb-9bfd-ba17eb70b737" />
Umap to visulaize the Subpatterns:
<img width="1500" height="1350" alt="pattern0_seqlet_umap" src="https://github.com/user-attachments/assets/9a4785c6-9254-430c-86cc-feffccdec3df" />
To understand what drives the heterogeneity in a Pattern, I prepare two kinds of plot, each answer different questions.
For the first kind
The plot on the left ask, what is the sequence difference within a Pattern?
The plot on the right ask, what is the attribution difference within a Pattern?
<img width="2400" height="1050" alt="design1_seq_vs_attr_umap" src="https://github.com/user-attachments/assets/60138751-0b99-4fce-a2c5-dfd8071defee" />
For the second kind
The plot ask similar question as the the first plot, but within one plot
The X axis is sequence similarity to cluster center, which measure how similar a seqlet's sequence is to the cluster center
The y axis is the attribution similarity to cluster center, which measure how similar a seqlet's attribution is to the cluster center
<img width="1350" height="1350" alt="design2_seq_vs_attr_axes" src="https://github.com/user-attachments/assets/42362391-d0bc-405a-9f8c-1e4117ff712b" />
