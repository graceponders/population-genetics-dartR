Population genetics of Cherax robustus (sand yabbies) at four locations using dartSeq. 
Pleopod muscle tissue was collected non-lethally from two mainland and two K’gari sites, genotyped by Diversity Array Technology (DArT) using the Illumina HiSeq2500 platform.

All data was filtered and processed in R using dartRverse v1.0.6, which still has a few bugs (e.g. 'gl.report.hwe' from 'ggtern' for Hardy-Wienberg equilibrium has not been carried over yet). Script sections are: 
1. filtering
2. pca visualisation
3. sex-linked markers
4. HWE (using 'pegas')
5. Population diversity metrics (using 'hierfstat')
6. Pairwise genetic distance/time since divergence
7. Ancestry (using 'LEA')
8. Neighbour joining tree (using 'ape')
