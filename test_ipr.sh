##!/usr/bash
alias interpro=/mlab/data/software/interproscan/interproscan.sh

## generating interproscan results
interpro -i ../output/non_redundant_zinc.fa -f tsv -iprlookup -pa -goterms -o ../output/non_redundant_zinc.ipr.tsv
echo "interproscan done"
date


# do functional characterization
./functional_characterization.R
echo "functional distances done"
date

## Calculate structural distances and all four measures vs. k for determining optimal k 
./measures.vs.k.R
echo "comparing functional and structural distances done"
date

## Graph the four measure on the same plot
./cluster_metric_graphs.R
echo "graphing clustering done"
date

## Print the results of cluster centers, counts, and average probabilities
./cluster.centers.prob.R
echo "data out"
date

## Graph hierarchical clustering of structural and functional cluster distances.
./hierarchical.clustering.R
echo "graphing other stuff"
date

## EC and IPR enrichment
./ec_ipr_enrichment.R
echo "doing ec and ipr enrichment"
date

echo "All done! Check ../output for files!"
