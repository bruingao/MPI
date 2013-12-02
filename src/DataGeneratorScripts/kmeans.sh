#!/bin/sh

#Number of Points
b=10
#Number of Points Clusters
k=2

#Number of DNA Strands
d=10
#Number of DNA Clusters
c=2
#Length of DNA Strands
l=20

echo ********GENERATE $b INPUT POINTS EACH IN $k CLUSTERS
python ./randomclustergen/generaterawdata.py -c $k -p $b -o input/cluster.csv

echo ********GENERATE $d INPUT DNA STRANDS OF LENGTH $l EACH IN $c CLUSTERS
python ./randomclustergen/generatednadata.py -c $c -d $d -l $l -o input/dnacluster.csv