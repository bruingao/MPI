#!/bin/sh

#Number of Points
b=1000000
#Number of Points Clusters
k=5

#Number of DNA Strands
d=1000000
#Number of DNA Clusters
c=5
#Length of DNA Strands
l=20

echo ********GENERATE $b INPUT POINTS EACH IN $k CLUSTERS
python ./DataGeneratorScripts/randomclustergen/generaterawdata.py -c $k -p $b -o ./input/cluster.csv

echo ********GENERATE $d INPUT DNA STRANDS OF LENGTH $l EACH IN $c CLUSTERS
python ./DataGeneratorScripts/randomclustergen/generatednadata.py -c $c -d $d -l $l -o ./input/dnacluster.csv
