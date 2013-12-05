#!/bin/sh


echo "k-means clustering of 2D points"
javac Cluster2D.java
java Cluster2D ../DataGeneratorScripts/input/cluster.csv 2

echo "k-means clustering of DNA strands"
javac ClusterDNA.java
java ClusterDNA ../DataGeneratorScripts/input/dnacluster.csv 2
