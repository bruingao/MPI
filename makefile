make:
	mpicc mpiKmeans/*.c -o mpiKmeans/mpiKmeans -Wl,-rpath /usr/local/lib/openmpi/lib/
	javac SeqCluster/*.java
clean:
	rm mpiKmeans/mpiKmeans
	rm SeqCluster/*.class
	rm input/*
