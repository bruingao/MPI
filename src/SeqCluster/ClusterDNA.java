import java.io.*;
import java.util.*;
import java.lang.Math;

public class ClusterDNA {

    public static void main(String[] argv) throws Exception {
        
        String inputFile=argv[0];
        int k_param=Integer.parseInt(argv[1]);
        
        // scan input file and store the DNA strands
        ArrayList<ArrayList<String>> DNAs = new ArrayList<ArrayList<String>>();
        Scanner aScanner = new Scanner(new File(inputFile));
        String aLine = null;
        do {
            aLine = aScanner.nextLine();
            String[] parts = aLine.split(",");
            ArrayList<String> aDNA=new ArrayList<String>();
            for (int i=0;i<parts.length;i++){
                aDNA.add(parts[i]);
            }
            DNAs.add(aDNA);
        } while (aScanner.hasNext());
        int lenDNA=DNAs.get(0).size();
        
        
        // choose initial centroids
        // use either one of the two strategies below
        
        // 1
        int maxSim=(int)(0.8f*lenDNA); // set it greater than maxSim in generatednadata.py
        ArrayList<ArrayList<String>> centroids = new ArrayList<ArrayList<String>>();
        int randInd=(int)(Math.random() * DNAs.size());
        //System.out.println(randInd);
        centroids.add(DNAs.get(randInd));
        for (int i=0;i<k_param;i++){
            while(tooSimilar(DNAs.get(randInd),centroids,maxSim)){
                randInd=(int)(Math.random() * DNAs.size());
                //System.out.println(randInd);
            }
            centroids.add(DNAs.get(randInd));
        }
        //System.out.println("ok");
        
        /*
        // 2
        ArrayList<ArrayList<String>> centroids = new ArrayList<ArrayList<String>>();
        for (int i=0;i<k_param;i++){
            centroids.add(DNAs.get(i));
        }
        */
        
        // iterate until centroids are stable
        boolean converged=false;
        while(!converged){
            
            int[][] tagCountA = new int[k_param][lenDNA];
            int[][] tagCountC = new int[k_param][lenDNA];
            int[][] tagCountG = new int[k_param][lenDNA];
            int[][] tagCountT = new int[k_param][lenDNA];
            
            ArrayList<Integer> tags=new ArrayList<Integer>();
            for (int i=0;i<DNAs.size();i++){
                ArrayList<Integer> sims=new ArrayList<Integer>();
                for (int j=0;j<k_param;j++){
                    sims.add(simDNA(DNAs.get(i),centroids.get(j)));
                }
                int tag=sims.indexOf(Collections.max(sims));
                tags.add(tag);
                for(int j=0;j<lenDNA;j++){
                    if(DNAs.get(i).get(j).equals("A")){
                        tagCountA[tag][j] += 1;
                    }else if(DNAs.get(i).get(j).equals("C")){
                        tagCountC[tag][j] += 1;
                    }else if(DNAs.get(i).get(j).equals("G")){
                        tagCountG[tag][j] += 1;
                    }else if(DNAs.get(i).get(j).equals("T")){
                        tagCountT[tag][j] += 1;
                    }
                }
            }
            
            converged=true;
            for(int i=0;i<k_param;i++){
                ArrayList<String> updatedCentroid=new ArrayList<String>();
                for(int j=0;j<lenDNA;j++){
                    if(tagCountA[i][j]>=tagCountC[i][j] && tagCountA[i][j]>=tagCountG[i][j] && tagCountA[i][j]>=tagCountT[i][j]){
                        updatedCentroid.add("A");
                    }else if(tagCountC[i][j]>=tagCountA[i][j] && tagCountC[i][j]>=tagCountG[i][j] && tagCountC[i][j]>=tagCountT[i][j]){
                        updatedCentroid.add("C");
                    }else if(tagCountG[i][j]>=tagCountA[i][j] && tagCountG[i][j]>=tagCountC[i][j] && tagCountG[i][j]>=tagCountT[i][j]){
                        updatedCentroid.add("G");
                    }else if(tagCountT[i][j]>=tagCountA[i][j] && tagCountT[i][j]>=tagCountC[i][j] && tagCountT[i][j]>=tagCountG[i][j]){
                        updatedCentroid.add("T");
                    }
                    if(! updatedCentroid.get(j).equals(centroids.get(i).get(j))){
                        converged=false;
                    }
                }
                centroids.set(i,updatedCentroid);
            }
            
            if(converged==true){
                for (int i=0;i<k_param;i++){
                    System.out.println("centroid "+i+":");
                    System.out.println(centroids.get(i));
                    System.out.println("points:");
                    for (int j=0;j<DNAs.size();j++){
                        if(tags.get(j)==i){
                            System.out.println(DNAs.get(j));
                        }
                    }
                }
            }
            
            //System.out.println(tags);
        }
        
    }
    
    
    // compute similarity between two DNA strands
    public static Integer simDNA(ArrayList<String> p1, ArrayList<String> p2) {
        int dlen=p1.size();
        int sim=dlen;
        for(int ii=0;ii<dlen;ii++){
            if(!p1.get(ii).equals(p2.get(ii))){
                sim--;
            }
        }
        return sim;
    }
    
    
    // decide if one DNA strand is too similar to any of the DNA strands in a list
    public static boolean tooSimilar(ArrayList<String> one, ArrayList<ArrayList<String>> list, Integer maxSim){
        for (int iter=0;iter<list.size();iter++){
            if(simDNA(one,list.get(iter)) > maxSim){
                return true;
            }
        }
        return false;
    }
    
}

