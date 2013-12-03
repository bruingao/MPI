import java.io.*;
import java.util.*;
import java.lang.Math;

public class Cluster2D {

    public static void main(String[] argv) throws Exception {
        
        String inputFile=argv[0];
        int k_param=Integer.parseInt(argv[1]);
        
        // scan input file and store the 2D points
        ArrayList<ArrayList<Float>> points = new ArrayList<ArrayList<Float>>();
        Scanner aScanner = new Scanner(new File(inputFile));
        String aLine = null;
        do {
            aLine = aScanner.nextLine();
            String[] parts = aLine.split(",");
            ArrayList<Float> xy=new ArrayList<Float>();
            xy.add(Float.parseFloat(parts[0]));
            xy.add(Float.parseFloat(parts[1]));
            points.add(xy);
        } while (aScanner.hasNext());
        
        
        // choose initial centroids
        // use either one of the two strategies below
        
        // 1
        float minDist=0.5f; // set it smaller than minDistance in generaterawdata.py
        ArrayList<ArrayList<Float>> centroids = new ArrayList<ArrayList<Float>>();
        int randInd=(int)(Math.random() * points.size());
        centroids.add(points.get(randInd));
        for (int i=0;i<k_param;i++){
            while(tooClose(points.get(randInd),centroids,minDist)){
                randInd=(int)(Math.random() * points.size());
            }
            centroids.add(points.get(randInd));
        }
        
        /*
        // 2
        ArrayList<ArrayList<Float>> centroids = new ArrayList<ArrayList<Float>>();
        for (int i=0;i<k_param;i++){
            centroids.add(points.get(i));
        }
        */
        
        // iterate until centroids are stable
        double threshold=1.0e-6;
        double incre=1.0e+6;
        while(incre>threshold){
            
            /*
            ArrayList<Float> tagX=new ArrayList<Float>();
            ArrayList<Float> tagY=new ArrayList<Float>();
            ArrayList<Integer> tagSize=new ArrayList<Integer>();
            for (int i=0;i<k_param;i++){
                tagX.add(0.0f);
                tagY.add(0.0f);
                tagSize.add(0);
            }
            */
            float[][] tagXY = new float[k_param][2];
            int[] tagSize = new int[k_param];
            ArrayList<Integer> tags=new ArrayList<Integer>();
            
            for (int i=0;i<points.size();i++){
                ArrayList<Float> dists=new ArrayList<Float>();
                for (int j=0;j<k_param;j++){
                    dists.add(dist(points.get(i),centroids.get(j)));
                }
                int tag=dists.indexOf(Collections.min(dists));
                tags.add(tag);
                //float sumX=tagX.get(tag)+points.get(i).get(0);
                //float sumY=tagY.get(tag)+points.get(i).get(1);
                //int sumSize=tagSize.get(tag)+1;
                //tagX.set(tag,sumX);
                //tagY.set(tag,sumY);
                //tagSize.set(tag,sumSize);
                tagXY[tag][0] += points.get(i).get(0);
                tagXY[tag][1] += points.get(i).get(1);
                tagSize[tag] += 1;
            }
            incre=0.0;
            for(int i=0;i<k_param;i++){
                //float meanX=tagX.get(i)/tagSize.get(i);
                //float meanY=tagY.get(i)/tagSize.get(i);
                //tagX.set(i,meanX);
                //tagY.set(i,meanY);
                tagXY[i][0] /= tagSize[i];
                tagXY[i][1] /= tagSize[i];
                double tempX2 = Math.pow((double)(tagXY[i][0]-centroids.get(i).get(0)),2);
                double tempY2 = Math.pow((double)(tagXY[i][1]-centroids.get(i).get(1)),2);
                incre = incre + Math.sqrt(tempX2+tempY2);
                ArrayList<Float> updatedCentroid=new ArrayList<Float>();
                updatedCentroid.add(tagXY[i][0]);
                updatedCentroid.add(tagXY[i][1]);
                centroids.set(i,updatedCentroid);
            }
            System.out.println(tags);
        }
        
    }
    
    
    // compute Euclidean distance between two points
    public static Float dist(ArrayList<Float> p1, ArrayList<Float> p2) {
        double xDiff=p1.get(0)-p2.get(0);
        double yDiff=p1.get(1)-p2.get(1);
        float dis=(float)Math.sqrt(Math.pow(xDiff,2)+Math.pow(yDiff,2));
        return dis;
    }
    
    
    // decide if one point is too close to any of the points in a list
    public static boolean tooClose(ArrayList<Float> one, ArrayList<ArrayList<Float>> list, Float minDist){
        for (int iter=0;iter<list.size();iter++){
            if(dist(one,list.get(iter))<minDist){
                return true;
            }
        }
        return false;
    }
    
}

