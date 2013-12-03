import java.io.*;
import java.util.*;
import java.lang.Math;

public class Cluster2D {

    public static void main(String[] argv) throws Exception {
        
        String inputFile=argv[0];
        int k_param=Integer.parseInt(argv[1]);
        
        ArrayList<ArrayList> points = new ArrayList<ArrayList>();
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
        
        ArrayList<ArrayList> centroids = new ArrayList<ArrayList>();
        for (int i=0;i<k_param;i++){
            centroids.add(points.get(i));
        }
        
        
        double threshold=1.0e-6;
        double incre=1.0e+6;
        while(incre>threshold){
        
            ArrayList<Float> tagX=new ArrayList<Float>();
            ArrayList<Float> tagY=new ArrayList<Float>();
            ArrayList<Integer> tagSize=new ArrayList<Integer>();
            for (int i=0;i<k_param;i++){
                tagX.add(0.0f);
                tagY.add(0.0f);
                tagSize.add(0);
            }
            ArrayList<Integer> tags=new ArrayList<Integer>();
            for (int i=0;i<points.size();i++){
                ArrayList<Float> dists=new ArrayList<Float>();
                for (int j=0;j<k_param;j++){
                    dists.add(dist(points.get(i),centroids.get(j)));
                }
                int tag=dists.indexOf(Collections.min(dists));
                tags.add(tag);
                float sumX=tagX.get(tag)+(float)points.get(i).get(0);
                float sumY=tagY.get(tag)+(float)points.get(i).get(1);
                int sumSize=tagSize.get(tag)+1;
                tagX.set(tag,sumX);
                tagY.set(tag,sumY);
                tagSize.set(tag,sumSize);
            }
            incre=0.0;
            for(int i=0;i<k_param;i++){
                float meanX=tagX.get(i)/tagSize.get(i);
                float meanY=tagY.get(i)/tagSize.get(i);
                tagX.set(i,meanX);
                tagY.set(i,meanY);
                double tempX2 = Math.pow((double)(meanX-(float)centroids.get(i).get(0)),2);
                double tempY2 = Math.pow((double)(meanY-(float)centroids.get(i).get(1)),2);
                incre = incre + Math.sqrt(tempX2+tempY2);
                ArrayList<Float> updatedCentroid=new ArrayList<Float>();
                updatedCentroid.add(meanX);
                updatedCentroid.add(meanY);
                centroids.set(i,updatedCentroid);
            }
            System.out.println(tags);
        }
        
    }
    
    public static Float dist(ArrayList<Float> p1, ArrayList<Float> p2) {
        double xDiff=p1.get(0)-p2.get(0);
        double yDiff=p1.get(1)-p2.get(1);
        float dis=(float)Math.sqrt(Math.pow(xDiff,2)+Math.pow(yDiff,2));
        return dis;
    }
    
}

