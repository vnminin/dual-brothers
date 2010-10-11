import java.io.*;


/**
* Class Partition provides basic structure for a partition/segment between change-points
* in a multiple change point model
*/

public class Partition  implements Serializable {

    /**
     * left end-point of a partition/segment - CHANGE-POINT
     */
    public BreakPoint  left;


    /**
     * right end-point of a partition/segment - NOT A CHANGE-POINT
     */
    public int right;



    /**
     * Sufficient statistic (l) of observed nucleotide patterns for a partition/segment.  l indexes pattern.
     */
    public int[] counts;


    /**
     * Observed nucleotide patterns ([l][i]) for a partition/segment.  l indexes pattern. i indexes taxon.
     */
    public int[][] data;


    /**
     * Current log likelihood of a partition.
     */
    public double cPartialLogLikelihood;

    
    /**
     * Current expected divergence \mu in a partition/segment.
     */
    public double cHyperParameter; 
    

    /**
     * Current log prior density of \mu  in a partition/segment.
     */
    public double cPartialLogHyperParameterPrior;


    /**
     * Current log prior density of \alhpa in a partition/segment.
     */
    public double log_alpha_prior;

  
    /**
     * Current tree number of a partition/segment.
     */
    public int cTree; // Holds the current model.  Should change tree to an abstract non-specific class
	

    /**
     * Current evolutionary matrix  of a partition/segment.
     */
    public QMatrix cMatrix; // Holds the current transition matrix.
   

    Partition(int l, int r){

	left = new BreakPoint(l, true, true);
	right = r;
    }


    Partition(int l, int r, boolean t, boolean p){

	left = new BreakPoint(l,t,p);
	right = r;
    }


    public boolean IsTopologyChange(){
	
	if (left.top_change)
	    return true;
	else 
	    return false;
    }




    public boolean IsParameterChange(){
	
	if (left.par_change)
	    return true;
	else
	    return false;
    }
    
    
    
}



 



