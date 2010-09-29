import java.util.*;
import java.io.*;
import corejava.*;


/**
* Topology break-points update kernel that can incorporate prior information
* Redefines the Update routine of the original TopChangeUpdate class
*/

public class InfTopChange extends TopChangeUpdate{

    public InfTopChange (DCPSampler s){

	super(s);
    }

    
    /**
     * Updates the locations of all topology break-points
     */

    protected void Update() {

	int lowerBound, upperBound , current, proposed,  landing_partition;
	Partition walk_part;
	double[] proposed_current_likelihood = {0,0}, prop_curr_small_like = {0,0};
	double logRatio =  0;

	

	/***************** start with the first partition ******************/
	
	Partition curr_part = (Partition)sampler.partition_list.elementAt(0);
	Partition second_part = curr_part;


	
	/************** find the second topology change point **************/
	
	int k = 1;
	while((!second_part.IsTopologyChange() || (k == 1)) && (k < sampler.partition_list.size())){
	    second_part = (Partition)sampler.partition_list.elementAt(k);
	    k++;
	}

	/************* if there is one topology change point *********
	 ************* start updating                        *********/ 

	boolean topology_change = second_part.IsTopologyChange() && (k != 1);

	
	if (topology_change){
	    
	    Partition prev_prev_part, prev_part;
	    
	    int prev_ch_pnt = k - 1, prev_prev_ch_pnt = 0, list_size = 0;

	    prev_prev_part = (Partition)sampler.partition_list.elementAt(prev_prev_ch_pnt);
	    prev_part = (Partition)sampler.partition_list.elementAt(prev_ch_pnt);
	
	    for(int t = k; t < sampler.partition_list.size(); t++) {

		double[] pPartialLogLikelihood = new double[sampler.partition_list.size()];

		list_size = sampler.partition_list.size();
		curr_part = (Partition)sampler.partition_list.elementAt(t);

		if (curr_part.IsTopologyChange()){
		    
		    current = prev_part.left.nuc_position;
		    
		    /***** determine the bounds for the proposed value *****/

		    lowerBound = prev_prev_part.left.nuc_position + 1;
		    upperBound = curr_part.left.nuc_position - 1;
		    
		    if (upperBound - lowerBound > 0) {

			proposed = sampler.ProposeNewXi(current, lowerBound, upperBound);
		
			/*****  Now proposed is symmetric and reflected and does not equal the original value *****/

			int testcounts[] = new int[curr_part.counts.length];

		
			if (proposed > current) {

			    /****** find the segment where proposed point landed *******/
			
			    landing_partition = sampler.LandingPart(proposed, prev_ch_pnt, t - 1);


			    /****** sum intermediate likelihoods of parameter partitions under both models *****/
			    
			    IntermediateTopLikelihoods(prev_ch_pnt, landing_partition, prev_prev_part.cTree, 
						    pPartialLogLikelihood, proposed_current_likelihood);

			    /***** add sites to classes before change point ******/
			    
			    walk_part = (Partition)sampler.partition_list.elementAt(landing_partition);

			    if (walk_part.left.nuc_position != proposed){
				
				UpdateTopCounts(walk_part, walk_part.left.nuc_position, proposed, prev_prev_part.cTree,
					     prev_part.cTree, testcounts, proposed_current_likelihood,
					     prop_curr_small_like);
			
			    }
			    
			    logRatio = proposed_current_likelihood[0] - proposed_current_likelihood[1]
				+ sampler.prior_inf.recombination[proposed]
				- sampler.prior_inf.recombination[current];

			    sampler.tries[2]++;

			    if (sampler.LogMHAccept(logRatio, sampler.set.unif.nextDouble())){
				prev_ch_pnt = MoveTopXiToRight(prev_ch_pnt, prev_prev_ch_pnt, landing_partition, 
							       testcounts, proposed, pPartialLogLikelihood, 
							       prop_curr_small_like);				    
				sampler.acceptancerate[2]++;
			    }
			}
			    
			else{// proposed < current
			    
			    /****** find the segment where proposed point landed *******/
			    
			    landing_partition = sampler.LandingPart(proposed, prev_prev_ch_pnt, prev_ch_pnt);
			    
			    /****** sum intermediate likelihoods of parameter partitions under both models *****/
			    
			    IntermediateTopLikelihoods(landing_partition + 1, prev_ch_pnt, prev_part.cTree, 
						       pPartialLogLikelihood, proposed_current_likelihood);
			    
			    /***** add sites to classes before change point ******/
			    
			    walk_part = (Partition)sampler.partition_list.elementAt(landing_partition);
			    

			    UpdateTopCounts(walk_part, proposed, walk_part.right + 1, prev_part.cTree,
					    prev_prev_part.cTree, testcounts, proposed_current_likelihood,
					    prop_curr_small_like);

			    
			    logRatio = proposed_current_likelihood[0] - proposed_current_likelihood[1]
				+ sampler.prior_inf.recombination[proposed]
				- sampler.prior_inf.recombination[current];
			    
			    sampler.tries[2]++;

			    if (sampler.LogMHAccept(logRatio, sampler.set.unif.nextDouble())){
				prev_ch_pnt = MoveTopXiToLeft(prev_ch_pnt, prev_prev_ch_pnt, landing_partition, 
							      testcounts,  proposed, pPartialLogLikelihood, 
							      prop_curr_small_like);
				sampler.acceptancerate[2]++;
				
			    }
			    
			}
			
		    }
		    
		
		
		    /***** reset topology change points *******/
		    prev_prev_ch_pnt = prev_ch_pnt;
		    t += sampler.partition_list.size() - list_size;
		    prev_ch_pnt = t;
		    prev_prev_part = (Partition)sampler.partition_list.elementAt(prev_prev_ch_pnt);
		    prev_part = (Partition)sampler.partition_list.elementAt(prev_ch_pnt);
		    

		}
	    }
	    
	    
	    /***** finish updating on the right end *******/
	    
	    current = prev_part.left.nuc_position;
	    
	    /***** determine the bounds for the proposed value *****/
	    
	    lowerBound = prev_prev_part.left.nuc_position + 1;
	    upperBound = sampler.lenSeq - 1;
	    
	    if (upperBound - lowerBound > 0) {

		double[] pPartialLogLikelihood = new double[sampler.partition_list.size()];
		
		proposed = sampler.ProposeNewXi(current, lowerBound, upperBound);
		
		/*****  Now proposed is symmetric and reflected and does not equal the original value *****/
		
		int testcounts[] = new int[curr_part.counts.length];
		
		
		if (proposed > current) {
		    
		    /****** find the segment where proposed point landed *******/
		    
		    landing_partition = sampler.LandingPart(proposed, prev_ch_pnt, sampler.partition_list.size() - 1);
		    
		    /****** sum intermediate likelihoods of parameter partitions under both models *****/
		
		    IntermediateTopLikelihoods(prev_ch_pnt, landing_partition, prev_prev_part.cTree, 
					       pPartialLogLikelihood, proposed_current_likelihood);

		    
		    /***** add sites to classes before change point ******/
		    
		    walk_part = (Partition)sampler.partition_list.elementAt(landing_partition);
		    
		    if (walk_part.left.nuc_position != proposed){
			
			UpdateTopCounts(walk_part, walk_part.left.nuc_position, proposed, prev_prev_part.cTree,
					prev_part.cTree, testcounts, proposed_current_likelihood,
					prop_curr_small_like);
		    }
		    
		    logRatio = proposed_current_likelihood[0] - proposed_current_likelihood[1]
			+ sampler.prior_inf.recombination[proposed]
			- sampler.prior_inf.recombination[current];

		    sampler.tries[2]++;
		    
		    if (sampler.LogMHAccept(logRatio, sampler.set.unif.nextDouble())){
			MoveTopXiToRight(prev_ch_pnt, prev_prev_ch_pnt, landing_partition, testcounts, 
					 proposed, pPartialLogLikelihood, prop_curr_small_like);

			sampler.acceptancerate[2]++;
		    }
		}
		
		else{// proposed < current
		    
		    /****** find the segment where proposed point landed *******/
		    
		    landing_partition = sampler.LandingPart(proposed, prev_prev_ch_pnt, prev_ch_pnt);
		    
		    /****** sum intermediate likelihoods of parameter partitions under both models *****/
		    
		    IntermediateTopLikelihoods(landing_partition + 1, prev_ch_pnt, prev_part.cTree, 
					       pPartialLogLikelihood, proposed_current_likelihood);
		    
		    
		    /***** add sites to classes before change point ******/
		    
		    walk_part = (Partition)sampler.partition_list.elementAt(landing_partition);
		    
		    UpdateTopCounts(walk_part, proposed, walk_part.right + 1, prev_part.cTree,
				    prev_prev_part.cTree, testcounts, proposed_current_likelihood,
				    prop_curr_small_like);
		    
		    logRatio = proposed_current_likelihood[0] - proposed_current_likelihood[1]
			+ sampler.prior_inf.recombination[proposed]
			- sampler.prior_inf.recombination[current];
		    
		    sampler.tries[2]++;

		    if (sampler.LogMHAccept(logRatio, sampler.set.unif.nextDouble())){
			MoveTopXiToLeft(prev_ch_pnt, prev_prev_ch_pnt, landing_partition, testcounts, 
					proposed, pPartialLogLikelihood, prop_curr_small_like);
			sampler.acceptancerate[2]++;
		    }
		    
		}
	    
	    }
	    
	}
	
    }

}


