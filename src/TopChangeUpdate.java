import java.util.*;
import java.io.*;
import corejava.*;


/**
* Topology break-points update kernel
*/

public class TopChangeUpdate {

    DCPSampler sampler;

    public TopChangeUpdate (DCPSampler s){

	sampler = s;
    }
    
    /**
     * Updates the locations of all topology break-points
     */

    protected void Update() {

	
	int lowerBound, upperBound , current, proposed,  landing_partition;
	Partition walk_part;

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

		double[] proposed_current_likelihood = {0,0}, prop_curr_small_like = {0,0};
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
			    
			    logRatio = proposed_current_likelihood[0] - proposed_current_likelihood[1];

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

			    
			    logRatio = proposed_current_likelihood[0] - proposed_current_likelihood[1];
			    
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
		double[] proposed_current_likelihood = {0,0}, prop_curr_small_like = {0,0};
			
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
		    
		    logRatio = proposed_current_likelihood[0] - proposed_current_likelihood[1];

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
		    
		    logRatio = proposed_current_likelihood[0] - proposed_current_likelihood[1];
		    
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




    protected void IntermediateTopLikelihoods(int start, int end, int p_tree, double[] part_likelihood, double[] p_c_like){
	Partition temp_part;
	
	for (int i = start; i < end; i++){
	    
	    temp_part = (Partition)sampler.partition_list.elementAt(i);
	    
	    part_likelihood[i] = sampler.PostTree[p_tree].LogLikelihood0(temp_part.cMatrix, temp_part.cHyperParameter, temp_part.counts, temp_part.data, temp_part.cMatrix.pi, sampler.handleGaps, 0, 0);
	    p_c_like[0] += part_likelihood[i];
	    
	    p_c_like[1] += temp_part.cPartialLogLikelihood;
	}
	
    }


    protected void UpdateTopCounts(Partition temp_part, int start, int end, int prop_tree, int curr_tree, 
				 int[] counts, double [] p_c_like, double[] p_c_small_like){
	
	for (int j = start; j < end; j++)
	    counts[sampler.indexSeq[j]]++;
	
	p_c_small_like[0] = sampler.PostTree[prop_tree].LogLikelihood0(temp_part.cMatrix, temp_part.cHyperParameter, counts, temp_part.data, temp_part.cMatrix.pi, sampler.handleGaps, 0, 0);
	p_c_like[0] += p_c_small_like[0];
	
	p_c_small_like[1] = sampler.PostTree[curr_tree].LogLikelihood0(temp_part.cMatrix, temp_part.cHyperParameter, counts, temp_part.data, temp_part.cMatrix.pi, sampler.handleGaps, 0, 0);
	p_c_like[1] += p_c_small_like[1];
	
    }


    protected int MoveTopXiToRight(int p_change, int p_p_change, int land_index, int[] counts,
				 int prop, double[] partial_likelihood, double[] p_c_small_like){
	
	Partition prev_part = (Partition)sampler.partition_list.elementAt(p_change);
	Partition prev_prev_part = (Partition)sampler.partition_list.elementAt(p_p_change);
	Partition land_part = (Partition)sampler.partition_list.elementAt(land_index);
	Partition walk_part;
       	
	int prev_tree = prev_part.cTree, new_p_change;

	/****** update trees and likelihoods on intermediate parameter partitions ********/

	// System.out.println("proposed: " + prop);
	
	for (int i = p_change; i < land_index; i++){
	    walk_part = (Partition)sampler.partition_list.elementAt(i);
	    walk_part.cTree = prev_prev_part.cTree;
	    walk_part.cPartialLogLikelihood = partial_likelihood[i];
	}
		
	if (prev_part.IsParameterChange()){/*** if moved point is topology and parameter change point ***/
	    
	   
	    /***** insert new topology point **********/
	    
	    if (land_part.left.nuc_position == prop){
		land_part.left.top_change = true;
		new_p_change = land_index;
	    }
	    else{
		
		sampler.InsertPartitionRight(land_part, prop, land_index, true, counts);
	    
		Partition new_part = (Partition)sampler.partition_list.elementAt(land_index + 1);
		
		new_part.cPartialLogLikelihood = land_part.cPartialLogLikelihood - p_c_small_like[1];
		land_part.cPartialLogLikelihood = p_c_small_like[0];
		
		new_part.cHyperParameter = land_part.cHyperParameter;
		new_part.cPartialLogHyperParameterPrior = land_part.cPartialLogHyperParameterPrior;
		new_part.log_alpha_prior = land_part.log_alpha_prior;
		new_part.cMatrix = land_part.cMatrix;
		
		new_part.cTree = prev_tree;
		land_part.cTree = prev_prev_part.cTree;
		new_p_change = land_index + 1;
	    }

	    /***** remove old topology point *****/    	    
	    
	    prev_part.left.top_change = false;

	}
	else{ /*** if moved point is topology, but not parameter change point ***/

	    Partition left_part = (Partition)sampler.partition_list.elementAt(p_change - 1);

	    if (p_change == land_index){
		
		/*** update change points ***/

		prev_part.left.nuc_position = prop;
		left_part.right = prop - 1;

		/*** update counts ***/

		for (int i = 0; i < counts.length; i++){
		    if (counts[i] > 0){
			prev_part.counts[i] -=  counts[i];
			left_part.counts[i] += counts[i];
		    }
		}

		/*** update likelihoods ***/

		prev_part.cPartialLogLikelihood -= p_c_small_like[1];
		left_part.cPartialLogLikelihood += p_c_small_like[0];


		new_p_change = p_change;

	    }
	    else{// p_change != land_index

		 /***** insert new topology point **********/
	    
		if (land_part.left.nuc_position == prop){
		    land_part.left.top_change = true;
		    new_p_change = land_index - 1;
		}
		else{
		    
		    sampler.InsertPartitionRight(land_part, prop, land_index, true, counts);
		    
		    Partition new_part = (Partition)sampler.partition_list.elementAt(land_index + 1);
		    
		    new_part.cPartialLogLikelihood = land_part.cPartialLogLikelihood - p_c_small_like[1];
		    land_part.cPartialLogLikelihood = p_c_small_like[0];
		    
		    new_part.cHyperParameter = land_part.cHyperParameter;
		    new_part.cPartialLogHyperParameterPrior = land_part.cPartialLogHyperParameterPrior;
		    new_part.log_alpha_prior = land_part.log_alpha_prior;
		    new_part.cMatrix = land_part.cMatrix;
		    
		    new_part.cTree = prev_tree;
		    land_part.cTree = prev_prev_part.cTree;

		    new_p_change = land_index;
		}
		
				
		left_part.right = prev_part.right;
		
		/*** update counts ***/

		for (int i = 0; i < left_part.counts.length; i++){
		    if (prev_part.counts[i] > 0){
			left_part.counts[i] += prev_part.counts[i];
		    }
		}

		/*** update likelihood ***/

		left_part.cPartialLogLikelihood += partial_likelihood[p_change];
	    

		/***** remove old topology point *****/ 

		sampler.partition_list.removeElementAt(p_change);
		
	    }
	    	    
	}

	return new_p_change;
    }

    


    protected int MoveTopXiToLeft(int p_change, int p_p_change, int land_index, int[] counts,
			  int prop, double[] partial_likelihood, double[] p_c_small_like){

	Partition prev_part = (Partition)sampler.partition_list.elementAt(p_change);
	Partition prev_prev_part = (Partition)sampler.partition_list.elementAt(p_p_change);
	Partition land_part = (Partition)sampler.partition_list.elementAt(land_index);
	Partition walk_part;
	
	int new_p_change;

	// System.out.println("proposed: " + prop);

	/****** update trees and likelihoods on intermediate parameter partitions ********/
	
	for (int i = land_index + 1; i < p_change; i++){
	    walk_part = (Partition)sampler.partition_list.elementAt(i);
	    walk_part.cTree = prev_part.cTree;
	    walk_part.cPartialLogLikelihood = partial_likelihood[i];
	}
       
	if (prev_part.IsParameterChange()){/*** if moved point is topology and parameter change point ***/
	    
	   
	    /***** insert new topology point **********/
	    
	    if (land_part.left.nuc_position == prop){
		land_part.left.top_change = true;
		land_part.cTree = prev_part.cTree;
		land_part.cPartialLogLikelihood = p_c_small_like[0];
		new_p_change = land_index;
	    }
	    else{
		
		sampler.InsertPartitionLeft(land_part, prop, land_index, true, counts);
	    
		Partition new_part = (Partition)sampler.partition_list.elementAt(land_index + 1);
		
		new_part.cPartialLogLikelihood = p_c_small_like[0];
		land_part.cPartialLogLikelihood = land_part.cPartialLogLikelihood - p_c_small_like[1];
		
		new_part.cHyperParameter = land_part.cHyperParameter;
		new_part.cPartialLogHyperParameterPrior = land_part.cPartialLogHyperParameterPrior;
		new_part.log_alpha_prior = land_part.log_alpha_prior;
		new_part.cMatrix = land_part.cMatrix;
		
		new_part.cTree = prev_part.cTree;
		
		new_p_change = land_index + 1;
	    }

	    /***** remove old topology point *****/    	    
	    
	    prev_part.left.top_change = false;

	}
	else{ /*** if moved point is topology, but not parameter change point ***/

	    Partition left_part = (Partition)sampler.partition_list.elementAt(p_change - 1);

	    if (p_change - 1 == land_index){
		
		if (left_part.left.nuc_position == prop){
		    left_part.left.top_change = true;
		    
		    left_part.right = prev_part.right;
		
		    /*** update counts ***/
		    
		    for (int i = 0; i < left_part.counts.length; i++){
			if (prev_part.counts[i] > 0){
			    left_part.counts[i] += prev_part.counts[i];
			}
		    }

		    /*** update likelihood ***/

		    left_part.cPartialLogLikelihood = p_c_small_like[0] + prev_part.cPartialLogLikelihood;

		    /*** update tree ***/
	       
		    left_part.cTree = prev_part.cTree;
		    

		    /***** remove old topology point *****/ 

		    sampler.partition_list.removeElementAt(p_change);

		    
		    new_p_change = p_change - 1;
		}
		else{
		    /*** update change points ***/
		    
		    prev_part.left.nuc_position = prop;
		    left_part.right = prop - 1;
		    
		    /*** update counts ***/
		    
		    for (int i = 0; i < counts.length; i++){
			if (counts[i] > 0){
			    prev_part.counts[i] += counts[i];
			    left_part.counts[i] -= counts[i];
			}
		    }

		    /*** update likelihoods ***/

		    prev_part.cPartialLogLikelihood += p_c_small_like[0];
		    left_part.cPartialLogLikelihood -= p_c_small_like[1];
		    
		    new_p_change = p_change;
		}

	    }
	    else{// p_change != land_index

		int removal_index;
		 /***** insert new topology point **********/

		if (land_part.left.nuc_position == prop){
		    land_part.left.top_change = true;
		    land_part.cTree = prev_part.cTree;
		    land_part.cPartialLogLikelihood = p_c_small_like[0];
		    removal_index = p_change;

		    new_p_change = land_index;
		}
		else{
		    
		    sampler.InsertPartitionLeft(land_part, prop, land_index, true, counts);
		    
		    Partition new_part = (Partition)sampler.partition_list.elementAt(land_index + 1);
		    
		    new_part.cPartialLogLikelihood = p_c_small_like[0];
		    land_part.cPartialLogLikelihood = land_part.cPartialLogLikelihood - p_c_small_like[1];
		    
		    new_part.cHyperParameter = land_part.cHyperParameter;
		    new_part.cPartialLogHyperParameterPrior = land_part.cPartialLogHyperParameterPrior;
		    new_part.log_alpha_prior = land_part.log_alpha_prior;
		    new_part.cMatrix = land_part.cMatrix;
		    
		    new_part.cTree = prev_part.cTree;
		    removal_index = p_change + 1;
		    
		    new_p_change = land_index + 1;
		}
		
				
		left_part.right = prev_part.right;
		
		/*** update counts ***/

		for (int i = 0; i < left_part.counts.length; i++){
		    if (prev_part.counts[i] > 0){
			left_part.counts[i] += prev_part.counts[i];
		    }
		}

		/*** update likelihood ***/

		left_part.cPartialLogLikelihood = partial_likelihood[p_change - 1] + prev_part.cPartialLogLikelihood;

		/*** update tree ***/
	       
		left_part.cTree = prev_part.cTree;
	    

		/***** remove old topology point *****/ 

		sampler.partition_list.removeElementAt(removal_index);

		
	    }

	    	    
	}
	
	return new_p_change;
	
    }
    


}
