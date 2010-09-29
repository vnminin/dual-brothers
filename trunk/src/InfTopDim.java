import java.util.*;
import java.io.*;
import corejava.*;


/**
* Adding and deleting new topology break-points kernel using informative priors
* Redefines AddOne; AddTwo and Delete functions
*/

public class InfTopDim extends TopDimChange{

    
    public InfTopDim (DCPSampler s){

	super(s);
    }
 
    /**
     * Adds a new topology break-point to the parameter space
     */

    protected void AddOne(){
	
	int proposed; // candidate for new topology change point
	int land_index, left_top_change, right_top_change;
	

	proposed = ProposeTopChangePointToInsert();
	
	
	/**** find partition, where new point landed ****/

	land_index = sampler.LandingPart(proposed, 0, sampler.partition_list.size() - 1);
	
	/**** find left and right parameter change point ****/

	left_top_change = FindLeftTopChange(land_index);
	right_top_change = FindRightTopChange(land_index);

	Partition land_part = (Partition)sampler.partition_list.elementAt(land_index);
	Partition left_part = land_part;

	if (left_top_change > 0){
	    left_part = (Partition)sampler.partition_list.elementAt(left_top_change - 1);
	}

	Partition right_part = land_part;

	if (right_top_change < sampler.partition_list.size()){
	    right_part = (Partition)sampler.partition_list.elementAt(right_top_change);
	}

	/**** compute counts on both sides of proposed new change point ****/

	int[] left_counts =  new int[land_part.counts.length];
	int[] right_counts =  new int[land_part.counts.length];

	if (land_part.left.nuc_position != proposed){
	    for (int i = land_part.left.nuc_position; i < proposed; i++)
		left_counts[sampler.indexSeq[i]]++;
	    
	    for (int i = 0; i < right_counts.length; i++)
		right_counts[i] = land_part.counts[i] - left_counts[i];
	}
	else{
	    right_counts = land_part.counts;
	}


	/**** propose new topologies ****/

	int left_tree;
	int right_tree;

	double[] pPartialLogLikelihood = new double[sampler.partition_list.size()];
	double pLogLikelihood = 0, cLogLikelihood = 0, left_likelihood = 0, right_likelihood = 0;
	double proposal_prob;
	boolean keep_on_left;

	Partition curr_part = left_part;


	if (sampler.set.unif.nextDouble() < 0.5){// keep current tree on the left
	    left_tree = land_part.cTree;
	    
	    if ((left_top_change == 0) || (right_top_change == sampler.partition_list.size())){
		right_tree = sampler.ProposeOneTree(left_tree);
		proposal_prob = 1.0/(double)(sampler.TreeSpan - 1);
	    }
	    else{
		right_tree = sampler.ProposeEndTree(right_part.cTree, left_tree);
		proposal_prob = 1.0/(double)(sampler.TreeSpan - 2);
	    }
	    
	    /**** compute current and proposal likelihood  on the right ****/
	    
	    for (int i = land_index + 1; i < right_top_change; i++){
		curr_part = (Partition)sampler.partition_list.elementAt(i);
		pPartialLogLikelihood[i] = 
		    sampler.PostTree[right_tree].LogLikelihood0(curr_part.cMatrix, curr_part.cHyperParameter, 
								     curr_part.counts, curr_part.data,
								     curr_part.cMatrix.pi, sampler.handleGaps, 0, 0);
		
		pLogLikelihood += pPartialLogLikelihood[i];
		cLogLikelihood += curr_part.cPartialLogLikelihood;
		
	    }

	    keep_on_left = true;
	    

	}
	else{// keep current tree on the right
	    right_tree = land_part.cTree;
	    
	    if ((left_top_change == 0) || (right_top_change == sampler.partition_list.size())){
		left_tree = sampler.ProposeOneTree(right_tree);
		proposal_prob = 1.0/(double)(sampler.TreeSpan - 1);
	    }
	    else{
		left_tree = sampler.ProposeEndTree(left_part.cTree, right_tree);
		proposal_prob = 1.0/(double)(sampler.TreeSpan - 2);
	    }

	    /**** compute current and proposal likelihood  on the left ****/

	    for (int i = left_top_change; i < land_index; i++){
		curr_part = (Partition)sampler.partition_list.elementAt(i);
		pPartialLogLikelihood[i] = 
		    sampler.PostTree[left_tree].LogLikelihood0(curr_part.cMatrix, curr_part.cHyperParameter, 
							       curr_part.counts, curr_part.data,
							       curr_part.cMatrix.pi, sampler.handleGaps, 0, 0);

		pLogLikelihood += pPartialLogLikelihood[i];
		cLogLikelihood += curr_part.cPartialLogLikelihood;

	    }
	    
	    keep_on_left = false;

	    
	}
	
	
	/**** compute current and proposal likelihood  on the landing partition ****/

	if (land_part.left.nuc_position != proposed){
	    left_likelihood = 
		sampler.PostTree[left_tree].LogLikelihood0(land_part.cMatrix, land_part.cHyperParameter, 
							   left_counts, land_part.data,
							   land_part.cMatrix.pi, sampler.handleGaps, 0, 0);

	    right_likelihood = 
		sampler.PostTree[right_tree].LogLikelihood0(land_part.cMatrix, land_part.cHyperParameter, 
							    right_counts, land_part.data,
							    land_part.cMatrix.pi, sampler.handleGaps, 0, 0);
	}
	else{
	    left_likelihood = 0;

	    if (keep_on_left){
		right_likelihood = 
		    sampler.PostTree[right_tree].LogLikelihood0(land_part.cMatrix, land_part.cHyperParameter, 
								right_counts, land_part.data,
								land_part.cMatrix.pi, sampler.handleGaps, 0, 0);
	    }
	    else{

		right_likelihood = land_part.cPartialLogLikelihood;
	    }
		
	}


	    
	 	
	pLogLikelihood += left_likelihood + right_likelihood;
	cLogLikelihood += land_part.cPartialLogLikelihood;

	/**** compute acceptance ratio ****/
	
	double logRatio = pLogLikelihood - cLogLikelihood
	    + Math.log(sampler.top_dkp1)
	    - Math.log(sampler.top_one_bk)
	    + sampler.prior_inf.recombination[proposed]
	    + Math.log(sampler.lenSeq - sampler.topology_changes)	
	    - Math.log(sampler.TreeSpan)
	    - Math.log(proposal_prob)
	    - Math.log(sampler.topology_changes);
		

	sampler.tries[4]++;
	
	if( sampler.LogMHAccept(logRatio, sampler.set.unif.nextDouble())){
	    
	    sampler.acceptancerate[4]++;
	    
	    /**** update intermediate likelihoods and topologies ****/

	    if (keep_on_left){
		for (int i = land_index + 1; i < right_top_change; i++){
		    curr_part = (Partition)sampler.partition_list.elementAt(i);
		    
		    curr_part.cPartialLogLikelihood = pPartialLogLikelihood[i];
		    curr_part.cTree = right_tree;
		}

	    }
	    else{
		for (int i = left_top_change; i < land_index; i++){
		    curr_part = (Partition)sampler.partition_list.elementAt(i);
		    
		    curr_part.cPartialLogLikelihood = pPartialLogLikelihood[i];
		    curr_part.cTree = left_tree;
		}

	    }
	    
	    /**** insert new topology change point ****/
	    
	    if (land_part.left.nuc_position == proposed){
		
		land_part.left.top_change = true;

		if (keep_on_left){
		    land_part.cPartialLogLikelihood = right_likelihood;
		    land_part.cTree = right_tree;
		}
	    }
	    else{//land_part.left.nuc_position != proposed
		Partition new_part = new Partition(proposed, land_part.right, true, false);
		
		
		/**** update left side ****/
		land_part.counts = left_counts;
		land_part.right = proposed - 1;
		land_part.cPartialLogLikelihood = left_likelihood;
		land_part.cTree = left_tree;
			
		
		/**** update right side ****/
		new_part.counts = right_counts;
		new_part.cPartialLogLikelihood = right_likelihood;
		new_part.cMatrix = land_part.cMatrix;
		new_part.cHyperParameter = land_part.cHyperParameter;
		new_part.cPartialLogHyperParameterPrior = land_part.cPartialLogHyperParameterPrior;
		new_part.cTree = right_tree;
		new_part.data = land_part.data;
		
		sampler.partition_list.insertElementAt(new_part, land_index + 1);
	    }

	    sampler.topology_changes++;

	}
	
    }



    /**
     * Adds two new topology break-points to the parameter space
     */

    
    protected void AddTwo(){
	
	int proposed = 0, left_proposed = 0, right_proposed = 0; // candidate for new topology change point
	int left_land_index = 0, right_land_index = 0, left_top_change = 0, right_top_change = 0;
	boolean enough_room = false;
	int proposing_second_reciprocal = 0;

	proposed = ProposeTopChangePointToInsert();

	if (sampler.set.unif.nextDouble() < 0.5){// draw the second change point on the left

	    right_land_index = sampler.LandingPart(proposed, 0, sampler.partition_list.size() - 1);
	    left_top_change = FindLeftTopChange(right_land_index);
	    
	    Partition left_part = (Partition)sampler.partition_list.elementAt(left_top_change);
	    
	    if (proposed - left_part.left.nuc_position > 1){
		enough_room = true;

		/*** propose the second topology change point ***/

		left_proposed = left_part.left.nuc_position + 1 +
		    (int)(sampler.set.unif.nextDouble()*(proposed - left_part.left.nuc_position - 1));

		proposing_second_reciprocal = proposed - left_part.left.nuc_position - 1;

		right_proposed = proposed;
		left_land_index = sampler.LandingPart(left_proposed, left_top_change, right_land_index);
		
		
	    }
	    
	}
	else{// draw the second change point on the right

	    left_land_index = sampler.LandingPart(proposed, 0, sampler.partition_list.size() - 1);
	    right_top_change = FindRightTopChange(left_land_index);
	    
	    int end = sampler.lenSeq;
	    
	    if (right_top_change < sampler.partition_list.size()){
		Partition right_part = (Partition)sampler.partition_list.elementAt(right_top_change);
		end = right_part.left.nuc_position;
	    }
	    
	    if (end - proposed > 1){
		enough_room = true;

		/*** proposed the second topology change point ***/
		
		right_proposed = proposed + 1 +
		    (int)(sampler.set.unif.nextDouble()*(end - proposed - 1));

		proposing_second_reciprocal = end - proposed - 1;

		left_proposed = proposed;
		right_land_index = sampler.LandingPart(right_proposed, left_land_index, right_top_change - 1);
		
		
	    }

	}



	if (enough_room){

	    if (left_land_index == right_land_index){// if both topology change points are in the same partition

		Partition land_part = (Partition)sampler.partition_list.elementAt(left_land_index);

		/**** compute counts on three new partitions formed by two new change points ****/
		
		int[] left_counts =  new int[land_part.counts.length];
		int[] middle_counts =  new int[land_part.counts.length];
		int[] right_counts =  new int[land_part.counts.length];
		
		for (int i = land_part.left.nuc_position; i < left_proposed; i++)
		    left_counts[sampler.indexSeq[i]]++;
		
		for (int i = left_proposed; i < right_proposed; i++)
		    middle_counts[sampler.indexSeq[i]]++;

		for (int i = 0; i < right_counts.length; i++)
		    right_counts[i] = land_part.counts[i] - left_counts[i] - middle_counts[i];


		/**** propose new topology in the middle ****/
		
		int middle_tree = sampler.ProposeOneTree(land_part.cTree);

		/**** compute current and proposal likelihods ****/
		
		double pLogLikelihood = 
		    sampler.PostTree[middle_tree].LogLikelihood0(land_part.cMatrix, land_part.cHyperParameter, 
								 middle_counts, land_part.data,
								 land_part.cMatrix.pi, sampler.handleGaps, 0, 0);
		double cLogLikelihood = 
		    sampler.PostTree[land_part.cTree].LogLikelihood0(land_part.cMatrix, land_part.cHyperParameter, 
								     middle_counts, land_part.data,
								     land_part.cMatrix.pi, sampler.handleGaps, 0, 0);
		
		double proposal_prob = 1.0/(double)(sampler.TreeSpan - 1);

		/**** compute acceptance ratio ****/
		
		double logRatio = pLogLikelihood - cLogLikelihood
		    + Math.log(sampler.top_dkp2)
		    - Math.log(sampler.top_two_bk)
		    - Math.log(sampler.topology_changes + 1)
		    + Math.log(sampler.TreeSpan - 1)
		    - Math.log(sampler.TreeSpan)
		    + Math.log(proposing_second_reciprocal)
		    + Math.log(sampler.lenSeq - sampler.topology_changes)
		    + sampler.prior_inf.recombination[left_proposed]
		    + sampler.prior_inf.recombination[right_proposed];


		sampler.tries[5]++;
		
		if( sampler.LogMHAccept(logRatio, sampler.set.unif.nextDouble())){
		    
		    sampler.acceptancerate[5]++;
		    
		    /**** insert two new parameter change points ****/
		    
		    if (land_part.left.nuc_position == left_proposed){
			
			Partition new_part = new Partition(right_proposed, land_part.right, true, false);

			/**** update right side ****/
			new_part.counts = right_counts;
			new_part.cPartialLogLikelihood = land_part.cPartialLogLikelihood - cLogLikelihood;
			new_part.cMatrix = land_part.cMatrix;
			new_part.cHyperParameter = land_part.cHyperParameter;
			new_part.cPartialLogHyperParameterPrior = land_part.cPartialLogHyperParameterPrior;
			new_part.cTree = land_part.cTree;
			new_part.data = land_part.data;
			
			/**** update left side ****/
			    
			land_part.right = right_proposed - 1;
			land_part.left.top_change = true;
			land_part.counts = middle_counts;
			land_part.cPartialLogLikelihood = pLogLikelihood;
			land_part.cTree = middle_tree;

			sampler.partition_list.insertElementAt(new_part, left_land_index + 1);


		    }
		    else{//land_part.left.nuc_position != left_proposed
			
			Partition right_new_part = new Partition(right_proposed, land_part.right, true, false);
			Partition middle_new_part = new Partition(left_proposed, right_proposed - 1, true, false);
			
			/**** update right side ****/
			right_new_part.counts = right_counts;
			right_new_part.cPartialLogLikelihood = 
			    sampler.PostTree[land_part.cTree].LogLikelihood0(land_part.cMatrix, land_part.cHyperParameter,
									 right_counts, land_part.data,
									 land_part.cMatrix.pi, sampler.handleGaps, 0, 0);;
			right_new_part.cMatrix = land_part.cMatrix;
			right_new_part.cHyperParameter = land_part.cHyperParameter;
			right_new_part.cPartialLogHyperParameterPrior = land_part.cPartialLogHyperParameterPrior;
			right_new_part.cTree = land_part.cTree;
			right_new_part.data = land_part.data;

			/**** update middle chunk ****/
			
			middle_new_part.counts = middle_counts;
			middle_new_part.cPartialLogLikelihood = pLogLikelihood;
			middle_new_part.cMatrix = land_part.cMatrix;
			middle_new_part.cHyperParameter = land_part.cHyperParameter;
			middle_new_part.cPartialLogHyperParameterPrior = land_part.cPartialLogHyperParameterPrior;
			middle_new_part.cTree = middle_tree;
			middle_new_part.data = land_part.data;

			
			/**** update left side ****/
			land_part.counts = left_counts;
			land_part.right = left_proposed - 1;
			land_part.cPartialLogLikelihood = land_part.cPartialLogLikelihood - 
			    cLogLikelihood - right_new_part.cPartialLogLikelihood;

			
			sampler.partition_list.insertElementAt(middle_new_part, left_land_index + 1);
			sampler.partition_list.insertElementAt(right_new_part, left_land_index + 2);
		    }
		    
		    sampler.topology_changes +=2;

		}

	    }
	    else{// if two topology change points are in different partitions

		Partition land_part_left = (Partition)sampler.partition_list.elementAt(left_land_index);
		Partition land_part_right = (Partition)sampler.partition_list.elementAt(right_land_index);

		/**** compute counts on the left and right of both new topology change points ****/
		
		int[] left_counts1 =  new int[land_part_left.counts.length];
		int[] right_counts1 =  new int[land_part_left.counts.length];
		int[] left_counts2 =  new int[land_part_left.counts.length];
		int[] right_counts2 =  new int[land_part_left.counts.length];
		
		for (int i = land_part_left.left.nuc_position; i < left_proposed; i++)
		    left_counts1[sampler.indexSeq[i]]++;
		
		for (int i = 0; i < right_counts1.length; i++)
		    right_counts1[i] = land_part_left.counts[i] - left_counts1[i];

		for (int i = land_part_right.left.nuc_position; i < right_proposed; i++)
		    left_counts2[sampler.indexSeq[i]]++;
		
		for (int i = 0; i < right_counts2.length; i++)
		    right_counts2[i] = land_part_right.counts[i] - left_counts2[i];

		/**** propose new topology in the middle ****/
		
		int middle_tree = sampler.ProposeOneTree(land_part_left.cTree);

		/**** compute current and proposal likelihods on the left and right ****/
		
		double pLogLikelihood = 0, cLogLikelihood = 0;

		double curr_right_likelihood1 = 
		    sampler.PostTree[land_part_left.cTree].LogLikelihood0(land_part_left.cMatrix, 
									  land_part_left.cHyperParameter, 
									  right_counts1, land_part_left.data,
									  land_part_left.cMatrix.pi, 
									  sampler.handleGaps, 0, 0);
		
		double prop_right_likelihood1 = 
		    sampler.PostTree[middle_tree].LogLikelihood0(land_part_left.cMatrix, land_part_left.cHyperParameter, 
								 right_counts1, land_part_left.data,
								 land_part_left.cMatrix.pi, sampler.handleGaps, 0, 0);
		cLogLikelihood += curr_right_likelihood1;
		pLogLikelihood += prop_right_likelihood1;


		double curr_left_likelihood2 = 0, prop_left_likelihood2 = 0;
		
		if (land_part_right.left.nuc_position != right_proposed){
		    curr_left_likelihood2 = 
			sampler.PostTree[land_part_left.cTree].LogLikelihood0(land_part_right.cMatrix, 
									      land_part_right.cHyperParameter, 
									      left_counts2, land_part_right.data,
									      land_part_right.cMatrix.pi, 
									      sampler.handleGaps, 0, 0);
		
		    prop_left_likelihood2 = 
			sampler.PostTree[middle_tree].LogLikelihood0(land_part_right.cMatrix, 
								     land_part_right.cHyperParameter,
								     left_counts2, land_part_right.data,
								     land_part_right.cMatrix.pi, 
								     sampler.handleGaps, 0, 0);
		

		    cLogLikelihood += curr_left_likelihood2;
		    pLogLikelihood += prop_left_likelihood2;

		}

		/**** compute current and proposal intermediate likelihoods ****/
		
		double[] pPartialLogLikelihood = new double[sampler.partition_list.size()];
		Partition curr_part = land_part_left;

		for (int i = left_land_index + 1; i < right_land_index; i++){
		    curr_part = (Partition)sampler.partition_list.elementAt(i);
		    pPartialLogLikelihood[i] = 
			sampler.PostTree[middle_tree].LogLikelihood0(curr_part.cMatrix, curr_part.cHyperParameter, 
								     curr_part.counts, curr_part.data,
								     curr_part.cMatrix.pi, sampler.handleGaps, 0, 0);
		    pLogLikelihood += pPartialLogLikelihood[i];
		    cLogLikelihood += curr_part.cPartialLogLikelihood;
		}

		double proposal_prob = 1.0/(double)(sampler.TreeSpan - 1);

		/**** compute acceptance ratio ****/
		
		double logRatio = pLogLikelihood - cLogLikelihood
		    + Math.log(sampler.top_dkp2)
		    - Math.log(sampler.top_two_bk)
		    - Math.log(sampler.topology_changes + 1)
		    + Math.log(sampler.TreeSpan - 1)
		    - Math.log(sampler.TreeSpan)
		    + Math.log(proposing_second_reciprocal)
		    + Math.log(sampler.lenSeq - sampler.topology_changes)
		    + sampler.prior_inf.recombination[left_proposed]
		    + sampler.prior_inf.recombination[right_proposed];

		sampler.tries[5]++;
		
		if( sampler.LogMHAccept(logRatio, sampler.set.unif.nextDouble())){
		    
		    sampler.acceptancerate[5]++;

		    /**** update intermediate likelihoods and topologies ****/

		    for (int i = left_land_index + 1; i < right_land_index; i++){
			curr_part = (Partition)sampler.partition_list.elementAt(i);
			
			curr_part.cPartialLogLikelihood = pPartialLogLikelihood[i];
			curr_part.cTree = middle_tree;

		    }

		    boolean add_first = false;

		    /**** insert the first new topology change point ****/
		    
		    if (land_part_left.left.nuc_position == left_proposed){
			
			land_part_left.left.top_change = true;
			
			land_part_left.cPartialLogLikelihood = prop_right_likelihood1;
			land_part_left.cTree = middle_tree;
		    }
		    else{//land_part.left.nuc_position != proposed
			Partition new_part = new Partition(left_proposed, land_part_left.right, true, false);
			
		
			/**** update left side ****/
			land_part_left.counts = left_counts1;
			land_part_left.right = left_proposed - 1;
			land_part_left.cPartialLogLikelihood -= curr_right_likelihood1;
						
			
			/**** update right side ****/
			new_part.counts = right_counts1;
			new_part.cPartialLogLikelihood = prop_right_likelihood1;
			new_part.cMatrix = land_part_left.cMatrix;
			new_part.cHyperParameter = land_part_left.cHyperParameter;
			new_part.cPartialLogHyperParameterPrior = land_part_left.cPartialLogHyperParameterPrior;
			new_part.cTree = middle_tree;
			new_part.data = land_part_left.data;
			
			sampler.partition_list.insertElementAt(new_part, left_land_index + 1);
			add_first = true;
		    }


		    /**** insert the second new topology change point ****/
		    
		    if (land_part_right.left.nuc_position == right_proposed){
			
			land_part_right.left.top_change = true;
		    }
		    else{//land_part.left.nuc_position != proposed
			Partition new_part = new Partition(right_proposed, land_part_right.right, true, false);

		
			/**** update right side ****/
			new_part.counts = right_counts2;
			new_part.cPartialLogLikelihood = land_part_right.cPartialLogLikelihood - curr_left_likelihood2;
			new_part.cMatrix = land_part_right.cMatrix;
			new_part.cHyperParameter = land_part_right.cHyperParameter;
			new_part.cPartialLogHyperParameterPrior = land_part_right.cPartialLogHyperParameterPrior;
			new_part.cTree = land_part_right.cTree;
			new_part.data = land_part_right.data;
			
		
			/**** update left side ****/
			land_part_right.counts = left_counts2;
			land_part_right.right = right_proposed - 1;
			land_part_right.cPartialLogLikelihood = prop_left_likelihood2;
			land_part_right.cTree = middle_tree;

			if (add_first){
			    sampler.partition_list.insertElementAt(new_part, right_land_index + 2);
			}
			else{
			    sampler.partition_list.insertElementAt(new_part, right_land_index + 1);
			}

		    }
		    
		    
		    sampler.topology_changes += 2;
		    
		}
		
	    }
	    
		
	}
	
    }


    /**
     * Removes a topology break-point (or two break-points) from the parameter space
     */


    protected void Delete(){
	
	int propose_to_delete;
	propose_to_delete = ProposeTopChangePointToDelete();


	/**** find left and right parameter change point ****/
	
	int left_top_change = FindLeftTopChange(propose_to_delete - 1);
	int right_top_change = FindRightTopChange(propose_to_delete);

	Partition proposed_part = (Partition)sampler.partition_list.elementAt(propose_to_delete);
	Partition left_part = (Partition)sampler.partition_list.elementAt(propose_to_delete - 1);
	Partition left_change_part = null, right_change_part = null;

	/**** decide which tree to keep after deleting the chosen point ****/

	int after_collapse_tree = 0;
	boolean inconsistency = false, keep_on_left = false;
	int left_to_delete = 0, right_to_delete = 0;
	int proposing_second_reciprocal = 0;
	
	if (sampler.set.unif.nextDouble() < 0.5){
	    after_collapse_tree = left_part.cTree;
	    keep_on_left = true;

	    /**** check if such a choise produced inconsistency ****/
	    
	    if (right_top_change != sampler.partition_list.size()){
		right_change_part = (Partition)sampler.partition_list.elementAt(right_top_change);

		if (right_change_part.cTree == after_collapse_tree){
		    inconsistency = true;
		    left_to_delete = propose_to_delete;
		    right_to_delete = right_top_change;
		    right_top_change = FindRightTopChange(right_to_delete);

		    if (right_top_change == sampler.partition_list.size()){
			proposing_second_reciprocal = sampler.lenSeq - left_to_delete - 1;
		    }
		    else{
			Partition very_right = (Partition)sampler.partition_list.elementAt(right_top_change);
			proposing_second_reciprocal = very_right.left.nuc_position - left_to_delete - 1;
		    }
		}
	    }

	}
	else{
	    after_collapse_tree = proposed_part.cTree;
	    keep_on_left = false;

	    /**** check if such a choise produced inconsistency ****/
	    
	    if (left_top_change != 0){
		left_change_part = (Partition)sampler.partition_list.elementAt(left_top_change - 1);

		if (left_change_part.cTree == after_collapse_tree){
		    inconsistency = true;
		    left_to_delete = left_top_change;
		    right_to_delete = propose_to_delete;
		    left_top_change = FindLeftTopChange(left_to_delete - 1);

		    Partition very_left = (Partition)sampler.partition_list.elementAt(left_top_change);
		    proposing_second_reciprocal = right_to_delete - very_left.left.nuc_position;
		}
	    }
	    
	}

	if (!inconsistency){// delete one topology change point
	    

	    double[] pPartialLogLikelihood = new double[sampler.partition_list.size()];
	    double pLogLikelihood = 0, cLogLikelihood = 0;

	    Partition curr_part = left_part;
	    
	    if (keep_on_left){

		/**** compute current and proposal likelihood on the left and on the right ****/
	
		for (int i = propose_to_delete; i < right_top_change; i++){
		    curr_part = (Partition)sampler.partition_list.elementAt(i);
		    pPartialLogLikelihood[i] = 
			sampler.PostTree[after_collapse_tree].LogLikelihood0(curr_part.cMatrix, curr_part.cHyperParameter,
									     curr_part.counts, curr_part.data,
									     curr_part.cMatrix.pi, 
									     sampler.handleGaps, 0, 0);
		    
		    pLogLikelihood += pPartialLogLikelihood[i];
		    cLogLikelihood += curr_part.cPartialLogLikelihood;
		    
		}
	    }
	    else{

		/**** compute current and proposal likelihood on the left and on the left ****/

		for (int i = left_top_change; i < propose_to_delete; i++){
		    curr_part = (Partition)sampler.partition_list.elementAt(i);
		    pPartialLogLikelihood[i] = 
			sampler.PostTree[after_collapse_tree].LogLikelihood0(curr_part.cMatrix, curr_part.cHyperParameter,
									     curr_part.counts, curr_part.data,
									     curr_part.cMatrix.pi, 
									     sampler.handleGaps, 0, 0);
		    
		    pLogLikelihood += pPartialLogLikelihood[i];
		    cLogLikelihood += curr_part.cPartialLogLikelihood;
		    
		}
		
	    }

	    /**** compute acceptance ratio ****/
	    
	    double logRatio = pLogLikelihood - cLogLikelihood
		- Math.log(sampler.top_dk)
		+ Math.log(sampler.top_one_bkm1)
		+ Math.log(sampler.topology_changes)
		+ Math.log(sampler.TreeSpan)
		- Math.log(sampler.TreeSpan - 2)
		- Math.log(sampler.lenSeq - sampler.topology_changes +1)
		- sampler.prior_inf.recombination[proposed_part.left.nuc_position];


	    sampler.tries[7]++;


	    if( sampler.LogMHAccept(logRatio, sampler.set.unif.nextDouble())){
		
		sampler.acceptancerate[7]++;
		
		
		if (keep_on_left){

		    /**** update intermediate likelihoods and topology on the right****/

		    for (int i = propose_to_delete + 1; i < right_top_change; i++){
			curr_part = (Partition)sampler.partition_list.elementAt(i);
			
			curr_part.cPartialLogLikelihood = pPartialLogLikelihood[i];
			curr_part.cTree = after_collapse_tree;
		    }

		    /**** delete topology change point ****/
		    
		    if (proposed_part.IsParameterChange()){
			
			proposed_part.left.top_change = false;
			proposed_part.cPartialLogLikelihood = pPartialLogLikelihood[propose_to_delete];
			proposed_part.cTree = after_collapse_tree;
		    }
		    else{//land_part.left.nuc_position != proposed

			/**** compute counts for new partition without a change point ****/
			
			int[] new_counts = new int[proposed_part.counts.length];
			for (int i = 0; i < new_counts.length; i++)
			    new_counts[i] = left_part.counts[i] + proposed_part.counts[i];
			
			/**** extend left partition ****/
			left_part.counts = new_counts;
			left_part.cPartialLogLikelihood += pPartialLogLikelihood[propose_to_delete];
			left_part.cTree = after_collapse_tree;
			left_part.right = proposed_part.right;
			
			
			sampler.partition_list.removeElementAt(propose_to_delete);
		    }

		   
		}
		else{
		    
		    /**** update intermediate likelihoods and topology on the left****/

		    for (int i = left_top_change; i < propose_to_delete; i++){
			curr_part = (Partition)sampler.partition_list.elementAt(i);
			
			curr_part.cPartialLogLikelihood = pPartialLogLikelihood[i];
			curr_part.cTree = after_collapse_tree;
		    }

		    
		    /**** delete topology change point ****/
		    
		    if (proposed_part.IsParameterChange()){
			
			proposed_part.left.top_change = false;
		    }
		    else{//land_part.left.nuc_position != proposed

			/**** compute counts for new partition without a change point ****/
		
			int[] new_counts = new int[proposed_part.counts.length];
			for (int i = 0; i < new_counts.length; i++)
			    new_counts[i] = left_part.counts[i] + proposed_part.counts[i];
			
			/**** extend left partition ****/
			left_part.counts = new_counts;
			left_part.cPartialLogLikelihood += proposed_part.cPartialLogLikelihood;
			left_part.right = proposed_part.right;
			
			
			sampler.partition_list.removeElementAt(propose_to_delete);
		    }

		}

		sampler.topology_changes--;
		
	    }
	    
	    

	}
	else{// delete two topology change points

	    double[] pPartialLogLikelihood = new double[sampler.partition_list.size()];
	    double pLogLikelihood = 0, cLogLikelihood = 0;

	    Partition curr_part = left_part;
	    
	   
	    /**** compute current and proposal likelihood  ****/
	    
	    for (int i = left_to_delete; i < right_to_delete; i++){
		curr_part = (Partition)sampler.partition_list.elementAt(i);
		pPartialLogLikelihood[i] = 
		    sampler.PostTree[after_collapse_tree].LogLikelihood0(curr_part.cMatrix, curr_part.cHyperParameter,
									 curr_part.counts, curr_part.data,
									 curr_part.cMatrix.pi, 
									 sampler.handleGaps, 0, 0);
		
		pLogLikelihood += pPartialLogLikelihood[i];
		cLogLikelihood += curr_part.cPartialLogLikelihood;
		
	    }

	    Partition left = (Partition)sampler.partition_list.elementAt(left_to_delete);
	    Partition right = (Partition)sampler.partition_list.elementAt(right_to_delete);

	  
	    /**** compute acceptance ratio ****/
	    
	    double logRatio = pLogLikelihood - cLogLikelihood
		+ Math.log(sampler.top_two_bkm2)
		- Math.log(sampler.top_dk)
		+ Math.log(sampler.topology_changes - 1)
		+ Math.log(sampler.TreeSpan)
		- Math.log(sampler.TreeSpan - 1)
		- Math.log(proposing_second_reciprocal)
		- Math.log(sampler.lenSeq - sampler.topology_changes + 1)
		- sampler.prior_inf.recombination[left.left.nuc_position]
		- sampler.prior_inf.recombination[right.left.nuc_position];

	    sampler.tries[8]++;


	    if( sampler.LogMHAccept(logRatio, sampler.set.unif.nextDouble())){
		
		sampler.acceptancerate[8]++;

		Partition first_proposed_part = (Partition)sampler.partition_list.elementAt(left_to_delete);
		Partition first_left_part = (Partition)sampler.partition_list.elementAt(left_to_delete - 1);
		Partition second_proposed_part = (Partition)sampler.partition_list.elementAt(right_to_delete);
		Partition second_left_part = (Partition)sampler.partition_list.elementAt(right_to_delete - 1);
		

		/**** update intermediate likelihoods and topology****/
		
		for (int i = left_to_delete; i < right_to_delete; i++){
		    curr_part = (Partition)sampler.partition_list.elementAt(i);
		    
		    curr_part.cPartialLogLikelihood = pPartialLogLikelihood[i];
		    curr_part.cTree = after_collapse_tree;
		}

		if (left_to_delete == right_to_delete - 1){
		
		    if (second_proposed_part.IsParameterChange()){
			second_proposed_part.left.top_change = false;

			if (first_proposed_part.IsParameterChange()){
			    first_proposed_part.left.top_change = false;
			    first_proposed_part.cTree = after_collapse_tree;
			}
			else{
		     
			    int[] new_counts = new int[first_proposed_part.counts.length];
			    for (int i = 0; i < new_counts.length; i++)
				new_counts[i] = first_left_part.counts[i] + first_proposed_part.counts[i];
			    
			    first_left_part.counts = new_counts;
			    first_left_part.cPartialLogLikelihood += first_proposed_part.cPartialLogLikelihood;
			    first_left_part.right = first_proposed_part.right;
			    
			    sampler.partition_list.removeElementAt(left_to_delete);
			}

		    }
		    else {

			if (first_proposed_part.IsParameterChange()){

			    first_proposed_part.left.top_change = false;

			    int[] new_counts = new int[first_proposed_part.counts.length];
			    for (int i = 0; i < new_counts.length; i++)
				new_counts[i] = first_proposed_part.counts[i] + second_proposed_part.counts[i];
			    
			    first_proposed_part.counts = new_counts;
			    first_proposed_part.cPartialLogLikelihood += second_proposed_part.cPartialLogLikelihood;
			    first_proposed_part.right = second_proposed_part.right;
			    first_proposed_part.cTree = after_collapse_tree;
			    
			    sampler.partition_list.removeElementAt(right_to_delete);

			}
			else{

			    int[] new_counts = new int[first_proposed_part.counts.length];
			    for (int i = 0; i < new_counts.length; i++)
				new_counts[i] = first_left_part.counts[i] + first_proposed_part.counts[i]
				    + second_proposed_part.counts[i];
			    
			    first_left_part.counts = new_counts;
			    first_left_part.cPartialLogLikelihood += first_proposed_part.cPartialLogLikelihood 
				+ second_proposed_part.cPartialLogLikelihood;
			    first_left_part.right = second_proposed_part.right;
			    
			    sampler.partition_list.removeElementAt(right_to_delete);
			    sampler.partition_list.removeElementAt(left_to_delete);
			}

		    }

		}
		else{

		    /**** delete RIGHT topology change point ****/
		    
		    if (second_proposed_part.IsParameterChange()){
			
			second_proposed_part.left.top_change = false;
			
		    }
		    else{//land_part.left.nuc_position != proposed
			
			/**** compute counts for new partitions without a change points ****/
			
			int[] new_counts = new int[second_proposed_part.counts.length];
			for (int i = 0; i < new_counts.length; i++)
			    new_counts[i] = second_left_part.counts[i] + second_proposed_part.counts[i];
			
			/**** extend left partition ****/
			second_left_part.counts = new_counts;
			second_left_part.cPartialLogLikelihood += second_proposed_part.cPartialLogLikelihood;
			second_left_part.right = second_proposed_part.right;
			
			
			sampler.partition_list.removeElementAt(right_to_delete);
		    }
		    

		   		    
		    /**** delete LEFT topology change point ****/
		    
		    if (first_proposed_part.IsParameterChange()){

			
			
			first_proposed_part.left.top_change = false;
			first_proposed_part.cPartialLogLikelihood = pPartialLogLikelihood[left_to_delete];
			first_proposed_part.cTree = after_collapse_tree;
		    }
		    else{//land_part.left.nuc_position != proposed
			
			/**** compute counts for new partitions without a change points ****/
			
			int[] new_counts = new int[first_proposed_part.counts.length];
			
			for (int i = 0; i < new_counts.length; i++)
			    new_counts[i] = first_left_part.counts[i] + first_proposed_part.counts[i];
			
			/**** extend left partition ****/
			first_left_part.counts = new_counts;
			
			if (right_to_delete - left_to_delete == 1){
			first_left_part.cPartialLogLikelihood += first_proposed_part.cPartialLogLikelihood;
			}
			else{
			    first_left_part.cPartialLogLikelihood += pPartialLogLikelihood[left_to_delete];
			}
			first_left_part.cTree = after_collapse_tree;
			first_left_part.right = first_proposed_part.right;
			
			
			sampler.partition_list.removeElementAt(left_to_delete);
		    }
		    
		}    
		
		sampler.topology_changes -= 2;
		
		
	    }   
	    

	}
	    

	
    }
    

}
