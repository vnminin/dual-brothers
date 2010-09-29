import java.util.*;
import java.io.*;
import corejava.*;


/**
* Adding and deleting new evolutionary change-points kernel
*/

public class ParDimChange {

    DCPSampler sampler;
    
    public ParDimChange (DCPSampler s){
	
	sampler = s;
    }
    
    
    /**
     * Adds a new evolutionary change-point to the parameter space
     */
    
    protected void Add(){
	
	int proposed; // candidate new parameter change point
	int land_index, left_par_change, right_par_change;
	

	proposed = ProposeParChangePointToInsert();


	/**** find partition, where new point landed ****/

	land_index = sampler.LandingPart(proposed, 0, sampler.partition_list.size() - 1);

	/**** find left and right parameter change point ****/

	left_par_change = FindLeftParChange(land_index);
	right_par_change = FindRightParChange(land_index);


	/**** compute start and end of the constant parameter interval ****/

	Partition left_change_part = (Partition)sampler.partition_list.elementAt(left_par_change);


	int start = left_change_part.left.nuc_position;
	int end = sampler.lenSeq;


	if (right_par_change != sampler.partition_list.size()){
	    Partition right_change_part = (Partition)sampler.partition_list.elementAt(right_par_change);
	    end = right_change_part.left.nuc_position;
	}


	/**** compute counts on both sides of proposed new change point ****/

	Partition land_part = (Partition)sampler.partition_list.elementAt(land_index);
	
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

	/**** compute relative lengths of the new segments ****/

	double length1 = proposed - start;
	double length2 = end - proposed;
	double weight1 = length1/(length1 + length2);
	double weight2 = length2/(length1 + length2);

	/**** propose new evolutionary parameters ****/

	double zAlpha = sampler.set.norm01.nextDouble(); 
	double zMu    = sampler.set.norm01.nextDouble(); 
	
	double pHyperParameter1 = sampler.DrawnMu1(land_part.cHyperParameter, weight2, sampler.sigmaM, zMu);
	double pHyperParameter2 = sampler.DrawnMu2(land_part.cHyperParameter, weight1, sampler.sigmaM, zMu);

	double cEP  = land_part.cMatrix.v[0];
	double pEP1 = sampler.DrawnAlpha1(cEP, weight2, sampler.sigmaA, zAlpha);
	double pEP2 = sampler.DrawnAlpha2(cEP, weight1, sampler.sigmaA, zAlpha);
	
	double[] ppEP1 = new double[2];
	double[] ppEP2 = new double[2];

	ppEP1[0] = pEP1;
	ppEP1[1] = pEP1;
	ppEP2[0] = pEP2;
	ppEP2[1] = pEP2;
	
	QMatrix pMatrix1 = new iHKYNoBoundFixPiMatrix(ppEP1, land_part.cMatrix.pi);
	QMatrix pMatrix2 = new iHKYNoBoundFixPiMatrix(ppEP2, land_part.cMatrix.pi);
		

	/**** compute current and proposal likelihood on the left and on the right ****/

	double[] pPartialLogLikelihood = new double[sampler.partition_list.size()];
	double pLogLikelihood = 0, cLogLikelihood = 0, left_likelihood = 0, right_likelihood = 0;

	Partition curr_part = left_change_part;

	for (int i = left_par_change; i < land_index; i++){
	    curr_part = (Partition)sampler.partition_list.elementAt(i);
	    pPartialLogLikelihood[i] = 
		sampler.PostTree[curr_part.cTree].LogLikelihood0(pMatrix1, pHyperParameter1, 
							 curr_part.counts, curr_part.data,
							 pMatrix1.pi, sampler.handleGaps, 0, 0);

	    pLogLikelihood += pPartialLogLikelihood[i];
	    cLogLikelihood += curr_part.cPartialLogLikelihood;

	}

	for (int i = land_index + 1; i < right_par_change; i++){
	    curr_part = (Partition)sampler.partition_list.elementAt(i);
	    pPartialLogLikelihood[i] = 
		sampler.PostTree[curr_part.cTree].LogLikelihood0(pMatrix2, pHyperParameter2, 
							 curr_part.counts, curr_part.data,
							 pMatrix2.pi, sampler.handleGaps, 0, 0);

	    pLogLikelihood += pPartialLogLikelihood[i];
	    cLogLikelihood += curr_part.cPartialLogLikelihood;

	}

	if (land_part.left.nuc_position != proposed){
	    left_likelihood = 
		sampler.PostTree[land_part.cTree].LogLikelihood0(pMatrix1, pHyperParameter1, 
							 left_counts, land_part.data,
							 pMatrix1.pi, sampler.handleGaps, 0, 0);
	}
	else{
	    left_likelihood = 0;
	}

	right_likelihood = 
	    sampler.PostTree[land_part.cTree].LogLikelihood0(pMatrix2, pHyperParameter2, 
							 right_counts, curr_part.data,
							 pMatrix2.pi, sampler.handleGaps, 0, 0);
	    
	 	
	pLogLikelihood += left_likelihood + right_likelihood;
	cLogLikelihood += land_part.cPartialLogLikelihood;


	/**** compute priors ****/
	double pPartialLogHyperParameterPrior1 = sampler.prior_inf.Prior_on_mu(pHyperParameter1);
	double pPartialLogHyperParameterPrior2 = sampler.prior_inf.Prior_on_mu(pHyperParameter2);
	double pLog_alpha_prior1 = sampler.prior_inf.Prior_on_alpha(pEP1);
	double pLog_alpha_prior2 = sampler.prior_inf.Prior_on_alpha(pEP2);
	

	/**** compute acceptance ratio ****/
	
	double logRatio = pLogLikelihood - cLogLikelihood
	    + pPartialLogHyperParameterPrior1 + pPartialLogHyperParameterPrior2
	    + pLog_alpha_prior1
	    + pLog_alpha_prior2
	    - land_part.log_alpha_prior
	    - land_part.cPartialLogHyperParameterPrior
	    //- Math.log(sampler.parameter_changes)
	    - sampler.logStandardNormalDensity(zAlpha)
	    - sampler.logStandardNormalDensity(zMu)
	    + sampler.logJacobian(cEP, pEP1, pEP2, sampler.sigmaA,
			  land_part.cHyperParameter, pHyperParameter1, pHyperParameter2, sampler.sigmaM);
	
	sampler.tries[6]++;


	if( sampler.LogMHAccept(logRatio, sampler.set.unif.nextDouble())){

	    sampler.acceptancerate[6]++;
	    
	    /**** update intermediate likelihoods and evolutionary parameters ****/
	    
	    for (int i = left_par_change; i < land_index; i++){
		curr_part = (Partition)sampler.partition_list.elementAt(i);
		
		curr_part.cPartialLogLikelihood = pPartialLogLikelihood[i];
		curr_part.cMatrix = pMatrix1;
		curr_part.cHyperParameter = pHyperParameter1;
		curr_part.cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior1;
		curr_part.log_alpha_prior = pLog_alpha_prior1;
	    }
	    
	    for (int i = land_index + 1; i < right_par_change; i++){
		curr_part = (Partition)sampler.partition_list.elementAt(i);
		
		curr_part.cPartialLogLikelihood = pPartialLogLikelihood[i];
		curr_part.cMatrix = pMatrix2;
		curr_part.cHyperParameter = pHyperParameter2;
		curr_part.cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior2;
		curr_part.log_alpha_prior = pLog_alpha_prior2;
	    }
	
	    
	    /**** insert new parameter change point ****/
	    
	    if (land_part.left.nuc_position == proposed){
		
		land_part.left.par_change = true;
		land_part.cPartialLogLikelihood = right_likelihood;
		land_part.cMatrix = pMatrix2;
		land_part.cHyperParameter = pHyperParameter2;
		land_part.cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior2;
		land_part.log_alpha_prior = pLog_alpha_prior2;

	    }
	    else{//land_part.left.nuc_position != proposed
		Partition new_part = new Partition(proposed, land_part.right, false, true);
		
		/**** update left side ****/
		land_part.counts = left_counts;
		land_part.cPartialLogLikelihood = left_likelihood;
		land_part.cMatrix = pMatrix1;
		land_part.cHyperParameter = pHyperParameter1;
		land_part.cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior1;
		land_part.log_alpha_prior = pLog_alpha_prior1;
		land_part.right = proposed - 1;


		/**** update right side ****/
		new_part.counts = right_counts;
		new_part.cPartialLogLikelihood = right_likelihood;
		new_part.cMatrix = pMatrix2;
		new_part.cHyperParameter = pHyperParameter2;
		new_part.cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior2;
		new_part.log_alpha_prior = pLog_alpha_prior2;
		new_part.cTree = land_part.cTree;
		new_part.data = land_part.data;
		
		sampler.partition_list.insertElementAt(new_part, land_index + 1);
	    }

	    sampler.parameter_changes++;
	    
	}   


    }


 
    /**
     * Removes an evolutionary change-point from the parameter space
     */


    protected void Delete(){
    
	int propose_to_delete;

	propose_to_delete = ProposeParChangePointToDelete();


	/**** find left and right parameter change point ****/

	int left_par_change = FindLeftParChange(propose_to_delete - 1);
	int right_par_change = FindRightParChange(propose_to_delete);


	/**** compute start and end of the constant parameter interval ****/

	Partition left_change_part = (Partition)sampler.partition_list.elementAt(left_par_change);


	int start = left_change_part.left.nuc_position;
	int end = sampler.lenSeq;

	if (right_par_change != sampler.partition_list.size()){
	    Partition right_change_part = (Partition)sampler.partition_list.elementAt(right_par_change);
	    end = right_change_part.left.nuc_position;
	}


	

	
	/**** compute counts for new partition without a change point ****/
	
	Partition proposed_part = (Partition)sampler.partition_list.elementAt(propose_to_delete);
	Partition left_part = (Partition)sampler.partition_list.elementAt(propose_to_delete - 1);

	int[] new_counts = new int[proposed_part.counts.length];

	for (int i = 0; i < new_counts.length; i++)
	    new_counts[i] = left_part.counts[i] + proposed_part.counts[i];

	/**** compute relative lengths of the new segments ****/

	double length1 = proposed_part.left.nuc_position - start;
	double length2 = end - proposed_part.left.nuc_position;
	double weight1 = length1/(length1 + length2);
	double weight2 = length2/(length1 + length2);

	/**** propose new evolutionary parameters ****/
	
	double cEP1;
	double cEP2;
	double pEP;
	double cHyperParameter1;
	double cHyperParameter2;
	double ppHyperParameter;
	
	double cPartialLogHyperParameterPrior1;
	double cPartialLogHyperParameterPrior2;
	double pPartialLogHyperParameterPrior;
	double cLog_alpha_prior1;
	double cLog_alpha_prior2;
	double pLog_alpha_prior;


	cHyperParameter1 = left_part.cHyperParameter;
	cHyperParameter2 = proposed_part.cHyperParameter;
	cEP1 = left_part.cMatrix.v[0];
	cEP2 = proposed_part.cMatrix.v[0];
	cPartialLogHyperParameterPrior1 = left_part.cPartialLogHyperParameterPrior;
	cPartialLogHyperParameterPrior2 = proposed_part.cPartialLogHyperParameterPrior;
	cLog_alpha_prior1 = left_part.log_alpha_prior;
	cLog_alpha_prior2 = proposed_part.log_alpha_prior;

	
	pEP = sampler.CondensedAlpha( cEP1, weight1, cEP2, weight2 );
	ppHyperParameter = sampler.CondensedMu( cHyperParameter1, weight1, cHyperParameter2, weight2 );
	
	double zAlpha = sampler.InverseZAlpha( cEP1, pEP, weight2, sampler.sigmaA );
	double zMu = sampler.InverseZMu( cHyperParameter1, ppHyperParameter, weight2,  sampler.sigmaM );
	
	double[] ppEP = new double[2];
	ppEP[0] = pEP;
	ppEP[1] = pEP;
	
	QMatrix pMatrix = new iHKYNoBoundFixPiMatrix(ppEP, proposed_part.cMatrix.pi);
	

	/**** compute current and proposal likelihood on the left and on the right ****/
	
	double[] pPartialLogLikelihood = new double[sampler.partition_list.size()];
	double pLogLikelihood = 0, cLogLikelihood = 0, new_likelihood = 0;
	
	Partition curr_part = left_part;

	for (int i = left_par_change; i < propose_to_delete; i++){
	    curr_part = (Partition)sampler.partition_list.elementAt(i);
	    pPartialLogLikelihood[i] = 
		sampler.PostTree[curr_part.cTree].LogLikelihood0(pMatrix, ppHyperParameter, 
							 curr_part.counts, curr_part.data,
							 pMatrix.pi, sampler.handleGaps, 0, 0);
	    
	    pLogLikelihood += pPartialLogLikelihood[i];
	    cLogLikelihood += curr_part.cPartialLogLikelihood;
	    
	}
	
	for (int i = propose_to_delete + 1; i < right_par_change; i++){
	    curr_part = (Partition)sampler.partition_list.elementAt(i);
	    pPartialLogLikelihood[i] = 
		sampler.PostTree[curr_part.cTree].LogLikelihood0(pMatrix, ppHyperParameter, 
							 curr_part.counts, curr_part.data,
							 pMatrix.pi, sampler.handleGaps, 0, 0);
	    
	    pLogLikelihood += pPartialLogLikelihood[i];
	    cLogLikelihood += curr_part.cPartialLogLikelihood;
	    
	}
	
	
	new_likelihood = 
	    sampler.PostTree[proposed_part.cTree].LogLikelihood0(pMatrix, ppHyperParameter, 
						     proposed_part.counts, proposed_part.data,
						     pMatrix.pi, sampler.handleGaps, 0, 0);
	
	 	
	pLogLikelihood += new_likelihood;
	cLogLikelihood += proposed_part.cPartialLogLikelihood;


	/**** compute priors ****/

	pPartialLogHyperParameterPrior = sampler.prior_inf.Prior_on_mu(ppHyperParameter);
	pLog_alpha_prior = sampler.prior_inf.Prior_on_alpha(pEP);

	/**** compute acceptance ratio ****/

	double logRatio = pLogLikelihood - cLogLikelihood
	    + pPartialLogHyperParameterPrior
	    + pLog_alpha_prior
	    - cLog_alpha_prior1
	    - cLog_alpha_prior2
	    - cPartialLogHyperParameterPrior1 - cPartialLogHyperParameterPrior2
	    //+ Math.log(sampler.parameter_changes - 1)
	    + sampler.logStandardNormalDensity(zAlpha)
	    + sampler.logStandardNormalDensity(zMu)
	    - sampler.logJacobian(pEP, cEP1, cEP2, sampler.sigmaA,
				  ppHyperParameter, cHyperParameter1, cHyperParameter2, sampler.sigmaM);
	
	sampler.tries[9]++;


	if( sampler.LogMHAccept(logRatio, sampler.set.unif.nextDouble())){
	    
	    sampler.acceptancerate[9]++;

	    /**** update intermediate likelihoods and evolutionary parameters ****/
	    
	    for (int i = left_par_change; i < propose_to_delete; i++){
		curr_part = (Partition)sampler.partition_list.elementAt(i);
		
		curr_part.cPartialLogLikelihood = pPartialLogLikelihood[i];
		curr_part.cMatrix = pMatrix;
		curr_part.cHyperParameter = ppHyperParameter;
		curr_part.cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior;
		curr_part.log_alpha_prior = pLog_alpha_prior;
	    }
	    
	    for (int i = propose_to_delete + 1; i < right_par_change; i++){
		curr_part = (Partition)sampler.partition_list.elementAt(i);
		
		curr_part.cPartialLogLikelihood = pPartialLogLikelihood[i];
		curr_part.cMatrix = pMatrix;
		curr_part.cHyperParameter = ppHyperParameter;
		curr_part.cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior;
		curr_part.log_alpha_prior = pLog_alpha_prior;
	    }

	    /**** delete parameter change point ****/
	    
	    if (proposed_part.IsTopologyChange()){
		
		proposed_part.left.par_change = false;
		proposed_part.cPartialLogLikelihood = new_likelihood;
		proposed_part.cMatrix = pMatrix;
		proposed_part.cHyperParameter = ppHyperParameter;
		proposed_part.cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior;
		proposed_part.log_alpha_prior = pLog_alpha_prior;

	    }
	    else{//land_part.left.nuc_position != proposed
			
		/**** extend left partition ****/
		left_part.counts = new_counts;
		left_part.cPartialLogLikelihood += new_likelihood;
		left_part.cMatrix = pMatrix;
		left_part.cHyperParameter = ppHyperParameter;
		left_part.cPartialLogHyperParameterPrior = pPartialLogHyperParameterPrior;
		left_part.log_alpha_prior = pLog_alpha_prior;
		left_part.right = proposed_part.right;
		
		
		sampler.partition_list.removeElementAt(propose_to_delete);
	    }

	    sampler.parameter_changes--;
	    
	} 
	

	
    }



    /** 
     * Returns a randomly chosen evolutionary change-point to remove
     */
    

    private int ProposeParChangePointToDelete(){
	
	int[] par_ch_index = new int[sampler.parameter_changes - 1];
	
	int j = 0;
	Partition curr_part;
	for (int i = 1; i < sampler.partition_list.size(); i++){
	    curr_part = (Partition)sampler.partition_list.elementAt(i);

	    if (curr_part.IsParameterChange()){
		par_ch_index[j] = i;
		j++;
	    }
	}
	
	int rnd = (int)(sampler.set.unif.nextDouble()*(sampler.parameter_changes - 1));
	

	return par_ch_index[rnd];
	
    }


    /** 
     * Returns a randomly chosen location of a new evolutionary change-point
     */

    private int ProposeParChangePointToInsert(){
	boolean good_point = false;
	Partition curr_part = null;
	int new_point = 0;
	
	while (!good_point){
	    new_point = (int)(sampler.set.unif.nextDouble()*sampler.lenSeq);
	    
	    good_point = true;
	    for (int i = 0; i < sampler.partition_list.size(); i++){
		curr_part = (Partition)sampler.partition_list.elementAt(i);
		if (curr_part.IsParameterChange() && new_point == curr_part.left.nuc_position)
		    good_point = false;
	    }
	    
	    //if (new_point == curr_part.right)
	    //good_point = false;
	    
	}
	
	return new_point;
    }


    
    /**
     * Returns the closest evolutionary change-point to the left
     * of the reference change-point
     * @param landing reference change-point
     */

    private int FindLeftParChange(int landing){
	Partition land_part = (Partition)sampler.partition_list.elementAt(landing);
	boolean par_change = false;
	
	if (land_part.IsParameterChange())
	    par_change = true;
    
	int i = landing - 1;
	Partition curr_part;
	while((!par_change) && (i >= 0)){
	    curr_part = (Partition)sampler.partition_list.elementAt(i);
	    if (curr_part.IsParameterChange())
		par_change = true;

	    i--;
	}
	
	return i+1;
	
    }


    /**
     * Returns the closest evolutionary change-point to the right
     * of the reference change-point
     * @param landing reference change-point
     */

    private int FindRightParChange(int landing){
	Partition land_part = (Partition)sampler.partition_list.elementAt(landing);
	boolean par_change = false;
	int right_change;

	if (landing == sampler.partition_list.size() - 1){
	    par_change = true;
	    right_change = sampler.partition_list.size();
	}
	else {
	    int i = landing + 1;
	    Partition curr_part;
	    while((!par_change) && (i < sampler.partition_list.size())){
		curr_part = (Partition)sampler.partition_list.elementAt(i);
		if (curr_part.IsParameterChange())
		    par_change = true;
		
		i++;
	    }
	    
	    if (!par_change)
		right_change = sampler.partition_list.size();
	    else
		right_change = i - 1;
	}

	return right_change;
	    
    }

}
