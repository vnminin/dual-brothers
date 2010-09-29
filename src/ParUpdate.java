import java.util.*;
import java.io.*;
import corejava.*;


/**
* Parameter update kernel
*/

public class ParUpdate {

    DCPSampler sampler;
    
    public ParUpdate (DCPSampler s){

	sampler = s;
    }

 /**
  * Sequentially updates evolutionary parameters alpha and mu at all partitions of the alignment
  * together with hyperparameters for alpha
  */

    
    protected void Update() {
	double pLogHyperParameterPrior = 0;
	double pLog_alpha_prior = 0;
	double pHyperParameter = 0;
	double pLogLikelihood = 0, cLogLikelihood = 0, logRatio = 0;
      	double[] pPartialLogLikelihood = new double[sampler.partition_list.size()];
	Partition temp_part;

	int parameter_start = 0;



	
	/***************** start with the first partition ******************/
	
	temp_part = (Partition)sampler.partition_list.elementAt(0);
	
	QMatrix pMatrix = temp_part.cMatrix.Proposal(sampler.set);
	if (pMatrix != null){
	    pHyperParameter = ProposeHyperParameter(temp_part.cHyperParameter,1);
	    pLogHyperParameterPrior = sampler.prior_inf.Prior_on_mu(pHyperParameter);
	    pLog_alpha_prior = sampler.prior_inf.Prior_on_alpha(pMatrix.v[0]);
	    pPartialLogLikelihood[0] = sampler.PostTree[temp_part.cTree].LogLikelihood0(pMatrix, pHyperParameter, temp_part.counts, temp_part.data, pMatrix.pi, sampler.handleGaps, 0, 0);
	    pLogLikelihood = pPartialLogLikelihood[0];
	    cLogLikelihood = temp_part.cPartialLogLikelihood;
	}
	
	
	for(int p = 1; p < sampler.partition_list.size(); p++) {

	    temp_part = (Partition)sampler.partition_list.elementAt(p);
	    
	    
	    /******** if proposed matrix is empty, skip the partition  ******/
	    
	    if (pMatrix == null){
		
		/****** if this is a parameter change point propose 
		 ****** new evolutionary parameters ***************/
		
		if (temp_part.IsParameterChange()){
		    parameter_start = p;
		    pMatrix = temp_part.cMatrix.Proposal(sampler.set);
		    if (pMatrix != null){
			pHyperParameter = ProposeHyperParameter(temp_part.cHyperParameter,1);
			pLogHyperParameterPrior = sampler.prior_inf.Prior_on_mu(pHyperParameter);
			pLog_alpha_prior = sampler.prior_inf.Prior_on_alpha(pMatrix.v[0]);
			pPartialLogLikelihood[p] = sampler.PostTree[temp_part.cTree].LogLikelihood0(pMatrix, pHyperParameter, temp_part.counts, temp_part.data, pMatrix.pi, sampler.handleGaps, 0, 0);
			pLogLikelihood = pPartialLogLikelihood[p];
			cLogLikelihood = temp_part.cPartialLogLikelihood;
		    }
		}
	    }
	    else{
		
		
		/********* if this is a parameter change point do M-H step *********/

		if (temp_part.IsParameterChange()){
		    
		    /******* step back to get current parameters ********/
		    
		    Partition previous_part = (Partition)sampler.partition_list.elementAt(p-1);

		    double cEP = previous_part.cMatrix.v[0];
		    double pEP = pMatrix.v[0];

		    logRatio = pLogLikelihood - cLogLikelihood
			+ Math.log(pHyperParameter) - Math.log(previous_part.cHyperParameter)
			+ Math.log(pEP) - Math.log(cEP)
			+ pLog_alpha_prior  
			- previous_part.log_alpha_prior
			+ pLogHyperParameterPrior - previous_part.cPartialLogHyperParameterPrior;
		    
		    sampler.tries[1]++;
		    
		    /******** update all previous partitions since 
		     ******** last parameter change point *********/

		    if( sampler.LogMHAccept( logRatio, sampler.set.unif.nextDouble() ) ) {

			sampler.acceptancerate[1]++;
			
			/******** update all previous partitions since 
			 ******** last parameter change point *********/
			
			
			Partition update_part;
			for (int j = parameter_start; j < p; j++){
			    update_part = (Partition)sampler.partition_list.elementAt(j);
			    
			    update_part.cMatrix = pMatrix;
			    update_part.cHyperParameter = pHyperParameter; 
			    update_part.cPartialLogHyperParameterPrior = pLogHyperParameterPrior;
			    update_part.log_alpha_prior = pLog_alpha_prior;
			    update_part.cPartialLogLikelihood = pPartialLogLikelihood[j];
			    
			}
			
		    }
		    
		    /********** Start new chain of equal paremater segments ********/
		    
		    parameter_start = p;
		    pMatrix = temp_part.cMatrix.Proposal(sampler.set);
		    
		    if (pMatrix != null){
			pHyperParameter = ProposeHyperParameter(temp_part.cHyperParameter,1);
			pLogHyperParameterPrior = sampler.prior_inf.Prior_on_mu(pHyperParameter);
			pLog_alpha_prior = sampler.prior_inf.Prior_on_alpha(pMatrix.v[0]);
			pPartialLogLikelihood[p] = sampler.PostTree[temp_part.cTree].LogLikelihood0(pMatrix, pHyperParameter, temp_part.counts, temp_part.data, pMatrix.pi, sampler.handleGaps, 0, 0);
			pLogLikelihood = pPartialLogLikelihood[p];
			cLogLikelihood = temp_part.cPartialLogLikelihood;
		    }

		}
		
		/******** if this is not a parameter change point, add current 
		 ******** and proposal likelihood to corresponding 
		 ******** comultative likelihoods****************************/
		
		else{
		    
		    pPartialLogLikelihood[p] = sampler.PostTree[temp_part.cTree].LogLikelihood0(pMatrix, pHyperParameter, temp_part.counts, temp_part.data, pMatrix.pi, sampler.handleGaps, 0, 0);
		    
		    pLogLikelihood += pPartialLogLikelihood[p];
		    cLogLikelihood += temp_part.cPartialLogLikelihood;
		}
		
	    }
	}
	
	if (pMatrix != null){

	    /********** Finish update on the right end ***************/

	    double cEP = temp_part.cMatrix.v[0];
	    double pEP = pMatrix.v[0];
	    	
	    logRatio = pLogLikelihood - cLogLikelihood
		+ Math.log(pHyperParameter) - Math.log(temp_part.cHyperParameter)
		+ Math.log(pEP) - Math.log(cEP)
		+ pLog_alpha_prior
		- temp_part.log_alpha_prior
		+ pLogHyperParameterPrior - temp_part.cPartialLogHyperParameterPrior;


	    //System.out.println("prop likelihood: " + pLogLikelihood);

	    sampler.tries[1]++;
	    
	    if( sampler.LogMHAccept( logRatio, sampler.set.unif.nextDouble() ) ) {
		
		sampler.acceptancerate[1]++;
		
		/******** update all previous partitions since 
		 ******** last parameter change point *********/
		
		Partition update_part;
		for (int j = parameter_start; j < sampler.partition_list.size(); j++){

		    update_part = (Partition)sampler.partition_list.elementAt(j);
		    
		    
		    update_part.cMatrix = pMatrix;
		    update_part.cHyperParameter = pHyperParameter; 
		    update_part.cPartialLogHyperParameterPrior = pLogHyperParameterPrior;
		    update_part.log_alpha_prior = pLog_alpha_prior;
		    update_part.cPartialLogLikelihood = pPartialLogLikelihood[j];
		}
		
	    }
	}
	
    }

    double ProposeHyperParameter(double current_value, double lambda){
	
	return current_value*Math.exp(lambda*(sampler.set.unif.nextDouble() - 0.5));
	
    }
    

}
    
