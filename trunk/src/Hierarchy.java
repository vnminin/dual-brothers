import java.util.*;
import java.io.*;
import corejava.*;


/**
* Updating hyperparameters for evolutionary parameters
*/

public class Hierarchy {

    DCPSampler sampler;

    public Hierarchy(DCPSampler s){
	sampler  = s;
    }

    
    protected void Update(){


	/****************** calculate alpha/mu mean statistics ***********************/

	Partition temp_part;

 	double sum_log_alpha = 0;
	double sum_log_mu = 0;

	for (int i = 0; i < sampler.partition_list.size(); i++){
	    temp_part = (Partition)sampler.partition_list.elementAt(i);
	    if (temp_part.IsParameterChange()){
		sum_log_alpha += Math.log(temp_part.cMatrix.v[0]);
		sum_log_mu += Math.log(temp_part.cHyperParameter);
	    }
	}

	/**************** update hyperparameter alpha_mean  **************************/
	

	double alpha_sample_mean = (sampler.prior_inf.alpha_mean_mean*sampler.prior_inf.alpha_variance + 
			      sampler.prior_inf.alpha_mean_variance*sum_log_alpha)/
	    (sampler.prior_inf.alpha_variance + 
	     sampler.parameter_changes*sampler.prior_inf.alpha_mean_variance);

	double alpha_sample_variance = 1.0/(sampler.parameter_changes/sampler.prior_inf.alpha_variance 
				      + 1.0/sampler.prior_inf.alpha_mean_variance);


	sampler.prior_inf.alpha_mean = Math.sqrt(alpha_sample_variance)*sampler.set.norm01.nextDouble() + alpha_sample_mean;


	/**************** update hyperparameter mu_mean  **************************/
	
	
	double mu_sample_mean = (sampler.prior_inf.mu_mean_mean*sampler.prior_inf.mu_variance + 
			      sampler.prior_inf.mu_mean_variance*sum_log_mu)/
	    (sampler.prior_inf.mu_variance + 
	     sampler.parameter_changes*sampler.prior_inf.mu_mean_variance);
	
	double mu_sample_variance = 1.0/(sampler.parameter_changes/sampler.prior_inf.mu_variance 
					 + 1.0/sampler.prior_inf.mu_mean_variance);
	
	
	sampler.prior_inf.mu_mean = Math.sqrt(mu_sample_variance)*sampler.set.norm01.nextDouble() + mu_sample_mean;
	

	/****************** calculate alpha/mu variance statistics ***********************/


	double dev_log_alpha = 0;
	double dev_log_mu = 0;

	for (int i = 0; i < sampler.partition_list.size(); i++){
	    temp_part = (Partition)sampler.partition_list.elementAt(i);
	    if (temp_part.IsParameterChange()){
		dev_log_alpha += (Math.log(temp_part.cMatrix.v[0]) - sampler.prior_inf.alpha_mean)*
		    (Math.log(temp_part.cMatrix.v[0]) - sampler.prior_inf.alpha_mean);

		dev_log_mu += (Math.log(temp_part.cHyperParameter) - sampler.prior_inf.mu_mean)*
		    (Math.log(temp_part.cHyperParameter) - sampler.prior_inf.mu_mean);
		
	    }
	}

	/***************** update hyperparameter alpha_variance ************************/

	
	double alpha_sample_shape = sampler.prior_inf.alpha_precision_shape + sampler.parameter_changes/2.0;
 
	double alpha_sample_scale = sampler.prior_inf.alpha_precision_scale + 0.5*dev_log_alpha;

	
	sampler.prior_inf.alpha_variance = 1.0/sampler.set.gammaV.nextDouble(alpha_sample_shape, alpha_sample_scale);
	
	
	


	/***************** update hyperparameter mu_variance ************************/

		
	double mu_sample_shape = sampler.prior_inf.mu_precision_shape + sampler.parameter_changes/2.0;
	
	double mu_sample_scale = sampler.prior_inf.mu_precision_scale + 0.5*dev_log_mu;

	
	sampler.prior_inf.mu_variance = 1.0/sampler.set.gammaV.nextDouble(mu_sample_shape, mu_sample_scale);


	/**************** update saved priors ********************/

	for (int j = 0; j < sampler.partition_list.size(); j ++){
	    Partition my_part = (Partition)sampler.partition_list.elementAt(j);
	    my_part.log_alpha_prior = sampler.prior_inf.Prior_on_alpha(my_part.cMatrix.v[0]);
	    my_part.cPartialLogHyperParameterPrior = sampler.prior_inf.Prior_on_mu(my_part.cHyperParameter);
	}
    }
}
