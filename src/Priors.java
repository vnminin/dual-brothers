import java.util.*;
import java.io.*;
import corejava.*;


/**
* Provides methods to compute and store prior probabilities
* for parameters of dual Multiple Change-Point model
*/

public class Priors implements PriorConst{

    /**
     * prior probability of recombination locations
     */
    public double recombination[];

    /**
     * indicator of informative prior
     */
    
    public boolean inf_recomb;
    
    /**
     * prior mean of alpha [log(alpha) ~ Normal(alpha_mean, alpha_variance)]
     */
    protected double alpha_mean;
    
    
    /** prior variance of alpha [log(alpha) ~ Normal(alpha_mean, alpha_variance)]
     *
     */
    protected double alpha_variance;


    /**
     * prior mean of mu [mu ~ Normal(mu_mean, mu_variance)]
     */
    protected double mu_mean;
    
    
    /** prior variance of mu [mu ~ Normal(alpha_mean, alpha_variance)]
     *
     */
    protected double mu_variance;

    


    public Priors(){
	recombination = new double[1];
    }



    public Priors(int sequence_length, double lambda){
	recombination = new double[sequence_length];

	for (int i =0; i < sequence_length; i++){
	    recombination[i] = lambda/(1.0*sequence_length);
	}
    }


    public Priors(int sequence_length, String prior_file){
	recombination = new double[sequence_length];
	
	BufferedReader input = null;
	try {
	    input = new BufferedReader( new FileReader(prior_file) );
	} catch (FileNotFoundException e) {
	    System.err.println("Error opening priorfile: " + prior_file);
	    System.exit(-1);
	}
	
	String my_string = null;
	int my_site = 0;

	try{
	    while( (my_string = input.readLine()) != null && my_site < sequence_length){
		StringTokenizer words = new StringTokenizer(my_string);
		String site = words.nextToken();
		double recomb_prior = Double.valueOf(words.nextToken()).doubleValue();
		
	

		if (recomb_prior <= 1){
		    recombination[my_site] = recomb_prior;
		    my_site++;
		}
		else{
		    System.err.println("Prior probability of recombination at line" + 
				       (my_site+1) + "must be between 0 and 1");
		}
		
	    }

	} catch (Exception e) {
	    System.err.println("Unable to parse line" + my_string);
	    System.exit(-1);
	    
	}

	
	
    }


    public Priors(Settings settings){

	if (settings.hyper_update){
	    alpha_mean = 0;
	    alpha_variance = 10;
	    mu_mean = 0;
	    mu_variance = 10;
	}
	else{
	    alpha_mean = settings.subst_hyper_mean;
	    alpha_variance = settings.subst_hyper_variance;
	    mu_mean = settings.diver_hyper_mean;
	    mu_variance =settings.diver_hyper_variance;
	}

        recombination = null;
	inf_recomb = false;
	
    }



    public Priors(Settings settings, String prior_file, int sequence_length){
    
	if (settings.hyper_update){
	    alpha_mean = 0;
	    alpha_variance = 10;
	    mu_mean = 0;
	    mu_variance = 10;
	}
	else{
	    alpha_mean = settings.subst_hyper_mean;
	    alpha_variance = settings.subst_hyper_variance;
	    mu_mean = settings.diver_hyper_mean;
	    mu_variance =settings.diver_hyper_variance;
	}
	
	recombination = new double[sequence_length];

	inf_recomb = true;
	
	BufferedReader input = null;
	try {
	    input = new BufferedReader( new FileReader(prior_file) );
	} catch (FileNotFoundException e) {
	    System.err.println("Error opening priorfile: " + prior_file);
	    System.exit(-1);
	}
	
	String my_string = null;
	int my_site = 0;

	try{
	    while( (my_string = input.readLine()) != null && my_site < sequence_length){
		StringTokenizer words = new StringTokenizer(my_string);
		//String site = words.nextToken();
		double recomb_prior = Double.valueOf(words.nextToken()).doubleValue();
		
		if (recomb_prior <= 1){
		    recombination[my_site] = recomb_prior;
		    my_site++;
		}
		else{
		    System.err.println("Prior probability of recombination at line" + 
				       (my_site+1) + "must be between 0 and 1");
		}
		
	    }

	} catch (Exception e) {
	    System.err.println("Unable to parse line" + my_string);
	    System.exit(-1);
	    
	}
	
    }
 


    public double Prior_on_mu(double mu){
	return lnLogNormalDensity(mu, mu_mean, mu_variance);
    }
    


    public double Prior_on_alpha(double alpha){
	
	return lnLogNormalDensity(alpha, alpha_mean, alpha_variance);
    }
    
 
    
    public static double lnExpDensity(double x) {
	return -x; // Exp(1)
    }
    
    
    public static double lnLogNormalDensity(double x, double mean, double sigma){

	double log_of_x = Math.log(x);
	return -1.0*(log_of_x - mean)*(log_of_x - mean)/(2.0*sigma) - log_of_x
	    - 0.5*Math.log(sigma) - normal_const;
    }
    
    
    public static double lnNormalDensity(double x, double mean, double sigma){
	
	return -0.5*Math.log(sigma) - (x - mean)*(x - mean)/(2.0*sigma) - normal_const;
	
    }
    
    
    


}
