/**
* Provides basic functions for Metropolis-Hasting sampling.
*/
public class Sampler {
	
	/**
	* Determines Metropolis-Hasting acceptance 
	* based on the natural log ratio of posterior probabilities.
	*/
    protected static boolean LogMHAccept(double logRatio, double random) {
	if( (logRatio > 0) || (random < Math.exp(logRatio)) )
	    return true;
	else
	    return false;
    }
	
	protected static double oneOverSqrtTwoPi = 0.399;
	protected static double logOneOverSqrtTwoPi = -0.9189;
	
	/**
	* Calculates the natural log probability of moving from->to in a Gaussian field.
	*/
	protected static double logAbsNormalDriverDensity(double from, double to, double variance) {
		double tau = 1.0 / variance;
		double sqrtTau = Math.sqrt(tau);
		return Math.log( oneOverSqrtTwoPi*sqrtTau*( Math.exp(-0.5*tau*(from-to)*(from-to))
			+ Math.exp(-0.5*tau*(from+to)*(from+to))
			) );
	}
	
	/**
	 * Calculates the natural log density at x under Normal distribution.
	 */
    protected static double logStandardNormalDensity(double x) {
		return logOneOverSqrtTwoPi - 0.5 * x * x;
		//return 0.0;
	}
	
	/**
	* Computes the logit transform of x.
	*/
	protected static double Logit(double x) {
		return Math.log( x / (1.0 - x) );
	}

	/**
	* Computes the inverse of the logit transform of x.
	*/
	protected static double InvLogit(double x) {
		double tmp = Math.exp(x);
		return tmp / (1.0 + tmp);
	}

}
