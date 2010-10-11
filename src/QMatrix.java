import java.io.Serializable;

import cern.jet.random.*;

/**
* Generic evolutionary matrix.
*
*/
public abstract class QMatrix  implements Serializable {

	/**
	* Evolutionary parameters.
	*/
	public double[] v;
	
	/**
	* Nucleotide stationary distribution.
	*/
	public double[] pi;

	/**
	* Returns the finite time transition probability of 
	* nucleotide i mutating to nucleotide j in time t.
	*/
	public abstract double Pt(int i, int j, double t);
	
	/**
	* Calculates the finite time transition probabilities of 
	* nucleotide i mutating to any state.  Replaces entries in mat 
	* with updated probabilities.
	*/
	public abstract void PtVector(int i, double t, double[] mat);

	/**
	* NOT IMPLEMENTED.
	*/
	public QMatrix() {
	}

	/**
	* NOT IMPLEMENTED.
	*/
	public double LogPrior() {
		return 0.0;
	}

	/**
	* Returns a [0,1) RV based on mean.  Draws a normally distributed
	* RV centered at mean from norm and reflects the result about zero and one.
	*/
	public static double gennor01(double mean, AbstractDistribution norm) {
		double r = mean + norm.nextDouble();
		while( (r > 1.0) || (r < 0.0) ) {
			if( r > 1.0 )
				r = 2.0 - r;
			else
				r = -r;
		}
		return r;
	}

	/**
	* Returns the 4x4 finite time transition matrix
	* evaluated at time t.
	*/
	public String toString(double t) {
		String rtn = new String("\tA\tG\tC\tT\r\n");
		rtn += "A";
		for(int i=0; i<4; i++)
			rtn += "\t" + Pt(0,i,t);
		rtn += "\r\n";
		rtn += "G";
		for(int i=0; i<4; i++)
			rtn += "\t" + Pt(1,i,t);
		rtn += "\r\n";
		rtn += "C";
		for(int i=0; i<4; i++)
			rtn += "\t" + Pt(2,i,t);
		rtn += "\r\n";
		rtn += "T";
		for(int i=0; i<4; i++)
			rtn += "\t" + Pt(3,i,t);
		rtn += "\r\n";
		return rtn;
	}

	/**
	* Returns a randomly drawn evolutionary matrix based on this
	* matrix.  Used in MH update proposals.
	*/
	public abstract QMatrix Proposal(Settings set);

}
