//import JSci.maths.statistics.*;

public class iHKYFixPiMatrix extends QMatrix {
	private double pR;
	private double pY;
	private double b;
	private double c;
	private double d;

	public iHKYFixPiMatrix() {
		// Starts at the JC69 constraint
		v = new double[2];
		v[0] = 1.0 / 3.0;
		v[1] = 1.0 / 3.0;
		//v[1] = v[0] = 1.0 / 2.0;
		pi = new double[4];
		pi[0] = 0.25;
		pi[1] = 0.25;
		pi[2] = 0.25;
		pi[3] = 0.25;
		pR = pi[0] + pi[1]; // pi_A + pi_G
		pY = pi[2] + pi[3]; // pi_C + pi_C
		b  = v[1];
		c  = v[0]*pR + b*pY;
		d  = b*pR    + v[0]*pY;
	}

	public iHKYFixPiMatrix(double[] inV, double[] inPi) {
		//System.err.println("Yo.");
		//System.exit(-1);
		v = inV;
		pi = inPi;
		pR = pi[0] + pi[1]; // pi_A + pi_G
		pY = pi[2] + pi[3]; // pi_C + pi_T
		b  = v[1];
		c  = v[0]*pR + b*pY;
		d  = b*pR    + v[0]*pY;
	}

	public final double Pt(int i, int j, double t) {
		int sw = i*4 + j;
		switch( sw ) {
			case  0 : return pi[0] + pi[0]*pY/pR/(1+b*t) + pi[1]/pR/(1+c*t);
			case  1 : return pi[1] + pi[1]*pY/pR/(1+b*t) - pi[1]/pR/(1+c*t);
			case  2 : return pi[2] - pi[2]/(1+b*t);
			case  3 : return pi[3] - pi[3]/(1+b*t);
			case  4 : return pi[0] + pi[0]*pY/pR/(1+b*t) - pi[0]/pR/(1+c*t);
			case  5 : return pi[1] + pi[1]*pY/pR/(1+b*t) + pi[0]/pR/(1+c*t);
			case  6 : return pi[2] - pi[2]/(1+b*t);
			case  7 : return pi[3] - pi[3]/(1+b*t);
			case  8 : return pi[0] - pi[0]/(1+b*t);
			case  9 : return pi[1] - pi[1]/(1+b*t);
			case 10 : return pi[2] + pi[2]*pR/pY/(1+b*t) + pi[3]/pY/(1+d*t);
			case 11 : return pi[3] + pi[3]*pR/pY/(1+b*t) - pi[3]/pY/(1+d*t);
			case 12 : return pi[0] - pi[0]/(1+b*t);
			case 13 : return pi[1] - pi[1]/(1+b*t);
			case 14 : return pi[2] + pi[2]*pR/pY/(1+b*t) - pi[2]/pY/(1+d*t);
			case 15 : return pi[3] + pi[3]*pR/pY/(1+b*t) + pi[2]/pY/(1+d*t);
			default :
				System.err.println("Invalid mutation: "+i+" -> "+j);
				System.exit(-1);
		}
		return 0;
	}

	public final void PtVector(int i, double t, double[] mat) {
		switch( i ) {
			case 0 :
				 mat[0] = pi[0] + pi[0]*pY/pR/(1+b*t) + pi[1]/pR/(1+c*t);
				 mat[1] = pi[1] + pi[1]*pY/pR/(1+b*t) - pi[1]/pR/(1+c*t);
				 mat[2] = pi[2] - pi[2]/(1+b*t);
				 mat[3] = pi[3] - pi[3]/(1+b*t);
				 return;
			case 1 :
				 mat[0] = pi[0] + pi[0]*pY/pR/(1+b*t) - pi[0]/pR/(1+c*t);
				 mat[1] = pi[1] + pi[1]*pY/pR/(1+b*t) + pi[0]/pR/(1+c*t);
				 mat[2] = pi[2] - pi[2]/(1+b*t);
				 mat[3] = pi[3] - pi[3]/(1+b*t);
				 return;
			case 2 :
				 mat[0] = pi[0] - pi[0]/(1+b*t);
				 mat[1] = pi[1] - pi[1]/(1+b*t);
				 mat[2] = pi[2] + pi[2]*pR/pY/(1+b*t) + pi[3]/pY/(1+d*t);
				 mat[3] = pi[3] + pi[3]*pR/pY/(1+b*t) - pi[3]/pY/(1+d*t);
				 return;
			case 3 :
				 mat[0] = pi[0] - pi[0]/(1+b*t);
				 mat[1] = pi[1] - pi[1]/(1+b*t);
				 mat[2] = pi[2] + pi[2]*pR/pY/(1+b*t) - pi[2]/pY/(1+d*t);
				 mat[3] = pi[3] + pi[3]*pR/pY/(1+b*t) + pi[2]/pY/(1+d*t);
		}
		return;
	}

/*	public QMatrix ProposedEPJump(NormalDistribution normal) {
		double[] p = new double[2];
		System.err.println("Deprecated.");
		System.exit(-1);
		return null; // All other models do not yet have proposals.
	}

	public QMatrix ProposedRhoJump(NormalDistribution normal) {
		System.err.println("Deprecated.");
		System.exit(-1);
		return null;
	}
*/
	public QMatrix Proposal(Settings set) {
		//System.err.println("Here");
		double[] p = new double[2];
	//	double[] npi = new double[4];
		p[0] = gennor01(v[0],set.normEP);
		p[1] = (1.0 - p[0]) / 2.0;
		return new iHKYFixPiMatrix(p,pi);
	}
/*
	public QMatrix Proposal(Settings set) {
		double[] p = new double[2];
		double[] npi = new double[4];
		p[0] = gennor01(v[0],set.normEP);
		p[1] = (1.0 - p[0]) / 2.0;

		// Use a tri-variate normal driver reflected about 0 and 1
		npi[3] = 1.0;
		for(int i=0; i<3; i++) {
			npi[i] = gennor01(pi[i],set.normPi);
			npi[3] -= npi[i];
		}
		if( npi[3] < 0.0 )
			return null;
		return new iHKYFixPiMatrix(p,npi);
	}
*/
	public static void main(String args[]) {
		double[] p  = { 0.5, 0.25 };
		p[0] = 0.5;
		p[1] = (1.0 - p[0]) / 2.0;
		double[] pi = { 0.25, 0.25, 0.25, 0.25 };
		QMatrix q = new iHKYFixPiMatrix(p,pi);
		System.out.println(q.toString(1));
		//Settings set = new Settings(1);
		//q = q.Proposal(set);
		//System.out.println(q.toString(1.0));
		//q = q.Proposal(set);
		//System.out.println(q.toString(1.0));
	}
}
