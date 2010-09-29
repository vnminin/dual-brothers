//import JSci.maths.statistics.*;

/**
* Hasegawa, Kishino and Yano (1985) evolutionary matrix.
*/
public class HKYMatrix extends QMatrix {
	private double pR;
	private double pY;
	private double b;
	private double c;
	private double d;

	public HKYMatrix() {
		// Starts at the JC69 constraint
		v = new double[2];
		v[0] = 1.0 / 3.0;
		v[1] = 1.0 / 3.0;
		pi = new double[4];
		pi[0] = 0.25;
		pi[1] = 0.25;
		pi[2] = 0.25;
		pi[3] = 0.25;
		preprocess();
	}

	public HKYMatrix(double[] inV, double[] inPi) {
		v = inV;
		pi = inPi;
		preprocess();
	}

	private void preprocess() {
		pR = pi[0] + pi[1]; // pi_A + pi_G
		pY = pi[2] + pi[3]; // pi_C + pi_T
		b  = v[1];
		c  = v[0]*pR + b*pY;
		d  = b*pR    + v[0]*pY;
	}

	public final double Pt(int i, int j, double t) {
		int sw = i*4 + j;
		switch( sw ) {
			case  0 : return pi[0] + Math.exp(-b*t)*pi[0]*pY/pR + Math.exp(-c*t)*pi[1]/pR;
			case  1 : return pi[1] + Math.exp(-b*t)*pi[1]*pY/pR - Math.exp(-c*t)*pi[1]/pR;
			case  2 : return pi[2] - Math.exp(-b*t)*pi[2];
			case  3 : return pi[3] - Math.exp(-b*t)*pi[3];
			case  4 : return pi[0] + Math.exp(-b*t)*pi[0]*pY/pR - Math.exp(-c*t)*pi[0]/pR;
			case  5 : return pi[1] + Math.exp(-b*t)*pi[1]*pY/pR + Math.exp(-c*t)*pi[0]/pR;
			case  6 : return pi[2] - Math.exp(-b*t)*pi[2];
			case  7 : return pi[3] - Math.exp(-b*t)*pi[3];
			case  8 : return pi[0] - Math.exp(-b*t)*pi[0];
			case  9 : return pi[1] - Math.exp(-b*t)*pi[1];
			case 10 : return pi[2] + Math.exp(-b*t)*pi[2]*pR/pY + Math.exp(-d*t)*pi[3]/pY;
			case 11 : return pi[3] + Math.exp(-b*t)*pi[3]*pR/pY - Math.exp(-d*t)*pi[3]/pY;
			case 12 : return pi[0] - Math.exp(-b*t)*pi[0];
			case 13 : return pi[1] - Math.exp(-b*t)*pi[1];
			case 14 : return pi[2] + Math.exp(-b*t)*pi[2]*pR/pY - Math.exp(-d*t)*pi[2]/pY;
			case 15 : return pi[3] + Math.exp(-b*t)*pi[3]*pR/pY + Math.exp(-d*t)*pi[2]/pY;
			default :
				System.err.println("Invalid mutation: "+i+" -> "+j);
				System.exit(-1);
		}
		return 0;
	}

	public final void PtVector(int i, double t, double[] mat) {
		switch( i ) {
			case 0 : mat[0] = pi[0] + Math.exp(-b*t)*pi[0]*pY/pR + Math.exp(-c*t)*pi[1]/pR;
				 mat[1] = pi[1] + Math.exp(-b*t)*pi[1]*pY/pR - Math.exp(-c*t)*pi[1]/pR;
				 mat[2] = pi[2] - Math.exp(-b*t)*pi[2];
				 mat[3] = pi[3] - Math.exp(-b*t)*pi[3];
				 return;
			case 1 : mat[0] = pi[0] + Math.exp(-b*t)*pi[0]*pY/pR - Math.exp(-c*t)*pi[0]/pR;
				 mat[1] = pi[1] + Math.exp(-b*t)*pi[1]*pY/pR + Math.exp(-c*t)*pi[0]/pR;
				 mat[2] = pi[2] - Math.exp(-b*t)*pi[2];
				 mat[3] = pi[3] - Math.exp(-b*t)*pi[3];
				 return;
			case 2 : mat[0] = pi[0] - Math.exp(-b*t)*pi[0];
				 mat[1] = pi[1] - Math.exp(-b*t)*pi[1];
				 mat[2] = pi[2] + Math.exp(-b*t)*pi[2]*pR/pY + Math.exp(-d*t)*pi[3]/pY;
				 mat[3] = pi[3] + Math.exp(-b*t)*pi[3]*pR/pY - Math.exp(-d*t)*pi[3]/pY;
				 return;
			case 3 : mat[0] = pi[0] - Math.exp(-b*t)*pi[0];
				 mat[1] = pi[1] - Math.exp(-b*t)*pi[1];
				 mat[2] = pi[2] + Math.exp(-b*t)*pi[2]*pR/pY - Math.exp(-d*t)*pi[2]/pY;
				 mat[3] = pi[3] + Math.exp(-b*t)*pi[3]*pR/pY + Math.exp(-d*t)*pi[2]/pY;
		}
		return;
	}

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
		return new HKYMatrix(p,npi);
	}

	public static void main(String args[]) {
		double[] p  = { 0.4, 0.3 };
		p[0] = //5.0 / 6.0;
		1.0 / 2.0;
		p[1] = (1.0 - p[0]) / 2.0;
		
		double[] pi = { 0.25, 0.25, 0.25, 0.25 };
		QMatrix q = new HKYMatrix(p,pi);
		System.out.println(q.toString(1.0));
		//Settings set = new Settings(1);
		//q = q.Proposal(set);
		//System.out.println(q.toString(1.0));
		//q = q.Proposal(set);
		//System.out.println(q.toString(1.0));
	}
}
