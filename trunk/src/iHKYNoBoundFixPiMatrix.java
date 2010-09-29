/**
* Hasegawa, Kishino and Yano (1985) evolutionary matrix with [0,\infty)
limits, fixed frequencies and integrated integrated branch lengths.
*/


public class iHKYNoBoundFixPiMatrix extends QMatrix {
    private double pR;
    private double pY;
    private double b;
    private double c;
    private double d;
    
    public iHKYNoBoundFixPiMatrix() {
	// Starts at the JC69 constraint
	v = new double[3];
	v[0] = 1.0 / 3.0;
	v[1] = 1.0 / 3.0;
	v[2] = 1.0 / 3.0;
	//v[1] = v[0] = 1.0 / 2.0;
	pi = new double[4];
	pi[0] = 0.25;
	pi[1] = 0.25;
	pi[2] = 0.25;
	pi[3] = 0.25;
	pR = pi[0] + pi[1]; // pi_A + pi_G
	pY = pi[2] + pi[3]; // pi_C + pi_C
	b  = 0.5 / (v[0]*pi[0]*pi[1] + v[1]*pi[2]*pi[3] + pR*pY);
	c  = (pY + v[0]*pR)*b;
	d  = (pR + v[1]*pY)*b;
    }
    
    public iHKYNoBoundFixPiMatrix(double[] inV, double[] inPi) {
	//System.err.println("HKY");
	v = inV;
	pi = inPi;
	pR = pi[0] + pi[1]; // pi_A + pi_G
	pY = pi[2] + pi[3]; // pi_C + pi_T
	b  = 0.5 / (v[0]*pi[0]*pi[1] + v[0]*pi[2]*pi[3] + pR*pY);
	c  = (pY + v[0]*pR)*b;
	d  = (pR + v[0]*pY)*b;
	//System.err.println(v.length);
	//System.exit(-1);
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
        

    public QMatrix Proposal(Settings set) {
	double lambda = 1;
	double[] p = new double[2];
	p[0] = v[0]*Math.exp(lambda*(set.unif.nextDouble() - 0.5));
	p[1] = p[0];
	
	return new iHKYNoBoundFixPiMatrix(p,pi);
    }
    
    public static void main(String args[]) {
	double[] p  = { 0.000001 };
	double[] pi = { 0.25, 0.25, 0.25, 0.25 };
	QMatrix q = new iHKYNoBoundFixPiMatrix(p,pi);
	System.out.println(q.v[0]);
	System.out.println(q.toString(1));
	//Settings set = new Settings(1);
	//q = q.Proposal(set);
	//System.out.println(q.toString(1.0));
	//q = q.Proposal(set);
	//System.out.println(q.toString(1.0));
    }
}
