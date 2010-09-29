import java.io.*;

/**
* Provides a commandline interface for DCPSampler.
* Runs the model with two change-point processes.
*/
public class DualBrothers {

    public static void main(String[] args) {
	int seed = 0;
	String cname = null;
	String dname = null;
	String oname = null;
	String pname = null;
	boolean prior = false;

	if (args.length == 5){
	    prior = true;
	}
	
	try {
	    seed = Integer.parseInt(args[0]);
	    cname = args[1];
	    dname = args[2];
	    oname = args[3];
	    if (prior){
		pname = args[4];
	    }
	} catch (Exception e) {
	    System.err.println();
	    System.err.println("java DualChange <seed integer> <commandfile name> <datafile name> <outfile name>");
	    System.err.println();
	    Settings temp_set = new Settings();
	    System.exit(-1);
	}
	
	Settings set = new Settings(seed,cname);
	CountStatistic cs = new CountStatistic(dname);
	PrintWriter output = OpenOutput(oname);
	Priors priors;

	prior = false;

	if (prior){
	    priors = new Priors(set, pname, cs.lenseq);
	}else{
	    priors = new Priors(set);
	}

	DCPSampler m = new DCPSampler(cs, set, output, priors);
	m.run();
		
	
	
    }

    /** Opens an output file to write to
     * @param oname output file name
     */
    
    public static PrintWriter OpenOutput(String oname) {
	PrintWriter output = null;
	try {
	    output = new PrintWriter(new BufferedWriter(new FileWriter(oname)));
	} catch (IOException e) {
	    System.err.println("Unable to open output file: "+oname);
	    System.exit(-1);
	}
	return output;
    }
    

    
}
