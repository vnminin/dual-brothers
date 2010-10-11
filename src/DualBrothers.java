import java.io.*;

/**
* Provides a commandline interface for DCPSampler.
* Runs the model with two change-point processes.
*/
public class DualBrothers {

	public static void printUsage()
	{
	    System.err.println();
	    System.err.println("java DualChange <seed integer> <commandfile name> <datafile name> <outfile name>");
	    System.err.println("or:");
	    System.err.println("java DualChange <checkpoint file>");
	    System.err.println();
	    new Settings();
	}

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
	
	File ckpnt_file = null;
	if( args.length == 1)
	{
		// assume the user is trying to resume a checkpointed run
		try{
			ckpnt_file = new File(args[0]);
			if(!ckpnt_file.exists())
				throw new RuntimeException("The checkpoint file \""+args[0]+"\" does not exist");
		}catch(Exception e){
			printUsage();
			System.err.println("Error opening the checkpoint file \""+args[0]+"\"");
		    System.exit(-1);
		}
	}else{
	try {
	    seed = Integer.parseInt(args[0]);
	    cname = args[1];
	    dname = args[2];
	    oname = args[3];
	    if (prior){
		pname = args[4];
	    }
	} catch (Exception e) {
		printUsage();
	    System.exit(-1);
	}
	}
	
	Settings set = null;
	DCPSampler m = null;
	
	// try to restore from a checkpoint
	if(ckpnt_file == null)
	{
		set = new Settings(seed,cname);
		try{
			ckpnt_file = new File(set.checkpoint_filename);
		}catch(Exception e){}
	}
	if(ckpnt_file != null && ckpnt_file.exists())
	{
    	try
    	{
    		FileInputStream fis = new FileInputStream(ckpnt_file);
    		ObjectInputStream ois = new ObjectInputStream(fis);
    		m = (DCPSampler)ois.readObject();
    		if(m.JumpNumber >= m.set.length)
    		{
    			System.err.println("\nThis run has completed.");
    			System.err.println("To start a new run delete the checkpoint file " + ckpnt_file);
    			System.err.println("Exiting now.");
    			System.exit(-1);
    		}
    		System.err.println("Attempting to resume a previous run from checkpoint file " + ckpnt_file);
    		PrintWriter output = resumeOutput(oname, m);
    		m.resumeState(output);
    		if(args.length > 1)
    			System.err.println("Delete the checkpoint file to prevent this behavior");
    	}
    	catch(FileNotFoundException fnfe){}
    	catch(IOException ioe){}
    	catch(ClassNotFoundException cnfe){}
	}

	// if we couldn't read a checkpoint then start a new run
	if(m == null)
	{

	CountStatistic cs = new CountStatistic(dname);
	PrintWriter output = OpenOutput(oname);
	Priors priors;

	prior = false;

	if (prior){
	    priors = new Priors(set, pname, cs.lenseq);
	}else{
	    priors = new Priors(set);
	}

		m = new DCPSampler(cs, set, output, priors);
	}
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
    
    /** 
     * Attempts to open the output file and seek to the previous checkpoint state
     * @param oname output file name
     * @param s		The sampler being resumed
     */
    
    public static PrintWriter resumeOutput(String oname, DCPSampler s) {
	PrintWriter output = null;
	try {
		RandomAccessFile raf = new RandomAccessFile(oname, "rw");
		// read up to the last checkpointed sample
		long ckpnt_samplecount = (s.JumpNumber - s.set.burnin) / s.set.subsample;
		ckpnt_samplecount = ckpnt_samplecount < 0 ? 0 : ckpnt_samplecount;
		long lineI = 0;
		try{
			for(; lineI <= ckpnt_samplecount; lineI++)
				raf.readLine();
		}catch(IOException ioe){}
		if(lineI == 0)
		{
			System.err.println("Warning: starting new posterior sample file.");
			System.err.println("Concatenate this file with any previous output before using TopoProfile or EPProfile");
		}
		else if(lineI <= ckpnt_samplecount)
		{
			System.err.println("\nError resuming.  Part of \"" + oname + "\" is missing or corrupt.");
			System.exit(-1);
		}
		// truncate here
		raf.setLength(raf.getFilePointer());
		raf.close();
		// open in append mode
	    output = new PrintWriter(new BufferedWriter(new FileWriter(oname, true)));
	} catch (IOException e) {
	    System.err.println("Unable to open output file: "+oname);
	    System.exit(-1);
	}
	return output;
    }

    
}
