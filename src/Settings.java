import cern.jet.random.engine.*;
import cern.jet.random.*;
import java.io.*;
import java.util.*;

/**
* Reads in configuration file and initializes pseudo-random number generators.
*/
public class Settings {

    /**
     * number of MCMC iterations
     */
    public int length;
    
    /**
     * length of the burn-in period.
     */
    public int burnin;

    /** 
     * subsampling step 
     */
    public int subsample;

    /**
     * max length of the jump during a change(break)-point update
     */
    public int lenWindow;
    
    /**
     * mean of the prior Poisson distribution of number of topology break-points
     */
    public double top_lambda;    

    /**
     * mean of the prior Poisson distribution of number of substitution change-points
     */
    public double par_lambda;

    /**
     * tunable constant for proposing new transition:transversion ratios
     * when adding a new substitution change-point
     */

    public double sigmaAlpha;

    /**
     * tunable constant for proposing new average divergences
     * when adding a new substitution change-point
     */

    public double sigmaMu;

    /**
     * tunable constant for rescaling probabilities of selecting
     * a particular move during rjMCMC
     */
    public double C;

    /**
     * Mean of log(kappa) ~ Normal
     */
    public double subst_hyper_mean;

    /**
     * Variance of log(kappa) ~ Normal
     */
    public double subst_hyper_variance;

    /**
     * Mean of log(mu) ~ Normal
     */
    public double diver_hyper_mean;

    /**
     * Variance of log(mu) ~ Normal
     */
    public double diver_hyper_variance;


    /**
     * Indicator of whether topology break-points are to be updated
     */

    public boolean top_breaks = true;

    /**
     * Indicator of whether substitution process change-points are to be updated
     */

    public boolean par_breaks = true;


    /**
     * Indicator of whether hyper parameters should be updated of fixed
     */
    
    public boolean hyper_update;

    /**
     * Starting tree topology
     */
    public Tree sTree;


    /**
     * List of trees to hold tree list supplied by user
     */

    String[] tree_list;

    /**
     * Indicator of "start_tree" option presence
     */

    public boolean start_tree_present;

    /**
     *Indicator of "tree_file" option presence
     */

    public boolean tree_file_present;


    /**
     * temporary integer to hold indicators
     */

    public int temp_int;

    public boolean hyper1 = false;    
    public boolean hyper2 = false;
    public boolean hyper3 = false;
    public boolean hyper4 = false;
   
    

    // distributions need during computations
    
    public AbstractDistribution normT1;
    public AbstractDistribution normT2;
    public AbstractDistribution normEP;
    public AbstractDistribution normPi;
    public AbstractDistribution normMu;
    public AbstractDistribution norm01;
    public AbstractDistribution normHyperEP;
    public AbstractDistribution normUV;
    public AbstractDistribution normP;
    public AbstractDistribution normNP;
    public AbstractDistribution unif;
    public AbstractDistribution dispersed_unif;
    public Gamma gamma;
    public Gamma gammaV;
    public Beta beta;
    public ChiSquare chisquare;
    public Normal normal;
    public RandomEngine rengine;
    
    
     public Settings() {
	 PrintOptions();
     }


    
    /** Loads settings fromthe command file and a seed
     * @param seed seed for the random number generator
     * @param cname name of the command file
     */

    public Settings(int seed, String cname) {
	rengine = new MersenneTwister(seed);

	// Fill in default parameters
	length    = 1100000;
	burnin    = 100000;
	subsample = 50;
	top_lambda = 1;
	par_lambda = 1;
	C = 0.45;
	lenWindow = 5;
	
	sigmaAlpha = 0.75;
	sigmaMu = 0.75;

	sTree = null;
	tree_list = null;
	
	// by defaul update break- and change-points

	top_breaks = true;
	par_breaks = true;

	// by default estimate hyperparameters

	hyper_update = true;

	hyper1 = false;    
	hyper2 = false;
	hyper3 = false;
	hyper4 = false;

	start_tree_present = false;
	tree_file_present = false;
	
	ReadCmdfile(cname);


	// check if options for updating hyper parameters were 
	// specified correctly

	if (hyper1 && hyper2 && hyper3 && hyper4){
	    hyper_update = false;
	    
	}else{
	    if (hyper1 || hyper2 || hyper3 || hyper4){
		PrintOptions();
		System.err.println("\n\nERROR! \nYou must either specify ALL hyper parameters");
		System.err.println("in command file or ommit them ALL\n");
		System.exit(-1);
	    }
	} 


	// check if strarting tree or tree list were supplied, but not both

	if ((start_tree_present && tree_file_present) || 
	    ((!start_tree_present) && (!tree_file_present))){
	    PrintOptions();
	    System.err.println("\n\nERROR! \nYou must supply either \"start_tree\" or \"tree_file\" but NOT BOTH\n");
	    System.exit(-1);
	}


	GenerateDistributions();

       // 	if(sTree == null ) {
// 	    System.err.println("No starting tree.");
// 	    System.exit(-1);
// 	}
	System.err.println("Loaded command file: "+cname);
    }


    

    private void CheckBoundsGr0(double x, String vname) {
	if( x <= 0 ) {
	    System.err.println("Parameter "+vname+" must be > 0.");
	    System.exit(-1);
	}
    }
    
    private void CheckBoundsGrEqual0(double x, String vname) {
	if( x < 0 ) {
	    System.err.println("Parameter "+vname+" must be >= 0.");
	    System.exit(-1);
	}
    }
    
    private void CheckBounds01(double x, String vname) {
	if( (x < 0) || (x > 1) ) {
	    System.err.println("Parameter "+vname+" must be >= 0 and <= 1.");
	    System.exit(-1);
	}
    }

    private void CheckBoolean(int x, String vname){
	if( (x != 0) && (x != 1) ){
	    System.err.println("Parameter "+vname+" must be = 0  or = 1.");
	    System.exit(-1);
	}

    }

    public String [] ReadTreesFromFile(String file_name){

	BufferedReader input = null;
	try {
	    input = new BufferedReader( new FileReader(file_name) );
	} catch (FileNotFoundException e) {
	    System.err.println("Error opening tree file: " + file_name);
	    System.exit(-1);
	}
	
	String my_string = null;
	int tree_index = 0;


	/* count trees in a file */

	try{
	    while( (my_string = input.readLine()) != null){
		tree_index++;
	    }
	    
	} catch (Exception e) {
	    System.err.println("Unable to parse line" + my_string);
	    System.exit(-1); 
	}

	if (tree_index <= 2){
	    System.err.println("\n Error \n Tree list must have at least 3 trees");
	    System.exit(-1);
	}


	/* read in individual tree strings into array */
	
	String[] tree_strings = new String[tree_index];

	try {
	    input = new BufferedReader( new FileReader(file_name) );
	} catch (FileNotFoundException e) {
	    System.err.println("Error opening tree file: " + file_name);
	    System.exit(-1);
	}



	try{
	    int j = 0;
	    while( (my_string = input.readLine()) != null){
		tree_strings[j] = my_string;
		j++;
	    }
	    
	} catch (Exception e) {
	    System.err.println("Unable to parse line" + my_string);
	    System.exit(-1); 
	}


	System.out.println("\n" + tree_index + " trees supplied:");
	for (int i = 0; i < tree_strings.length; i ++){
	    System.out.println(tree_strings[i]);
	}
	System.out.println();

	return tree_strings;
    }

    public Tree[] LoadTrees(){
	Tree[] trees = new Tree[tree_list.length];

	for (int i = 0; i < tree_list.length; i++){
	    trees[i] = new iTree(tree_list[i]);
	}

	return trees;
    }
    
    
    public void ReadCmdfile(String fname) {
	BufferedReader br = null;
	try {
	    br = new BufferedReader( new FileReader(fname) );
	} catch (FileNotFoundException e) {
	    System.err.println("Error opening cmdfile: "+fname);
	    System.exit(-1);
	}
	String s = null;
	try {
	    while( (s = br.readLine()) != null ) {
		StringTokenizer st = new StringTokenizer(s);
		String cmd = st.nextToken();
		String option = st.nextToken();
		System.out.println(cmd);
		System.out.println(option);

		// number of MCMC iterations
		if( cmd.compareTo("length:") == 0 ) {
		    length = Integer.parseInt(option);
		    CheckBoundsGr0(length,"length");
		    
		    // length of the burnin period
		} else if( cmd.compareTo("burnin:") == 0 ) {
		    burnin = Integer.parseInt(option);
		    CheckBoundsGrEqual0(burnin,"burnin");
		    
		    // subsampling step
		} else if( cmd.compareTo("subsample:") == 0 ) {
		    subsample = Integer.parseInt(option);
		    CheckBoundsGr0(subsample,"subsample");

		    // max length of the jump during a change(break)-point update
		} else if( cmd.compareTo("window_length:") == 0 ) {
		    lenWindow = Integer.parseInt(option);
		    CheckBoundsGr0(lenWindow,"window_length");

		    
		    // mean of the prior Poisson distribution of number of topology break-points
		} else if( cmd.compareTo("top_lambda:") == 0) {
		    top_lambda = Double.valueOf(option).doubleValue();
		    CheckBoundsGr0(top_lambda,"top_lambda");

		    // mean of the prior Poisson distribution of number of substitution change-points
		} else if( cmd.compareTo("par_lambda:") == 0) {
		    par_lambda = Double.valueOf(option).doubleValue();
		    CheckBoundsGr0(par_lambda,"par_lambda");

		    // tunable constant for proposing new transition:transversion ratios
		    // when adding a new substitution change-point

		} else if( cmd.compareTo("sigma_alpha:") == 0) {
		    sigmaAlpha = Double.valueOf(option).doubleValue();
		    CheckBoundsGr0(sigmaAlpha,"sigma_alpha");

		    // tunable constant for proposing new average divergence
		    // when adding a new substitution change-point

		} else if( cmd.compareTo("sigma_mu:") == 0) {
		    sigmaMu    = Double.valueOf(option).doubleValue();
		    CheckBoundsGr0(sigmaMu,"sigma_mu");

		    // tunable constant for rescaling probabilities of selecting
		    // a particular move during rjMCMC

		} else if( cmd.compareTo("C:") == 0) {
		    C = Double.valueOf(option).doubleValue();
		    CheckBounds01(C,"C");
		    		    
		    // Mean of log(kappa) ~ Normal
		    
		} else if( cmd.compareTo("subst_hyper_mean:") == 0) {
		    subst_hyper_mean = Double.valueOf(option).doubleValue();
		    hyper1 = true;
		    
		    // Variance of log(kappa) ~ Normal
		    
		} else if( cmd.compareTo("subst_hyper_variance:") == 0) {
		    subst_hyper_variance = Double.valueOf(option).doubleValue();
		    CheckBoundsGr0(subst_hyper_variance,"subst_hyper_variance");
		    hyper2 = true;
		    
		    // Mean of log(mu) ~ Normal
		    
		} else if( cmd.compareTo("diver_hyper_mean:") == 0) {
		    diver_hyper_mean = Double.valueOf(option).doubleValue();
		    hyper3 = true;
		    
		    // Variance of log(mu) ~ Normal
		    
		} else if( cmd.compareTo("diver_hyper_variance:") == 0) {
		    diver_hyper_variance = Double.valueOf(option).doubleValue();
		    CheckBoundsGr0(subst_hyper_variance,"diver_hyper_variance");
		    hyper4 = true;
		    
		    // Starting tree topology
		    
		} else if( cmd.compareTo("start_tree:") == 0 ) {
		    sTree = new iTree(option);
		    start_tree_present = true;

		    // Name of the file with tree list

		} else if (cmd.compareTo("tree_file:") == 0){
		    tree_file_present = true;
		    tree_list = ReadTreesFromFile(option);

		    // Indicator of whether topology break-points are to be updated

		} else if( cmd.compareTo("top_breaks:") == 0 ) {
		    temp_int = Integer.parseInt(option);
		    CheckBoolean(temp_int,"top_breaks");
		    if (temp_int == 0) {
			top_breaks = false;
		    }
		    else {
			top_breaks = true;
		    }

		} else if( cmd.compareTo("par_breaks:") == 0 ) {
		    temp_int = Integer.parseInt(option);
		    CheckBoolean(temp_int,"par_breaks");
		    if (temp_int == 0){
			par_breaks = false;
		    }
		    else{
			par_breaks = true;
		    }
		    

		} else { 
		    PrintOptions();
		    System.err.println("\nUnable to parse line: "+s);
		    System.exit(-1);
		}
		

	    }
	} catch (Exception e) {
	    PrintOptions();
	    System.err.println("\nUnable to parse line huy: "+s);
	    System.exit(-1);
	}
	
    }


    public void PrintOptions(){
	System.err.println("\nList of all possible command file options");
	System.err.println("--------------------------------------------");

	System.err.println("length: integer  - length of the Markov chain");
	System.err.println("burnin: integer - number of initial iterations to be discarded");
	System.err.println("subsample: integer");
	System.err.println("window_length: integer");
	System.err.println("top_lambda: real > 0");
	System.err.println("par_lambda: real > 0");
	System.err.println("sigma_alpha: real > 0");
	System.err.println("sigma_mu: real > 0");
	System.err.println("C: real > 0");
	System.err.println("subst_hyper_mean: real ");
	System.err.println("subst_hyper_variance: real > 0");
	System.err.println("diver_hyper_mean: real");
	System.err.println("diver_hyper_variance: real > 0");
	System.err.println("start_tree: tree in Newick format");
	System.err.println("tree_file: file name with trees in Newick format");
	System.err.println("top_breaks: 1 = update, 0 = do not update (default: 1)");
	System.err.println("par_breaks: 1 = update, 0 = do not update (default: 1)");

	System.err.println("\nCAUTION: You must supply either \"start_tree\" or \"tree_file\" but NOT BOTH\n");

	System.err.println("-------------------------------------------");


    }
    
    public void GenerateDistributions() {
	norm01  = new Normal(0,1,rengine);
	gamma   = new Gamma(1,1,rengine);
	gammaV  = new Gamma(2.1,1.1,rengine);
	beta    = new Beta(1,1,rengine);
	chisquare = new ChiSquare(1,rengine);
	normal = new Normal(0,1,rengine);
	unif    = new Uniform(rengine);
	dispersed_unif = new Uniform(0,100,rengine);
    }

    public static void main(String[] args) {
	
	String cname = null;
	
	
	try {
	    cname = args[0];

	} catch (Exception e) {
	    System.err.println("Error parsing the command line.");
	    System.exit(-1);
	}
	
	Settings set = new Settings(100,cname);
    }
}
