import java.util.*;
import java.io.*;

import corejava.*;


/**
* Multiple change-point sampler extended to handle two types of change points 
* (topology break-points and evolutionary change-points)
*/

public class DCPSampler extends Sampler  implements Serializable {

    /**
     * List of current partitions/segments
     * @see Partition
     */
    protected Vector partition_list;


    /**
     * Mapping from site location along alignment to pattern location in counts/data.
     */
    protected int[] indexSeq;


    /**
     * Length of alignment.
     */
    protected int lenSeq;



    /**
     * Current number of topology change points
     */
    protected int topology_changes;

    /**
     * Current number of parameter change points
     */
    protected int parameter_changes;


    /**
     * Holds the observed data 
     * @see CountStatistic
     */
    protected CountStatistic mData[];


	
    /**
     * Number of change-points with the same topology on both sides.
     */
    protected int cSameT;

    
    /**
     * Current chain step.
     */
    protected int JumpNumber = 1;


    /**
     * Current total log likelihood.
     */
    protected double cLogLikelihood; // Holds the current loglikelihood


    /**
     * Acceptance counts [i] of MH proposal steps i.
     */
    protected double acceptancerate[];


    /**
     * Trial counts [i] of MH proposal steps i.
     */
    protected double tries[];


    /**
     * Names of proposal steps
     */

    protected String acc_names[] = {"acceptance probability of topology_update   ", 
				"acceptance probability of parameter_update  ", 
				"acceptance probability of top_xi_update     ", 
				"acceptance probability of par_xi_update     ", 
				"acceptance probability of add_one_top_xi    ",
				"acceptance probability of add_two_top_xi    ",
				"acceptance probability of add_par_xi        ", 
				"acceptance probability of remove_one_top_xi ",
				"acceptance probability of remove_two_top_xi ",
				"acceptance probability of remove_par_xi     ", 
				"acceptance probability of fixed_sampler_prob"
                                };

    protected String try_names[] = {"total number of tries of topology_update   ", 
				    "total number of tries of parameter_update  ", 
				    "total number of tries of top_xi_update     ", 
				    "total number of tries of par_xi_update     ", 
				    "total number of tries of add_one_top_xi    ",
				    "total number of tries of add_two_top_xi    ",
				    "total number of tries of add_par_xi        ", 
				    "total number of tries of remove_one_top_xi ",
				    "total number of tries of remove_two_top_xi ",
				    "total number of tries of remove_par_xi     ", 
				    "total number of tries of fixed_sampler_prob"
    };


    /**
     * Length of uniform window proposal for updating change-point locations.
     */
    protected int lenWindow;


    /**
     * Log of 1/2.
     */
    private double logWeight = -0.693;
    

    /**
     * One topology change point birth mixing probabilitiy given topology_changes.
     */
    double top_one_bk;


    /**
     * Two topology change points birth mixing probabilitiy given topology_changes.
     */
    double top_two_bk;

     /**
     * Birth mixing probabilitiy given parameter_changes.
     */
    double par_bk;

    
    /**
     * Death mixing probability given topology_changes.
     */
    double top_dk;


    /**
     * Death mixing probability given parameter_changes.
     */
    double par_dk;


    /**
     * One topology change point birth mixing probabilitiy given topology_changes - 1.
     */
    double top_one_bkm1;


    /**
     * Two topology change points birth mixing probabilitiy given topology_changes - 1.
     */
    double top_two_bkm2;


    /**
     * Birth mixing probability given parameter_changes - 1.
     */
    double par_bkm1;

    
    /**
     * Death mixing probability given topology_changes + 1.
     */
    double top_dkp1;


    /**
     * Death mixing probability given topology_changes + 1.
     */
    double top_dkp2;
    

    /**
     * Death mixing probability given parameter_changes + 1.
     */
    double par_dkp1;
    
    
    /**
     * Possible trees.
     */
    protected Tree[] PostTree = null;


    /**
     * Current sampler settings.
     */
    protected Settings set;


    /**
     * Results file.
     */
    protected PrintWriter output;

    /**
     * Priors
     */
    protected Priors prior_inf;

    /**
     * Indicates if sites with gaps should be included in likelihood (true) or ignored (false).
     */
    protected boolean handleGaps = true;

    

    /**
     * Hyperprior mean on parameter_changes
     */
    protected double par_lambda;


    /**
     * Hyperprior mean on topology_changes
     */
    protected double top_lambda;

    
    double logLength;
    double logParLambda;
    double logTopLambda;
    double logTwo;
    double logNTaxa;
    double logTreeSpan;
    double oneOverTreeSpan;

    /**
     * Length of possible trees in PostTree
     */
    int TreeSpan;
    double sigmaA;
    double sigmaM;
    double birthMixture;
    double C;

    public Tree[] ApproxPost;


   
    

    protected double[] ChooseStartingPi(CountStatistic cs) {
	double[] rtn = new double[4];
	int length = cs.sequences.length;
	for(int i=0; i<length; i++) {
	    rtn[0] += (cs.sequences[i]).Composition(0);
	    rtn[1] += (cs.sequences[i]).Composition(1);
	    rtn[2] += (cs.sequences[i]).Composition(2);
	    rtn[3] += (cs.sequences[i]).Composition(3);
	}
	rtn[0] /= length;
	rtn[1] /= length;
	rtn[2] /= length;
	rtn[3] /= length;
	return rtn;
    }
 

    protected double ChooseStartingAlpha() {
	//double rtn = 0.5;
	return 100*set.unif.nextDouble();
    }

    
    protected double[] ChooseStartingEP() {
	double[] rtn = new double[2];
	rtn[0] = ChooseStartingAlpha();
	rtn[1] = rtn[0];
	return rtn;
    }


    protected double ChooseStartingMu() {
	return set.unif.nextDouble(); 
    }

    
    protected int ChooseStartingTree(int i) {
	int new_tree = 0;

	if ( i == 0 )
	    new_tree = (int)(set.unif.nextDouble() * PostTree.length);
	else{
	    boolean good_tree = false;
	    
	    Partition prev_part = (Partition)partition_list.elementAt(i-1);
	    
	    while(!good_tree) {
		new_tree = (int)(set.unif.nextDouble()*PostTree.length);

		if ( new_tree != prev_part.cTree)
		    good_tree = true;
	    }
	}
	
	return new_tree;
    }

    
    protected QMatrix ChooseStartingQ(int i) {
	return new iHKYNoBoundFixPiMatrix(ChooseStartingEP(), ChooseStartingPi(mData[i]));
    }
    
    
    /** Ininializes parameters of the model
     *  @param part_list list of labeled change-points
     */
    

    protected void InitialState(Vector part_list) {
      	Partition curr_segment, last_segment = new Partition(0,0);
	cLogLikelihood = 0;

	

	int j = 0;
	for(int i = 0; i < partition_list.size(); i++){
	    curr_segment = (Partition)partition_list.elementAt(i);

	    
	    if (curr_segment.IsTopologyChange()){
		curr_segment.cTree = ChooseStartingTree(i);
		j++;
	    }
	    else
		curr_segment.cTree = last_segment.cTree;

 
	    if (curr_segment.IsParameterChange()){
		curr_segment.cHyperParameter = ChooseStartingMu();
		curr_segment.cMatrix = ChooseStartingQ(0);
	    }else {
		curr_segment.cHyperParameter = last_segment.cHyperParameter;
		curr_segment.cMatrix = last_segment.cMatrix;
	    }

  
	    curr_segment.cPartialLogHyperParameterPrior = 
		prior_inf.Prior_on_mu(curr_segment.cHyperParameter);

	    //Calculate initial likelihoods
	    curr_segment.cPartialLogLikelihood = PostTree[curr_segment.cTree].LogLikelihood0(curr_segment.cMatrix, curr_segment.cHyperParameter,curr_segment.counts, curr_segment.data, curr_segment.cMatrix.pi, handleGaps,0,0);

	    cLogLikelihood += curr_segment.cPartialLogLikelihood;
	    last_segment = curr_segment;
	}
		
	System.err.println("Starting Likelihood: "+cLogLikelihood);
	//cSameT = determineSameT(partition_list);
		
    }
    

    protected void SetParameters() {
	// Move the remainder to an init routine later
	lenWindow = set.lenWindow;
	par_lambda = set.par_lambda;
	top_lambda = set.top_lambda;
	
	logLength = Math.log(lenSeq);
	logParLambda = Math.log(par_lambda);
	logTopLambda = Math.log(top_lambda);
	logTwo = Math.log(2.0);
	Partition temp_segment = (Partition)partition_list.firstElement();
	logNTaxa = Math.log( PostTree[temp_segment.cTree].numLeaves );
	TreeSpan = PostTree.length;
	logTreeSpan = Math.log(TreeSpan);
	oneOverTreeSpan = 1.0 / (double) TreeSpan;
	System.err.println("1/Span = "+oneOverTreeSpan);

	sigmaA = set.sigmaAlpha;
	sigmaM = set.sigmaMu;

	C = set.C;
    }


    /**
     * Forms tree space either from starting tree by removing the last taxon 
     * and reattaching it to all possible branches of "parental tree" or 
     * by loading tree list, supplied by user.
     */

    protected void SetPostTree() {
	if (set.start_tree_present){
	    PostTree = set.sTree.EnumerateLastTaxon();
	}
	else{
	    try{
		PostTree = set.LoadTrees();
	    } catch (Exception e) {
		System.err.println("\n Error! \n Something is wrong with the tree list");
		System.exit(-1);
	    }
	}
    }


    /** Setting up the model for the 
     *  first round of sampling
     *  @param ics contains the data and data statistics
     *  @param iset contains settings from the command file
     *  @param ioutput name of the output file
     */


    public DCPSampler(CountStatistic ics, Settings iset, PrintWriter ioutput, Priors ipriors) {
	set = iset;
	output = ioutput;
	prior_inf = ipriors;

	/****** Set up arrays to hold acceptance statistics *******/
	
	acceptancerate 	= new double[acc_names.length];
	tries      	= new double[acc_names.length];
       
	lenSeq = ics.lenseq;
	partition_list = new Vector();
	Partition segment = new Partition(0,lenSeq-1);
	parameter_changes = 1;
	topology_changes = 1;

	//Straighted out data partitions
	mData = new CountStatistic[1];
	mData[0] = ics;
	segment.counts = ics.countsMatrix();
	segment.data   = ics.dataMatrix();

	partition_list.addElement(segment);
      
	indexSeq = ics.siteMap();
	int nPatterns = ics.lenunique;
	
	// Add change-points to the initial state if needed
	
// 	AddChangePoint(80, false, true);
// 	AddChangePoint(160, false, true);
// 	AddChangePoint(240, false, true);
// 	AddChangePoint(320, false, true);
// 	AddChangePoint(400, false, true);
// 	AddChangePoint(480, true, false);
// 	AddChangePoint(560, true, false);
// 	AddChangePoint(640, true, false);
// 	AddChangePoint(720, false, true);
// 	AddChangePoint(800, false, true);
	
	//AddChangePoint(300, true, false);


	SetPostTree();
	InitialState(partition_list);
	SetParameters();
	
	JumpNumber = 1;
	if( JumpNumber > set.burnin)
	    saveEstimates();
    }

    
    public DCPSampler() {
	set = null;
    }

    protected double DrawnAlpha1(double alpha, double weight, double sigma, double z) {
	return alpha * Math.exp(weight * sigma * z);
    }


    protected double DrawnAlpha2(double alpha, double weight, double sigma, double z) {
	return alpha * Math.exp(- weight * sigma * z);
    }


    protected double DrawnMu1(double mu, double weight, double sigma, double z) {
	return mu * Math.exp(weight * sigma * z);
    }


    protected double DrawnMu2(double mu, double weight, double sigma, double z) {
	return mu * Math.exp(- weight * sigma * z);
    }


    protected double CondensedAlpha(double alpha1, double weight1, double alpha2, double weight2 ) {
	return Math.exp( weight1 * Math.log(alpha1) + weight2 * Math.log(alpha2) );
    }


    protected double CondensedMu(double mu1, double weight1, double mu2, double weight2) {
	return Math.exp( weight1 * Math.log(mu1) + weight2 * Math.log(mu2) );
    }


    protected double InverseZMu(double mu1, double mu0, double weight1, double sigma) {
	return Math.log( mu1 / mu0 ) / weight1 / sigma;
    }

	
    protected double InverseZAlpha( double alpha1, double alpha0, double weight1, double sigma) {
	return Math.log(alpha1 / alpha0 ) / weight1 / sigma;
    }


    protected double logJacobian(double alpha, double alpha1, double alpha2, double sigmaa,
			       double mu, double mu1, double mu2, double sigmam) {
	return Math.log( alpha1 )  + Math.log( alpha2 )
	    + Math.log( mu1 ) + Math.log( mu2 )
	    + Math.log( sigmaa ) + Math.log( sigmam )
	    - Math.log( alpha ) - Math.log( mu );
    }
    


    String delimiter = new String("\t");
    static Format fmtIter  = new Format("%-6d");
    Format fmtN     = new Format(" %2d");
    Format fmtLogD  = new Format(" %10.2f");
    Format fmtProp  = new Format(" %.4f");
    Format fmtD     = new Format(" %7.4f");
    Format fmtDBig  = new Format(" %8.3f");
    static Format fmtDE	   = new Format(" %12.4E");
    
    
    protected StringBuffer OutputLine() {
	StringBuffer sb = new StringBuffer();
	sb.append(fmtIter.form(JumpNumber));
	sb.append(fmtN.form(partition_list.size()));
	sb.append(fmtN.form(cSameT));
	Partition temp_part;

	for(int i = 0; i < partition_list.size(); i++) {
	    temp_part = (Partition)partition_list.elementAt(i);

	    sb.append(" ");
	    sb.append(PostTree[temp_part.cTree].Root.toString());
	    sb.append(fmtLogD.form(temp_part.cPartialLogLikelihood));
	    sb.append(fmtProp.form(temp_part.cMatrix.v[0]));
	    sb.append(fmtProp.form(temp_part.cMatrix.pi[0]));
	    sb.append(fmtProp.form(temp_part.cMatrix.pi[1]));
	    sb.append(fmtProp.form(temp_part.cMatrix.pi[2]));
	    sb.append(fmtProp.form(temp_part.cMatrix.pi[3]));
	    sb.append(fmtD.form(temp_part.cHyperParameter));
	}

	for(int i = 0; i < partition_list.size(); i++) {
	    temp_part = (Partition)partition_list.elementAt(i);
	    sb.append(" ");
	    sb.append(temp_part.left.nuc_position);
	}
	
	return sb;
    }



    
    /**
     * Saves estimates of a MCMC update step
     * also checks if the number of steps exceeds 
     * user defined chain length and exits the program
     * when it happens. At the end outputs date and
     * acceptance rates of M-H steps.
     **/

    public final void saveEstimates() {
	output.println(OutputLine());
	output.flush(); //May speed up output to remove this line.
	// Test to see if we are finished with this chain
	if( JumpNumber >= set.length ) {
	    Date now = new Date();
	    System.out.println("# Finished @ "+now);
	    for(int i = 0; i < acceptancerate.length; i++ ) {
		if( tries[i] > 0 ){
		    System.out.println(try_names[i] + " = " + tries[i]);
		    System.out.println(acc_names[i] + " = "+((double) acceptancerate[i] / (double)(tries[i]) ));
		}
		else
		    System.out.println(acc_names[i] + " not used");
	    }
	    output.close();
	    System.exit(0);
	}
    }
    


    TopUpdate topology = new TopUpdate(this);
    ParUpdate parameter = new ParUpdate(this);
    ParChangeUpdate parameter_change_points = new ParChangeUpdate(this);
    ParDimChange parameter_dim_change = new ParDimChange(this);
    Hierarchy hierarchical_parameters = new Hierarchy(this);
    
    TopChangeUpdate topology_change_points;
    TopDimChange topology_dim_change;
    

    /** Runs the sampler
     *
     */
    
    public void run() {
	int sincePrint = 0;
	double alpha;
	if (TreeSpan ==2){
	    alpha = 0.5;
	}else{
	    alpha = 1.0/(double)(TreeSpan - 1);
	}

 	if(prior_inf.inf_recomb){
 	    topology_change_points = new InfTopChange(this);
 	    topology_dim_change = new InfTopDim(this);
	}else{
 	    topology_change_points = new TopChangeUpdate(this);
	    topology_dim_change = new TopDimChange(this);
 	}
	

 	while( true ) {

	    top_one_bk = (1-alpha)*(double)top_lambda / (double)(topology_changes);
	    if(top_one_bk > 1.0 ) top_one_bk = 1.0;
	    
	    top_two_bk = alpha*(double)top_lambda * (double)top_lambda/ 
		((double)(topology_changes)*(double)(topology_changes + 1));
	    if(top_two_bk > 1.0 ) top_two_bk = 1.0;

	    top_dk = (double)(topology_changes - 1) / (double)top_lambda;
	    if( top_dk > 1.0 ) top_dk = 1.0;

	    par_bk = (double)par_lambda / (double)(parameter_changes);
	    if(par_bk > 1.0 ) par_bk = 1.0;

	    par_dk = (double)(parameter_changes - 1) / (double)par_lambda;
	    if( par_dk > 1.0 ) par_dk = 1.0;

	    top_one_bkm1 = (1-alpha)*(double)top_lambda / (double)(topology_changes - 1);
	    if( top_one_bkm1 > 1.0 ) top_one_bkm1 = 1.0;

	    top_two_bkm2 = alpha*(double)top_lambda * (double)top_lambda/ 
		((double)(topology_changes - 2)*(double)(topology_changes - 1));
	    if(top_two_bkm2 > 1.0 ) top_two_bkm2 = 1.0;

	    top_dkp1 = (double)(topology_changes) / (double)top_lambda;
	    if( top_dkp1 > 1.0 ) top_dkp1 = 1.0;

	    top_dkp2 = (double)(topology_changes + 1) / (double)top_lambda;
	    if( top_dkp2 > 1.0 ) top_dkp2 = 1.0;
	    

	    if (set.top_breaks){
		top_one_bk *= C;
		top_two_bk *= C;
		top_dk *= C;
	    }else{
		top_one_bk = 0;
		top_two_bk = 0;
		top_dk = 0;
	    }

	    if (set.top_breaks){
		par_bk *= C;
		par_dk *= C;
	    }else{
		par_bk = 0;
		par_dk = 0;
	    }
		
	    top_one_bkm1 *= C;
	    top_two_bkm2 *= C;
	    top_dkp1 *=C;
	    
	    double BirthOrDeathOrMove = 0;


	    BirthOrDeathOrMove = set.unif.nextDouble();
		
	    tries[10]++;

	    if (BirthOrDeathOrMove < top_one_bk){

		if  (topology_changes < lenSeq)
		    topology_dim_change.AddOne();
		


	    }else if (BirthOrDeathOrMove < top_one_bk + top_two_bk){

		if (topology_changes < lenSeq - 1)
		    topology_dim_change.AddTwo();

		    
	    }else if (BirthOrDeathOrMove < top_one_bk + top_two_bk + top_dk){

		if (topology_changes > 1)
		    topology_dim_change.Delete();


	    }else if (BirthOrDeathOrMove < top_one_bk + top_two_bk + top_dk + par_bk){

		if (parameter_changes < lenSeq)
		    parameter_dim_change.Add();


	    }else if (BirthOrDeathOrMove < top_one_bk + top_two_bk + top_dk + par_bk + par_dk){
		
		if (parameter_changes > 1)
		    parameter_dim_change.Delete();



	    }else {
		acceptancerate[10]++;

		parameter.Update();
 		topology.Update();
		if (set.hyper_update){
		    hierarchical_parameters.Update();
		}
		
		if(parameter_changes > 1)
 		    parameter_change_points.Update();
		
		
 		if (topology_changes > 1)

 		    topology_change_points.Update();

	    }


	    JumpNumber++;

	    
	    sincePrint++;
	    if( (sincePrint >= set.subsample) ) {
		sincePrint = 0;
		System.out.println(JumpNumber);// + " " + prior_inf.alpha_mean + " " + prior_inf.alpha_variance
		//+ " " + prior_inf.mu_mean + " " + prior_inf.mu_variance);
		if( JumpNumber > set.burnin ){
		    saveEstimates();
		   
		}
 	    }
	    
 	}

    }



    /** debugging routine that checks likelihood calculations
     * @param label this can be used in order to find out after what step
     * an inconsistency in likelihood occured
     */

    void Debug(String label){
	
	for (int j = 0; j < partition_list.size(); j ++){
	    
	    Partition my_part = (Partition)partition_list.elementAt(j);
	    int[] my_counts = new int[my_part.counts.length];
	    
	    for(int k = my_part.left.nuc_position; k <= my_part.right; k++)
		my_counts[indexSeq[k]]++;
	
	    // double[] EP = new double [2];
// 	    EP[0] = 1000;
// 	    EP[1] = 1000;
// 	    QMatrix pMatrix = new iHKYNoBoundFixPiMatrix(EP, my_part.cMatrix.pi);

	    double my_like = PostTree[my_part.cTree].LogLikelihood0(my_part.cMatrix, my_part.cHyperParameter, my_counts, my_part.data, my_part.cMatrix.pi, handleGaps, 0, 0);
	    
	    if (Math.abs(my_like - my_part.cPartialLogLikelihood) > 0.000001){
		System.out.println(label + " " + JumpNumber);
		System.exit(1);
	    }

	    // System.out.println("likelihood: " + my_like + " alpha: " + my_part.cMatrix.v[0]);
	}
	
    }
    

    /**
     * Returns the number of the proposed tree. The tree is drawn uniformly 
     * from the all possible trees so, that it is not equal to the current tree,
     * left tree or right tree. If there are only 3 possible trees and
     * left and right trees are not equal, returning value is -1, indicating
     * that current value can not be updated.
     * @param curr_tree  # of the current tree 
     * @param left_tree  # of the tree at the left constant topology segment
     * @param right_tree  # of the tree at the right constatnt topology segment
     */
    

    protected int ProposeMiddleTree(int right_tree, int curr_tree, int left_tree){

	boolean good_tree = false;
	int new_tree = 0;
	

	if ((left_tree != right_tree) && (TreeSpan == 3))
	    return -1;

	while( !good_tree ) {

	    new_tree = (int)(set.unif.nextDouble()*TreeSpan);
	    
	    if((new_tree != right_tree) && (new_tree != left_tree) && (new_tree != curr_tree))
		good_tree = true;
	    
	}
	
	return new_tree;
    }
    
    /**
     * Returns the number of the proposed tree. The tree is drawn uniformly 
     * from the all possible trees so, that it is not equal to the current tree
     * and to the tree next to it.
     * @param curr_tree - # of the current tree 
     * @param next_tree - # of the tree next to the current (left or right)
     */

    
    protected int ProposeEndTree(int next_tree, int update_tree){
	boolean good_tree = false;
	int new_tree = 0;
	

	while( !good_tree ) {
	    new_tree = (int)(set.unif.nextDouble()*TreeSpan);

	    if((new_tree != next_tree) && (new_tree != update_tree))
		good_tree = true;
	    
	}
	
	return new_tree;
    }


    /**
     * Returns the number of the proposed tree. The tree is drawn uniformly 
     * from the all possible trees so, that it is not equal to the current tree
     * and to the tree next to it.
     * @param curr_tree - # of the current tree 
     */


    protected int ProposeOneTree(int update_tree){
	boolean good_tree = false;
	int new_tree = 0;

	while( !good_tree ) {

	    new_tree = (int)(set.unif.nextDouble()*TreeSpan);

	    if(new_tree != update_tree)
		good_tree = true;
	    
	}
	
	return new_tree;
    }
    



    /**
     * Generates a new location for a change-point randomly
     */

    protected int ProposeNewXi(int curr, int lower_bound, int upper_bound){
	
	int prop = curr;
	
    	while (prop == curr) {
	    
	    int prnd = (int)(set.unif.nextDouble()*lenWindow) + 1;
	    
	    if (set.unif.nextDouble() < 0.5)
		prnd *= -1;
	    
	    prop = curr + prnd;
	    
	    while ((prop < lower_bound) || (prop > upper_bound)){
		if (prop < lower_bound)
		    prop = 2*lower_bound - prop;
		else
		    prop = 2*upper_bound - prop;
		
	    }
	    
	}

	return prop;
    }

    
    /** 
     * Inserts a new partition to the left of a given partition
     */

    
    protected void InsertPartitionRight(Partition part, int new_point, int land_index, boolean top, int[] counts){
	
	Partition new_part = new Partition(new_point, part.right, top, !top);
	new_part.counts = new int[counts.length];

	part.right = new_point - 1;
	
	new_part.data = part.data;

	for (int i = 0; i < counts.length; i++){
		new_part.counts[i] = part.counts[i] - counts[i];
		part.counts[i] = counts[i];
	}
	
	partition_list.insertElementAt(new_part, land_index + 1);
	
    }

    
    
    /** 
     * Inserts a new partition to the right of a given partition
     */

    
    protected void InsertPartitionLeft(Partition part, int new_point, int land_index, boolean top, int[] counts){
	
	Partition new_part = new Partition(new_point, part.right, top, !top);
	new_part.counts = new int[counts.length];

	part.right = new_point - 1;
	
	new_part.data = part.data;
	
	for (int i = 0; i < counts.length; i++){
		new_part.counts[i] = counts[i];
		part.counts[i] -= counts[i];
	}

	partition_list.insertElementAt(new_part, land_index + 1);
	
    }
    


    /** Method that allows mannually add change-points to change the starting 
     *  state of the parameter space. Frequently used during debugging to 
     *  add a change point and block the steps that add new ones and possibly change their location
     *  @param new_point location of the new change-point
     *  @param t true = topology break-point
     *  @param p true = evolutionary change-point
     */ 
    
    public void AddChangePoint(int new_point, boolean t, boolean p) {
	
	int landing = 0;
	landing = LandingPart(new_point, 0, partition_list.size() - 1);

	Partition part_split = (Partition)partition_list.elementAt(landing);
	
	if (part_split.left.nuc_position == new_point){

	    if((part_split.IsTopologyChange() && t) || (part_split.IsParameterChange() && p)){
		System.out.println("Error");
		System.exit(1);
	    }
	    else{
		if (t)
		    part_split.left.top_change = true;
		if (p)
		    part_split.left.par_change = true;
	    }
	}
	else{
	    
	    Partition new_partition;
	    
	    new_partition = new Partition(new_point, part_split.right, t, p);
	    
	    part_split.right = new_point - 1;
	    
	    // Update counts for the new partition

	    new_partition.counts = new int[part_split.counts.length];

	    for (int i = new_point; i <= new_partition.right; i++)
		new_partition.counts[indexSeq[i]]++;

	    for (int i = 0; i < part_split.counts.length; i++)
		part_split.counts[i] -= new_partition.counts[i];
	    
	    // Update data for the new partition
	    
	    new_partition.data = part_split.data;

	    // Finaly insert new partition into the list
	    
	    partition_list.insertElementAt(new_partition, landing + 1);

	}

	if(t)
	    topology_changes++;
	if(p)
	    parameter_changes++;
	
    }

 
  
    /** Returns a partition to which a new point should belong
     *  @param new_point location of a new change-point
     *  @param start the first partition the search will be started from
     *  @param end the last parition search will check
     */

    public int LandingPart(int new_point, int start, int end){
	
	int i = start;
	
	Partition temp = (Partition)partition_list.elementAt(start);
	
	while (!((new_point >= temp.left.nuc_position) && (new_point <= temp.right)) && i <= end){
	    i++;
	    if (i <= end)
		temp = (Partition)partition_list.elementAt(i);
	    else
		System.err.println("Something is wrong");
	}
	


	return i;
    }

    
    /** Print partitions together with parameters associated with them
     *  and likelihoods and counts usefull for debugging
     */

    public void print_partitions(){
	for(int i = 0; i < partition_list.size(); i++){
	    Partition temp = (Partition)partition_list.elementAt(i);
	    System.out.print("(" + temp.left.nuc_position);
	    System.out.print(", " + temp.right+")");
	    if (temp.IsTopologyChange())
		if(temp.IsParameterChange())
		    System.out.print("- t+p");
		else
		    System.out.print("- t");
	    else
		if (temp.IsParameterChange())
		    System.out.print("- p");
		else
		    System.out.print("something wrong");

	    System.out.println(" tree: " + temp.cTree);
	    //System.out.print(" Mu: " + temp.cHyperParameter);
	    //System.out.print(" Alpha: " + temp.cMatrix.v[0]);
	    //System.out.print(" likel: " + temp.cPartialLogLikelihood);
	    //System.out.print(" counts: ");
	    //for(int j = 0; j < 20; j++)
	    //System.out.print(temp.counts[j] + "  ");
	    
 	    System.out.println("partition number: " + i);
	    
	    
	}
	
    }
    
   
} 
