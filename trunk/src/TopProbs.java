import java.util.*;
import java.io.*;
import corejava.*;
//import lava.clib.stdio.*;

/**
* Constructs the site-by-site marginal topology estimates.
* Input comes from DCPSampler
*/



public class TopProbs {

    public Vector tree_names;
    public Vector tree_counts;


    //PrintfFormatString pfs2 = new PrintfFormatString(" %5.3f");
    
    Format fmtInt  = new Format("%-6d");
    Format fmtDouble     = new Format(" %5.3f");
    
    public TopProbs(String fname, String oname, int seqlen) {
	BufferedReader br = null;
	
	tree_counts = new Vector();
	tree_names = new Vector();
	String temp_token;
	int num_changes;
	
	//	int temp[] = new int[10000];
	//temp[4] = 23;
	//tree_counts.addElement(temp);
	//int temp2[] = (int[])tree_counts.elementAt(0);

	
	
	try {
	    br = new BufferedReader(new FileReader(fname));
	} catch (IOException e) {
	    System.err.println("Unable to open file: "+fname);
	    System.exit(-1);
	}
	PrintWriter output = null;
	try {
	    output = new PrintWriter(new BufferedWriter(new FileWriter(oname)));
	} catch (IOException e) {
	    System.err.println("Unable to open output file: "+oname);
	    System.exit(-1);
	}
	


	int line = 0;


	try {
	    while( br.ready() ) {
		StringTokenizer st = new StringTokenizer(br.readLine());
		temp_token = st.nextToken(); //skip number of iterations
		num_changes = Integer.parseInt(st.nextToken());
		temp_token = st.nextToken(); //skip zero

		System.out.println("num_changes: " + num_changes);

		String[] trees = new String[num_changes];
		int[] tree_numbers = new int[num_changes];
		int cur_index;

		for (int i = 0; i < num_changes; i++){
		    trees[i] = st.nextToken(); // read in a tree
		    cur_index = TreeIndex(trees[i]);
		    if (cur_index == -1){
			//tree_names.addElement(trees[i]);
			int[] temp = new int[seqlen];
			tree_counts.addElement(temp);			    
		    }else{

		    }
		    temp_token = st.nextToken(); //skip likelihood
		    temp_token = st.nextToken(); //skip \kappa
		    temp_token = st.nextToken(); //skip frequencies
		    temp_token = st.nextToken(); //skip frequencies
		    temp_token = st.nextToken(); //skip frequencies
		    temp_token = st.nextToken(); //skip frequencies
		    temp_token = st.nextToken(); //skip \mu		     
		}
	     // 				for(int i=0; i<end; i++) {
// 					f[i] = st.nextToken();
// 				}
// 				len -= headerlen;
// 				len /= segmentlen;
// 				for(int i=0; i<len; i++) {
// 					TreeCount probe = new TreeCount(f[headerlen+(i*(segmentlen-1))]);
// 					int loc = alltrees.indexOf(probe);
// 					if( loc == -1 ) {
// 						alltrees.add(probe);
// 					} else {
// 						((TreeCount)alltrees.get(loc)).count++;
// 					}
// 				}
// 				lines++;
// 			}
// 			br.close();
// 		} catch (IOException e) {
// 			System.err.println("Unable to parse line #"+(lines+1));
// 			System.exit(-1);
// 		}


		line++;
		output.close();
	    }
	} catch (IOException e) {
	    System.err.println("Unable to parse line #"+(line+1));
	    System.exit(-1);
	}



    }


    public int TreeIndex(String tree){
	boolean in_the_list = false;
	int index = -1;
	int i = 0;

	while (!in_the_list && i < tree_names.size()){ 
	    if (tree.equals((String)tree_names.elementAt(i))){
		index = i;
		in_the_list = true;
	    }
	    
	    i++;
	}

	if (in_the_list){
	    return index;
	}else{
	    return -1;
	}

	    
    }
    

    public static void main(String[] arg) {
	
	
	try {
	    TopProbs tp = new TopProbs(arg[0],arg[1],Integer.parseInt(arg[2]));

	    int myindex = tp.TreeIndex("fdfd");

	    System.out.println("tree list size " + tp.tree_names.size());
	} catch (Exception e) {
	    System.err.println("Unable to parse commandline.");
	    System.err.println("java TopProbs <posterior filename> <output filename> <sequence length>");
	    System.exit(-1);
	}
	

	

    }
}
