import java.util.*;
import java.io.*;
import corejava.*;
//import lava.clib.stdio.*;

/**
* Constructs the site-by-site marginal topology estimates.
* Input comes from CPSampler
*/
public class TopologyProfile {

	int headerlen = 3;
	int segmentlen = 9;
	//PrintfFormatString pfs1 = new PrintfFormatString("%-6d");
	//PrintfFormatString pfs2 = new PrintfFormatString(" %5.3f");

	Format fmtInt  = new Format("%-6d");
	Format fmtDouble     = new Format(" %5.3f");

	public TopologyProfile(String fname, String oname, int seqlen) {
		BufferedReader br = null;
		ArrayList alltrees = new ArrayList();
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
		int lines = 0;
		try {
			while( br.ready() ) {
				StringTokenizer st = new StringTokenizer(br.readLine());
				int len = st.countTokens();
				int end = len;
				String[] f = new String[end];
				for(int i=0; i<end; i++) {
					f[i] = st.nextToken();
				}
				len -= headerlen;
				len /= segmentlen;
				for(int i=0; i<len; i++) {
					TreeCount probe = new TreeCount(f[headerlen+(i*(segmentlen-1))]);
					int loc = alltrees.indexOf(probe);
					if( loc == -1 ) {
						alltrees.add(probe);
					} else {
						((TreeCount)alltrees.get(loc)).count++;
					}
				}
				lines++;
			}
			br.close();
		} catch (IOException e) {
			System.err.println("Unable to parse line #"+(lines+1));
			System.exit(-1);
		}
		Collections.sort(alltrees);
		System.out.println("Read in "+lines+" lines.");
		System.out.println("Found "+alltrees.size()+" unique trees.");
		System.out.print("Tree list:");
		for(int i=0; i<alltrees.size(); i++) {
			System.out.print(" "+((TreeCount)alltrees.get(i)).tree );
		}
		System.out.println();
		int[][] post = new int[alltrees.size()][seqlen];
		// Reopen and reread post file
		try {
			br = new BufferedReader(new FileReader(fname));
		} catch (IOException e) {
			System.err.println("Unable to open file: "+fname);
			System.exit(-1);
		}
		try {
			for(int l=0; l<lines; l++) {
				StringTokenizer st = new StringTokenizer(br.readLine());
				int len = st.countTokens();
				int end = len;
				String[] f = new String[end];
				for(int i=0; i<end; i++) {
					f[i] = st.nextToken();
				}
				len -= headerlen;
				len /= segmentlen;
				int[] chpt = new int[len];
				for(int i=0; i<len-1; i++) {
					chpt[len-i-1] = Integer.parseInt(f[end-1-i]);
				}
				for(int i=0; i<len; i++) {
					TreeCount probe = new TreeCount(f[headerlen+(i*(segmentlen-1))]);
					int loc = alltrees.indexOf(probe);
					int start = chpt[i];
					int stop = seqlen;
					if( (i+1) < len )
						stop = chpt[i+1];
					for(int j=start; j<stop; j++)
						post[loc][j]++;
				}
			}
		} catch (IOException e) {
			System.err.println("Unable to parse line #"+(lines+1));
			System.exit(-1);
		}
		// Output results
		for(int i=0; i<seqlen; i++) {
		/*	StringBuffer sb = new StringBuffer(
				//pfs1.sprintf(i+1));
				Printf.format( pfs1, new Object[] {
				 	new Integer(i+1)
				} ));
		*/
			StringBuffer sb = new StringBuffer(fmtInt.form(i+1));
			for(int j=0; j<alltrees.size(); j++) {
			/*	sb.append(// pfs2.sprintf(
					Printf.format(pfs2, new Object[] {
						new Double( post[j][i]/(double)lines)
					} ) );*/
				sb.append(fmtDouble.form( (post[j][i]/(double)lines) ) );
			}
			output.println(sb.toString());
		}
		output.close();
	}

	public static void main(String[] arg) {
		try {
			new TopologyProfile(arg[0],arg[1],Integer.parseInt(arg[2]));
		} catch (Exception e) {
			System.err.println("Unable to parse commandline.");
			System.err.println("java TopologyProfile <posterior filename> <output filename> <sequence length>");
			System.exit(-1);
		}
	}
}
