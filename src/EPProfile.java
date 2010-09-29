import java.io.*;
import java.util.*;
import corejava.*;
//import lava.clib.stdio.*;

/**
* Constructs the site-by-site marginal evolutionary parameter estimates.
* Input comes from CPSampler
*/
public class EPProfile {

	int stepsize = 100;
	int segmentlen = 9;
	int headlen = 3;
	int firstEP = 5;
	int firstMu = 10;
	int firstTree = 3;

//	static PrintfFormatString pfs = new PrintfFormatString(
//		"%-6d %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f" );

	Format fmtInt  = new Format("%-6d");
	Format fmtDouble     = new Format(" %6.4f");


	EPProfile(String fname, String oname, int seqlen) {

		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(fname));
		} catch (IOException e) {
			System.err.println("Unable to open input file: "+fname);
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
		ArrayList datalist = new ArrayList();
		try {
			while( br.ready() ) {
				datalist.add(br.readLine());
				lines++;
			}
			br.close();
		} catch (IOException e) {
			System.err.println("Unable to read line #"+(lines+1));
			System.exit(-1);
		}

	//	System.out.println("Read "+lines+" lines.");
	//	System.exit(-1);

		int start = 0;
		int stop  = stepsize - 1;
		if( stop >= seqlen )
			stop = seqlen - 1;

		while( start < seqlen ) {
			double[] COPs = new double[stop-start+1];
			double[] eps  = new double[stop-start+1];
			double[] mus  = new double[stop-start+1];
			//double[] rs  = new double[stop-start+1];
			double[][] dep = new double[stop-start+1][lines];
			double[][] dmu = new double[stop-start+1][lines];
			//double[][] dr  = new double[stop-start+1][lines];
			int median = (int) (0.500*lines);
			int q1     = (int) (0.025*lines);
			int q2     = (int) (0.975*lines);
			int line = 0;
			for(ListIterator li = datalist.listIterator(0) ; li.hasNext() ; line++) {
				StringTokenizer st = new StringTokenizer((String)li.next());
				int len = st.countTokens();
				int end = len;
				String[] f = new String[end];
				for(int i=0; i<end; i++) {
					f[i] = st.nextToken();
				}
				len -= headlen;
				len /= segmentlen;
				//System.err.println(len);
				int[] chpt = new int[len];
				for(int i=0; i<len-1; i++) {
					chpt[len-i-1] = Integer.parseInt(f[end-i-1]);
				}
				chpt[0] = 0;
				for(int i=0; i<len; i++) {
					double ep = (new Double(f[firstEP +(i*(segmentlen-1))])).doubleValue();
					double mu = (new Double(f[firstMu+(i*(segmentlen-1))])).doubleValue();
					//double r  =  2.0 * ep / (1.0 - ep);
					int end2 = 0;
					if( (i+1) == len )
						end2 = seqlen;
					else
						end2 = chpt[i+1];
					for(int j=chpt[i]; j<end2; j++) {
						if( (j >= start) && (j <= stop) ) {
							dep[j-start][line] = ep;
							dmu[j-start][line] = mu;
							//dr[j-start][line]  = r;
							eps[j-start] += ep;
							mus[j-start] += mu;
							//rs[j-start]  += r;
						}
					}
					if( (i>0) && (chpt[i] >= start) && (chpt[i] <= stop) ) {
						//System.err.print("Testing "+chpt[i]);
						String t1 = f[firstTree+((i-1)*(segmentlen-1))];
						String t2 = f[firstTree+((i  )*(segmentlen-1))];
						//System.err.print(" "+t1+" "+t2);
						if( t1.compareTo(t2) != 0 )
							//System.err.println(" delta");
							COPs[chpt[i]-start]++;
						//} else {
						//	System.err.println(" same");
						//}
					}
				}
				f = null;
			}

			for(int i=start; i<=stop; i++) {
				eps[i-start]  /= lines;
				mus[i-start]  /= lines;
				COPs[i-start] /= lines;
				//rs[i-start] /= lines;
				try {
					sort(dep[i-start]);
					sort(dmu[i-start]);
					//sort(dr[i-start]);
				} catch (Exception e) {
					System.err.println("Sorting error.");
					System.exit(-1);
				}
				StringBuffer sb = new StringBuffer();
			/*	sb.append( //pfs.sprintf( new Object[] {
			 	Printf.format( pfs, new Object[] {
					new Integer(i+1),
					new Double(eps[i-start]),
					new Double(mus[i-start]),
					new Double(dep[i-start][q1]),
					new Double(dep[i-start][median]),
					new Double(dep[i-start][q2]),
					new Double(dmu[i-start][q1]),
					new Double(dmu[i-start][median]),
					new Double(dmu[i-start][q2]),
					new Double(COPs[i-start])
					//new Double(dr[i-start][q1]),
					//new Double(dr[i-start][median]),
					//new Double(dr[i-start][q2])
				} ) );
			*/
				sb.append(fmtInt.form(i+1));
				sb.append(fmtDouble.form(eps[i-start]));
				sb.append(fmtDouble.form(mus[i-start]));
				sb.append(fmtDouble.form(dep[i-start][q1]));
				sb.append(fmtDouble.form(dep[i-start][median]));
				sb.append(fmtDouble.form(dep[i-start][q2]));
				sb.append(fmtDouble.form(dmu[i-start][q1]));
				sb.append(fmtDouble.form(dmu[i-start][median]));
				sb.append(fmtDouble.form(dmu[i-start][q2]));
				sb.append(fmtDouble.form(COPs[i-start]));
				output.println(sb);
			}
			start = stop + 1;
			stop  += stepsize;
			if( stop >= seqlen )
				stop = seqlen - 1;
		}
	//	try {
			output.close();
	//	} catch (IOException e) {
	//		System.err.println("Error closing "+oname);
	//		System.exit(-1);
	//	}
	}


	void QuickSort(double a[], int lo0, int hi0) throws Exception {
	  int lo = lo0;
	  int hi = hi0;
	  double mid;

	  // pause for redraw
	  //pause(lo, hi);
	  if ( hi0 > lo0)
	  {

             /* Arbitrarily establishing partition element as the midpoint of
              * the array.
              */
             mid = a[ ( lo0 + hi0 ) / 2 ];

             // loop through the array until indices cross
             while( lo <= hi )
             {
        	/* find the first element that is greater than or equal to
        	 * the partition element starting from the left Index.
        	 */
        	while( ( lo < hi0 ) && ( a[lo] < mid ) )
        	   ++lo;

        	/* find an element that is smaller than or equal to
        	 * the partition element starting from the right Index.
        	 */
        	while( ( hi > lo0 ) && ( a[hi] > mid ) )
        	   --hi;

        	// if the indexes have not crossed, swap
        	if( lo <= hi )
        	{
        	   swap(a, lo, hi);
        	   // pause
        	   //pause();

        	   ++lo;
        	   --hi;
        	}
             }

             /* If the right index has not reached the left side of array
              * must now sort the left partition.
              */
             if( lo0 < hi )
        	QuickSort( a, lo0, hi );

             /* If the left index has not reached the right side of array
              * must now sort the right partition.
              */
             if( lo < hi0 )
        	QuickSort( a, lo, hi0 );

	  }
       }

       private void swap(double a[], int i, int j)
       {
	  double T;
	  T = a[i];
	  a[i] = a[j];
	  a[j] = T;

       }

       public void sort(double a[]) throws Exception
       {
	  QuickSort(a, 0, a.length - 1);
       }


	public static void main(String arg[]) {
		int slen = 0;
		try {
			slen = Integer.parseInt(arg[2]);
		} catch (Exception e) {
			System.err.println("try:");
			System.err.println("java EPProfile <posterior filename> <output filename> <sequence len>");
			System.exit(-1);
		}
		new EPProfile(arg[0],arg[1],slen);
		System.exit(0);
	}
}
