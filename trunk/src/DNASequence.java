import java.util.*;
import java.io.*;

/**
* Input and output operations for nucleotide sequences.
*/
public class DNASequence extends Object {
	// All indices for nucleotide-types are array[4], where:
	//		0 = A
	//		1 = G
	//		2 = C
	//		3 = T
	//
	//	       -9 = gap....
	//
	private int numBases = 4;

	// Default nucleotide composition proportions from A. salina (Eukaryote).
	// Values from Schadt (1998).

	// Actually values are maximum likelihood estimates from same paper.
	//private static double Priors[][] = { {0.255, 0.279, 0.234, 0.232},{'A', 'G', 'C', 'T'} };
	public static double Priors[][] = { {0.32, 0.10, 0.31, 0.27},{'A', 'G', 'C', 'T'} };
	public static int baseLength = Priors[0].length;

	/**
	* Internal representation of sequence.
	* Translation: 0=A, 1=G, 2=C, 3=T, -9=all else.
	*/
	public int Strand[];

	public String name;

	public int span() {
		return numBases;
	}
	public DNASequence(DNASequence inSeq) {
		int l = inSeq.Strand.length;
		Strand = new int[l];
		for(int i=0; i<l; i++)
			Strand[i] = inSeq.Strand[i];
	}
	public DNASequence(int inLength) {
		// Initializer to create a random sequence of length inLength
		Strand = new int[inLength];

		// Calculate the total of the prior probabilities; this routine does not assume that
		// they sum => 1.
		double pSum = 0;
		for(int i=0; i<baseLength; i++)
			pSum += Priors[0][i];

		// Iterate thru each elements of the strand and randomly generate a base given the priors
		for(int i=0; i<Strand.length; i++)
			Strand[i] = Stat.RandomGivenPrior(Priors[0],pSum);
	}

	public DNASequence(int inLength, double[] rho) {
		// Initializer to create a random sequence of length inLength
		Strand = new int[inLength];

		// Calculate the total of the prior probabilities; this routine does not assume that
		// they sum => 1.
		double pSum = 0;
		for(int i=0; i<baseLength; i++)
			pSum += rho[i];

		// Iterate thru each elements of the strand and randomly generate a base given the priors
		for(int i=0; i<Strand.length; i++)
			Strand[i] = Stat.RandomGivenPrior(rho,pSum);
	}
	public DNASequence(String inStr) {
		// Initializer to allocate a sequence in memory from a string
		Strand = new int[inStr.length()];
		for(int i=0; i<inStr.length(); i++ )
			Strand[i] = Nucleotide2Int(inStr.charAt(i));
	}
	public String toString() {
		//return new String(Strand);
		char rtn[] = new char[Strand.length];
		for(int i=0; i<Strand.length; i++)
			rtn[i] = (char) Priors[1][ Strand[i] ];
		return new String(rtn);
	}

	public static char Int2Nucleotide(int i) {
		if( i == -9 )
			return '-';
		else
			return (char) Priors[1][i];
	}

	public static int Nucleotide2Int(char inChar) {
		int found = -9;
		for(int i=0; i<baseLength; i++) {
			if( inChar == (char) Priors[1][i] ) {
				found = i;
				break;
			}
		}
		// -9 should throw up some OutOfBoundsArrayExceptions somewhere to catch the error.
		return found;
	}

	public int Length() {
		return Strand.length;
	}
	public int getBase(int inBase) {
		return Strand[inBase];
	}
	public char getBaseChar(int inBase) {
		if( Strand[inBase] != -9 )
			return (char) Priors[1][ Strand[inBase] ];
		else
			return '-';
	}
	public double Composition(int in) {
		double rtn = 0;
		double len = 0;
		for(int i=0; i<Strand.length; i++) {
			if( Strand[i] == in )
				rtn++;
			if( Strand[i] != -9 )
				len++;
		}
		return (rtn/len);
	}
	public double Composition(char inChar) {
		double rtn = 0;
		double len = 0;
		int comp = Nucleotide2Int(inChar);
		for(int i=0; i<Strand.length; i++) {
			if( Strand[i] == comp )
				rtn++;
			if( Strand[i] != -9 )
				len++;
		}
		return (rtn/len);
	}
	//public double[][] CompBySite() {

	/**
	* Reads in a set of sequences in Phylip 4.0 format.
	*/
	public static DNASequence[] ReadPhyllip(String fname, boolean noGaps) {
		BufferedReader br = null;
		try {
			br = new BufferedReader( new FileReader(fname));
		} catch (FileNotFoundException e) {
			System.err.println("Error reading file: " + fname);
			System.exit(-1);
		}
		int numSeqs = 0;
		try {
			StringTokenizer st = new  StringTokenizer(br.readLine());
			numSeqs = Integer.valueOf(st.nextToken()).intValue();
		} catch (Exception e) {
			System.err.println("Unable to determine the number of sequences.");
			System.exit(-1);
		}
		String[] strSeq = new String[numSeqs];
		String[] allNames = new String[numSeqs];
		for(int i=0; i<numSeqs; i++)
			strSeq[i] = new String();
		String str;
		try {
			int cSeq = 0;
			while( (str = br.readLine())!=null ) {
				StringTokenizer st = new StringTokenizer(str);
				if( st.hasMoreTokens() ) {
					if( allNames[cSeq] == null ) {
						String name = st.nextToken(); // Just pass by the name for now.
						allNames[cSeq] = name;
					}
					while  (st.hasMoreTokens()) {
						String seg = st.nextToken();
						try {
							Integer i = Integer.valueOf(seg); // Ignore numbers at the end
						} catch (NumberFormatException e) {
							strSeq[cSeq] += seg;
						}
					}
					cSeq++;
					if( cSeq == numSeqs )
						cSeq = 0;
				} else
					cSeq = 0;
			}
		} catch (IOException e) {
			System.err.println("Error reading sequences data.");
			System.exit(-1);
		}
		// Now check for and remove gaps, assuming aligned sequences.
 		if( noGaps ) {
			int l=0;
			while(l < strSeq[0].length() ) {
			//	System.out.print(l+" ");
					boolean remove = false;
				for(int j=0; j<numSeqs; j++) {
					if( Nucleotide2Int(strSeq[j].charAt(l)) == -9 )
						remove = true;
				}
				if( remove ) {
					for(int j=0; j<numSeqs; j++) {
						strSeq[j] = strSeq[j].substring(0,l) + strSeq[j].substring(l+1);
					}
				} else
					l++;
			}
		} else {
			System.err.println("Loading data with gaps included.");
		}

		//for(int i=0; i<5; i++)
		//	System.out.println(strSeq[i]);
		// First we assume these are aligned sequences with NO GAPS
		System.err.println("Loaded data for taxa:");
		DNASequence[] all = new DNASequence[numSeqs];
		for(int i=0; i<numSeqs; i++) {
			all[i] = new DNASequence(strSeq[i]);
			all[i].name = allNames[i];
			System.err.println(all[i].name);
		}
		return all;
	}
	public static void PrintJanet(DNASequence[] inSeq, PrintStream out) {
		int numSeqs = inSeq.length;
		int idx = 0; // Start display at base 0.
		int[] counts = new int[numSeqs]; // These are initiated as 0.
		int[] lengths = new int[numSeqs];
		int basesRemaining = 0;
	//	System.out.println(numSeqs);
		for(int i=0; i<numSeqs; i++) {
			counts[i] = 0;
			lengths[i] = inSeq[i].Length();
			if( lengths[i] > basesRemaining )
				basesRemaining = lengths[i];   // Now basesRemaining equals the longest sequence
		}
		// Build concensus sequence
		char con[] = new char[basesRemaining];
		for(int i=0; i<basesRemaining; i++) { // Repeat through each base
			int c[] = new int[4]; // Fills with zeros
			for(int j=0; j<numSeqs; j++) { // Examine each sequences
				if( i < lengths[j] )
					c[inSeq[j].Strand[i]]++;
			}
			// How to determine which is greatest?
			// Start searching at c[] = numSeqs, and decrease by one until we find the first highest
			int k = numSeqs;
			int ch = 0;
			boolean found = false;
			while( !found ) {
				for(int l=0; l<4; l++) {
					if( c[l] == k ) {
						found = true;
						ch = l;
					}
				}
				k--;
			}
			con[i] = (char) Priors[1][ch];
		}
	//	System.out.println(String.valueOf(con));
		int conprint = 0;
		while( basesRemaining > 0 ) {
			// The following prints a single page.
			// First print the concensus
			int countme = 0;
			char segcon[] = new char[68];
			for(int i=0; i<68; i++) {
				if( (conprint+i) < con.length ) {
					segcon[i] = con[conprint+i];
					countme++;
				}
				else
					segcon[i] = ' ';
			}
						System.out.println("Concensus "+String.valueOf(segcon)+" "+(conprint+countme));
			conprint += 68;
			for(int i=0; i<numSeqs; i++) {
				int colRemaining = basesRemaining;
				String name = new String("T"+String.valueOf(i)); // Each taxa now has the name Tx
				//for(int j=name.length(); j<10; j++)
				//	name += " ";

				out.print(name+" "); // Let's first try to use a tab. If that doesn't work use above
				char[] seg = new char[68];
				for(int j=0; j<68; j++) {  // There are 68 bases per line
					//char[] seg = new char[68];
					if( counts[i] < lengths[i] ) {
						if( con[counts[i]] != inSeq[i].getBaseChar(counts[i]) )
							seg[j] = inSeq[i].getBaseChar(counts[i]++);
						else {
							seg[j] = '-';
							counts[i]++;
						}

					} else {
						if( colRemaining > 0 )
							seg[j] = '*';
						else
							seg[j] = ' ';
					}
					colRemaining--;
				}
				out.print(String.valueOf(seg));
				out.println(" "+counts[i]);
			}
			out.println(); // A blank line between pages.
			basesRemaining -= 68;
		}
	}
	public static void PrintPhyllip(DNASequence[] inSeq, PrintStream out) {
		int numSeqs = inSeq.length;
		out.println(numSeqs+"  "+inSeq[0].Length()); // First line is the total number of sequences and the total length of each sequence.
		int idx = 0; // Start display at base 0.
		int[] counts = new int[numSeqs]; // These are initiated as 0.
		int[] lengths = new int[numSeqs];
		//t[] indices = new int[numSeqs];
		int basesRemaining = 0;
		for(int i=0; i<numSeqs; i++) {
			counts[i] = 0;
			//indices[i] = 0;
			lengths[i] = inSeq[i].Length();
			if( lengths[i] > basesRemaining )
				basesRemaining = lengths[i];
		}
		//int colRemaining = basesRemaining;
		// loop pages  until end
		// not implemented yet.
		//basesRemaining *= numSeqs;
		while( basesRemaining > 0 ) {
			// The following prints a single page.
			for(int i=0; i<numSeqs; i++) {
				int colRemaining = basesRemaining;
				//String line = new String("                                                                                                     "); // an 80?? char line
				// Janet's program needs the sequence between 11 and 65 spaces from the left.
				// Columns 1-10 are the names.
				String name = new String("S"+String.valueOf(i)); // Each taxa now has the name Tx
				for(int j=name.length(); j<10; j++)
					name += " ";
				//System.out.println("!"+name+"!");
				//System.out.println("AAAAAAAAAAAAAAAAAAAAAAAAAAA");
				out.print(name); // Let's first try to use a tab.
				//for(int j=String.valueOf(i).length(); j<10; j++)
				//	out.print(" ");
				for(int j=0; j<5; j++) {
					char[] seg = new char[10];
					for(int k=0; k<10; k++) {
						//basesRemaining--;
						if( counts[i] < lengths[i] ) {
							//seg[k] = inSeq[i].getBaseChar(k);
							//counts[i]++;
							seg[k] = inSeq[i].getBaseChar(counts[i]++);
						} else {
							if( colRemaining > 0 )
								seg[k] = '-';
							else
								seg[k] = ' ';
						}
						colRemaining--;
						// check to see if the value is missing, if not then increase count and print
						// check is not yet implemented.
					}
					out.print(String.valueOf(seg));
					out.print(" "); // space after every ten bases
				}
				out.println("\t"+counts[i]);
			}
			out.println(); // A blank line between pages.
			basesRemaining -= 50;
		}
	}
}

