import java.io.Serializable;
import java.util.*;

/**
* Generates the sufficient statistics of the observed data.
*/
public class CountStatistic  implements Serializable{
	public int lenseq;
	public int lenunique;
	public int ntaxa;
	public DNASequence[] sequences;
	private int[][] data;
	private int[] counts;
	private int[] map;
	//private TreeSet sorted;
	//private Site comparor;
	public int[][] dataMatrix() { return data; }
	public int[] countsMatrix() { return counts; }
	public int[] siteMap() { return map; }

	public CountStatistic(DNASequence[] inSeq) {

	    sequences = inSeq;
	    lenseq = sequences[0].Length();
	    ntaxa = inSeq.length;
	    Site comparor = new Site();
	    TreeSet tsorted = new TreeSet(comparor);
	    for(int i=0; i<lenseq; i++) {
		int[] sdata = new int[ntaxa];
		for(int j=0; j<ntaxa; j++) {
		    sdata[j] = sequences[j].Strand[i];
		}
		tsorted.add(new Site(sdata));
	    }
	    lenunique = tsorted.size();
	    System.err.println("Number of unique sites: "+lenunique);
	    Vector lsorted = new Vector(tsorted);
	    map    = new int[lenseq];
	    counts = new int[lenunique];
	    for(int i=0; i<lenseq; i++) {
		int[] sdata = new int[ntaxa];
		for(int j=0; j<ntaxa; j++) {
		    sdata[j] = sequences[j].Strand[i];
		}
		Site probe = new Site(sdata);
		int loc = Collections.binarySearch(lsorted,probe,comparor);
		map[i] = loc;
			counts[loc]++;
	    }
	    data = new int[lenunique][];
	    for(int i=0; i<lenunique; i++) {
		data[i] = new int[ntaxa];
		for(int j=0; j<ntaxa; j++) {
		    data[i][j] = ((Site)(lsorted.get(i))).data[j];
		    System.err.print(DNASequence.Int2Nucleotide(data[i][j])+" ");
		}
		System.err.println(": "+counts[i]);
	    }
	}

    /** Loads data from the file with the name dname
     * @param dname name of the data file
     */


    public CountStatistic(String data_file_name) {
	
	sequences = DNASequence.ReadPhyllip(data_file_name,false);
	
	lenseq = sequences[0].Length();
	ntaxa = sequences.length;
	Site comparor = new Site();
	TreeSet tsorted = new TreeSet(comparor);
	for(int i=0; i<lenseq; i++) {
	    int[] sdata = new int[ntaxa];
	    for(int j=0; j<ntaxa; j++) {
		    sdata[j] = sequences[j].Strand[i];
	    }
	    tsorted.add(new Site(sdata));
	}
	lenunique = tsorted.size();
	System.err.println("Number of unique sites: "+lenunique);
	Vector lsorted = new Vector(tsorted);
	map    = new int[lenseq];
	counts = new int[lenunique];
	for(int i=0; i<lenseq; i++) {
	    int[] sdata = new int[ntaxa];
	    for(int j=0; j<ntaxa; j++) {
		sdata[j] = sequences[j].Strand[i];
	    }
	    Site probe = new Site(sdata);
	    int loc = Collections.binarySearch(lsorted,probe,comparor);
	    map[i] = loc;
			counts[loc]++;
	}
	data = new int[lenunique][];
	for(int i=0; i<lenunique; i++) {
	    data[i] = new int[ntaxa];
	    for(int j=0; j<ntaxa; j++) {
		data[i][j] = ((Site)(lsorted.get(i))).data[j];
		System.err.print(DNASequence.Int2Nucleotide(data[i][j])+" ");
	    }
	    System.err.println(": "+counts[i]);
	}
	
	System.err.println("Loaded data file: "+data_file_name);
    }
    
    public static void main(String[] arg) {
	DNASequence[] data = DNASequence.ReadPhyllip(arg[0],false);
	new CountStatistic(data);
    }
    
}
