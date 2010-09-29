import java.util.*;

public class TreeCount extends Object implements Comparable
{

	public String tree;
	public int count;
	
	public TreeCount(String s) {
		tree = s;
		count = 1;
	}
	
	public String toString() {
		return tree+" "+count;
	}
	
	public boolean equals(Object i1) {
		TreeCount tc1 = (TreeCount) i1;
		if( tree.compareTo(tc1.tree) == 0 )
			return true;
		else
			return false;
	}
	
	public int compareTo(Object i1) {
		TreeCount tc1 = (TreeCount)i1;
		int rtn = 0;
		if( count < tc1.count )
			rtn = 1;
		else if( count > tc1.count)
			rtn = -1;
		return rtn;
	}
}
