import jal.Object.*;
import java.util.*;

public class Site extends Object implements Comparator 
{

	public int[] data;
	
	public Site(int[] indata) {
		data = indata;
	}
	
	public Site(int l) {
		data = new int[l];
	}
	
	public Site() {
	}
	
	public int compare(Object i1, Object i2) {
		int rtn = 0;
		Site s1 = (Site) i1;
		Site s2 = (Site) i2;
		int len = s1.data.length;
		for(int i=0; i<len && rtn == 0; i++) {
			if( s1.data[i] < s2.data[i] )
				rtn--;
			else if( s1.data[i] > s2.data[i] )
				rtn++;
		}
		return rtn;
	}

}
