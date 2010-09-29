/**
* Class BreakPoint marks points as topology break-points, evolutionary  
* change-points or both
*/

public class BreakPoint {

    /**
     * position of a change-point
     */

    int nuc_position;
    
    /**
     * indicator of topology break-point
     */
    
    boolean top_change;
    
    /**
     * indicator of evolutionary change-point
     */
    
    boolean par_change;

    BreakPoint(){
	nuc_position = 0;
	top_change = true;
	par_change = true;
    }

    BreakPoint(int x, boolean t, boolean p){
	nuc_position = x;
	top_change = t;
	par_change = p;
    }


}
