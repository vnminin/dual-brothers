
public class Stat {
	public static int RandomGivenPrior(double[] inPrior) {
		// Returns the index of the random event given priors, with marginal = 1
		double random =  Math.random();  // Returns a 0 > uniform value <= pSum
		int j = -1;
		do 
			random -= inPrior[++j];
		while( random > 0 ); // Should this really be 'random >= 0'????
		return j;
	}
	public static int RandomGivenPrior(double[] inPrior, double pSum) {
		double random = pSum * Math.random();  // Returns a 0 > uniform value <= pSum
		int j = -1;
		do 
			random -= inPrior[++j];
		while( random > 0 ); // Should this really be 'random >= 0'????
		return j;
	}
	public static int RandomGivenPriorCalcMargin(double[] inPrior) {
		double pSum = 0;
		for(int i=0; i<inPrior.length; i++)
			pSum += inPrior[i];
		double random = pSum * Math.random();  // Returns a 0 > uniform value <= pSum
		int j = -1;
		do 
			random -= inPrior[++j];
		while( random > 0 ); // Should this really be 'random >= 0'????
		return j;
	}
}
