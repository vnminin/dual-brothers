import cern.jet.random.*;

public abstract class Tree {
	public Node[] Leaf;
	public Node[] Branch;
	public int numLeaves;
	public int numBranches;
	public Node Root;

	public double LogPrior() {
		return Math.log(1 / Span());
	}

	public double Span() {
		int n = numLeaves;
	//		n--; // Assumes unrooted!
		return (double)( factorial(2*n-5)/(Math.pow(2,n-3)*factorial(n-3)) );
	}

	public double LogSpan() {
		int n = numLeaves;
//			n--; // Assumes unrooted!!!
		return   logFactorial(2*n-5) - (n-3)*Math.log(2) - logFactorial(n-3);
	}

	private int factorial(int n) {
		int rtn = 1;
		for(int i=n; i>1; i--)
			rtn *= i;
		return rtn;
	}

	private double logFactorial(int n) {
		double rtn = 0;
		for(int i=n; i>1; i--)
			rtn += Math.log(i);
		return rtn;
	}

	public void SnapInBranch(Node current) {
		Branch[numBranches++] = current;
	}

	public void SnapInLeaf(Node current, int ID) {
		Leaf[ID] = current;
		numLeaves++;
	}

	public void BalanceTree() {
		Root.Right.Balance();
	}

	public boolean equals(Tree test) {
		boolean rtn = false;
		String testStr = test.toStringNoBranches();
		String thisStr = this.toStringNoBranches();
	//	System.err.println(testStr+" "+thisStr);
		if( thisStr.compareTo(testStr) == 0 )
			rtn = true;
		return rtn;
	}

	protected void SproutLeaf(Node current, Node neighbor, int ID, AbstractDistribution unif) {
		Node newN = new Node(neighbor.Up,current,neighbor,-1, /*inStock, inStock[0],*/ ID);
		// Determine if neighbor is a left or right child or root
		if( neighbor == Root ) {
			Root = neighbor.Up = current.Up = newN;
			Branch[numBranches] = newN;
		} else {
			Branch[numBranches]  = current.Up = newN;
			if( neighbor == (neighbor.Up).Left ) // true if left child
				(neighbor.Up).Left = newN;
			else // true if right child
				(neighbor.Up).Right = newN;
			neighbor.Up = newN;
		}
		// Randomly assign right/left for newN.
		if( unif.nextDouble() < 0.5 ) {
			Node temp = newN.Left;
			newN.Left = newN.Right;
			newN.Right = temp;
		}
	}
	protected Node CopyDownTree(Node inNode, Node inParent) {
		Node node = null;
		if( inNode.isBranch() ) { // This is a branching node
			node = new Node(inParent,null,null,-1,inNode.uID,inNode.isLikelihoodDone);
			SnapInBranch(node);
			node.Left = CopyDownTree(inNode.Left,node);
			node.Right = CopyDownTree(inNode.Right,node);
		} else { // This is a leaf
			node = new Node(inParent,null,null,inNode.ID,inNode.uID,inNode.isLikelihoodDone);
			Leaf[inNode.ID] = node;
			numLeaves++;
		}
		return node;
	}

	abstract Tree ProposedJump(int method, boolean fix, AbstractDistribution unif);
	abstract Tree JointBranchAndTopology(AbstractDistribution nglobal, AbstractDistribution nlocal, AbstractDistribution unif,
						double pmix, boolean fixRoot, int bnum,Boolean update);
	abstract double LogLikelihood0(QMatrix Q, double mean, int[] counts, int[][] data, double[] pi, boolean handleGaps, double K, int bnum);
	abstract Tree[] EnumerateLastTaxon();
	abstract double SumOfBranchLengths();
	public abstract String toString();
	public abstract String toStringNoBranches();
	//public abstract boolean equals(Tree test);
}
