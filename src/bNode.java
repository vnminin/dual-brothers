import java.util.*;
//import lava.clib.stdio.*;
import corejava.*;

/**
* Provides operations on nodes with branch lengths.
*/
public class bNode //extends Object implements Cloneable
{

	//static PrintfFormatString pfs1 = new PrintfFormatString("%d:%.6f");
	//static PrintfFormatString pfs2 = new PrintfFormatString("):%.6f");
	static Format fmtDouble = new Format("%.6f");

	public bNode Up;
	public bNode Right;
	public bNode Left;
	public int ID;
//	public int uID;
	public boolean isBranch;
	public double Stock = -9;
	//public double Stock2;
	public boolean isLikelihoodDone = false;
	//public double calcStock;
	public double[] ConditionalLikelihood;
	//public double[][] cl;


	public bNode(bNode inUp, bNode inRight, bNode inLeft, int inID)
	{
		Up = inUp;
		Right = inRight;
		Left = inLeft;
		ID = inID;
		if( ID == -1 ) // This is a new branch
			isBranch = true;
		else {
			isBranch = false;
			//uID = ID;
		}
		//Stock = null; // Undefined.
		//Stock2 = -9;
		//calcStock = -9;
		// Might want to comment out next line in all Node creators.
		ConditionalLikelihood = new double[4];
	}
	public bNode(bNode inUp, bNode inRight, bNode inLeft, int inID, double inStock)
	{
		Up = inUp;
		Right = inRight;
		Left = inLeft;
		ID = inID;
		//uID = inuID;
		if( ID == -1 ) // This is a new branch
			isBranch = true;
		else
			isBranch = false;
		Stock = inStock;
		//Stock2 = inStock2;
		//calcStock = -9;
		ConditionalLikelihood = new double[4];
	}
/*	public Node(Node inUp, Node inRight, Node inLeft, int inID, double inStock, int inuID,boolean yn)
	{
		Up = inUp;
		Right = inRight;
		Left = inLeft;
		ID = inID;
		uID = inuID;
		if( ID == -1 ) // This is a new branch
			isBranch = true;
		else
			isBranch = false;
		Stock = inStock;
		//Stock2 = inStock2;
		//calcStock = -9;
		//isLikelihoodDone = yn;
		ConditionalLikelihood = new double[4];
	}
*/	public bNode(String inStr, bTree inTree, bNode inParent, boolean isroot)
	{
		Up = inParent;
		//Stock = 1;
		//Stock = new double[nStocks];
		//for(int i=0; i<nStocks; i++)
		//   Stock[i] = 1;
		Stock = 1;
		ConditionalLikelihood = new double[4];
		if( inStr.charAt(0) == '(' )
		{	// This is a branch
			isBranch = true;
			ID = -1;
			if( isroot )
				inTree.Root = this;
			else
				inTree.SnapInBranch(this);
			// Check to see if stock length is given.
			if( inStr.charAt(inStr.length()-1) != ')' )
			{ // We probably have a stock length to parse
				int i = inStr.length() - 1;
				while( inStr.charAt(i) != ':' )
					i--;
			//	try
			//	{ // Conversion from String into Double may cause an error
					Double tmp = new Double(inStr.substring(i+1,inStr.length()));
					//for(int j=0; j<nStocks; j++)
					Stock = tmp.doubleValue();
					//Stock2 = Stock[0];
			//	} catch (Exception e)
			//	{
			//		Stock = -9;
			//		System.err.println("Error parsing branch stock length.");
			//	}
				inStr = inStr.substring(0,i);
			}
			// Find dividing ',' for this node
			int nLevel = 0;
			int i = 0;
			char current = inStr.charAt(i);
			while( (current != ',') || (nLevel != 1) )
			{
				if( current == '(' )
					nLevel++;
				if( current == ')' )
					nLevel--;
				current = inStr.charAt(++i);
			}
			// Recursivedly add on children
			Left = new bNode(inStr.substring(1,i),inTree,this,false);
			Right = new bNode(inStr.substring(i+1,inStr.length()-1),inTree,this,false);
		}
		else
		{
			// This is a new leaf
			// First check to see if a stock length is given.  But it's easier to do this here
			// as the form is always either #:LENGTH or #
			// First we break the string into segments based on ":"
			// When we try to get the second element, if it doesn't exist, an exception is called
			// which is caught and stock length is defaulted to -9
			StringTokenizer tok = new StringTokenizer(inStr,":");
			String intStr = tok.nextToken();
			try
			{
				Double temp = new Double(tok.nextToken());
				//for(int j=0; j<nStocks; j++)
				Stock = temp.doubleValue();
				//Stock2 = Stock[0];
			} catch (NoSuchElementException e)
			{
				//Stock = -9;  already defaulted above.
			}
			Integer temp = new Integer(intStr); // This is liable to cause an exception
			ID = temp.intValue();
			Left = Right = null;
			inTree.SnapInLeaf(this,ID);
		}
	}
/*	public void ForceDownReprocess() {
		this.isLikelihoodDone = false;
		if( this.Left != null ) {
			Left.ForceDownReprocess();
			Right.ForceDownReprocess();
		}
	}
*/
/*	public String toStringJustStocks() {
		if( isBranch ) {
			if( Stock == null )
				return Left.toStringJustStocks()+Right.toStringJustStocks();
			else
				return Left.toStringJustStocks()+Right.toStringJustStocks()+" "+Stock[i];
		} else {
			return " "+Stock;
		}
	}
*/

	public String toStringWithStocks()
	{
		if( isBranch )
		{
			if( Stock == -9 )
				return "("+Left.toStringWithStocks()+","+Right.toStringWithStocks()+")";
			else
				return "("+Left.toStringWithStocks()+","+Right.toStringWithStocks()+"):"+
				//	Printf.format(pfs2, new Object[] { new Double(Stock) } );
				//"):"+Stock;
				fmtDouble.form(Stock);
		}
		else
		{
			if( Stock == -9 )
				return String.valueOf(ID);
			else
				return //ID+":"+Stock;
					//Printf.format(pfs1, new Object[] { new Integer(ID), new Double(Stock) } );
					ID+":" + fmtDouble.form(Stock);
		}
	}
/*	public String toStringWithStocks2()
	{
		if( isBranch )
		{
			if( Stock2 == -9 )
				return "("+Left.toStringWithStocks2()+","+Right.toStringWithStocks2()+")";
			else
				return "("+Left.toStringWithStocks2()+","+Right.toStringWithStocks2()+
					Printf.format(pfs2, new Object[] { new Double(Stock2) } );
				//"):"+Stock2;
		}
		else
		{
			if( Stock2 == -9 )
				return String.valueOf(ID);
			else
				return //ID+":"+Stock2;
					Printf.format(pfs1, new Object[] { new Integer(ID), new Double(Stock2) } );
		}
	}
*/	public String toString()
	{
		if( isBranch )
			return "("+Left.toString()+","+Right.toString()+")";
		else
			return String.valueOf(ID);
	}

	public int Balance() {
		//Node temp;
		if( isBranch ) {
			int i = Left.Balance();
			int j = Right.Balance();
			if( i < j )
				return i;
			else {
				bNode temp = Left;
				Left = Right;
				Right = temp;
				return j;
			}
		} else
			return ID;
	}

	public void fillVectorWithIDs(Vector v) {
		if( isBranch ) {
			Left.fillVectorWithIDs(v);
			Right.fillVectorWithIDs(v);
		} else {
			v.addElement(new Integer(ID));
		}
		return;
	}
	public void drawIDFromVector(Vector v,bNode[] array) {
		if( isBranch ) {
			Left.drawIDFromVector(v,array);
			Right.drawIDFromVector(v,array);
		} else {
			int random = (int)(Math.random()*v.size());
			Integer i = (Integer) v.elementAt(random);
			v.removeElementAt(random);
			ID = i.intValue();
			array[ID] = this;
		}
		return;
	}
	public void setParent(bNode inUp)
	{
		Up = inUp;
	}
	public void setRight(bNode inRight)
	{
		Right = inRight;
	}
	public void setLeft(bNode inLeft)
	{
		Left = inLeft;
	}
	public bNode getParent()
	{
		return Up;
	}
	public bNode getRight()
	{
		return Right;
	}
	public bNode getLeft()
	{
		return Left;
	}
	public bNode getBrother() {
		if (this == (this.Up).Right)
			return (this.Up).Left;
		return (this.Up).Right;
	}
	public boolean isLeftChild() {
		if (this == Up.Left)
			return true;
		return false;
	}
	public bNode getSelf() {
		return this;
	}
	public double getHeight() {
		double rtn = 0;
		bNode current = this;
		while( current.getParent() != null ) {
			rtn += current.getStock();
			current = current.getParent();
		}
		return rtn;
	}

	public int getNumberNodes() {
		int rtn = 0;
		bNode current = this;
		while( current.Up != null ) {
			rtn++;
			current = current.Up;
		}
		return rtn;
	}

	public int getNumberChildren() {
		if( this.Left != null )
			return(1 + Left.getNumberChildren() + Right.getNumberChildren() );
		else
			return 0;
	}


	public boolean[] LowerPartition(int N) {
	    boolean[] result = new boolean[N];
	    String me = this.toString();
	    if( this.Right != null ) {
		for(int i=0; i<N; i++) {
		    if( (me.indexOf(i+",") != -1) || (me.indexOf(i+")") != -1) ) {
			result[i] = true;
		    }
		}
	    } else {
		result[this.ID] = true;
	    }
	    return result;
	}

	public boolean[] UpperPartition(int N) {
	    boolean[] result = new boolean[N];
	    if( this.Right == null ) {
		for(int i=0; i<N; i++) {
		    if( i == this.ID )
			result[i] = false;
		    else
			result[i] = true;
		}
	    } else {
		String me = this.toString();
		for(int i=0; i<N; i++) {
		    if( (me.indexOf(i+",") != -1) || (me.indexOf(i+")") != -1) ) {
			result[i] = false;
		    } else {
			result[i] = true;
		    }
		}
	    }
	    return result;
	}

/*
	public double getHeightBeforeNode(Node n) {
		double rtn = 0;
		Node current = this;
		while( current.getParent() != n ) {
			rtn += current.getStock();
			current = current.getParent();
		}
		return rtn;
	}
*/
	public int CountDown() {
		int rtn = 1;
		if( Right != null )
			rtn += Right.CountDown();
		if( Left != null )
			rtn += Left.CountDown();
		return rtn;
	}

	public bNode GetNthChild(int n) {
		if( n == 0 )
			return this;
		bNode result;
		int nleft = Left.CountDown();
		if( n <= nleft ) { // then the resultant is somewhere on the left
			result = Left.GetNthChild(n-1);
		} else { // resultant is somewhere on the right
			result = Right.GetNthChild(n-1-nleft);
		}
		return result;
	}

	public int getID()
	{
		return ID;
	}
	public void setID(int inID)
	{
		ID = inID;
	}
	public boolean isBranch()
	{
		return isBranch;
	}
	public double getStock()
	{
		return Stock;
	}
	public void removeStock()
	{
		Stock = -9;
		//Stock2 = -9;
	}
	public void setStock(double in)
	{
		Stock = in;
		//Stock2 = in[0];
	}

/*

  Since I made all of the data bits public, most of the above
  routines can be removed.

*/

}
