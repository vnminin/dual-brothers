/*

Started rewritting as an optimized copy.


*/

import java.util.*;
//import lava.clib.stdio.*;
import cern.jet.random.*;



/**
* Provides operations on tree with branch lengths.
*/
public class bTree extends Tree {

	public bNode[] Leaf;
	public bNode[] Branch;
	//public int numLeaves;
	//public int numBranches;
	public bNode Root;
	public double QQ = 1.0;


	public void BalanceTree() {
		Root.Right.Balance();
	}

	Tree[] EnumerateLastTaxon() {
		System.err.println("Enumeration not yet implemented for trees with branches.");
		System.exit(-1);
		return null;
	}


	public double LogLikelihood0(QMatrix Q, double mean, int[] counts, int[][] data, double[] pi, boolean handleGaps, double K, int bnum) {
		// Returns the log likelihood of the data in C, given the topology of this tree, Q and
		// a proposed distribution of bases at the root by summing over all possible states
		// for the internal nodes.

		int numSites = counts.length;
		double LL = 0;
		boolean[] repeats = new boolean[numLeaves]; // Automatically fills with false
		double[] tCL = new double[4]; // Automatically zeros all entries. for temporarily storage.
		bNode cNode = null;
		double stock = mean;
		for(int i=0; i<numSites; i++) {
			if( counts[i] != 0 ) {
			       for(int j=0; j<numLeaves; j++)
			       		Leaf[j].isLikelihoodDone = false;
			       for(int j=1; j<numBranches; j++)
			       		Branch[j].isLikelihoodDone = false;
			       for(int j=0; j<numLeaves; j++) { 	// Here, let's calculate the likelihood of a single word
					cNode = Leaf[j]; // Start at the bottom of each leaf
					if( data[i][cNode.ID] >= 0 ) {
						if( j != 0 )
							stock = cNode.Stock;
						else
							stock = Root.Right.Stock;
						cNode.ConditionalLikelihood[0] = Q.Pt(0,data[i][cNode.ID],stock);
						cNode.ConditionalLikelihood[1] = Q.Pt(1,data[i][cNode.ID],stock);
						cNode.ConditionalLikelihood[2] = Q.Pt(2,data[i][cNode.ID],stock);
						cNode.ConditionalLikelihood[3] = Q.Pt(3,data[i][cNode.ID],stock);
					} else {
						// Gaps.... use pi as prior on site.
						cNode.ConditionalLikelihood[0] = 1;
						cNode.ConditionalLikelihood[1] = 1;
						cNode.ConditionalLikelihood[2] = 1;
						cNode.ConditionalLikelihood[3] = 1;
					}
					cNode.isLikelihoodDone = true;
					while( cNode.getBrother().isLikelihoodDone && (cNode.Up.Up != Root) ) { // Propogate up tree until we hit the root or an unfinished section
						cNode = cNode.Up;
						bNode lChild = cNode.Left;
						bNode rChild = cNode.Right;
							stock = cNode.Stock;
							tCL[0] = lChild.ConditionalLikelihood[0] * rChild.ConditionalLikelihood[0];
							tCL[1] = lChild.ConditionalLikelihood[1] * rChild.ConditionalLikelihood[1];
							tCL[2] = lChild.ConditionalLikelihood[2] * rChild.ConditionalLikelihood[2];
							tCL[3] = lChild.ConditionalLikelihood[3] * rChild.ConditionalLikelihood[3];
							//cNode.ConditionalLikelihood[0] = 0;
							cNode.ConditionalLikelihood[0] = Q.Pt(0,0,stock) * tCL[0];
							cNode.ConditionalLikelihood[0] += Q.Pt(0,1,stock) * tCL[1];
							cNode.ConditionalLikelihood[0] += Q.Pt(0,2,stock) * tCL[2];
							cNode.ConditionalLikelihood[0] += Q.Pt(0,3,stock) * tCL[3];
							//cNode.ConditionalLikelihood[1] = 0;
							cNode.ConditionalLikelihood[1] = Q.Pt(1,0,stock) * tCL[0];
							cNode.ConditionalLikelihood[1] += Q.Pt(1,1,stock) * tCL[1];
							cNode.ConditionalLikelihood[1] += Q.Pt(1,2,stock) * tCL[2];
							cNode.ConditionalLikelihood[1] += Q.Pt(1,3,stock) * tCL[3];
							//cNode.ConditionalLikelihood[2] = 0;
							cNode.ConditionalLikelihood[2] = Q.Pt(2,0,stock) * tCL[0];
							cNode.ConditionalLikelihood[2] += Q.Pt(2,1,stock) * tCL[1];
							cNode.ConditionalLikelihood[2] += Q.Pt(2,2,stock) * tCL[2];
							cNode.ConditionalLikelihood[2] += Q.Pt(2,3,stock) * tCL[3];
							//cNode.ConditionalLikelihood[3] = 0;
							cNode.ConditionalLikelihood[3] = Q.Pt(3,0,stock) * tCL[0];
							cNode.ConditionalLikelihood[3] += Q.Pt(3,1,stock) * tCL[1];
							cNode.ConditionalLikelihood[3] += Q.Pt(3,2,stock) * tCL[2];
							cNode.ConditionalLikelihood[3] += Q.Pt(3,3,stock) * tCL[3];
						cNode.isLikelihoodDone = true;
					}
			       }
			       cNode = cNode.Up;
			       bNode lChild = cNode.Left;
			       bNode rChild = cNode.Right;
			       bNode fChild = Leaf[0];
			       double L;
			       L  = lChild.ConditionalLikelihood[0] * rChild.ConditionalLikelihood[0] * fChild.ConditionalLikelihood[0] * pi[0];
			       L += lChild.ConditionalLikelihood[1] * rChild.ConditionalLikelihood[1] * fChild.ConditionalLikelihood[1] * pi[1];
			       L += lChild.ConditionalLikelihood[2] * rChild.ConditionalLikelihood[2] * fChild.ConditionalLikelihood[2] * pi[2];
			       L += lChild.ConditionalLikelihood[3] * rChild.ConditionalLikelihood[3] * fChild.ConditionalLikelihood[3] * pi[3];

			       LL += Math.log(L) * counts[i];
			}
		}
		return LL;
	}

	public Tree Test() {
		return null;
	}

	public Tree JointBranchAndTopology(AbstractDistribution nglobal, AbstractDistribution nlocal, AbstractDistribution unif,
			double pmix, boolean fixRoot, int bnum,Boolean update) {
		bTree pTree = new bTree(this);
		//int span = Root.Right.Stock.length;
		if( unif.nextDouble() < pmix ) {
			for(int i=1; i<numLeaves; i++)
				pTree.Leaf[i].Stock = Math.abs(pTree.Leaf[i].Stock + nglobal.nextDouble());
			pTree.Branch[1].Stock = Math.abs(pTree.Branch[1].Stock + nglobal.nextDouble());
			for(int i=2; i<numBranches; i++) {
				bNode selected = pTree.Branch[i];
				selected.Stock = selected.Stock + nglobal.nextDouble();
				if( selected.Stock < 0 ) {
					selected.Stock = - selected.Stock;
					bNode parent = selected.Up;
					bNode moving = null;
					boolean SelIsLeftKid = false;
					if( selected == parent.Left )
						SelIsLeftKid = true;
					if( unif.nextDouble() < 0.5 ) {
						moving = selected.Left;
						if( SelIsLeftKid ) {
							// reattach as right
							bNode aunt = parent.Right;
							aunt.Up = selected;
							selected.Left = aunt;
							parent.Right = moving;
							moving.Up = parent;
						} else {
							// reattach as left;
							bNode aunt = parent.Left;
							aunt.Up = selected;
							selected.Left = aunt;
							parent.Left = moving;
							moving.Up = parent;
						}
					} else {
						moving = selected.Right;
						if( SelIsLeftKid ) {
							bNode aunt = parent.Right;
							aunt.Up = selected;
							selected.Right = aunt;
							parent.Right = moving;
							moving.Up = parent;
						} else {
							bNode aunt = parent.Left;
							aunt.Up = selected;
							selected.Right = aunt;
							parent.Left = moving;
							moving.Up = parent;
						}
					}
					pTree.BalanceTree();
				}
			}
		} else {
			boolean found = false;
			int r = 0;
			r = (int) (unif.nextDouble()*(2*numLeaves - 3));
			if( r == 0 )
				pTree.Root.Right.Stock = Math.abs(pTree.Root.Right.Stock + nlocal.nextDouble());
			else if( r < pTree.numLeaves )
				pTree.Leaf[r].Stock = Math.abs(pTree.Leaf[r].Stock + nlocal.nextDouble());
			else { // update an internal branch
				r -= pTree.numLeaves;
				r += 2;
				bNode selected = pTree.Branch[r];
				selected.Stock = selected.Stock + nlocal.nextDouble();
				if( selected.Stock < 0 ) {
					selected.Stock = - selected.Stock;
					bNode parent = selected.Up;
					bNode moving = null;
					boolean SelIsLeftKid = false;
					if( selected == parent.Left )
						SelIsLeftKid = true;
					if( unif.nextDouble() < 0.5 ) {
						moving = selected.Left;
						if( SelIsLeftKid ) {
							// reattach as right
							bNode aunt = parent.Right;
							aunt.Up = selected;
							selected.Left = aunt;
							parent.Right = moving;
							moving.Up = parent;
						} else {
							// reattach as left;
							bNode aunt = parent.Left;
							aunt.Up = selected;
							selected.Left = aunt;
							parent.Left = moving;
							moving.Up = parent;
						}
					} else {
						moving = selected.Right;
						if( SelIsLeftKid ) {
							bNode aunt = parent.Right;
							aunt.Up = selected;
							selected.Right = aunt;
							parent.Right = moving;
							moving.Up = parent;
						} else {
							bNode aunt = parent.Left;
							aunt.Up = selected;
							selected.Right = aunt;
							parent.Left = moving;
							moving.Up = parent;
						}
					}
					pTree.BalanceTree();
				}
			}
		}
		return pTree;
	}

	public double SumOfBranchLengths() {
		double rtn = 0;
		for(int i=1; i<numLeaves; i++)
			rtn += Leaf[i].Stock;
		for(int i=1; i<numBranches; i++)
			rtn += Branch[i].Stock;
		return rtn;
	}

	public Tree ProposedJump(int method, boolean fix, AbstractDistribution unif) {
		bTree pTree = null;
		int prandom;
		bNode pendant = null;
		switch( method ) {
		case 0:
			System.exit(-1);
			break;
		case 5: // Permute children on internal node, allows re-rooting.
			pTree = new bTree(this);
			// Select the parent node, can be any of n - 2 internal nodes or the root
			int r = (int)(unif.nextDouble()*numBranches);
			bNode parent = pTree.Branch[r];
			String current = parent.toString();
			String proposed = current;
			// Make a vector of current IDs from all offspring
			Vector IDs = new Vector();
			while( current.compareTo(proposed) == 0 ) {
				parent.fillVectorWithIDs(IDs);
				// Randomly sort permutation
				parent.drawIDFromVector(IDs,pTree.Leaf);
				proposed = parent.toString();
			}
			//QQ = new Double(1.0); // symmetric jump
			break;
		case 6: // Permute children on internal node, does NOT allow re-rooting.
			pTree = new bTree(this);
			// Select the parent node, can be any of n - 2 internal nodes
			while( (parent = pTree.Branch[(int)(unif.nextDouble()*numBranches)]) == pTree.Root )
				; // reselect
			//parent.ForceDownReprocess();
			// Make a vector of current IDs from all offspring
			IDs = new Vector();

			current = parent.toString();
			proposed = current;
			while( current.compareTo(proposed) == 0 ) {
				//System.err.println("GOT HERE AGAIN!!!");
				parent.fillVectorWithIDs(IDs);
				// Randomly sort permutation
				parent.drawIDFromVector(IDs,pTree.Leaf);
				proposed = parent.toString();
			}
			//QQ = new Double(1.0); // symmetric jump
			break;
/*		case 4: // Local Rearrangement based on LPD
			double oldLength = 0;
			pTree = new Tree(this);
			int prandom;
			Node pendant = null;
			// Choose a pendant leaf at random
			do {
				prandom = (int)(Math.random() * numLeaves);
				pendant = pTree.Leaf[prandom];
			} while( pendant.Up == pTree.Root );
			//System.err.println("Tree: "+pTree.Root.toString());
			if( pendant.Up != pTree.Root ) {
				//pendant.Up.ForceDownReprocess();
			//	System.err.println("Pendant: "+pendant.toString());
				Node pneighbor = pendant.getBrother();
				System.err.println("PNeighbor: "+pneighbor.toString());
				Node cut = pendant.Up;
				if( cut.isLeftChild() )
					(cut.Up).Left = pneighbor;
				else
					(cut.Up).Right = pneighbor;
				pneighbor.Up = cut.Up;
				double oldStock = pneighbor.Stock;
				//pneighbor.Stock = new double[oldStock.length];
				//for(int j=0; j<oldStock.length; j++) {
				pneighbor.Stock = oldStock + cut.Stock;
				//}
				System.err.println("Cut tree: "+pTree.Root.toString());
				// Pendant has now been removed from the tree.  pneighbor.Up points to
				// the internal node around which we reattach.
				// if center != Root, then there are two places to reattach beyond it
				// if center == Root, then there is only one place to reattach beyond it.
				// if pneight == Leaf, then there are no places to reattach beyond it
				// if pneighbor != lead, then there are two places to reattach beyond it
				// => 1,2,3, or 4 possible reattachment locations.

				if( pneighbor.Up == pTree.Root ) {
					// only one place beyond
					if( pneighbor.isBranch() ) {
						// two possible children
						prandom = (int) (Math.random() * 3);
						System.err.println("CASE A "+prandom);
					} else {
						prandom = 2;
						System.err.println("CASE C "+prandom);
					}
				} else {
					if( pneighbor.isBranch() ) {
						prandom = (int) (Math.random() * 4);
						System.err.println("CASE A "+prandom);
					} else {
						prandom = 2 + (int) (Math.random() * 2);
						System.err.println("CASE D "+prandom);
					}
				}
				// if the neighbor center is the root, then there are only
				// one place to reattach (near pneighbor), else one can also
				// reattach above the center.
				Node neighbor;
				if( prandom == 0 ) {
					// reattach between center and pneighbor
					neighbor = pneighbor.Left;
				} else if ( prandom == 1 ) {
					neighbor = pneighbor.Right;
				} else if ( prandom == 2 ) {
					neighbor = pneighbor.getBrother();
				} else {
					neighbor = pneighbor.Up;
				}

				// Let's now re-attach our pendant to it's neighbor.
				double rand = Math.random();
				int span = neighbor.Stock.length;
				double oldNeighborStock = neighbor.Stock;
				//neighbor.Stock = new double[span];
				//cut.Stock = new double[span];
				//for(int j=0; j<span; j++) {
				    //double split = neighbor.Stock * Math.random();
				cut.Stock = oldNeighborStock * rand;
				neighbor.Stock = oldNeighborStock * (1-rand);
				}
 				cut.Up = neighbor.Up;
				if( neighbor.isLeftChild() ) // true if left child
					(neighbor.Up).Left = cut;
				else // true if right child
					(neighbor.Up).Right = cut;
				// Add a little randomness into the local branch lengths.
				//NormalDistribution n = new NormalDistribution(0,0.1);
				//cut.Stock += n.inverse(Math.random());
				//neighbor.Stock += n.inverse(Math.random());
				//pendant.Stock += n.inverse(Math.random());
				//if( (cut.Stock <= 0) || (neighbor.Stock <= 0) || (pendant.Stock <= 0) )
				//	return null; // Bad jump

				if( pendant.isLeftChild() )
					cut.Right = neighbor;
				else
					cut.Left = neighbor;
				neighbor.Up = cut;
				//cut.ForceDownReprocess();
			}
			// Do nothing if we are around the root, as this would just be re-rooting.
			break;
	*/
		case 1: // Li Pearl and Doss's pendant rearrangement, does not re-root trees. supposedly.

			// 1. clone tree.
			pTree = new bTree(this);
			// 2. choose a pendant leaf at random
			//prandom = (int)(Math.random() * numLeaves);
			//pendant = pTree.Leaf[prandom];
			//if( pendant.Up != pTree.Root ) {
			do {
				pendant = pTree.Leaf[(int)(unif.nextDouble()*numLeaves)];
			} while (pendant.Up == pTree.Root);
				//penddant.Up.ForceDownReprocess();
				// No need to re-root.
				// 2b. clip out extra node.
			bNode pneighbor = pendant.getBrother();
			bNode cut = pendant.Up; // the branching node to be removed
			if( cut.isLeftChild() )  // cut is a left child
				(cut.Up).Left = pneighbor;
			else // cut is a right child
				(cut.Up).Right = pneighbor;
			pneighbor.Up = cut.Up;
			//int span = pneighbor.Stock.length;
			double oldPneighborStock = pneighbor.Stock;
			pneighbor.Stock = oldPneighborStock + cut.Stock; // Conserve stocks.
			// 3. choose a potential branch to which to reattach the pendant/
			// the pendant can re-attach as nearest-neighbor to the remaining n - 1
			// leaves or n - 3 internal nodes (does not include root).  This leaves
			// 2n - 4 locations.
			bNode neighbor;
			int random = (int)(unif.nextDouble() * (2*numLeaves - 4));
			if( random < (numLeaves - 1) ) { // next to a leaf, but skip pendant
				if( random < pendant.ID )
					neighbor = pTree.Leaf[random];
				else
					neighbor = pTree.Leaf[random+1];
				if( fix && (neighbor.Up == pTree.Root) ) {
					return null;
				}
			} else { // next to a branch, but skip cut and root
				random -= numLeaves;
				random++;
				int i = 0;
				for(int j=0; j<=random; j++) {
					if( pTree.Branch[j] == cut )
						i++;
					if( pTree.Branch[j] == pTree.Root )
						i++;
				}
				if( pTree.Branch[random+i] == pTree.Root )
					i++;
				if( pTree.Branch[random+i] == cut )
					i++;
				neighbor = pTree.Branch[random+i];
			}

			// Let's now re-attach our pendant to it's neighbor.
			// first calculate proposal density
			//System.err.println("new length = "+neighbor.Stock);
			//QQ = new Double( neighbor.Stock / oldLength );
			double oldNeighborStock = neighbor.Stock;
			//int span = neighbor.Stock.length;
			double rand = unif.nextDouble();
			//cut.Stock = new double[span];
			//neighbor.Stock = new double[span];
			pTree.QQ = 1;
			//for(int j=0; j<span; j++) {
			pTree.QQ *= neighbor.Stock / pneighbor.Stock;
			double split = oldNeighborStock * rand;
			cut.Stock = split;
			neighbor.Stock = oldNeighborStock - split;
			//}
 			cut.Up = neighbor.Up;
			if( neighbor.isLeftChild() ) // true if left child
				(neighbor.Up).Left = cut;
			else // true if right child
				(neighbor.Up).Right = cut;
			// Add a little randomness into the local branch lengths.
			if( pendant.isLeftChild() )
				cut.Right = neighbor;
			else
				cut.Left = neighbor;
			neighbor.Up = cut;
			break;
			//else { // Let's now deal with re-rooting.
		case 7: // cluster re-arrangement
			// 1. clone tree.
			pTree = new bTree(this);
			// 2. choose a pendant internal node at random
			//    exclude root and direct child of root
			do {
				prandom = (int)(unif.nextDouble() * numBranches);
				pendant = pTree.Branch[prandom];
			} while( (pendant == pTree.Root) || (pendant.Up == pTree.Root) );

			// now cut this internal node out of the tree.
			pneighbor = pendant.getBrother();
			cut = pendant.Up; // the branching node to be removed
			if( cut.isLeftChild() )  // cut is a left child
				(cut.Up).Left = pneighbor;
			else // cut is a right child
				(cut.Up).Right = pneighbor;
			pneighbor.Up = cut.Up;
			double oldPNeighborStock = pneighbor.Stock;
			//pneighbor.Stock = new double[oldPNeighborStock.length];
			//for(int j=0; j<oldPNeighborStock.length; j++ )
			 pneighbor.Stock = oldPNeighborStock + cut.Stock; // Conserve stocks.
			// 3. choose a potential branch/leaf to which to reattach the pendant/
			// the pendant can re-attach as nearest-neighbor to the remaining
			// leaves or internal nodes (does not include root).
			// how many are remaining on the tree?
			// try to count them recursively.
			int count = pTree.Root.CountDown();
			count--;
			if( fix ) // do not allow root (#0) or outgroup (#1)
				prandom = (int)(unif.nextDouble() * (count - 1)) + 2;
			else // do not all root (#0)
				prandom = (int)(unif.nextDouble() * count) + 1;
			neighbor = pTree.Root.GetNthChild(prandom);
			oldNeighborStock = neighbor.Stock;
			//span = neighbor.Stock.length;
			//neighbor.Stock = new double[span];
			//cut.Stock = new double[span];
			pTree.QQ = 1;
			rand = unif.nextDouble();
			//for(int j=0; j<span; j++) {
			    // Let's now re-attach our pendant to it's neighbor.
			split = oldNeighborStock * rand;
			//QQ = new Double( neighbor.Stock / oldLength );
			pTree.QQ *= oldNeighborStock / oldPNeighborStock;
			cut.Stock = split;
			neighbor.Stock = oldNeighborStock - split;
			//}
			cut.Up = neighbor.Up;
			if( neighbor.isLeftChild() ) // true if left child
				(neighbor.Up).Left = cut;
			else // true if right child
				(neighbor.Up).Right = cut;
			if( pendant.isLeftChild() )
				cut.Right = neighbor;
			else
				cut.Left = neighbor;
			neighbor.Up = cut;
			break;
		}
		if( pTree != null ) {
		//	System.err.println(pTree.toStringNoBranches());
			pTree.BalanceTree();
		}
		return pTree;
	}
	public void SnapInBranch(bNode current) {
		Branch[numBranches++] = current;
	}
	public void SnapInLeaf(bNode current, int ID) {
		Leaf[ID] = current;
		numLeaves++;
	}
	private void SproutLeaf(bNode current, bNode neighbor, double inStock, AbstractDistribution unif)
	{
		bNode newN = new bNode(neighbor.Up,current,neighbor,-1, inStock);
		// Determine if neighbor is a left or right child or root
		if( neighbor == Root )
		{
			neighbor.Stock = inStock;
			//neighbor.Stock2 = inStock;
			Root = neighbor.Up = current.Up = newN;
			Root.Stock = -9;
			Branch[numBranches] = newN;
		}
		else
		{
			Branch[numBranches]  = current.Up = newN;
			//boolean rooted = false;
			if( neighbor == (neighbor.Up).Left ) // true if left child
				(neighbor.Up).Left = newN;
			else // true if right child
				(neighbor.Up).Right = newN;
			neighbor.Up = newN;
		}
		// Randomly assign right/left for newN.
		if( unif.nextDouble() < 0.5 ) {
			bNode temp = newN.Left;
			newN.Left = newN.Right;
			newN.Right = temp;
		}
	}

	public void YankHardOnTree(bNode parent, bNode child, double connector){
		// This still needs to re-adjust stock lengths properly.
		bNode t1;
		bNode t2;
		bNode t3;
		if( child.Up == parent ) {
			child.Stock = connector;
			//child.Stock2 = connector2;
			// Then finished
			return;
		}
		double saveStock = child.Stock;
		//double saveStock2 = child.Stock2;
		t1 = child.Up;
		child.Up = parent;
		child.Stock = connector;
		//child.Stock2 = connector2;
		if( child.Right == parent )
			child.Right = t1;
		else
			child.Left = t1;
		if( child.Right != null ) {
			if( child.Right == t1 )
				YankHardOnTree(child,child.Right,saveStock);
			else
				YankHardOnTree(child, child.Right, child.Right.Stock);
		}
		if( child.Left != null ) {
			if( child.Left == t1 )
				YankHardOnTree(child,child.Left,saveStock);
			else
				YankHardOnTree(child, child.Left, child.Left.Stock);
		}
		//child.Stock = parent.Stock
		return;
	}

	public double logBranchesPrior(double hyper) {
		// Returns the log prior for all branches given hyper
		// each branch is distributed as exponential with mean = hyper
		//int count = 0;
		double total = 0;
		for(int i=0; i<numLeaves; i++) {
			total += Leaf[i].Stock;
		}
		for(int i=0; i<numBranches; i++) {
//			if( Branch[i] != Root ) {
				total += Branch[i].Stock;
//			}
		}
		//count--;  //don't over count the two child branches off the root.
		double rtn = 0;
		rtn -= (2*numLeaves - 3) * Math.log(hyper);
		rtn -= total / hyper;
		//System.err.println("branch prior: "+rtn);
		return rtn;
	}


	private bNode CopyDownTree(bNode inNode, bNode inParent) {
		bNode node = null;
		if( inNode.isBranch() ) { // This is a branching node
			node = new bNode(inParent,null,null,-1,inNode.Stock);
			SnapInBranch(node);
			node.Left = CopyDownTree(inNode.Left,node);
			node.Right = CopyDownTree(inNode.Right,node);
		} else { // This is a leaf
			node = new bNode(inParent,null,null,inNode.ID,inNode.Stock);
			Leaf[inNode.ID] = node;
			numLeaves++;
		}
		return node;
	}

	public bTree(bTree inT) {
		// Make a clone of the inT tree
		int n = inT.numLeaves;
		Leaf = new bNode[n];
		Branch = new bNode[n-1];
		numLeaves = 0;
		Root = new bNode(null,null,null,-1);
		Branch[0] = Root;
		numBranches = 1;
		Root.Left = CopyDownTree((inT.Root).Left,Root);
		Root.Right = CopyDownTree((inT.Root).Right,Root);
	}




	public bTree(int inLeaves, double mu, AbstractDistribution unif)
	{
		// Allocate necessary arrays for easy access to tree components.
		// Rooted trees have n - 1 internal nodes (including root)
		Leaf = new bNode[inLeaves];
		Branch = new bNode[inLeaves-1];
		numLeaves = numBranches = 0;
		int numB = 2*inLeaves-2;
		//double sum = 0;
		double[] branch = new double[numB];
		for(int i=0; i<numB; i++) {
			branch[i] = - mu * Math.log(unif.nextDouble());
			//sum += branch[i];
		}
		//sum /= numB;
		int b =0;
		Branch[numBranches] = Root = new bNode(null,null,null,-1); // The root needs a stock length in case we pop up higher than it.
		numBranches++;
		Leaf[numLeaves] = new bNode(Root,null,null,0,branch[b++]); numLeaves++; // base stock length == 1
		Leaf[numLeaves] = new bNode(Root,null,null,1,branch[b++]); numLeaves++;

		Root.Right = Leaf[1];
		Root.Left = Leaf[0];

		// Add each new leaf by:
		// 1. Adding as nearest-neighbor, the next node leaf, to a randomly selected node
		// Key. If we add nodes in numeric order (0,1,...,n) then we don't have to balance the tree and it is still randomly uniform
		while( numLeaves < inLeaves )
		{
			Leaf[numLeaves] = new bNode(null,null,null,numLeaves,branch[b++]);
			// There is equal probability of sprouting as nearest neighbor to
			// any existing leaf or branch (including popping up a new root).
			// Let n = # of existing nodes, then Prob = 1/n.
			int random = (int)(unif.nextDouble() * (numLeaves+numBranches));
			if( random < numLeaves )
				SproutLeaf(Leaf[numLeaves],Leaf[random],branch[b++], unif);

			else {
	/****/			//Sprout as nearest neighbor to a branch.
				SproutLeaf(Leaf[numLeaves],Branch[random-numLeaves],branch[b++], unif);
			}
			numBranches++;
			numLeaves++; // repeat through all leaves
		}
		// So that we are not biasing the root position.
		// Redraw them each from ibeta / 2;
		//gamma = new GammaDistribution(2);
		//Root.Right.Stock = ibeta*gamma.inverse(Math.random());
		//Root.Left.Stock = ibeta*gamma.inverse(Math.random());
		//for(int j=0; j<bnum; j++) {
		Root.Right.Stock = - (mu/2.0) * Math.log(unif.nextDouble());
		Root.Left.Stock = - (mu/2.0) * Math.log(unif.nextDouble());
		//for(int j=1; j<bnum; j++) {
		//    Root.Right.Stock[j] = Root.Right.Stock[0];
		//    Root.Left.Stock[j] = Root.Left.Stock[0];
		//}
		//}
		//Root.Right.Stock2 = - (mu/2.0) * Math.log(Math.random());
		//Root.Left.Stock2 = - (mu/2.0) * Math.log(Math.random());
		//Root.Right.Stock2 = Root.Right.Stock[0];
		//Root.Left.Stock2 = Root.Left.Stock[0];
	}
	public bTree(String inStr)
	{
		// First determine the size of the Leaf and Branch arrays necessary to hold the tree
		int nCount = 0;
		for(int i=0; i<inStr.length(); i++) {
			if( inStr.charAt(i) == '(')
				nCount++;
		}
		// We should also count ')' and make sure that this total equals the count of '('
		Leaf = new bNode[nCount+1];
		Branch = new bNode[nCount];
		numBranches = numLeaves = 0;
		//  recursively builds the whole tree
		Branch[numBranches++] = Root = new bNode(inStr,this,null,true);
		Root.Stock = -9;
	}

	public String toString() {
		return Root.toStringWithStocks();
	}

	public String toStringNoBranches() {
		return Root.toString();
	}

	public static void main(String a[]) {
		Tree t = new bTree("(0:0,(1:2,(2:2,3:2):2):2)");

		QMatrix q = new HKYMatrix();
		DNASequence[] mLeafSequence = DNASequence.ReadPhyllip("CAR.phy",false); // Read in data
		CountStatistic cs = new CountStatistic(mLeafSequence);
		int[] counts = cs.countsMatrix();
		int[][] data   = cs.dataMatrix();
		System.out.println(t.toString());
		System.out.println(
			t.LogLikelihood0(q, 1.0, counts, data,
			q.pi, true, 1.0 ,0)
			);
/*		System.out.println(
			((iTree)t).LogLikelihood1(q, 1.0, counts, data,
			q.pi, true, 1.0 ,0)
			);

		Tree[] te = t.EnumerateLastTaxon();
		for(int i=0; i<te.length; i++) {
			System.out.println(te[i].Root.toString());
		}
		System.exit(0); */
	}


}

