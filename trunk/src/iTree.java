import java.util.*;
import cern.jet.random.*;

public class iTree extends Tree {

	public Tree JointBranchAndTopology(AbstractDistribution nglobal, AbstractDistribution nlocal, AbstractDistribution unif,
						double pmix, boolean fixRoot, int bnum,Boolean update) {
		System.err.println("Joint branch proposal not implemented for integrated trees.");
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
		Node cNode = null;
		double stock;
		for(int i=0; i<numSites; i++) {
			if( counts[i] != 0 ) {
			       for(int j=0; j<numLeaves; j++) { 	// Here, let's calculate the likelihood of a single word
					cNode = Leaf[j]; // Start at the bottom of each leaf
					if( data[i][cNode.ID] >= 0 ) {
						if( j == 0 )
							stock = 0;
						else
							stock = mean;
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
					while( (cNode != Root) && cNode.getBrother().isLikelihoodDone ) { // Propogate up tree until we hit the root or an unfinished section
						cNode = cNode.Up;
						Node lChild = cNode.Left;
						Node rChild = cNode.Right;
						if (cNode != Root) {
							stock = mean;
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
						} else {
							cNode.ConditionalLikelihood[0] = lChild.ConditionalLikelihood[0] * rChild.ConditionalLikelihood[0];
							cNode.ConditionalLikelihood[1] = lChild.ConditionalLikelihood[1] * rChild.ConditionalLikelihood[1];
							cNode.ConditionalLikelihood[2] = lChild.ConditionalLikelihood[2] * rChild.ConditionalLikelihood[2];
							cNode.ConditionalLikelihood[3] = lChild.ConditionalLikelihood[3] * rChild.ConditionalLikelihood[3];
						}
						cNode.isLikelihoodDone = true;
					}
			       }
			       // By this point all conditional likelihoods are calculated and cNode = Root. The conditionals
			       // for the Root are still in tCL
			       // Posterior likelihood = SUM( conditional * PI_i ) of all i in span.
			       double L = 0;
			       L += pi[0]*cNode.ConditionalLikelihood[0];
			       L += pi[1]*cNode.ConditionalLikelihood[1];
			       L += pi[2]*cNode.ConditionalLikelihood[2];
			       L += pi[3]*cNode.ConditionalLikelihood[3];
			       LL += Math.log(L) * counts[i];
			}
		}
		return LL;
	}

	public double LogLikelihood1(QMatrix Q, double mean, int[] counts, int[][] data, double[] pi, boolean handleGaps, double K, int bnum) {
		// Returns the log likelihood of the data in C, given the topology of this tree, Q and
		// a proposed distribution of bases at the root by summing over all possible states
		// for the internal nodes.

		int numSites = counts.length;
		double LL = 0;
		boolean[] repeats = new boolean[numLeaves]; // Automatically fills with false
		double[] tCL = new double[4]; // Automatically zeros all entries. for temporarily storage.
		Node cNode = null;
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
					//	if( j == 0 )
					//		stock = 0;
					//	else
					//	stock = mean;
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
					//System.err.println("Trying: "+cNode.ID);
					//System.err.println(cNode.toString());
					//System.err.println(cNode.getBrother().toString());
					//if( cNode.getBrother().isLikelihoodDone )
					//	System.err.println("true");
					//else
					//	System.err.println("false");
					while( cNode.getBrother().isLikelihoodDone && (cNode.Up.Up != Root) ) { // Propogate up tree until we hit the root or an unfinished section
						//System.err.println("c");
						cNode = cNode.Up;
						Node lChild = cNode.Left;
						Node rChild = cNode.Right;
						//if (cNode != Root) {
							//stock = mean;
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
					//	} else {
					//		cNode.ConditionalLikelihood[0] = lChild.ConditionalLikelihood[0] * rChild.ConditionalLikelihood[0]
					//			* Q.Pt(0,data[i][Leaf[0].ID],stock);
					//		cNode.ConditionalLikelihood[1] = lChild.ConditionalLikelihood[1] * rChild.ConditionalLikelihood[1]
					//			* Q.Pt(1,data[i][Leaf[0].ID],stock);
					//		cNode.ConditionalLikelihood[2] = lChild.ConditionalLikelihood[2] * rChild.ConditionalLikelihood[2]
					//			* Q.Pt(2,data[i][Leaf[0].ID],stock);
					//		cNode.ConditionalLikelihood[3] = lChild.ConditionalLikelihood[3] * rChild.ConditionalLikelihood[3]
					//			* Q.Pt(3,data[i][Leaf[0].ID],stock);
					//	}
						cNode.isLikelihoodDone = true;
					}
			       }
			       cNode = cNode.Up;
			       Node lChild = cNode.Left;
			       Node rChild = cNode.Right;
			       Node fChild = Leaf[0];
			       double L;
			       L  = lChild.ConditionalLikelihood[0] * rChild.ConditionalLikelihood[0] * fChild.ConditionalLikelihood[0] * pi[0];
			       L += lChild.ConditionalLikelihood[1] * rChild.ConditionalLikelihood[1] * fChild.ConditionalLikelihood[1] * pi[1];
			       L += lChild.ConditionalLikelihood[2] * rChild.ConditionalLikelihood[2] * fChild.ConditionalLikelihood[2] * pi[2];
			       L += lChild.ConditionalLikelihood[3] * rChild.ConditionalLikelihood[3] * fChild.ConditionalLikelihood[3] * pi[3];

			       // By this point all conditional likelihoods are calculated and cNode = Root. The conditionals
			       // for the Root are still in tCL
			       // Posterior likelihood = SUM( conditional * PI_i ) of all i in span.
			    //   double L = 0;
			    //  L += pi[0]*cNode.ConditionalLikelihood[0];
			    //   L += pi[1]*cNode.ConditionalLikelihood[1];
			    //   L += pi[2]*cNode.ConditionalLikelihood[2];
			    //   L += pi[3]*cNode.ConditionalLikelihood[3];
			       LL += Math.log(L) * counts[i];
			}
		}
		return LL;
	}


	public Tree ProposedJump(int method, boolean fix, AbstractDistribution unif) {
		iTree pTree = null;
		switch( method ) {
		case 0:
		// Our first jump proposal is just another random tree.
			//pTree = new Tree(this.numLeaves);
			System.exit(-1);
			//; // symmetric jump
			break;
		case 5: // Permute children on internal node, allows re-rooting.
			pTree = new iTree(this);
			// Select the parent node, can be any of n - 2 internal nodes or the root
			int r = (int)(unif.nextDouble()*numBranches);
			Node parent = pTree.Branch[r];
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
			pTree = new iTree(this);
			// Select the parent node, can be any of n - 2 internal nodes
			//int r = (int)(unif.nextDouble()*numBranches);
			//Node parent;
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
				prandom = (int)(unif.nextDouble() * numLeaves);
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
				double[] oldStock = pneighbor.Stock;
				pneighbor.Stock = new double[oldStock.length];
				for(int j=0; j<oldStock.length; j++) {
				    pneighbor.Stock[j] = oldStock[j] + cut.Stock[j];
				}
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
						prandom = (int) (unif.nextDouble() * 3);
						System.err.println("CASE A "+prandom);
					} else {
						prandom = 2;
						System.err.println("CASE C "+prandom);
					}
				} else {
					if( pneighbor.isBranch() ) {
						prandom = (int) (unif.nextDouble() * 4);
						System.err.println("CASE A "+prandom);
					} else {
						prandom = 2 + (int) (unif.nextDouble() * 2);
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
				double rand = unif.nextDouble();
				int span = neighbor.Stock.length;
				double[] oldNeighborStock = neighbor.Stock;
				neighbor.Stock = new double[span];
				cut.Stock = new double[span];
				for(int j=0; j<span; j++) {
				    //double split = neighbor.Stock * unif.nextDouble();
				    cut.Stock[j] = oldNeighborStock[j] * rand;
				    neighbor.Stock[j] = oldNeighborStock[j] * (1-rand);
				}
 				cut.Up = neighbor.Up;
				if( neighbor.isLeftChild() ) // true if left child
					(neighbor.Up).Left = cut;
				else // true if right child
					(neighbor.Up).Right = cut;
				// Add a little randomness into the local branch lengths.
				//NormalDistribution n = new NormalDistribution(0,0.1);
				//cut.Stock += n.inverse(unif.nextDouble());
				//neighbor.Stock += n.inverse(unif.nextDouble());
				//pendant.Stock += n.inverse(unif.nextDouble());
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
			pTree = new iTree(this);
			//System.out.println(pTree.Root.toString());
			Node pendant = null;
			// 2. choose a pendant leaf at random
			//prandom = (int)(unif.nextDouble() * numLeaves);
			//pendant = pTree.Leaf[prandom];
			//if( pendant.Up != pTree.Root ) {
			//do {
			//	pendant = pTree.Leaf[(int)(unif.nextDouble()*numLeaves)];
			//} while (pendant.Up == pTree.Root);
				pendant = pTree.Leaf[(int)(unif.nextDouble()*(numLeaves-1))+1];
				//penddant.Up.ForceDownReprocess();
				// No need to re-root.
				// 2b. clip out extra node.
				Node pneighbor = pendant.getBrother();
				Node cut = pendant.Up; // the branching node to be removed
				if( cut.isLeftChild() )  // cut is a left child
					(cut.Up).Left = pneighbor;
				else // cut is a right child
					(cut.Up).Right = pneighbor;
				pneighbor.Up = cut.Up;
				//int span = pneighbor.Stock.length;
				//double[] oldPneighborStock = pneighbor.Stock;
				//pneighbor.Stock = new double[span];
				//for(int j=0; j<span; j++)
				//    pneighbor.Stock[j] = oldPneighborStock[j] + cut.Stock[j]; // Conserve stocks.
				//oldLength = pneighbor.Stock;
				//System.err.println("old length = "+oldLength);
				// 3. choose a potential branch to which to reattach the pendant/
				// the pendant can re-attach as nearest-neighbor to the remaining n - 1
				// leaves or n - 3 internal nodes (does not include root).  This leaves
				// 2n - 4 locations.
				Node neighbor;
				int random = (int)(unif.nextDouble() * (2*numLeaves - 4));
				if( random < (numLeaves - 1) ) { // next to a leaf, but skip pendant
					if( random < pendant.ID )
						neighbor = pTree.Leaf[random];
					else
						neighbor = pTree.Leaf[random+1];
					//System.out.println("new nearest neighbor is Leaf "+neighbor.ID);
					if( fix && (neighbor.Up == pTree.Root) ) {
						//System.err.println("Pendant Neighbor's parent is the root.");
						//System.exit(-1);
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
				if( neighbor == pTree.Root ) {
					System.err.println("I should never get here!");
					System.exit(-1);
				}
				// Let's now re-attach our pendant to it's neighbor.
				// first calculate proposal density
				//System.err.println("new length = "+neighbor.Stock);
				//QQ = new Double( neighbor.Stock / oldLength );
				//double[] oldNeighborStock = neighbor.Stock;
				//int span = neighbor.Stock.length;
				//double rand = unif.nextDouble();
				//cut.Stock = new double[span];
				//neighbor.Stock = new double[span];
				//pTree.QQ = 1;
				//pTree.RR = oldNeighborStock[0];
				//for(int j=0; j<span; j++) {
				//    pTree.QQ *= neighbor.Stock[j] / pneighbor.Stock[j];
				//    double split = oldNeighborStock[j] * rand;
				//    cut.Stock[j] = split;
				//    neighbor.Stock[j] = oldNeighborStock[j] - split;
				//}
 				cut.Up = neighbor.Up;
				if( neighbor.isLeftChild() ) // true if left child
					(neighbor.Up).Left = cut;
				else // true if right child
					(neighbor.Up).Right = cut;
				// Add a little randomness into the local branch lengths.
				//NormalDistribution n = new NormalDistribution(0,0.1);
				//cut.Stock += n.inverse(unif.nextDouble());
				//neighbor.Stock += n.inverse(unif.nextDouble());
				//pendant.Stock += n.inverse(unif.nextDouble());
				//if( (cut.Stock <= 0) || (neighbor.Stock <= 0) || (pendant.Stock <= 0) )
				//	return null; // Bad jump

				if( pendant.isLeftChild() )
					cut.Right = neighbor;
				else
					cut.Left = neighbor;
				neighbor.Up = cut;
				//cut.ForceDownReprocess();
				//break;
			//}  // That's the end of pendant branch...
			break;
			//else { // Let's now deal with re-rooting.
	//	case 7: break;
		case 7: // cluster re-arrangement
			// 1. clone tree.
			int prandom;
			if( numBranches <= 3 )
				return null;
			pTree = new iTree(this);
			// 2. choose a pendant internal node at random
			//    exclude root and direct child of root
			do {
				prandom = (int)(unif.nextDouble() * numBranches);
				pendant = pTree.Branch[prandom];
			} while( (pendant == pTree.Root) || (pendant.Up == pTree.Root) );
			//break;

			// now cut this internal node out of the tree.
				pneighbor = pendant.getBrother();
				cut = pendant.Up; // the branching node to be removed
				if( cut.isLeftChild() )  // cut is a left child
					(cut.Up).Left = pneighbor;
				else // cut is a right child
					(cut.Up).Right = pneighbor;
				pneighbor.Up = cut.Up;

		//	break; /*
				//double[] oldPNeighborStock = pneighbor.Stock;
				//pneighbor.Stock = new double[oldPNeighborStock.length];
				//for(int j=0; j<oldPNeighborStock.length; j++ )
				//    pneighbor.Stock[j] = oldPNeighborStock[j] + cut.Stock[j]; // Conserve stocks.
				//oldLength = pneighbor.Stock;
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
				//oldNeighborStock = neighbor.Stock;
				//span = neighbor.Stock.length;
				//neighbor.Stock = new double[span];
				//cut.Stock = new double[span];
				//pTree.QQ = 1;
				//rand = unif.nextDouble();
				//for(int j=0; j<span; j++) {
				    // Let's now re-attach our pendant to it's neighbor.
				//    double split = oldNeighborStock[j] * rand;
				//QQ = new Double( neighbor.Stock / oldLength );
				//    pTree.QQ *= oldNeighborStock[j] / oldPNeighborStock[j];
				//    cut.Stock[j] = split;
				//    neighbor.Stock[j] = oldNeighborStock[j]- split;
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
				return pTree;

		}
		return pTree;
	}

	public void YankHardOnTree(Node parent, Node child){
		Node t1;
		Node t2;
		Node t3;
		if( child.Up == parent ) {
			return;
		}
		t1 = child.Up;
		child.Up = parent;
		if( child.Right == parent )
			child.Right = t1;
		else
			child.Left = t1;
		if( child.Right != null ) {
			if( child.Right == t1 )
				YankHardOnTree(child,child.Right);
			else
				YankHardOnTree(child, child.Right);
		}
		if( child.Left != null ) {
			if( child.Left == t1 )
				YankHardOnTree(child,child.Left);
			else
				YankHardOnTree(child, child.Left);
		}
		return;
	}

	public Tree[] EnumerateLastTaxon() {
		BalanceTree();
		Tree[] rtn = new iTree[2*numLeaves-5];
		for(int i=0; i<rtn.length; i++)
			rtn[i] =  new iTree(this.Root.toString());
		// note last taxon is always the right child on a balanced tree
		// Nearest neighbor to taxa 0
		rtn[0] = new iTree(this.Root.toString());
		Node lt = rtn[0].Leaf[numLeaves-1];
		Node p = lt.Up;
		Node sib = p.Left;
		if( p.isLeftChild() )
			p.Up.Left = sib;
		else
			p.Up.Right = sib;
		sib.Up = p.Up;
		p.Up = rtn[0].Root;
		p.Left = rtn[0].Root.Right;
		p.Left.Up = p;
		rtn[0].Root.Right = p;
		rtn[0].BalanceTree();
		rtn[0] = new iTree(rtn[0].Root.toString());
		// Nearest neighbor to the remaining taxa
		for(int i=1; i<numLeaves-1; i++) {
			lt = rtn[i].Leaf[numLeaves-1];
			p = lt.Up;
			sib = p.Left;
			if( p.isLeftChild() )
				p.Up.Left = sib;
			else
				p.Up.Right = sib;
			sib.Up = p.Up;

			Node neighbor = rtn[i].Leaf[i];
			Node nup = neighbor.Up;
			if( neighbor.isLeftChild() )
				nup.Left = p;
			else
				nup.Right = p;
			p.Left = neighbor;
			neighbor.Up = p;
			rtn[i] = new iTree(rtn[i].Root.toString());
			rtn[i].BalanceTree();
		}
		// Nearest neighbor to remaining internal branches
		int c = numLeaves - 1;
		for(int i=0; i<numBranches; i++) {
			if( (Branch[i] == Root) ||
			    (Branch[i].Up == Root) ||
			    (Branch[i].Up.Right == Leaf[numLeaves-1]) ||
			    (Branch[i].Right == Leaf[numLeaves-1]) )
			    ;
			else {
				lt = rtn[c].Leaf[numLeaves-1];
				p = lt.Up;
				sib = p.Left;
				if( p.isLeftChild() )
					p.Up.Left = sib;
				else
					p.Up.Right = sib;
				sib.Up = p.Up;
				Node bottom = rtn[c].Branch[i];
				Node top = bottom.Up;
				if( bottom.isLeftChild() )
					top.Left = p;
				else
					top.Right = p;
				p.Up = top;
				p.Left = bottom;
				bottom.Up = p;
				rtn[c] = new iTree(rtn[c].Root.toString());
				rtn[c].BalanceTree();
				c++;
			}
		}
		return rtn;
	}

	public iTree(iTree inT) {
		// Make a clone of the inT tree
		int n = inT.numLeaves;
		Leaf = new Node[n];
		Branch = new Node[n-1];
		numLeaves = 0;
		Root = new Node(null,null,null,-1,inT.Root.uID,inT.Root.isLikelihoodDone);
		Branch[0] = Root;
		numBranches = 1;
		Root.Left = CopyDownTree((inT.Root).Left,Root);
		Root.Right = CopyDownTree((inT.Root).Right,Root);
	}

	public iTree(int inLeaves, AbstractDistribution unif) {
		// Allocate necessary arrays for easy access to tree components.
		// Rooted trees have n - 1 internal nodes (including root)
		Leaf = new Node[inLeaves];
		Branch = new Node[inLeaves-1];
		numLeaves = numBranches = 0;
		Branch[numBranches] = Root = new Node(null,null,null,-1,inLeaves+numBranches); // The root needs a stock length in case we pop up higher than it.
		numBranches++;
		Leaf[numLeaves] = new Node(Root,null,null,0,numLeaves++); // base stock length == 1
		Leaf[numLeaves] = new Node(Root,null,null,1,numLeaves++);

		Root.Right = Leaf[1];
		Root.Left = Leaf[0];

		// Add each new leaf by:
		// 1. Adding as nearest-neighbor, the next node leaf, to a randomly selected node
		// Key. If we add nodes in numeric order (0,1,...,n) then we don't have to balance the tree and it is still randomly uniform
		while( numLeaves < inLeaves )
		{
			Leaf[numLeaves] = new Node(null,null,null,numLeaves,numLeaves);
			// There is equal probability of sprouting as nearest neighbor to
			// any existing leaf or branch (including popping up a new root).
			// Let n = # of existing nodes, then Prob = 1/n.
			// This has been modified to return unrooted trees ONLY
			int random = (int)(unif.nextDouble() * (numLeaves+numBranches-2));
			if( random < (numLeaves-1) )
				SproutLeaf(Leaf[numLeaves],Leaf[random+1],inLeaves+numBranches, unif);

			else {
	/****/			//Sprout as nearest neighbor to a branch.
				SproutLeaf(Leaf[numLeaves],Branch[random-numLeaves+2],inLeaves+numBranches, unif);
			}
			numBranches++;
			numLeaves++; // repeat through all leaves
		}
	}

	public iTree(String inStr) {
		// First determine the size of the Leaf and Branch arrays necessary to hold the tree
		int nCount = 0;
		for(int i=0; i<inStr.length(); i++) {
			if( inStr.charAt(i) == '(')
				nCount++;
		}
		// We should also count ')' and make sure that this total equals the count of '('
		Leaf = new Node[nCount+1];
		Branch = new Node[nCount];
		numBranches = numLeaves = 0;
		//  recursively builds the whole tree
		Branch[numBranches++] = Root = new Node(inStr,this,null,true);
	}

	public String toString() {
		return Root.toString();
	}

	public String toStringNoBranches() {
		return Root.toString();
	}

	public double SumOfBranchLengths() {
		System.err.println("No branches lengths are defined.");
		System.exit(-1);
		return 0;
	}

	public static void main(String a[]) {
		Tree t = new iTree("(0,(1,(2,3)))");
		QMatrix q = new HKYMatrix();
		DNASequence[] mLeafSequence = DNASequence.ReadPhyllip("CAR.phy",false); // Read in data
		CountStatistic cs = new CountStatistic(mLeafSequence);
		int[] counts = cs.countsMatrix();
		int[][] data   = cs.dataMatrix();
		/* **** */
		System.out.println(
			t.LogLikelihood0(q, 2.0, counts, data,
			q.pi, true, 1.0 ,0)
			);
		System.out.println(
			((iTree)t).LogLikelihood1(q, 1.0, counts, data,
			q.pi, true, 1.0 ,0)
			);

		Tree[] te = t.EnumerateLastTaxon();
		for(int i=0; i<te.length; i++) {
			System.out.println(te[i].Root.toString());
		}
		System.exit(0);
	}
}

