import java.util.*;
import java.io.*;
import corejava.*;


/**
* Topology update kernel
*/

public class TopUpdate{

    DCPSampler sampler;
    
    public TopUpdate (DCPSampler s){

	sampler = s;
    }
    

     
    /**
     * Sequentially updates tree topology at all partitions of the alignment
     * 
     */


    protected void Update(){

	double pLogLikelihood = 0, cLogLikelihood = 0, logRatio = 0;
	double[] pPartialLogLikelihood = new double[sampler.partition_list.size()];
	int pTree;


	/***************** start with the first partition ******************/
	
	Partition curr_part = (Partition)sampler.partition_list.elementAt(0);
	Partition temp_part = null;
	Partition second_part = curr_part;


	
	/************** find the second topology change point - (k-1) **************/
	
	int k = 1;
	while((!second_part.IsTopologyChange() || (k == 1)) && (k < sampler.partition_list.size())){
	    second_part = (Partition)sampler.partition_list.elementAt(k);
	    k++;
	}

	

	/************* if there is no topology change point ***************
	 ************* propose one tree from uniform distribution *********/ 


	if (!second_part.IsTopologyChange() || (k == 1)){
	    
	    pTree = sampler.ProposeOneTree(curr_part.cTree);

	    pLogLikelihood = 0;
	    cLogLikelihood = 0;
	    
	    for(int i = 0; i < sampler.partition_list.size(); i++){
		curr_part = (Partition)sampler.partition_list.elementAt(i);
		
		pPartialLogLikelihood[i] = sampler.PostTree[pTree].LogLikelihood0(curr_part.cMatrix, curr_part.cHyperParameter, curr_part.counts, curr_part.data, curr_part.cMatrix.pi, sampler.handleGaps, 0, 0);
	    
	    pLogLikelihood += pPartialLogLikelihood[i];
	    cLogLikelihood += curr_part.cPartialLogLikelihood;
	    }

	    logRatio = pLogLikelihood - cLogLikelihood; // uniform prior
	    sampler.tries[0]++;
	    
	    if( sampler.LogMHAccept( logRatio, sampler.set.unif.nextDouble() ) ) {
		
		sampler.acceptancerate[0]++;
		
		/******** update all previous partitions since 
		 ******** last topology change point *********/
		
		Partition update_part;
		for (int j = 0; j < sampler.partition_list.size(); j++){
		    update_part = (Partition)sampler.partition_list.elementAt(j);
		    update_part.cTree = pTree;
		    update_part.cPartialLogLikelihood = pPartialLogLikelihood[j];
		}
		
	    } 
	}

	/************ if there is at least one topology change point ********
	 ************ seguantially update trees along partitions with *******
	 ************ constant topology rejecting topologies equal to *******
	 ************ those of neighbor tree constant segments        ********/

	else {

	    pTree = sampler.ProposeEndTree(second_part.cTree, curr_part.cTree);

	
	    for (int i = 0; i < k - 1; i++){
		temp_part = (Partition)sampler.partition_list.elementAt(i);
		pPartialLogLikelihood[i] = sampler.PostTree[pTree].LogLikelihood0(temp_part.cMatrix, temp_part.cHyperParameter, temp_part.counts, temp_part.data, temp_part.cMatrix.pi, sampler.handleGaps, 0, 0);
		pLogLikelihood += pPartialLogLikelihood[i];
		cLogLikelihood += temp_part.cPartialLogLikelihood;
	    }
	    
	    logRatio = pLogLikelihood - cLogLikelihood; // uniform prior
	    
	    sampler.tries[0]++;
	    
	    if( sampler.LogMHAccept( logRatio, sampler.set.unif.nextDouble() ) ) {
		
		sampler.acceptancerate[0]++;
		
		/******** update all previous partitions since 
		 ******** last topology change point *********/
		
		Partition update_part;
		for (int j = 0; j < k - 1; j++){
		    update_part = (Partition)sampler.partition_list.elementAt(j);
		    update_part.cTree = pTree;
		    update_part.cPartialLogLikelihood = pPartialLogLikelihood[j];
		}
		
	    }
	    
	
	    int  prev_ch_pnt = k - 1, prev_prev_ch_pnt = 0;
	

	    for(int t = k; t < sampler.partition_list.size(); t++) {


		curr_part = (Partition)sampler.partition_list.elementAt(t);
	    
		/********* if this is a tolology change point do M-H step *********/
	    
		if (curr_part.IsTopologyChange()){

		    Partition prev_prev_part = (Partition)sampler.partition_list.elementAt(prev_prev_ch_pnt);
		    Partition prev_part = (Partition)sampler.partition_list.elementAt(prev_ch_pnt);

		    pTree = sampler.ProposeMiddleTree(curr_part.cTree, prev_part.cTree, prev_prev_part.cTree);
		    if (pTree != -1){

			pLogLikelihood = 0;
			cLogLikelihood = 0;

			for (int i = prev_ch_pnt; i < t; i++){
			    temp_part = (Partition)sampler.partition_list.elementAt(i);

			    pPartialLogLikelihood[i] = sampler.PostTree[pTree].LogLikelihood0(temp_part.cMatrix, temp_part.cHyperParameter, temp_part.counts, temp_part.data, temp_part.cMatrix.pi, sampler.handleGaps, 0, 0);
			    
			    pLogLikelihood += pPartialLogLikelihood[i];
			    cLogLikelihood += temp_part.cPartialLogLikelihood;
			}
			
			logRatio = pLogLikelihood - cLogLikelihood; // uniform prior
			
			sampler.tries[0]++;
			
			if( sampler.LogMHAccept( logRatio, sampler.set.unif.nextDouble() ) ) {
			    
			    sampler.acceptancerate[0]++;
			    
			    /******** update all previous partitions since 
			     ******** last topology change point *********/
			    
			    Partition update_part;
			    for (int j = prev_ch_pnt; j < t; j++){
				update_part = (Partition)sampler.partition_list.elementAt(j);
				update_part.cTree = pTree;
				update_part.cPartialLogLikelihood = pPartialLogLikelihood[j];
			    }
			}
		    }

		    /******* Start new chain of equal topology partitions ****/
		    
		    prev_prev_ch_pnt = prev_ch_pnt;
		    prev_ch_pnt = t;
		}

	
	    }
    
	    /******** Finish update on the right end *****************/

	    Partition prev_part = (Partition)sampler.partition_list.elementAt(prev_ch_pnt);
	    Partition prev_prev_part = (Partition)sampler.partition_list.elementAt(prev_prev_ch_pnt);
	     
	    pTree = sampler.ProposeEndTree(prev_part.cTree, prev_prev_part.cTree);

	    
	    
	    pLogLikelihood = 0;
	    cLogLikelihood = 0;
	    
	    for (int i = prev_ch_pnt; i < sampler.partition_list.size(); i++){
		temp_part = (Partition)sampler.partition_list.elementAt(i);
		pPartialLogLikelihood[i] = sampler.PostTree[pTree].LogLikelihood0(temp_part.cMatrix, temp_part.cHyperParameter, temp_part.counts, temp_part.data, temp_part.cMatrix.pi, sampler.handleGaps, 0, 0);
		pLogLikelihood += pPartialLogLikelihood[i];
		cLogLikelihood += temp_part.cPartialLogLikelihood;
	    }
	    
	    
	    logRatio = pLogLikelihood - cLogLikelihood; // uniform prior
	    
	    sampler.tries[0]++;
	    
	    if( sampler.LogMHAccept( logRatio, sampler.set.unif.nextDouble() ) ) {
		
		sampler.acceptancerate[0]++;
		
		/******** update all previous partitions since 
		 ******** last topology change point *********/
		
		Partition update_part;
		for (int j = prev_ch_pnt; j < sampler.partition_list.size(); j++){
		    update_part = (Partition)sampler.partition_list.elementAt(j);
		    update_part.cTree = pTree;
		    update_part.cPartialLogLikelihood = pPartialLogLikelihood[j];
		}
		
	    }
	    
	    
	}
	
    }
    
    
}
