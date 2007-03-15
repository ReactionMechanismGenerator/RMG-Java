//!********************************************************************************
//!
//!    RMG: Reaction Mechanism Generator                                            
//!
//!    Copyright: Jing Song, MIT, 2002, all rights reserved
//!     
//!    Author's Contact: jingsong@mit.edu
//!
//!    Restrictions:
//!    (1) RMG is only for non-commercial distribution; commercial usage
//!        must require other written permission.
//!    (2) Redistributions of RMG must retain the above copyright
//!        notice, this list of conditions and the following disclaimer.
//!    (3) The end-user documentation included with the redistribution,
//!        if any, must include the following acknowledgment:
//!        "This product includes software RMG developed by Jing Song, MIT."
//!        Alternately, this acknowledgment may appear in the software itself,
//!        if and wherever such third-party acknowledgments normally appear.
//!  
//!    RMG IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED 
//!    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
//!    OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
//!    DISCLAIMED.  IN NO EVENT SHALL JING SONG BE LIABLE FOR  
//!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
//!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
//!    OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;  
//!    OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  
//!    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
//!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
//!    THE USE OF RMG, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//! 
//!******************************************************************************



package jing.rxnSys;


import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import jing.param.*;
import jing.rxn.NegativeRateException;
import jing.rxn.Reaction;
import jing.rxn.TemplateReaction;
import jing.chem.ChemGraph;
import jing.chem.Species;
import jing.chem.SpeciesDictionary;
//import RMG;
//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\RateBasedRME.java                                                                  
//----------------------------------------------------------------------------

//## class RateBasedRME 
public class RateBasedRME implements ReactionModelEnlarger {
    
	protected HashSet includeSpecies = null; //these species are included into the core even if they have very 
	  										 //low flux.
    // Constructors
    
    public  RateBasedRME() {
    }
    
	public void addIncludeSpecies(HashSet p_includeSpecies){
		if (includeSpecies == null) {
			includeSpecies = p_includeSpecies;
		}
		else {
			System.out.println("IncludeSpecies have already been added!!");
			System.exit(0);
		}
			
	}
	
    //## operation enlargeReactionModel(ReactionSystem) 
    public void enlargeReactionModel(ReactionSystem p_reactionSystem) {
        //#[ operation enlargeReactionModel(ReactionSystem) 
        ReactionModel rm = p_reactionSystem.getReactionModel();
        if (!(rm instanceof CoreEdgeReactionModel)) throw new InvalidReactionModelTypeException();
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)rm;
        
        PresentStatus ps = p_reactionSystem.getPresentStatus();
        String maxflux="";
		double startTime = System.currentTimeMillis();
        Species next = getNextCandidateSpecies(cerm,ps,maxflux);
		
		double findSpeciesTime = (System.currentTimeMillis()-startTime)/1000/60;
		
		Global.diagnosticInfo.append(next.getChemkinName()+"\t" + maxflux+"\t"+ ((RateBasedVT)p_reactionSystem.finishController.validityTester).Rmin+ "\t" + findSpeciesTime +"\t");
		
		
		
        System.out.print("\nAdd a new reacted Species:");
        System.out.println(next.getName());
        Temperature temp = new Temperature(298, "K");
        double H = next.calculateH(temp);
        double S = next.calculateS(temp);
        double G = next.calculateG(temp);
        double Cp = next.calculateCp(temp);
        System.out.println("Thermo\t" + String.valueOf(H) + '\t' + String.valueOf(S)+ '\t' + String.valueOf(G)+ '\t' + String.valueOf(Cp));
        
		
		startTime = System.currentTimeMillis();
		
        if (cerm.containsAsReactedSpecies(next)) 
        	throw new InvalidNextCandidateSpeciesException();
        else {
        	cerm.moveFromUnreactedToReactedSpecies(next);
        	cerm.moveFromUnreactedToReactedReaction();
        	
        }
		Global.moveUnreactedToReacted = (System.currentTimeMillis()-startTime)/1000/60; 
		
        // generate new reaction set
		startTime = System.currentTimeMillis();
		LinkedHashSet newReactionSet = p_reactionSystem.getReactionGenerator().react(cerm.getReactedSpeciesSet(),next);
		
		double enlargeTime = (System.currentTimeMillis()-startTime)/1000/60;
		
		
		
		startTime = System.currentTimeMillis();
		StringBuilder restartFileContent= new StringBuilder();
		try{
			File allReactions = new File ("Restart/allReactions.txt");
			FileWriter fw = new FileWriter(allReactions, true);
			//Species species = (Species) iter.next();
			for(Iterator iter=newReactionSet.iterator();iter.hasNext();){
				
				Reaction reaction = (Reaction) iter.next();
				if (cerm.categorizeReaction(reaction)==-1)
					restartFileContent.append(reaction.toRestartString(Global.temperature) + "\n");
								
			}
			
			//restartFileContent += "\nEND";
			fw.write(restartFileContent.toString());
			fw.close();
		}
		catch (IOException e){
			System.out.println("Could not write the added Reactions to the allReactions file");
        	System.exit(0);
		}
		
		double restartTime = (System.currentTimeMillis()-startTime)/1000/60;
        
		Global.diagnosticInfo.append(Global.moveUnreactedToReacted + "\t" +enlargeTime+"\t" + restartTime +"\t");
		
        // partition the reaction set into reacted reaction set and unreacted reaction set
        // update the corresponding core and edge model of CoreEdgeReactionModel
        cerm.addReactionSet(newReactionSet);
        //String return_string;
        return;
        
        
        //#]
    }
    
	public boolean presentInIncludedSpecies(Species p_species){
		Iterator iter = includeSpecies.iterator();
		while (iter.hasNext()){
			Species spe = (Species)iter.next();
			Iterator isomers = spe.getResonanceIsomers();
			while (isomers.hasNext()){
				ChemGraph cg = (ChemGraph)isomers.next();
				if (cg.equals(p_species.getChemGraph())) 
					return true;
			}
		}
		return false;
	}
	
    //## operation getNextCandidateSpecies(CoreEdgeReactionModel,PresentStatus) 
    public Species getNextCandidateSpecies(CoreEdgeReactionModel p_reactionModel, PresentStatus p_presentStatus, String maxflux) {
        //#[ operation getNextCandidateSpecies(CoreEdgeReactionModel,PresentStatus) 
    	LinkedHashSet unreactedSpecies = p_reactionModel.getUnreactedSpeciesSet();
         
        Species maxSpecies = null;
        double maxFlux = 0;
        Species maxIncludedSpecies = null;
		double maxIncludedFlux = 0;
		
        Iterator iter = unreactedSpecies.iterator();
        while (iter.hasNext()) {
        	Species us = (Species)iter.next();
        	//double thisFlux = Math.abs(p_presentStatus.getSpeciesStatus(us).getFlux());
			//System.out.println(p_presentStatus.unreactedSpeciesFlux[83]);
			//System.exit(0);
			double thisFlux = Math.abs(p_presentStatus.unreactedSpeciesFlux[us.getID()]);
			if (includeSpecies != null && includeSpecies.contains(us)) {
				if (thisFlux > maxIncludedFlux) {
	        		maxIncludedFlux = thisFlux;
	        		maxIncludedSpecies = us;
	        	}
			}
			else {
				if (thisFlux > maxFlux) {
	        		maxFlux = thisFlux;
	        		maxSpecies = us;
	        	}
			}
        	
        }
		
		if (maxIncludedSpecies != null){
			System.out.println("Instead of "+maxSpecies.toChemkinString()+" with flux "+ maxFlux + " "+ maxIncludedSpecies.toChemkinString() +" with flux " + maxIncludedFlux);
			maxFlux = maxIncludedFlux;
			maxSpecies = maxIncludedSpecies;
			includeSpecies.remove(maxIncludedSpecies);
		}
		
        maxflux = ""+maxFlux;
        if (maxSpecies == null) throw new NullPointerException();
        
		
		
        LinkedHashSet ur = p_reactionModel.getUnreactedReactionSet();
        
		HashMap significantReactions = new HashMap();
		int reactionWithSpecies = 0;
		
        for (Iterator iur = ur.iterator(); iur.hasNext();) {
        	Reaction r = (Reaction)iur.next();
        	double flux = 0;
			Temperature p_temperature = p_presentStatus.temperature;
			if (r.contains(maxSpecies)){
				reactionWithSpecies++;
				if (r instanceof TemplateReaction) {
	        		flux = ((TemplateReaction)r).calculateTotalPDepRate(p_temperature);
	        	}
	        	else {
	        	 	flux = r.calculateTotalRate(p_temperature);
	        	}
	        	if (flux > 0) {
	        		for (Iterator rIter=r.getReactants(); rIter.hasNext();) {
						Species spe = (Species)rIter.next();
	        		    
	        		    double conc = (p_presentStatus.getSpeciesStatus(spe)).getConcentration();
	        			if (conc<0)
	        				throw new NegativeConcentrationException(spe.getName() + ": " + String.valueOf(conc));
	        		    flux *= conc;
						
	        		}
					

	        	}
	        	else {
	        		throw new NegativeRateException(r.toChemkinString(Global.temperature) + ": " + String.valueOf(flux));
	        	}
				if (flux > 0.01 * maxFlux)
					significantReactions.put(r,new Double(flux));
			}
        	
        }
		
        System.out.print("Time: ");
        System.out.println(p_presentStatus.getTime());
        System.out.println("unreacted Spe with highest flux: " + String.valueOf(maxFlux));
		System.out.println("The total number of unreacted reactions with this species is "+reactionWithSpecies);
		Iterator reactionIter = significantReactions.keySet().iterator();
		while (reactionIter.hasNext()){
			Reaction r = (Reaction)reactionIter.next();
			System.out.println("1.\t"+r.getStructure().toChemkinString(r.hasReverseReaction())+"\t"+significantReactions.get(r));
		}
        
        return maxSpecies;
        //#]
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\RateBasedRME.java
*********************************************************************/

