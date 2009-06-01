////////////////////////////////////////////////////////////////////////////////
//
//	RMG - Reaction Mechanism Generator
//
//	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
//	RMG Team (rmg_dev@mit.edu)
//
//	Permission is hereby granted, free of charge, to any person obtaining a
//	copy of this software and associated documentation files (the "Software"),
//	to deal in the Software without restriction, including without limitation
//	the rights to use, copy, modify, merge, publish, distribute, sublicense,
//	and/or sell copies of the Software, and to permit persons to whom the
//	Software is furnished to do so, subject to the following conditions:
//
//	The above copyright notice and this permission notice shall be included in
//	all copies or substantial portions of the Software.
//
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//	DEALINGS IN THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////



package jing.rxnSys;


import jing.chem.*;
import java.util.*;

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\ConversionTT.java                                                                  
//----------------------------------------------------------------------------

/**
Using conversion to control the termination of reaction systerm.  ie, reaction terminates at a given conversion of  some important reactants.
*/
//## class ConversionTT 
public class ConversionTT implements TerminationTester {
    
    protected LinkedList speciesGoalConversionSet;
    
    // Constructors
    
    //## operation ConversionTT() 
    private  ConversionTT() {
        {
            speciesGoalConversionSet=new LinkedList();
        }
        //#[ operation ConversionTT() 
        //#]
    }
    //## operation ConversionTT(HashSet) 
    public  ConversionTT(LinkedList p_scs) {
        {
            speciesGoalConversionSet=new LinkedList();
        }
        //#[ operation ConversionTT(HashSet) 
        speciesGoalConversionSet = p_scs;
        //#]
    }
    
    //## operation deleteSpeciesConversionSet(SpeciesConversion) 
    public void deleteSpeciesConversionSet(SpeciesConversion p_SpeciesConversion) {
        //#[ operation deleteSpeciesConversionSet(SpeciesConversion) 
        speciesGoalConversionSet.remove(p_SpeciesConversion);
        p_SpeciesConversion=null;
        //#]
    }
    
    //## operation isReactionTerminated(InitialStatus,PresentStatus) 
    public boolean isReactionTerminated(InitialStatus p_initialStatus, PresentStatus p_presentStatus) {
        //#[ operation isReactionTerminated(InitialStatus,PresentStatus) 
        Iterator iter = speciesGoalConversionSet.iterator();
        while (iter.hasNext()) {
        	SpeciesConversion sc = (SpeciesConversion)iter.next();
        	Species spe = sc.getSpecies();
        	SpeciesStatus iss = p_initialStatus.getSpeciesStatus(spe);
        	SpeciesStatus pss = p_presentStatus.getSpeciesStatus(spe);
        	double init_conc = iss.getConcentration();
        	double pres_conc = pss.getConcentration();
        	/* MRH on 1-Jun-2009
        		Updated from pres_conc>init_conc to pres_conc>1.01*init_conc
        			For certain cases, ODESovler would fail on the first step
        			The ODESolver would change the limiting reactant concentration
        			from 5.xxxxxxx3e-7 to 5.xxxxxxx5e-7 (for example), thus causing
					RMG to throw an exception.
			*/
        	if (pres_conc>(1.01*init_conc) || init_conc == 0) throw new InvalidConversionException();
        	double conversion = (init_conc-pres_conc)/init_conc;
        	if (conversion < sc.getConversion()) return false;
        }
        return true;
        //#]
    }
    
    public Iterator getSpeciesGoalConversionSet() {
        Iterator iter=speciesGoalConversionSet.iterator();
        return iter;
    }
   
    //11/1/07 gmagoon: created alternate accessor to return LinkedList; UPDATE: not needed
    //public LinkedList getSpeciesGoalConversionSetList(){
    //    return speciesGoalConversionSet;
    //}
    
    public void addSpeciesGoalConversionSet(SpeciesConversion p_SpeciesConversion) {
        speciesGoalConversionSet.add(p_SpeciesConversion);
    }
    
    public void removeSpeciesGoalConversionSet(SpeciesConversion p_SpeciesConversion) {
        speciesGoalConversionSet.remove(p_SpeciesConversion);
    }
    
    public void clearSpeciesGoalConversionSet() {
        speciesGoalConversionSet.clear();
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\ConversionTT.java
*********************************************************************/

