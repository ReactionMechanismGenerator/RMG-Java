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


import java.util.*;

import jing.chem.Species;

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\RateBasedVT.java                                                                  
//----------------------------------------------------------------------------

/**
In rate-based model generation process, the criteria is:
Rj<e*Rchar
in which, Rj is the unreacted species flux, e is the user tolerance, and Rchar is the system characteristic flux calculated by L2 norm of reacted species fluxes.
*/
//## class RateBasedVT 
public class RateBasedVT implements ValidityTester {
    
    protected double Rmin = 0;		//## attribute Rmin 
    
    protected double tolerance;		//## attribute tolerance 
    
    
    // Constructors
    
    //## operation RateBasedVT() 
    protected  RateBasedVT() {
        //#[ operation RateBasedVT() 
        //#]
    }
    //## operation RateBasedVT(double) 
    public  RateBasedVT(double p_tolerance) {
        //#[ operation RateBasedVT(double) 
        tolerance = p_tolerance;
        //#]
    }
    
    //## operation calculateRchar(PresentStatus) 
    public double calculateRchar(PresentStatus p_presentStatus) {
        //#[ operation calculateRchar(PresentStatus) 
        double minflux = 0;
        Iterator iter = p_presentStatus.getSpeciesStatus();
        while (iter.hasNext()) {
        	SpeciesStatus ss = (SpeciesStatus)iter.next();
        	if (ss.isReactedSpecies()) {
        		double flux = ss.getFlux();
        		minflux += flux*flux;
        	}
        }
        // calculate L2 norm of reacted species
        return Math.sqrt(minflux);
        
        
        //#]
    }
    
    //## operation calculateRmin(PresentStatus) 
    public double calculateRmin(PresentStatus p_presentStatus) {
        //#[ operation calculateRmin(PresentStatus) 
        Rmin = tolerance*calculateRchar(p_presentStatus);
        return Rmin;
        //#]
    }
    
    //## operation isModelValid(ReactionSystem) 
    public boolean isModelValid(ReactionSystem p_reactionSystem) {
        //#[ operation isModelValid(ReactionSystem) 
        // check if all the unreacted species has their fluxes under the system min flux
        PresentStatus ps = p_reactionSystem.getPresentStatus();
        calculateRmin(ps);
        for (Iterator iter =((CoreEdgeReactionModel) p_reactionSystem.getReactionModel()).getUnreactedSpeciesSet().iterator(); iter.hasNext(); ) {
        	Species s = (Species)iter.next();
        	if (ps.unreactedSpeciesFlux[s.getID()] > Rmin) 
            {
                System.out.println("Exceeded largest permitted flux for convergence (tolerance="+tolerance+"): " + Rmin);
                return false;
            }
        	
        }
        //System.out.println("The minimum flux required is " + Rmin);
		
        return true;
        
        
        //#]
    }
    
    public double getTolerance() {
        return tolerance;
    }
    
    public void setTolerance(double p_tolerance) {
        tolerance = p_tolerance;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\RateBasedVT.java
*********************************************************************/

