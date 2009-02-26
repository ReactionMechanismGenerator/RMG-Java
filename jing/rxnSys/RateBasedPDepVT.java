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


import jing.rxn.*;
import jing.chem.*;
import java.util.*;
import jing.param.*;

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\RateBasedPDepVT.java                                                                  
//----------------------------------------------------------------------------

//## class RateBasedPDepVT 
public class RateBasedPDepVT extends RateBasedVT {
    
    
    // Constructors
    
    //## operation RateBasedPDepVT(double) 
    public  RateBasedPDepVT(double p_tolerance) {
        //#[ operation RateBasedPDepVT(double) 
        super(p_tolerance);
        //#]
    }
    public  RateBasedPDepVT() {
    }
    
    //## operation isModelValid(ReactionSystem) 
    public boolean isModelValid(ReactionSystem p_reactionSystem) {
        
		// Check if all the unreacted species has their fluxes under the system min flux
        if (!super.isModelValid(p_reactionSystem))
			return false;

		// check if all the networks has their leak fluxes under the system min flux
        PresentStatus ps = p_reactionSystem.getPresentStatus();
        calculateRmin(ps);

		for (Iterator iter = PDepNetwork.getNetworks().iterator(); iter.hasNext(); ) {
        	PDepNetwork pdn = (PDepNetwork )iter.next();
        	if (pdn.getLeakFlux(ps) > Rmin)
            {
                System.out.println("Exceeded largest permitted flux for convergence (tolerance="+tolerance+"): " + Rmin);
                return false;
            }

        }
        
        return true;


        //#]

    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\RateBasedPDepVT.java
*********************************************************************/

