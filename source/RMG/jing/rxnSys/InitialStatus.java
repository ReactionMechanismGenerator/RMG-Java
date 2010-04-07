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
import jing.param.Pressure;
import jing.param.Temperature;


public class InitialStatus extends SystemSnapshot {
    
    protected static double colliderLimit = 0.1;		//## attribute colliderLimit 
    
    // Constructors
    
    private  InitialStatus() {
    }
	
	
	public  InitialStatus(LinkedHashMap p_speciesStatus, Temperature p_temperature, Pressure p_pressure) {
		super(new ReactionTime(0,"S"),p_speciesStatus,p_temperature,p_pressure);
	}
    
	
    public HashMap identifyColliders() {
		
        HashMap result = new HashMap();
        double totalMole = getTotalMole();
        double adjTotalMole = 0;
        for (Iterator iter = getSpeciesStatus(); iter.hasNext(); ) {
        	SpeciesStatus ss = (SpeciesStatus)iter.next();
        	double conc = ss.getConcentration();
        	if (conc > colliderLimit*totalMole) {
        		adjTotalMole += conc;
        		result.put(ss.getSpecies(), new Double(conc)); 
        	}
        }
        for (Iterator iter = this.inertGas.keySet().iterator(); iter.hasNext(); ) {//gmagoon 6/23/09: changed to use this.inertGas... rather than just inertGas... to accommodate now non-static variable inertGas
        	Object key = iter.next();
        	double conc = ((Double)inertGas.get(key)).doubleValue(); 
        	if (conc < 0) {
				double aTol = ReactionModelGenerator.getAtol();
				//if (Math.abs(conc) < aTol) conc = 0;
				//else throw new NegativeConcentrationException("InertGas");
				if (conc < -100.0 * aTol)
					throw new NegativeConcentrationException("Inert Gas has negative concentration: " + String.valueOf(conc));
        	}
        	if (conc > colliderLimit*totalMole) {
        		adjTotalMole += conc;
        		result.put(key, new Double(conc)); 
        	}
        }
        
        if (result.isEmpty()) throw new InvalidSystemCompositionException();
        
        return result;
		
    }
    private static double getColliderLimit() {
        return colliderLimit;
    }
}
