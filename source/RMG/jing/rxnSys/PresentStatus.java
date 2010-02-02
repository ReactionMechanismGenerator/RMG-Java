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
import jing.param.Pressure;
import jing.param.Temperature;

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\PresentStatus.java                                                                  
//----------------------------------------------------------------------------

//## class PresentStatus 
public class PresentStatus extends SystemSnapshot {
    
    
    // Constructors
    
	   //## operation PresentStatus(ReactionTime,HashMap,Temperature,Pressure) 
    public  PresentStatus(ReactionTime p_reactionTime, LinkedHashMap p_speciesStatus, Temperature p_temperature, Pressure p_pressure) {
        //#[ operation PresentStatus(ReactionTime,HashMap,Temperature,Pressure) 
        super(p_reactionTime, p_speciesStatus, p_temperature, p_pressure);
        //#]
    }
    //## operation PresentStatus(SystemSnapshot) 
    public  PresentStatus(SystemSnapshot p_systemSnapshot) {
        //#[ operation PresentStatus(SystemSnapshot) 
        time = p_systemSnapshot.time;
        speciesStatus = p_systemSnapshot.speciesStatus;
        temperature = p_systemSnapshot.temperature;
        pressure = p_systemSnapshot.pressure;
		unreactedSpeciesFlux = p_systemSnapshot.unreactedSpeciesFlux;
		inertGas = p_systemSnapshot.inertGas;
        //#]
    }
    public  PresentStatus() {
    }
    
    //## operation getPresentTime() 
    public ReactionTime getPresentTime() {
        //#[ operation getPresentTime() 
        return super.getTime();
        //#]
    }
	
	/*public double getUnreactedSpeciesFlux(Species species) {
		return unreactedSpeciesFlux[species.getID()];
	}*/
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\PresentStatus.java
*********************************************************************/

