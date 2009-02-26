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

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\ConstantPM.java                                                                  
//----------------------------------------------------------------------------

//## class ConstantPM 
public class ConstantPM implements PressureModel {
    
    protected Pressure pressure;		//## attribute pressure 
    
    
    // Constructors
    
    //## operation ConstantPM(Pressure) 
    public  ConstantPM(Pressure p_pressure) {
        //#[ operation ConstantPM(Pressure) 
        pressure = p_pressure;
        //#]
    }
    //## operation ConstantPM(double,String) 
    public  ConstantPM(double p_pressure, String p_unit) {
        //#[ operation ConstantPM(double,String) 
        pressure = new Pressure(p_pressure, p_unit);
        
        
        //#]
    }
    public  ConstantPM() {
    }
    
    //## operation getPressure(ReactionTime) 
    public Pressure getPressure(ReactionTime p_reactionTime) {
        //#[ operation getPressure(ReactionTime) 
        return pressure;
        //#]
    }
    
    public Pressure getPressure() {
        return pressure;
    }
    
    public void setPressure(Pressure p_pressure) {
        pressure = p_pressure;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\ConstantPM.java
*********************************************************************/

