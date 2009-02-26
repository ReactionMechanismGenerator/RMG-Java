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



package jing.param;


import java.util.*;

//## package jing::param 

//----------------------------------------------------------------------------
// jing\param\GasConstant.java                                                                  
//----------------------------------------------------------------------------

//## class GasConstant 
public class GasConstant {
    
    
    // Constructors
    
    public  GasConstant() {
    }
    
    //## operation getCCAtmMolK() 
    public static double getCCAtmMolK() {
        //#[ operation getCCAtmMolK() 
        //return 82.059;
    	return 82.053;
        //#]
    }
    
    //## operation getCalMolK() 
    public static double getCalMolK() {
        //#[ operation getCalMolK() 
        return 1.987;
        //#]
    }
    
    //## operation getJMolK() 
    public static double getJMolK() {
        //#[ operation getJMolK() 
        return 8.314;
        //#]
    }
    
    //## operation getKcalMolK() 
    public static double getKcalMolK() {
        //#[ operation getKcalMolK() 
        return 0.001987;
        //#]
    }
    
    //## operation getStandard() 
    public static double getStandard() {
        //#[ operation getStandard() 
        return getJMolK();
        //#]
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\param\GasConstant.java
*********************************************************************/

