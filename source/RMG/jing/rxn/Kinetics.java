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



package jing.rxn;


import java.util.*;
import jing.param.Temperature;
import jing.mathTool.UncertainDouble;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\Kinetics.java                                                                  
//----------------------------------------------------------------------------

/**
This interface defines the functionality of rate calculation.
k=f(T)
Note, this k is k infinity, which is the high pressure limit rate of a reaction.  therefore, it is only a function of temperature. 
*/
//## class Kinetics 
public interface Kinetics {
    
    
    // Argument Temperaturetemperature : 
    /**
    The temperature that the rate is calculated.
    */
    //## operation calculateRate(Temperature,double) 
    double calculateRate(Temperature temperature, double Hrxn);
	
	double calculateRate(Temperature temperature);
    
    //## operation getA() 
    UncertainDouble getA();
    
    //## operation getAValue() 
    double getAValue();
    
    //## operation getComment() 
    String getComment();
    
    //## operation getE() 
    UncertainDouble getE();
    
    //## operation getEValue() 
    double getEValue();
    
    //## operation getN() 
    UncertainDouble getN();
    
    //## operation getNValue() 
    double getNValue();
    
    //## operation getRank() 
    int getRank();
    
	//## operation getTRange()
	String getTRange();
	
    //## operation getSource() 
    String getSource();
	
	void setSource(String p_string);
	
	void setComments(String p_string);
	
	void setFromPrimaryKineticLibrary(boolean p_boolean);
	
	boolean getFromPrimaryKineticLibrary();
    
    //## operation multiply(double) 
    Kinetics multiply(double p_multiple);
    
    //## operation repOk() 
    boolean repOk();
    
    //## operation toChemkinString() 
    String toChemkinString(double Hrxn, Temperature p_temperature, boolean includeComments);
    
    //## operation toString() 
    String toString();
	
	boolean equals(Kinetics p_k);
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\Kinetics.java
*********************************************************************/

