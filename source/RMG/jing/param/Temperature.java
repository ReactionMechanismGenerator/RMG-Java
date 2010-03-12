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
// jing\param\Temperature.java                                                                  
//----------------------------------------------------------------------------

//## class Temperature 
public class Temperature {
    
    /**
    It is the standard temperature in kelvin.
    */
    protected double standard;		//## attribute standard 
    
    
    // Constructors
    
    //## operation Temperature(double,String) 
    public  Temperature(double p_T, String p_unit) {
        //#[ operation Temperature(double,String) 
        if (p_unit.compareToIgnoreCase("K") == 0) {	
        	standard = p_T;
        }
        else if (p_unit.compareToIgnoreCase("C") == 0) {
        	standard = p_T + 273.15;
        }               
        else if (p_unit.compareToIgnoreCase("F") == 0) {
        	standard = (p_T - 32)*5/9 + 273.15;
        }
        else {
        	throw new InvalidUnitException();
        }
        //#]
    }
    public  Temperature() {
    }
    
    //## operation clone() 
    public Object clone() {
        //#[ operation clone() 
        return new Temperature(standard, getStandardUnit());
        //#]
    }
    
    //## operation equals(Object) 
    public boolean equals(Object p_temperature) {
        //#[ operation equals(Object) 
        if (!(p_temperature instanceof Temperature)) return false;
        Temperature t = (Temperature)p_temperature;
        
        return (t.getStandard()==getStandard());
        //#]
    }
    
    //## operation getC() 
    public double getC() {
        //#[ operation getC() 
        return (standard - 273.15);
        //#]
    }
    
    //## operation getF() 
    public double getF() {
        //#[ operation getF() 
        return ((standard-273.15)*9/5+32);
        //#]
    }
    
    //## operation getK() 
    public double getK() {
        //#[ operation getK() 
        return standard;
        //#]
    }
    
    //## operation getStandardUnit() 
    public static String getStandardUnit() {
        //#[ operation getStandardUnit() 
        return "K";
        //#]
    }
    
    //## operation hashCode() 
    public int hashCode() {
        //#[ operation hashCode() 
        return (int)standard;
        //#]
    }

    public String toString() {
        return String.format("%g K",this.getK() );
    }

    public double getStandard() {
        return standard;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\param\Temperature.java
*********************************************************************/

