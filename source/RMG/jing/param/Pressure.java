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
// jing\param\Pressure.java                                                                  
//----------------------------------------------------------------------------

//## class Pressure 
public class Pressure {
    
    protected double standard;		//## attribute standard 
    
    
    // Constructors
    
    //## operation Pressure(double,String) 
    public  Pressure(double p_pressure, String p_unit) {
        //#[ operation Pressure(double,String) 
        if (p_unit.compareToIgnoreCase("atm") == 0) {	
        	standard = p_pressure * 101325;
        }
        else if (p_unit.compareToIgnoreCase("psi") == 0) {
        	standard = p_pressure * 100000 / 14.5038;
        }               
        else if (p_unit.compareToIgnoreCase("bar") == 0) {
        	standard = p_pressure * 100000;
        }
        else if (p_unit.compareToIgnoreCase("pa") == 0) {
        	standard = p_pressure;
        }
        else if (p_unit.compareToIgnoreCase("torr") == 0) {
        	standard = p_pressure *100000 / 750.061;
        }
        else  {
        	throw new InvalidUnitException();
        }
        //#]
    }
    //## operation Pressure() 
    public  Pressure() {
        //#[ operation Pressure() 
        //#]
    }
    
    //## operation clone() 
    public Object clone() {
        //#[ operation clone() 
        return new Pressure(standard, getStandardUnit());
        //#]
    }
    
    //## operation equals(Object) 
    public boolean equals(Object p_pressure) {
        //#[ operation equals(Object) 
        if (!(p_pressure instanceof Pressure)) return false;
        Pressure p = (Pressure)p_pressure;
        
        return (p.getStandard()==getStandard());
        //#]
    }
    
    //## operation getAtm() 
    public double getAtm() {
        //#[ operation getAtm() 
        return standard/101325;
        //#]
    }
    
    //## operation getBar() 
    public double getBar() {
        //#[ operation getBar() 
        return standard/100000;
        //#]
    }
    
    //## operation getPa() 
    public double getPa() {
        //#[ operation getPa() 
        return standard;
        //#]
    }
    
    //## operation getPsi() 
    public double getPsi() {
        //#[ operation getPsi() 
        return standard*14.5038/100000;
        //#]
    }
    
    //## operation getStandardUnit() 
    public static String getStandardUnit() {
        //#[ operation getStandardUnit() 
        return "pa";
        //#]
    }
    
    //## operation hashCode() 
    public int hashCode() {
        //#[ operation hashCode() 
        return (int)standard;
        //#]
    }
    
    public String toString() {
        return String.format("%g Bar",this.getBar() );
    }
    
    protected double getStandard() {
        return standard;
    }
    
    protected void setStandard(double p_standard) {
        standard = p_standard;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\param\Pressure.java
*********************************************************************/

