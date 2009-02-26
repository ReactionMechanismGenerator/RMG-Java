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

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\ReactionTime.java                                                                  
//----------------------------------------------------------------------------

//## class ReactionTime 
public class ReactionTime {
    
    protected double time;		//## attribute time 
    
    protected String unit;		//## attribute unit 
    
    
    // Constructors
    
    //## operation ReactionTime(double,String) 
    public  ReactionTime(double p_time, String p_unit) {
        //#[ operation ReactionTime(double,String) 
        if (p_time < 0) throw new InvalidReactionTimeException();
        time = p_time;
        
        if (p_unit.compareToIgnoreCase("sec")==0 || p_unit.compareToIgnoreCase("second")==0 || p_unit.compareToIgnoreCase("s")==0)
        	unit = "SEC";
        else if (p_unit.compareToIgnoreCase("min")==0 || p_unit.compareToIgnoreCase("minute")==0 || p_unit.compareToIgnoreCase("m")==0)
        	unit = "MIN";
        else if (p_unit.compareToIgnoreCase("hr")==0 || p_unit.compareToIgnoreCase("hour")==0 || p_unit.compareToIgnoreCase("h")==0)
        	unit = "HR";
        else if (p_unit.compareToIgnoreCase("day")==0 || p_unit.compareToIgnoreCase("d")==0)
        	unit = "DAY";
        else throw new InvalidReactionTimeUnitException();
        
        	
        	
        //#]
    }
    public  ReactionTime() {
    }
    
    //## operation add(ReactionTime) 
    public ReactionTime add(ReactionTime p_reactionTime) {
        //#[ operation add(ReactionTime) 
        if (!(p_reactionTime.getUnit()).equals(getUnit())) {
        	double present = getStandardTime();
        	double add = p_reactionTime.getStandardTime();
        	
        	return new ReactionTime(present+add, getStandardUnit());
        }
        else {
        	double present = getTime();
        	double add = p_reactionTime.getTime();
        	
        	return new ReactionTime(present+add, getUnit());
        }
        
        
        //#]
    }
    
    //## operation equals(Object) 
    public boolean equals(Object p_time) {
        //#[ operation equals(Object) 
        if (!(p_time instanceof ReactionTime)) return false;
        
        return ((ReactionTime)p_time).getStandardTime() == getStandardTime();
        
        
        //#]
    }
    
    //## operation getStandardTime() 
    public double getStandardTime() {
        //#[ operation getStandardTime() 
        if (unit.equals("SEC")) return time;
        else if (unit.equals("MIN")) return 60*time;
        else if (unit.equals("HR")) return 3600*time;
        else if (unit.equals("DAY")) return 3600*24*time;
        else {
        	throw new InvalidReactionTimeUnitException();
        }
        //#]
    }
    
    //## operation getStandardUnit() 
    public static String getStandardUnit() {
        //#[ operation getStandardUnit() 
        return "S";
        //#]
    }
    
    //## operation hashCode() 
    public int hashCode() {
        //#[ operation hashCode() 
        return (int)getStandardTime();
        //#]
    }
    
    //## operation minus(ReactionTime) 
    public ReactionTime minus(ReactionTime p_reactionTime) {
        //#[ operation minus(ReactionTime) 
        if (!(p_reactionTime.getUnit()).equals(getUnit())) {
        	double present = getStandardTime();
        	double add = p_reactionTime.getStandardTime();
        	
        	if (present < add) throw new InvalidReactionTimeException("negative time");
        	
        	return new ReactionTime(present-add, getStandardUnit());
        }
        else {
        	double present = getTime();
        	double add = p_reactionTime.getTime();
        	
        	if (present < add) throw new InvalidReactionTimeException("negative time");
        	
        	return new ReactionTime(present-add, getUnit());                           
        }
        
        
        //#]
    }
    
    //## operation reach(ReactionTime) 
    public boolean reach(ReactionTime p_time) {
        //#[ operation reach(ReactionTime) 
        if (getStandardTime() >= p_time.getStandardTime()) return true;
        else return false;
        
        
        //#]
    }
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        return String.valueOf(getTime()) + getUnit();
        
        
        //#]
    }
    
    public double getTime() {
        return time;
    }
    
    public String getUnit() {
        return unit;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\ReactionTime.java
*********************************************************************/

