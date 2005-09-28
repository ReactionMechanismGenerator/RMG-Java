//!********************************************************************************
//!
//!    RMG: Reaction Mechanism Generator                                            
//!
//!    Copyright: Jing Song, MIT, 2002, all rights reserved
//!     
//!    Author's Contact: jingsong@mit.edu
//!
//!    Restrictions:
//!    (1) RMG is only for non-commercial distribution; commercial usage
//!        must require other written permission.
//!    (2) Redistributions of RMG must retain the above copyright
//!        notice, this list of conditions and the following disclaimer.
//!    (3) The end-user documentation included with the redistribution,
//!        if any, must include the following acknowledgment:
//!        "This product includes software RMG developed by Jing Song, MIT."
//!        Alternately, this acknowledgment may appear in the software itself,
//!        if and wherever such third-party acknowledgments normally appear.
//!  
//!    RMG IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED 
//!    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
//!    OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
//!    DISCLAIMED.  IN NO EVENT SHALL JING SONG BE LIABLE FOR  
//!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
//!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
//!    OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;  
//!    OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  
//!    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
//!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
//!    THE USE OF RMG, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//! 
//!******************************************************************************



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

