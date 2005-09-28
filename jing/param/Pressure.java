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

