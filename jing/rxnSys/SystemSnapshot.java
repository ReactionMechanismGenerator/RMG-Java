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
import jing.chem.Species;
import jing.param.Pressure;
import jing.param.Temperature;

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\SystemSnapshot.java                                                                  
//----------------------------------------------------------------------------

/**
The real status of the system at a specific reaction time.  For example: the T, P, composition, 
*/
//## class SystemSnapshot 
public class SystemSnapshot {
    
    protected static HashMap inertGas = new HashMap();		//## attribute inertGas 
    
    protected double totalMole = -1;		//## attribute totalMole 
	protected Pressure pressure;		//## attribute pressure 
	    
	protected Temperature temperature;		//## attribute temperature 
	       
    protected HashMap speciesStatus;
    protected ReactionTime time;
    
    // Constructors
    
    //## operation SystemSnapshot() 
    public  SystemSnapshot() {
        {
            speciesStatus=new HashMap();
        }
        initRelations();
        //#[ operation SystemSnapshot() 
        time = new ReactionTime(0,"S");
        
        
        //#]
    }
    //## operation SystemSnapshot(ReactionTime,HashMap) 
    public  SystemSnapshot(ReactionTime p_reactionTime, HashMap p_speciesStatus, Temperature p_temperature, Pressure p_pressure) {
        {
            speciesStatus=new HashMap();
        }
        initRelations();
        //#[ operation SystemSnapshot(ReactionTime,HashMap) 
        time = p_reactionTime;
        speciesStatus = p_speciesStatus;
	      temperature = p_temperature;
	        pressure = p_pressure;
	        //#]
	    }
    
    //## operation addSpeciesStatus(HashMap) 
    public void addSpeciesStatus(HashMap p_speciesStatus) {
        //#[ operation addSpeciesStatus(HashMap) 
        if (speciesStatus == null) speciesStatus = new HashMap();
        
        if (p_speciesStatus == null) return;
        
        speciesStatus.putAll(p_speciesStatus);
        
        return;
        //#]
    }
    
    //## operation getInertGas(String) 
    public double getInertGas(String p_name) {
        //#[ operation getInertGas(String) 
        Double c = (Double)inertGas.get(p_name);
        return c.doubleValue(); 
        //#]
    }
    
    //## operation getInertGas() 
    public Iterator getInertGas() {
        //#[ operation getInertGas() 
        return inertGas.keySet().iterator();
        //#]
    }
    
    //## operation getSpeciesStatus() 
    public Iterator getSpeciesStatus() {
        //#[ operation getSpeciesStatus() 
        Iterator iter=speciesStatus.values().iterator();
        return iter;
        //#]
    }
    
    //## operation getSpeciesStatus(Species) 
    public SpeciesStatus getSpeciesStatus(Species p_species) {
        //#[ operation getSpeciesStatus(Species) 
        return (SpeciesStatus)speciesStatus.get(p_species);
        //#]
    }
    
    //## operation getTotalInertGas() 
    public double getTotalInertGas() {
        //#[ operation getTotalInertGas() 
        double total = 0;
        for (Iterator iter = inertGas.values().iterator(); iter.hasNext(); ) {
        	double c = ((Double)iter.next()).doubleValue(); 
        	if (c < 0) throw new NegativeConcentrationException("InertGas");
        	total += c;
        }
        return total;
        
        
        //#]
    }
    
    //## operation getTotalMole() 
    public double getTotalMole() {
        //#[ operation getTotalMole() 
        if (totalMole<0) {
        	totalMole = getTotalInertGas();
        	for (Iterator iter = getSpeciesStatus(); iter.hasNext(); ) {
        		SpeciesStatus ss = (SpeciesStatus)iter.next();
        		if (ss.getConcentration()<0) throw new NegativeConcentrationException();
        		totalMole += ss.getConcentration();
        	}
        }
        return totalMole;
        //#]
    }
    
    //## operation putInertGas(String,double) 
    public static void putInertGas(String p_name, double p_concentration) {
        //#[ operation putInertGas(String,double) 
        inertGas.put(p_name, new Double(p_concentration));
        //#]
    }
    
    //## operation putSpeciesStatus(SpeciesStatus) 
    public void putSpeciesStatus(SpeciesStatus p_speciesStatus) {
        //#[ operation putSpeciesStatus(SpeciesStatus) 
        if (p_speciesStatus == null) throw new NullPointerException("SpeciesStatus");
        if (!p_speciesStatus.repOk()) throw new InvalidSpeciesStatusException();
        speciesStatus.put(p_speciesStatus.getSpecies(),p_speciesStatus);
        
        
        //#]
    }
    
    public void deleteSpeciesStatus(SpeciesStatus p_SpeciesStatus) {
        Iterator iter=speciesStatus.keySet().iterator();
        while(iter.hasNext()) {
          Object key = iter.next();
          if (speciesStatus.get(key).equals(p_SpeciesStatus)) {
          	speciesStatus.remove(key);
          	break;
          }
        };
        p_SpeciesStatus=null;
    }
    
    public ReactionTime getTime() {
        return time;
    }
    
    public ReactionTime newTime() {
        time = new ReactionTime();
        return time;
    }
    
	  public Pressure getPressure() {
	        return pressure;
	    }
	    
	    public void setPressure(Pressure p_pressure) {
	        pressure = p_pressure;
	    }
	    
	    public Temperature getTemperature() {
	        return temperature;
	    }
	    
	    public void setTemperature(Temperature p_temperature) {
	        temperature = p_temperature;
	    }
	    
	
    public void deleteTime() {
        time=null;
    }
    
    protected void initRelations() {
        time = newTime();
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\SystemSnapshot.java
*********************************************************************/

