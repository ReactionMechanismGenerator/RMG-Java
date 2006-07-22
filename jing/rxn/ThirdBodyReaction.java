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



package jing.rxn;


import jing.chem.*;
import java.util.*;
import jing.param.*;
import jing.rxnSys.*;
import jing.rxnSys.SystemSnapshot;

//## package jing::rxn 

//----------------------------------------------------------------------------
//jing\rxn\ThirdBodyReaction.java                                                                  
//----------------------------------------------------------------------------

//## class ThirdBodyReaction 
public class ThirdBodyReaction extends Reaction {
  
  protected HashMap weightMap = new HashMap();		//## attribute weightMap 
  
  
  // Constructors
  
  //## operation ThirdBodyReaction() 
  protected  ThirdBodyReaction() {
      //#[ operation ThirdBodyReaction() 
      //#]
  }
  
  //## operation calculateRate(SystemSnapshot) 
  public double calculateRate(SystemSnapshot p_systemSnapshot) {
      //#[ operation calculateRate(SystemSnapshot) 
      Temperature temp = p_systemSnapshot.getTemperature();
      //Kinetics k = getKinetics();
      //double rate = k.calculateRate(temp, calculateHrxn(temp));
      double rate = super.calculateRate(temp);
      rate *= calculateThirdBodyCoefficient(p_systemSnapshot);
      return rate;
      //#]
  }
  
  public HashMap getWeightMap(){
	  return weightMap;
  }
  //## operation calculateThirdBodyCoefficient(SystemSnapshot) 
  public double calculateThirdBodyCoefficient(SystemSnapshot p_presentStatus) {
      //#[ operation calculateThirdBodyCoefficient(SystemSnapshot) 
      double coef_total = 0;
	  Set colliders = weightMap.keySet();
	  
      for (Iterator iter = p_presentStatus.getSpeciesStatus(); iter.hasNext(); ) {
      	SpeciesStatus ss = (SpeciesStatus)iter.next();
      	double conc = ss.getConcentration();
      	double coef = 1;
      	
      	String name = ss.getSpecies().getName();
		Iterator colliderIter = colliders.iterator();
		while (colliderIter.hasNext()){
			String colliderName = (String)colliderIter.next();
			if (colliderName.compareToIgnoreCase(name) == 0){
				coef = ((Double)weightMap.get(colliderName)).doubleValue();
				break;
			}
				
		}
      	/*if (weightMap.containsKey(name)) {
      		coef = ((Double)weightMap.get(name)).doubleValue();
      	}*/
      	
      	coef_total += coef*conc;
      }
      
      for (Iterator iter = p_presentStatus.getInertGas(); iter.hasNext(); ) {
      	String name = (String)iter.next();
      	double conc = p_presentStatus.getInertGas(name);
      	double coef = 1;
      	
      	if (weightMap.containsKey(name)) {
      		coef = ((Double)weightMap.get(name)).doubleValue();
      	}
      	
      	coef_total += coef*conc;
      }
      
      return coef_total;
      //#]
  }
  
//## operation calculateThirdBodyCoefficient(SystemSnapshot) 
  public double calculateThirdBodyCoefficientForInerts(SystemSnapshot p_presentStatus) {
      //#[ operation calculateThirdBodyCoefficient(SystemSnapshot) 
      double coef_total = 0;
	  Set colliders = weightMap.keySet();
      
      for (Iterator iter = p_presentStatus.getInertGas(); iter.hasNext(); ) {
      	String name = (String)iter.next();
      	double conc = p_presentStatus.getInertGas(name);
      	double coef = 1;
      	
      	if (weightMap.containsKey(name)) {
      		coef = ((Double)weightMap.get(name)).doubleValue();
      	}
      	
      	coef_total += coef*conc;
      }
      
      return coef_total;
      //#]
  }
  
  //## operation formPDepSign(String) 
  public String formPDepSign(String p_string) {
      //#[ operation formPDepSign(String) 
      StringTokenizer st = new StringTokenizer(p_string, "=");
      String s1 = st.nextToken();
      s1 += "+m=";
      String s2 = st.nextToken();
      s2 += "+m";
      return (s1+s2);
      
      //#]
  }
  
  //## operation generateReverseReaction() 
  public void generateReverseReaction() {
      //#[ operation generateReverseReaction() 
      ThirdBodyReaction r = new ThirdBodyReaction();
      r.structure = getStructure().generateReverseStructure();
      r.kinetics = getKinetics();
      r.comments = "Reverse reaction";
      r.weightMap = weightMap;
      	
      r.setReverseReaction(this);
      this.setReverseReaction(r);
      
      //this.setReverseReaction(null);
      
      return;
      //#]
  }
  
  //## operation make(Reaction,HashMap) 
  public static ThirdBodyReaction make(Reaction p_reaction, HashMap p_thirdBodyList) {
      //#[ operation make(Reaction,HashMap) 
      ThirdBodyReaction tbr = new ThirdBodyReaction();
      tbr.structure = p_reaction.getStructure();
      tbr.kinetics = p_reaction.getKinetics();
      tbr.comments = p_reaction.getComments();
      tbr.weightMap = p_thirdBodyList;
      tbr.generateReverseReaction(); 
      
      return tbr;
      
      
      //#]
  }
  
  //## operation putThirdBodyCoefficient(String,double) 
  public void putThirdBodyCoefficient(String p_name, double p_coefficient) {
      //#[ operation putThirdBodyCoefficient(String,double) 
      weightMap.put(p_name,new Double(p_coefficient));
      
      
      //#]
  }
  
  //## operation toChemkinString() 
  public String toChemkinString(Temperature p_temperature) {
      //#[ operation toChemkinString() 
      String s = getStructure().toChemkinString(true);
      s = formPDepSign(s);
      s = s + '\t' + getKinetics().toChemkinString(calculateHrxn(p_temperature),p_temperature) + '\n';
      
      String tbr = "";
      // write 3rd-body efficiencies
      for(Iterator iter = weightMap.keySet().iterator(); iter.hasNext();) {
      	Object key = iter.next();
      	double tbe = ((Double)weightMap.get(key)).doubleValue();
      	String name = key.toString();
      	Species spe = SpeciesDictionary.getInstance().getSpeciesFromName(name);
      	if (spe!=null) name = spe.getChemkinName();
      	if (name.equals("AR")) name = "Ar";
      	String next = name + "/" + tbe + "/ ";
      	if ((tbr.length()+next.length())>=80) {
      		s += tbr + '\n';
      		tbr = "";
      	}
      	tbr += next;
      
      }
      
      s += tbr;
      return s;
      
      //#]
  }
  
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ThirdBodyReaction.java
*********************************************************************/

