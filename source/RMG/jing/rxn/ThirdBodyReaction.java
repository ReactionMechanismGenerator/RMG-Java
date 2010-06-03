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
      double rate = super.calculateTotalRate(temp);
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
  public StringBuilder formPDepSign(StringBuilder p_string) {
      //#[ operation formPDepSign(String) 
      StringTokenizer st = new StringTokenizer(p_string.toString(), "=");
      StringBuilder s1 = new StringBuilder(st.nextToken());
      if (this instanceof TROEReaction || this instanceof LindemannReaction)
    	  s1.append("(+m)=");
      else
    	  s1.append("+m=");
      s1.append(st.nextToken());
      if (this instanceof TROEReaction || this instanceof LindemannReaction)
    	  s1.append("(+m)");
      else
    	  s1.append("+m");
      return s1;
      
      //#]
  }
  
  public StringBuilder formPDepSignForRestart(StringBuilder p_string) {
      StringTokenizer st = new StringTokenizer(p_string.toString(), "=");
      StringBuilder s1 = new StringBuilder(st.nextToken());
      if (this instanceof TROEReaction || this instanceof LindemannReaction)
    	  s1.append("(+m) = ");
      else
    	  s1.append("+m = ");
      s1.append(st.nextToken());
      if (this instanceof TROEReaction || this instanceof LindemannReaction)
    	  s1.append("(+m)");
      else
    	  s1.append("+m");
      return s1;
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
      StringBuilder s = getStructure().toChemkinString(true);
      s = formPDepSign(s);
      for (int i=0; i<getKinetics().length; i++) {
    	  s.append("\t" + getKinetics()[i].toChemkinString(calculateHrxn(p_temperature),p_temperature, true) + '\n');
      }
      
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
      		s.append(tbr + '\n');
      		tbr = "";
      	}
      	tbr += next;
      
      }
      
      s.append(tbr)	;
      return s.toString();
      
      //#]
  }
  
  public String toChemkinString(Temperature p_temperature, Pressure p_pressure) {
	  return toChemkinString(p_temperature) + "\n";
  }
  
  public String toRestartString(Temperature t) {
      StringBuilder s = getStructure().toRestartString(true);
      s = formPDepSignForRestart(s);
      for (int i=0; i<getKinetics().length; i++) {
    	  s.append("\t" + getKinetics()[i].toChemkinString(calculateHrxn(t),t,false) + " 0.0 0.0 0.0\n");
      }
      
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
//      	if ((tbr.length()+next.length())>=80) {
//      		s.append(tbr + '\n');
//      		tbr = "";
//      	}
      	tbr += next;
      
      }
      
      s.append(tbr)	;
      return s.toString();
  }
  
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ThirdBodyReaction.java
*********************************************************************/

