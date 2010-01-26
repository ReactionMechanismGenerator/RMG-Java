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

import jing.mathTool.MathTool;
import jing.param.*;
import jing.chem.Species;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.rxn.PDepReaction;
import jing.rxn.Reaction;

//## package jing::rxnSys

//----------------------------------------------------------------------------
// jing\rxnSys\SystemSnapshot.java
//----------------------------------------------------------------------------

/**
The real status of the system at a specific reaction time.  For example: the T, P, composition,
*/
//## class SystemSnapshot
public class SystemSnapshot {

    protected LinkedHashMap inertGas = new LinkedHashMap();		//## attribute inertGas //6/23/09 gmagoon: it seems that for multiple T,P cases, there will be multiple system snapshots that have separate concentrations, so I think this should not be static

    protected Pressure pressure;		//## attribute pressure

    protected Temperature temperature;		//## attribute temperature

    protected double totalMole = -1;		//## attribute totalMole

    protected LinkedHashMap speciesStatus;
    protected double [] senStatus;//svp
    protected ReactionTime time;
    protected LinkedList reactionList;//svp
    protected LinkedHashMap IDTranslator;//svp
	protected double[] unreactedSpeciesFlux = null;
	protected double[] reactionFlux ;

    // Constructors

    //## operation SystemSnapshot()
    public  SystemSnapshot() {
        {
            speciesStatus = new LinkedHashMap();
        }
        initRelations();
        //#[ operation SystemSnapshot()
        time = new ReactionTime(0,"S");


        //#]
    }
    //## operation SystemSnapshot(ReactionTime,HashMap,Temperature,Pressure)
    public  SystemSnapshot(ReactionTime p_reactionTime, LinkedHashMap p_speciesStatus, Temperature p_temperature, Pressure p_pressure) {
        {
            speciesStatus = new LinkedHashMap();
        }
        initRelations();
        //#[ operation SystemSnapshot(ReactionTime,HashMap,Temperature,Pressure)
        time = p_reactionTime;
        speciesStatus = p_speciesStatus;
        temperature = p_temperature;
        pressure = p_pressure;
        //#]
    }

    //## operation SystemSnapshot(ReactionTime,HashMap,HashMap)
    //svp
    public SystemSnapshot(ReactionTime p_reactionTime, LinkedHashMap p_speciesStatus, double [] p_sensitivityStatus, Temperature p_temperature, Pressure p_pressure) {
      {
        speciesStatus = new LinkedHashMap();
      }
      initRelations();
      time = p_reactionTime;
      speciesStatus = p_speciesStatus;
      senStatus = p_sensitivityStatus;
      temperature = p_temperature;
        pressure = p_pressure;
    }

    //## operation addSensitivity(LinkedList)
    //svp
    public void addSensitivity(double [] p_senStatus) {
      //#[ operation addSensitivity(LinkedList)
      senStatus = p_senStatus;
      //#]
    }



    //## operation addSpeciesStatus(HashMap)
    public void addSpeciesStatus(LinkedHashMap p_speciesStatus) {
        //#[ operation addSpeciesStatus(HashMap)
        if (speciesStatus == null) speciesStatus = new LinkedHashMap();

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

    //## operation getSensitivityStatus()
    //svp
      public double [] getSensitivityStatus(){
        //#[operation getSensitivityStatus()
        
        return senStatus;
      }

      //## operation getSensitivityStatus(int)
      //svp
        public double getSensitivityStatus(int p_int){
          //#[operation getSensitivityStatus(int)
          return senStatus[p_int];
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
                if (c < 0) {
    				double aTol = ReactionModelGenerator.getAtol();
    				//if (Math.abs(c) < aTol) c = 0;
    				//else throw new NegativeConcentrationException("InertGas");
					if (c < -100.0 * aTol)
						throw new NegativeConcentrationException("Total inert gas has negative concentration: " + String.valueOf(c));
                }
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
                        if (ss.getConcentration()<0) {
            				double aTol = ReactionModelGenerator.getAtol();
            				//if (Math.abs(ss.getConcentration()) < aTol) ss.setConcentration(0.0);
            				//else throw new NegativeConcentrationException();
							if ( ss.getConcentration() < -100.0 * aTol )
								throw new NegativeConcentrationException("Total inert gas has negative concentration: " + String.valueOf(ss.getConcentration()));
                        }
                        totalMole += ss.getConcentration();
                }
        }
        return totalMole;
        //#]
    }

    //## operation isTPCConsistent()
    public boolean isTPCConsistent(DynamicSimulator ds, FinishController fc) {
        //#[ operation isTPCConsistent()
        double T = temperature.getK();
        double P = pressure.getPa();
        double R = GasConstant.getStandard();

        double realConc = P/(R*T*Math.pow(10,6.0));
        double conc = getTotalMole();

        double error = Math.abs(realConc-conc)/realConc;

        if (error < 0.00001) {
        	
        	return true;
        }
        else {
//        	renormalizing the concentration
        	HashSet hs= new HashSet(inertGas.keySet());
        	for (Iterator iter = hs.iterator(); iter.hasNext(); ) {
                String name = (String)iter.next();
        		double c = (Double)inertGas.get(name);
                c = c*realConc/conc;
                //System.out.println(MathTool.formatDouble(c, 15,6));
                inertGas.remove(name);
                inertGas.put(name, c);
        	}
        	for (Iterator iter = getSpeciesStatus(); iter.hasNext(); ) {
                SpeciesStatus ss = (SpeciesStatus)iter.next();
                //System.out.println(String.format("%1.6e",ss.getConcentration()*realConc/conc));
                ss.setConcentration(ss.getConcentration()*realConc/conc);
        	}
        	if (fc.terminationTester instanceof ConversionTT) {
        		double [] conversion = ds.getConversion();
        		for (int j=0; j< conversion.length; j++){
        			conversion[j] = conversion[j]*realConc/conc;
        		}
        	}
        	totalMole = realConc;
        	return false;
        }


        //#]
    }

    //## operation putInertGas(String,double)
    public void putInertGas(String p_name, double p_concentration) {
        //#[ operation putInertGas(String,double)
        inertGas.put(p_name, new Double(p_concentration));
        //#]
    }

    //## operation putSensitivityStatus(int)
    //svp
      public void putSensitivityStatus(int index, double p_senStatus) {
        //#[ operation putSensitivityStatus(int)
        senStatus[index] = p_senStatus;
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

    public Pressure getPressure() {
        return pressure;
    }

    //svp
    public LinkedList getReactionList(){
      return reactionList;
    }


	public LinkedList getUniqueReactionList() {
		LinkedList uniqueList = new LinkedList();
		for (Iterator iter = reactionList.iterator(); iter.hasNext(); ) {
			Reaction rxn = (Reaction) iter.next();
			if (rxn instanceof PDepReaction) {
				// PDepReaction reactions are already handled correctly
				uniqueList.add(rxn);
			}
			else {
				if (!uniqueList.contains(rxn) && !uniqueList.contains(rxn.getReverseReaction()))
					uniqueList.add(rxn);
			}
		}
		return uniqueList;
	}

    //svp
    //## operation getRealID(Species)
      public int getRealID(Species p_species) {
          //#[ operation getRealID(Species)
          Integer id = (Integer)IDTranslator.get(p_species);
          return id.intValue();
          //#]
      }

      //## operation setIDTranslator(HashMap)
    //svp
      public void setIDTranslator(LinkedHashMap p_hashMap){
        //#[ operation setIDTranslator(HashMap)
        IDTranslator = p_hashMap;
        //#]
      }



    public void setPressure(Pressure p_pressure) {
        pressure = p_pressure;
    }

    //## operation setReactionList(LinkedList)
    //svp
      public void setReactionList(LinkedList p_list){
        //#[ operation setReactionList(LinkedList)
        reactionList = p_list;
        //#]
      }

      public void setReactionFlux(double [] p_reactionFlux) {
    	  reactionFlux = p_reactionFlux;
      }

    public Temperature getTemperature() {
        return temperature;
    }

    public void setTemperature(Temperature p_temperature) {
        temperature = p_temperature;
    }

    public void deleteSpeciesStatus(SpeciesStatus p_SpeciesStatus) {
        Iterator iter = speciesStatus.keySet().iterator();
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

    public void deleteTime() {
        time=null;
    }

    protected void initRelations() {
        time = newTime();
    }
	
	public double getUnreactedSpeciesFlux(Species species) {
		return unreactedSpeciesFlux[species.getID()];
	}

	public double[] getUnreactedSpeciesFlux() {
		return unreactedSpeciesFlux;
	}

}
/*********************************************************************
        File Path	: RMG\RMG\jing\rxnSys\SystemSnapshot.java
*********************************************************************/


