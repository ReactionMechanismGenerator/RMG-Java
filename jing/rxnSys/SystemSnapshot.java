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
import jing.param.*;
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

    protected Pressure pressure;		//## attribute pressure

    protected Temperature temperature;		//## attribute temperature

    protected double totalMole = -1;		//## attribute totalMole

    protected HashMap speciesStatus;
    protected LinkedList senStatus;//svp
    protected ReactionTime time;
    protected LinkedList reactionList;//svp
    protected HashMap IDTranslator;//svp

    // Constructors

    //## operation SystemSnapshot()
    public  SystemSnapshot() {
        {
            speciesStatus = new HashMap();
        }
        initRelations();
        //#[ operation SystemSnapshot()
        time = new ReactionTime(0,"S");


        //#]
    }
    //## operation SystemSnapshot(ReactionTime,HashMap,Temperature,Pressure)
    public  SystemSnapshot(ReactionTime p_reactionTime, HashMap p_speciesStatus, Temperature p_temperature, Pressure p_pressure) {
        {
            speciesStatus = new HashMap();
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
    public SystemSnapshot(ReactionTime p_reactionTime, HashMap p_speciesStatus, LinkedList p_sensitivityStatus, Temperature p_temperature, Pressure p_pressure) {
      {
        speciesStatus = new HashMap();
        senStatus = new LinkedList();
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
    public void addSensitivity(LinkedList p_senStatus) {
      //#[ operation addSensitivity(LinkedList)
      senStatus = p_senStatus;
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

    //## operation getSensitivityStatus()
    //svp
      public Iterator getSensitivityStatus(){
        //#[operation getSensitivityStatus()
        Iterator iter = senStatus.iterator();
        return iter;
      }

      //## operation getSensitivityStatus(int)
      //svp
        public SensitivityStatus getSensitivityStatus(int p_int){
          //#[operation getSensitivityStatus(int)
          return (SensitivityStatus)senStatus.get(p_int);
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

    //## operation isTPCConsistent()
    public boolean isTPCConsistent() {
        //#[ operation isTPCConsistent()
        double T = temperature.getK();
        double P = pressure.getPa();
        double R = GasConstant.getStandard();

        double realConc = P/(R*T*Math.pow(10,6.0));
        double conc = getTotalMole();

        double error = Math.abs(realConc-conc)/realConc;

        if (error < 0.01) return true;
        else return false;


        //#]
    }

    //## operation putInertGas(String,double)
    public static void putInertGas(String p_name, double p_concentration) {
        //#[ operation putInertGas(String,double)
        inertGas.put(p_name, new Double(p_concentration));
        //#]
    }

    //## operation putSensitivityStatus(int)
    //svp
      public void putSensitivityStatus(int index, SensitivityStatus p_senStatus) {
        //#[ operation putSensitivityStatus(int)
        int spe_num = p_senStatus.getSpeciesNumber();
        int rxn_num = p_senStatus.getReactionNumber();
        senStatus.add(index, p_senStatus);
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
      public void setIDTranslator(HashMap p_hashMap){
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

}
/*********************************************************************
        File Path	: RMG\RMG\jing\rxnSys\SystemSnapshot.java
*********************************************************************/

