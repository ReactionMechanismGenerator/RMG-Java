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


import jing.chem.*;
import java.util.*;

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\ConversionTT.java                                                                  
//----------------------------------------------------------------------------

/**
Using conversion to control the termination of reaction systerm.  ie, reaction terminates at a given conversion of  some important reactants.
*/
//## class ConversionTT 
public class ConversionTT implements TerminationTester {
    
    protected LinkedList speciesGoalConversionSet;
    
    // Constructors
    
    //## operation ConversionTT() 
    private  ConversionTT() {
        {
            speciesGoalConversionSet=new LinkedList();
        }
        //#[ operation ConversionTT() 
        //#]
    }
    //## operation ConversionTT(HashSet) 
    public  ConversionTT(LinkedList p_scs) {
        {
            speciesGoalConversionSet=new LinkedList();
        }
        //#[ operation ConversionTT(HashSet) 
        speciesGoalConversionSet = p_scs;
        //#]
    }
    
    //## operation deleteSpeciesConversionSet(SpeciesConversion) 
    public void deleteSpeciesConversionSet(SpeciesConversion p_SpeciesConversion) {
        //#[ operation deleteSpeciesConversionSet(SpeciesConversion) 
        speciesGoalConversionSet.remove(p_SpeciesConversion);
        p_SpeciesConversion=null;
        //#]
    }
    
    //## operation isReactionTerminated(InitialStatus,PresentStatus) 
    public boolean isReactionTerminated(InitialStatus p_initialStatus, PresentStatus p_presentStatus) {
        //#[ operation isReactionTerminated(InitialStatus,PresentStatus) 
        Iterator iter = speciesGoalConversionSet.iterator();
        while (iter.hasNext()) {
        	SpeciesConversion sc = (SpeciesConversion)iter.next();
        	Species spe = sc.getSpecies();
        	SpeciesStatus iss = p_initialStatus.getSpeciesStatus(spe);
        	SpeciesStatus pss = p_presentStatus.getSpeciesStatus(spe);
        	double init_conc = iss.getConcentration();
        	double pres_conc = pss.getConcentration();
        	if (pres_conc>init_conc || init_conc == 0) throw new InvalidConversionException();
        	double conversion = (init_conc-pres_conc)/init_conc;
        	if (conversion < sc.getConversion()) return false;
        }
        return true;
        //#]
    }
    
    public Iterator getSpeciesGoalConversionSet() {
        Iterator iter=speciesGoalConversionSet.iterator();
        return iter;
    }
   
    //11/1/07 gmagoon: created alternate accessor to return LinkedList; UPDATE: not needed
    //public LinkedList getSpeciesGoalConversionSetList(){
    //    return speciesGoalConversionSet;
    //}
    
    public void addSpeciesGoalConversionSet(SpeciesConversion p_SpeciesConversion) {
        speciesGoalConversionSet.add(p_SpeciesConversion);
    }
    
    public void removeSpeciesGoalConversionSet(SpeciesConversion p_SpeciesConversion) {
        speciesGoalConversionSet.remove(p_SpeciesConversion);
    }
    
    public void clearSpeciesGoalConversionSet() {
        speciesGoalConversionSet.clear();
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\ConversionTT.java
*********************************************************************/

