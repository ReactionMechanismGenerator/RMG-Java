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

//## package jing::rxnSys 

//----------------------------------------------------------------------------
// jing\rxnSys\SpeciesStatus.java                                                                  
//----------------------------------------------------------------------------

//## class SpeciesStatus 
public class SpeciesStatus {
    
    protected double concentration;		//## attribute concentration 
    
    protected double flux;		//## attribute flux 
    
    protected Species species;		//## attribute species 
    
    protected int type = -1;		//## attribute type 
    
    
    // Constructors
    
    //## operation SpeciesStatus() 
    private  SpeciesStatus() {
        //#[ operation SpeciesStatus() 
        //#]
    }
    // Argument intp_type : 
    /**
    if species is reacted (core) species, type =1;
    if species is unreacted (edge) species, type = 0;
    */
    //## operation SpeciesStatus(Species,int,double,double) 
    public  SpeciesStatus(Species p_species, int p_type, double p_concentration, double p_flux) {
        //#[ operation SpeciesStatus(Species,int,double,double) 
        species = p_species;
        type = p_type;
        concentration = p_concentration;
        flux = p_flux;
        
        
        //#]
    }
    
    //## operation isReactedSpecies() 
    public boolean isReactedSpecies() {
        //#[ operation isReactedSpecies() 
        return (type == 1);
        //#]
    }
    
    //## operation isUnreactedSpecies() 
    public boolean isUnreactedSpecies() {
        //#[ operation isUnreactedSpecies() 
        return (type == 0);
        //#]
    }
    
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        if (species == null || !species.repOk()) return false;
        
        if (concentration < 0) return false;
        
        if (type == 0) {
        	// unreacted species will have a non-negative formation rate
        	if (flux < 0) return false;
        }
        return true;
        //#]
    }
    
    public double getConcentration() {
        return concentration;
    }
    
    public double getFlux() {
        return flux;
    }
    
    public void setFlux(double p_flux) {
    	flux = p_flux;
    }
    
    public void setConcentration(double p_conc) {
    	concentration = p_conc;
    }
    
    public Species getSpecies() {
        return species;
    }

	//11/1/07 gmagoon: added accessor for speciesType
    public int getSpeciesType() {
        return type;
    }   
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\SpeciesStatus.java
*********************************************************************/

