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

