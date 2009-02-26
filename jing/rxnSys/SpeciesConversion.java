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
// jing\rxnSys\SpeciesConversion.java                                                                  
//----------------------------------------------------------------------------

//## class SpeciesConversion 
public class SpeciesConversion {
    
    protected double conversion;		//## attribute conversion 
    
    protected Species species;		//## attribute species 
    
    
    // Constructors
    
    //## operation SpeciesConversion(Species,double) 
    public  SpeciesConversion(Species p_species, double p_conversion) {
        //#[ operation SpeciesConversion(Species,double) 
        species = p_species;
        conversion = p_conversion;
        //#]
    }
    public  SpeciesConversion() {
    }
    
    public double getConversion() {
        return conversion;
    }
    
    public void setConversion(double p_conversion) {
        conversion = p_conversion;
    }
    
    public Species getSpecies() {
        return species;
    }
    
    public void setSpecies(Species p_species) {
        species = p_species;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\SpeciesConversion.java
*********************************************************************/

