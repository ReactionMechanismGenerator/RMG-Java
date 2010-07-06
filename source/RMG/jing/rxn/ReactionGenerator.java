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


import java.util.*;
import jing.chem.Species;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\ReactionGenerator.java                                                                  
//----------------------------------------------------------------------------

/**
This interface defines a standard to react a list of species as reactants to form a reaction system.  Any reaction generator should implement this interface.
*/
//## class ReactionGenerator 
public interface ReactionGenerator {
    
    
    /**
    This abstract operation defines the way a reaction generator generates a reaction system from a species list.
    The returned value should be a reaction system.
    */
    // Argument HashSetp_speciesseed : 
    /**
    a set of reactants.
    */
    //## operation react(HashSet) 
	LinkedHashSet react(LinkedHashSet p_speciesseed);
    
    //## operation react(HashSet,Species) 
	LinkedHashSet react(LinkedHashSet p_speciesSet, Species p_species, String p_rxnFamilyName);
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ReactionGenerator.java
*********************************************************************/

