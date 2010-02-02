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



package jing.chem;


import java.util.*;

//## package jing::chem 

//----------------------------------------------------------------------------
// jing\chem\ChemNodeElement.java                                                                  
//----------------------------------------------------------------------------

/**
Atom represents the stuff that could be put in the node site in a chem graph.  Atom is designed as an interface inherited by ChemElement, Radical, and Ion.  There are three types of Atom:  
(1) ChemElement: no free radical, no charge on the central chemical element
(2) Radical: with free radical, no chare on the central chemical element
(3) Ion: no free radical, with charge on teh central chemical element. (not implement in this package, but
users can extend to define Ion if they want.)
*/
//## class ChemNodeElement 
public interface ChemNodeElement {
    
    
    //## operation changeRadical(int,String) 
    ChemNodeElement changeRadical(int p_radical, String p_spin);
    
    //## operation getFreeElectron() 
    FreeElectron getFreeElectron();
    
    //## operation getName() 
    String getName();
    
    //## operation isAny() 
    boolean isAny();
    
    //## operation isCarbon() 
    boolean isCarbon();
    
    //## operation isHydrogen() 
    boolean isHydrogen();
    
    //## operation isNonH() 
    boolean isNonH();
    
    //## operation isOxygen() 
    boolean isOxygen();
    
    //## operation isRadical() 
    boolean isRadical();
    
    // Added by MRH on 18-Jun-2009
    //	Hardcoding Si and S into RMG-java
    boolean isSilicon();
    
    boolean isSulfur();
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\ChemNodeElement.java
*********************************************************************/

