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
// jing\chem\FGAtom.java                                                                  
//----------------------------------------------------------------------------

//## class FGAtom 
public class FGAtom implements ChemNodeElement {
    
    protected static HashMap dictionary = new HashMap();		//## attribute dictionary 
    
    protected String name;		//## attribute name 
    
    protected FGElement fgElement;
    protected FreeElectron freeElectron;
    
    // Constructors
    
    //## operation FGAtom(FGElement,FreeElectron) 
    public  FGAtom(FGElement p_fgElement, FreeElectron p_freeElectron) throws InvalidChemNodeElementException, InvalidFreeElectronException {
        //#[ operation FGAtom(FGElement,FreeElectron) 
        fgElement = p_fgElement;
        freeElectron = p_freeElectron; 
        
        name = fgElement.getName();
        if (freeElectron != null) {
        	int order = freeElectron.getOrder();
        	if (order == 0) name = name;
        	else if (order == 1) name = name + " mono-radical";
         	else if (order ==2) name = name + " bi-radical " + freeElectron.getSpin();
         	else if (order ==3) name = name + " tri-radical";
         	else throw new InvalidChemNodeElementException();
        }
        else throw new InvalidFreeElectronException();
        //#]
    }
    public  FGAtom() {
    }
    
    //## operation changeRadical(int,String) 
    public ChemNodeElement changeRadical(int p_radical, String p_spin) {
        //#[ operation changeRadical(int,String) 
        FGAtom newAtom = null;
        FreeElectron fe = getFreeElectron();
        if (fe == null) {
        	throw new InvalidFreeElectronException();
        }
        else {
        	int order = fe.getOrder()+p_radical;
        	if (order < 0) throw new InvalidFreeElectronException();
        	String name = String.valueOf(order);
        	if (order==2 && p_spin != null) {
        		name = name + p_spin;
        	}
        	FreeElectron newfe = FreeElectron.make(name);
        	newAtom = FGAtom.make(getFgElement(),newfe);
        }
        return newAtom;
        
        
        //#]
    }
    
    //## operation getFreeElectron() 
    public FreeElectron getFreeElectron() {
        //#[ operation getFreeElectron() 
        return freeElectron;
        //#]
    }
    
    //## operation getName() 
    public String getName() {
        //#[ operation getName() 
        return name;
        //#]
    }
    
    //## operation getType() 
    public String getType() {
        //#[ operation getType() 
        return fgElement.getType();
        //#]
    }
    
    //## operation isAny() 
    public boolean isAny() {
        //#[ operation isAny() 
        return (getFgElement().isAny());
        //#]
    }
    
    //## operation isCarbon() 
    public boolean isCarbon() {
        //#[ operation isCarbon() 
        return (getFgElement().isCarbon());
        //#]
    }
    
    //## operation isHydrogen() 
    public boolean isHydrogen() {
        //#[ operation isHydrogen() 
        return (getFgElement().isHydrogen());
        //#]
    }
    
    //## operation isNonH() 
    public boolean isNonH() {
        //#[ operation isNonH() 
        return (getFgElement().isNonH());
        //#]
    }
    
    //## operation isOxygen() 
    public boolean isOxygen() {
        //#[ operation isOxygen() 
        return (getFgElement().isOxygen());
        //#]
    }
    
    //## operation isRadical() 
    public boolean isRadical() {
        //#[ operation isRadical() 
        if (freeElectron == null) throw new InvalidFreeElectronException();
        if (freeElectron.getOrder() == 0) return false;
        else return true;
        
        
        //#]
    }
    
    public boolean isSilicon() {
        return (getFgElement().isSilicon());
    }
    
    public boolean isSulfur() {
        return (getFgElement().isSulfur());
    }
    
    //## operation make(FGElement,FreeElectron) 
    public static FGAtom make(FGElement p_fgElement, FreeElectron p_freeElectron) {
        //#[ operation make(FGElement,FreeElectron) 
        FGAtom fgatom = new FGAtom(p_fgElement, p_freeElectron);
        FGAtom old = (FGAtom)(dictionary.get(fgatom.getName()));
        if (old == null) {
        	dictionary.put(fgatom.getName(), fgatom);
        	return fgatom;
        }
        else {
        	fgatom = null;
        	return old;
        }
        
        
        //#]
    }
    
    //## operation setFgElement(FGElement) 
    public void setFgElement(FGElement p_fgElement) {
        //#[ operation setFgElement(FGElement) 
        fgElement = p_fgElement;
        //#]
    }
    
    public static HashMap getDictionary() {
        return dictionary;
    }
    
    public static void setDictionary(HashMap p_dictionary) {
        dictionary = p_dictionary;
    }
    
    public void setName(String p_name) {
        name = p_name;
    }
    
    public FGElement getFgElement() {
        return fgElement;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\FGAtom.java
*********************************************************************/

