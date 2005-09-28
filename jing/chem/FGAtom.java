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

