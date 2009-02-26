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


import jing.chem.*;
import java.util.*;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\KineticsTemplate.java                                                                  
//----------------------------------------------------------------------------

//## class KineticsTemplate 
public class KineticsTemplate {
    
    protected HashSet key;		//## attribute key 
    
    protected Kinetics kinetics;
    
    // Constructors
    
    //## operation KineticsTemplate(HashSet,Kinetics) 
    protected  KineticsTemplate(HashSet p_key, Kinetics p_kinetics) {
        //#[ operation KineticsTemplate(HashSet,Kinetics) 
        key = p_key;
        kinetics = p_kinetics;
        
        
        
        //#]
    }
    public  KineticsTemplate() {
    }
    
    //## operation printKey() 
    public String printKey() throws InvalidKineticsTemplateException, InvalidFunctionalGroupException {
        //#[ operation printKey() 
        String s = "";
        for (Iterator iter = key.iterator(); iter.hasNext(); ) {
        	Object fg = iter.next();
        	if (fg instanceof String) {
        		s = s + (String)fg;
        	}
        	else if (fg instanceof Matchable) {
        		s = s + ((Matchable)fg).getName();
        	}
        	else {
        		throw new InvalidKineticsKeyException();
        	}
        	if (iter.hasNext()) s = s + " + ";
        }
        if (s.equals("")) throw new InvalidKineticsTemplateException();
        
        return s;
        	
        	
        //#]
    }
    
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        return kinetics.repOk();
        //#]
    }
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        return "Use kinetics for: "+ printKey() + '\n' + kinetics.toString();
        //#]
    }
    
    public HashSet getKey() {
        return key;
    }
    
    public void setKey(HashSet p_key) {
        key = p_key;
    }
    
    public Kinetics getKinetics() {
        return kinetics;
    }
    
    public void setKinetics(Kinetics p_Kinetics) {
        kinetics = p_Kinetics;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\KineticsTemplate.java
*********************************************************************/

