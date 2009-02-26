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


import java.io.*;
import jing.chem.*;
import java.util.*;
import jing.mathTool.*;
import jing.chemUtil.*;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\KineticsTemplateLibrary.java                                                                  
//----------------------------------------------------------------------------

//## class KineticsTemplateLibrary 
public class KineticsTemplateLibrary {
    
    protected LinkedHashMap kineticsTemplate;
    
    // Constructors
    
    public  KineticsTemplateLibrary() {
        {
            kineticsTemplate=new LinkedHashMap();
        }
    }
    
    //## operation addKinetics(HashSet,Kinetics) 
    public KineticsTemplate addKinetics(LinkedHashSet p_fgc, Kinetics p_kinetics) {
        //#[ operation addKinetics(HashSet,Kinetics) 
        KineticsTemplate old = getKineticsTemplate(p_fgc);
        // if there is already a number in the library, and it is not the root node, output information
        // otherwise, add the kt number into library, or replace the root number
        if (old != null && old.getKinetics().getRank() != 0) {
        	// if rate contant already exist, output a warning message, keep the original number, don't replace
        	String s = "";
        	for (Iterator iter = p_fgc.iterator(); iter.hasNext(); ) {
        		Object fg = iter.next();
        		if (fg instanceof String) {
        			s = s + fg;
        		}
        		else if (fg instanceof Matchable) {
        			s = s + ((Matchable)fg).getName();
        		}
        		else {
        			throw new InvalidKineticsKeyException();
        		}
        		if (iter.hasNext()) s = s + " + ";
        	}
					System.out.println(" Multiple values found for: " + s + ". Ignoring " + p_kinetics.getTRange() +"K data with A="+p_kinetics.getAValue());
        	return old;
        }
        else {
        	KineticsTemplate kt = new KineticsTemplate(p_fgc,p_kinetics);
        	addKineticsTemplate(kt);
        	return kt;
        }
        
        
        
        //#]
    }
    
    //## operation addKineticsTemplate(KineticsTemplate) 
    private Object addKineticsTemplate(KineticsTemplate p_kineticsTemplate) {
        //#[ operation addKineticsTemplate(KineticsTemplate) 
        return kineticsTemplate.put(p_kineticsTemplate.getKey(),p_kineticsTemplate);
        //#]
    }
    
    //## operation getKinetics(HashSet) 
    public Kinetics getKinetics(LinkedHashSet p_key) {
        //#[ operation getKinetics(HashSet) 
        KineticsTemplate kt = getKineticsTemplate(p_key);
        if (kt==null) return null;
        else return kt.getKinetics();
        //#]
    }
    
    //## operation getKineticsTemplate(HashSet) 
    public KineticsTemplate getKineticsTemplate(LinkedHashSet p_key) {
        //#[ operation getKineticsTemplate(HashSet) 
        return (KineticsTemplate)(kineticsTemplate.get(p_key));
        //#]
    }
    
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        // check if each kinetics is fine
        Iterator iter = getKineticsTemplate();
        while (iter.hasNext()) {
        	KineticsTemplate kt = getKineticsTemplate((LinkedHashSet)iter.next());
        	if (!kt.repOk()) return false;
        }
        return true;
        
        
        
        
        //#]
    }
    
    //## operation size() 
    public int size() {
        //#[ operation size() 
        return kineticsTemplate.size();
        //#]
    }
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        String s ="";
        int index = 0;
        
        Iterator key_iter = getKineticsTemplate();
        while (key_iter.hasNext()) {
        	index++;
        	LinkedHashSet key = (LinkedHashSet)key_iter.next();
        	KineticsTemplate kt = getKineticsTemplate(key);
        	s = s + kt.toString() + '\n';
        }
        
        return s;
        
        
        
        //#]
    }
    
    public Iterator getKineticsTemplate() {
        Iterator iter=kineticsTemplate.keySet().iterator();
        return iter;
    }
    
    public void clearKineticsTemplate() {
        kineticsTemplate.clear();
    }
    
    public void removeKineticsTemplate(KineticsTemplate p_KineticsTemplate) {
        Iterator iter=kineticsTemplate.keySet().iterator();
        while(iter.hasNext()) {
          Object key = iter.next();
          if (kineticsTemplate.get(key).equals(p_KineticsTemplate)) {
          	kineticsTemplate.remove(key);
          	break;
          }
        };
    }
    
    public void removeKineticsTemplate(LinkedHashSet key) {
        KineticsTemplate p_KineticsTemplate = getKineticsTemplate(key);
        kineticsTemplate.remove(key);
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\KineticsTemplateLibrary.java
*********************************************************************/

