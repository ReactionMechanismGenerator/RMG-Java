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
import jing.mathTool.*;
import jing.param.Global;
import jing.param.Temperature;

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\RateConstant.java                                                                  
//----------------------------------------------------------------------------

//## class RateConstant 
public class RateConstant {
    
    protected String comments = "";		//## attribute comments 
    
    protected int distance = -1;		//## attribute distance 
    
    protected KineticsTemplate kineticsTemplate;
    protected HashSet kineticsTemplateFromLibrary;
    
    // Constructors
    
    //## operation RateConstant(KineticsTemplate,int) 
    public  RateConstant(KineticsTemplate p_kineticsTemplate, int p_distance) {
        {
            kineticsTemplateFromLibrary=new HashSet();
        }
        //#[ operation RateConstant(KineticsTemplate,int) 
        kineticsTemplate = p_kineticsTemplate;
        distance = p_distance;
        
        
        //#]
    }
    //## operation RateConstant(KineticsTemplate,HashSet,int) 
    public  RateConstant(KineticsTemplate p_kineticsTemplate, HashSet p_kineticsTemplateFromLibrary, int p_distance) {
        {
            kineticsTemplateFromLibrary=new HashSet();
        }
        //#[ operation RateConstant(KineticsTemplate,HashSet,int) 
        kineticsTemplate = p_kineticsTemplate;
        kineticsTemplateFromLibrary = p_kineticsTemplateFromLibrary;
        distance = p_distance;
        
        
        //#]
    }
    public  RateConstant() {
        {
            kineticsTemplateFromLibrary=new HashSet();
        }
    }
    
    //## operation calculateRate(Temperature,double) 
    public double calculateRate(Temperature p_temperature, double p_Hrxn) {
        //#[ operation calculateRate(Temperature,double) 
        return kineticsTemplate.getKinetics().calculateRate(p_temperature,p_Hrxn);
        //#]
    }
    
    //## operation getKinetics() 
    public Kinetics getKinetics() {
        //#[ operation getKinetics() 
        return kineticsTemplate.getKinetics();
        //#]
    }
    
    //## operation getRank() 
    public int getRank() {
        //#[ operation getRank() 
        return getKinetics().getRank();
        //#]
    }
    
    //## operation printKey() 
    public String printKey() {
        //#[ operation printKey() 
        return kineticsTemplate.printKey();
        
        
        //#]
    }
    
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        return kineticsTemplate.repOk();
        
        
        //#]
    }
    
    //## operation toChemkinString() 
    /*public String toChemkinString() {
        //#[ operation toChemkinString() 
        return getKineticsTemplate().getKinetics().toChemkinString(Global.temperature);
        //#]
    }*/
    
    //## operation toString() 
    public String toString() {
        //#[ operation toString() 
        String s = "Kinetics: " + kineticsTemplate.toString() + " Distance = " + String.valueOf(distance) + '\n';
        if (kineticsTemplateFromLibrary != null && kineticsTemplateFromLibrary.size()>0) {
        	s = s + "Use average of: " + '\n';
        	Iterator iter = kineticsTemplateFromLibrary.iterator();
        	while (iter.hasNext()) {
        		KineticsTemplate kt = (KineticsTemplate)iter.next();
        		s = s + kt.toString() + '\n';
        	}
        }
        
        return s;
        		
        //#]
    }
    
    public String getComments() {
        return comments;
    }
    
    public void setComments(String p_comments) {
        comments = p_comments;
    }
    
    public int getDistance() {
        return distance;
    }
    
    public KineticsTemplate getKineticsTemplate() {
        return kineticsTemplate;
    }
    
    public Iterator getKineticsTemplateFromLibrary() {
        Iterator iter=kineticsTemplateFromLibrary.iterator();
        return iter;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\RateConstant.java
*********************************************************************/

