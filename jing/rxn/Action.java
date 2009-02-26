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

//## package jing::rxn 

//----------------------------------------------------------------------------
// jing\rxn\Action.java                                                                  
//----------------------------------------------------------------------------

//## class Action 
public class Action {
    
    final protected static int BREAK_BOND = 1;		//## attribute BREAK_BOND 
    
    final protected static int CHANGE_BOND = 0;		//## attribute CHANGE_BOND 
    
    final protected static int FORM_BOND = 2;		//## attribute FORM_BOND 
    
    final protected static int GAIN_RADICAL = 3;		//## attribute GAIN_RADICAL 
    
    final protected static int LOSE_RADICAL = 4;		//## attribute LOSE_RADICAL 
    
    protected Object element;		//## attribute element 
    
    protected LinkedList site;		//## attribute site 
    
    protected int type;		//## attribute type 
    
    
    // Constructors
    
    //## operation Action(String,LinkedList,Object) 
    public  Action(String p_type, LinkedList p_site, Object p_element) throws InvalidActionException {
        //#[ operation Action(String,LinkedList,Object) 
        if (p_type == null) throw new InvalidActionException();
        
        if (p_type.compareToIgnoreCase("CHANGE_BOND")==0) {
        	type = CHANGE_BOND;
        }
        else if (p_type.compareToIgnoreCase("BREAK_BOND")==0) {
        	type = BREAK_BOND;
        }
        else if (p_type.compareToIgnoreCase("FORM_BOND")==0) {
        	type = FORM_BOND;
        }
        else if (p_type.compareToIgnoreCase("GAIN_RADICAL")==0) {
        	type = GAIN_RADICAL;
        }
        else if (p_type.compareToIgnoreCase("LOSE_RADICAL")==0) {
        	type = LOSE_RADICAL;
        }
        else {
        	throw new InvalidActionException();
        }
        
        site = p_site;
        element = p_element;
        
        
        //#]
    }
    public  Action() {
    }
    
    //## operation generateReverse() 
    public Action generateReverse() throws InvalidActionException {
        //#[ operation generateReverse() 
        Action reverse;
        int type = getType();
        switch (type) {
        	case Action.CHANGE_BOND:
        		{
        		int change = ((Integer)getElement()).intValue();
        		Integer newElement = new Integer(-change);
        		reverse = new Action("CHANGE_BOND",site,newElement);
        	    break;
        	    }
        	case Action.BREAK_BOND:
        		{
        		reverse = new Action("FORM_BOND",site,getElement());
        	    break;
        	    }
        	case Action.FORM_BOND:
        		{
        		reverse = new Action("BREAK_BOND",site,getElement());
        	    break;
        	    }
        	case Action.GAIN_RADICAL:
        		{
        		reverse = new Action("LOSE_RADICAL",site,getElement());
        	    break;
        	    }
        	case Action.LOSE_RADICAL:
        		{
        		reverse = new Action("GAIN_RADICAL",site,getElement());
        	    break;
        	    }
           	default:
           		throw new InvalidActionException("unknown action");
        }
        
        return reverse;
        
        
        //#]
    }
    
    //## operation getElement() 
    public Object getElement() {
        //#[ operation getElement() 
        return element;
        //#]
    }
    
    //## operation getSite() 
    public Iterator getSite() {
        //#[ operation getSite() 
        Iterator iter = site.iterator();
        return iter;
        //#]
    }
    
    public int getType() {
        return type;
    }
    
    public void setType(int p_type) {
        type = p_type;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\Action.java
*********************************************************************/

