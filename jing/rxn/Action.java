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

