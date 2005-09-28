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
// jing\chem\FunctionalGroupCollection.java                                                                  
//----------------------------------------------------------------------------

//## class FunctionalGroupCollection 
public class FunctionalGroupCollection implements Matchable {
    
    protected String name;		//## attribute name 
    
    protected HashSet functionalGroups;
    
    // Constructors
    
    //## operation FunctionalGroupCollection(String) 
    public  FunctionalGroupCollection(String p_name) {
        {
            functionalGroups=new HashSet();
        }
        //#[ operation FunctionalGroupCollection(String) 
        name = p_name;
        //#]
    }
    public  FunctionalGroupCollection() {
        {
            functionalGroups=new HashSet();
        }
    }
    
    //## operation addFunctionalGroups(FunctionalGroup) 
    public void addFunctionalGroups(FunctionalGroup p_FunctionalGroup) {
        //#[ operation addFunctionalGroups(FunctionalGroup) 
        functionalGroups.add(p_FunctionalGroup);
        //#]
    }
    
    //## operation addFunctionalGroups(Matchable) 
    public void addFunctionalGroups(Matchable p_functionalGroup) throws InvalidFunctionalGroupException {
        //#[ operation addFunctionalGroups(Matchable) 
        if (p_functionalGroup instanceof FunctionalGroup) {
        	FunctionalGroup fg = (FunctionalGroup)p_functionalGroup;
        	addFunctionalGroups(fg);
        }
        else if (p_functionalGroup instanceof FunctionalGroupCollection) {
        	FunctionalGroupCollection fgc = (FunctionalGroupCollection)p_functionalGroup;
        	addFunctionalGroups(fgc);
        }
        else throw new InvalidFunctionalGroupException();
        
        
        //#]
    }
    
    //## operation addFunctionalGroups(FunctionalGroupCollection) 
    public void addFunctionalGroups(FunctionalGroupCollection p_fgc) {
        //#[ operation addFunctionalGroups(FunctionalGroupCollection) 
        Iterator iter = p_fgc.getFunctionalGroups();
        while (iter.hasNext()) {
        	FunctionalGroup fg = (FunctionalGroup)iter.next();
        	addFunctionalGroups(fg);
        }
        
        
        //#]
    }
    
    //## operation contains(FunctionalGroup) 
    public boolean contains(FunctionalGroup p_fg) {
        //#[ operation contains(FunctionalGroup) 
        if (functionalGroups == null) return false;
        else return functionalGroups.contains(p_fg);
        //#]
    }
    
    //## operation isSubAtCentralNodes(Matchable) 
    public boolean isSubAtCentralNodes(Matchable p_functional) {
        //#[ operation isSubAtCentralNodes(Matchable) 
        if (this == p_functional) return false;
        if (!(p_functional instanceof FunctionalGroupCollection)) return false;
        
        Collection c1 = functionalGroups;
        Collection c2 = ((FunctionalGroupCollection)p_functional).functionalGroups;
        
        if (c2.size()==c1.size() && c2.containsAll(c1)) return false;
        
        boolean found = true;
        Iterator iter1 = c1.iterator();
        while (iter1.hasNext()) {
        	found = false;
        	FunctionalGroup fg1 = (FunctionalGroup)iter1.next();
        	Iterator iter2 = c2.iterator();
        	while (iter2.hasNext()) {
        		FunctionalGroup fg2 = (FunctionalGroup)iter2.next();
        		if (fg1.equals(fg2) || fg1.isSubAtCentralNodes(fg2)) {
        			found = true;
        			break;
        		}
        	}
        	if (!found) return false;
        }
        return true;
        
        
        
        //#]
    }
    
    public String getName() {
        return name;
    }
    
    public void setName(String p_name) {
        name = p_name;
    }
    
    public Iterator getFunctionalGroups() {
        Iterator iter=functionalGroups.iterator();
        return iter;
    }
    
    public void removeFunctionalGroups(FunctionalGroup p_FunctionalGroup) {
        functionalGroups.remove(p_FunctionalGroup);
    }
    
    public void clearFunctionalGroups() {
        functionalGroups.clear();
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\FunctionalGroupCollection.java
*********************************************************************/

