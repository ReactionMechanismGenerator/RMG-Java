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



package jing.chemUtil;


import java.util.*;

//## package jing::chemUtil 

//----------------------------------------------------------------------------
// jing\chemUtil\MatchedSite.java                                                                  
//----------------------------------------------------------------------------

//## class MatchedSite 
public class MatchedSite {
    
    protected HashMap center = new HashMap();		//## attribute center 
    
    protected HashSet periphery = new HashSet();		//## attribute periphery 
    
    protected int redundancy = 1;		//## attribute redundancy 
    
    
    // Constructors
    
    //## operation MatchedSite(Collection) 
    public  MatchedSite(Collection p_gcSet) {
        //#[ operation MatchedSite(Collection) 
        for (Iterator iter = p_gcSet.iterator(); iter.hasNext() ; ) {
        	GraphComponent gc = (GraphComponent)iter.next();
        	if (gc instanceof Node) {
        		Node n = (Node)gc;
        		if (n.isCentralNode()) {
        			Integer cid = n.getCentralID();
        			center.put(cid, n);
        		}
        		else {
        			periphery.add(n);
        		}
        	}
        	else if (gc instanceof Arc) {
        	}
        	else {
        		throw new InvalidGraphComponentException();
        	}
        }
        //#]
    }
    public  MatchedSite() {
    }
    
    //## operation addPeriphery(Node) 
    public boolean addPeriphery(Node p_node) {
        //#[ operation addPeriphery(Node) 
        if (p_node == null) throw new NullPointerException("Node");
        if (contains(p_node)) return false;
        else {
        	periphery.add(p_node);                    
        	return true;
        }
        //#]
    }
    
    //## operation contains(Node) 
    public boolean contains(Node p_node) {
        //#[ operation contains(Node) 
        if (p_node == null) throw new NullPointerException("Node");
        return (containsAsCenter(p_node) || containsAsPeriphery(p_node));
        //#]
    }
    
    //## operation containsAsCenter(Node) 
    public boolean containsAsCenter(Node p_node) {
        //#[ operation containsAsCenter(Node) 
        if (p_node == null) throw new NullPointerException("Node");
        return center.containsValue(p_node);
        //#]
    }
    
    //## operation containsAsPeriphery(Node) 
    public boolean containsAsPeriphery(Node p_node) {
        //#[ operation containsAsPeriphery(Node) 
        if (p_node == null) throw new NullPointerException("Node");
        return periphery.contains(p_node);
        //#]
    }
    
    //## operation equals(Object) 
    public boolean equals(Object p_matchedSite) {
        //#[ operation equals(Object) 
        if (this == p_matchedSite) return true;
        
        if (!(p_matchedSite instanceof MatchedSite)) return false;
        
        MatchedSite ms = (MatchedSite)p_matchedSite;
        
        return (center.equals(ms.center) && periphery.equals(ms.periphery));
        
        
        //#]
    }
    
    //## operation hashCode() 
    public int hashCode() {
        //#[ operation hashCode() 
        return center.hashCode() + periphery.hashCode();
        //#]
    }
    
    //## operation merge(MatchedSite,MatchedSite) 
    public static MatchedSite merge(MatchedSite p_ms1, MatchedSite p_ms2) {
        //#[ operation merge(MatchedSite,MatchedSite) 
        MatchedSite merged = new MatchedSite();
        
        HashMap c1 = p_ms1.getCenter();
        HashMap c2 = p_ms2.getCenter();
        HashMap cm = merged.getCenter();
        HashSet p1 = p_ms1.getPeriphery();
        HashSet p2 = p_ms2.getPeriphery();
        HashSet pm = merged.getPeriphery();
        
        cm.putAll(c1);
        pm.addAll(p1);
        
        for (Iterator cIter = c2.keySet().iterator(); cIter.hasNext();) {
        	Integer centralID = (Integer)cIter.next();
        	Node n = (Node)c2.get(centralID);
        	if (cm.containsKey(centralID) || cm.containsValue(n)) {
        		return null;
        	}
        	else if (p1.contains(n) || p2.contains(n)) {
        		return null;
        	}
        	else {
        		cm.put(centralID,n);
        	}
        }
        
        for (Iterator pIter = p2.iterator(); pIter.hasNext();) {
        	Node n = (Node)pIter.next();
        	if (pm.contains(n) || cm.containsValue(n)) {
        		return null;
        	}
        	else {
        		pm.add(n);
        	}
        }
        
        return merged;
        
        
        
        
        
        //#]
    }
    
    //## operation putCenter(Integer,Node) 
    public boolean putCenter(Integer p_centralID, Node p_node) {
        //#[ operation putCenter(Integer,Node) 
        if (p_node == null) throw new NullPointerException("Node");
        if (contains(p_node)) return false;
        else {
        	center.put(p_centralID, p_node);
        	return true;
        }
        
        
        //#]
    }
    
    //## operation union(Collection) 
    public static MatchedSite union(Collection p_matchedSites) {
        //#[ operation union(Collection) 
        if (p_matchedSites.isEmpty()) return null;
        
        try {
        	Iterator iter = p_matchedSites.iterator();
        	MatchedSite result = (MatchedSite)iter.next();
        	while (iter.hasNext()) {
        		MatchedSite ms = (MatchedSite)iter.next();
        		result = merge(result, ms);
        		if (result==null) return null;
        	}
        	
        	return result;
        }
        catch (ClassCastException e) {
        	return null;
        }
        
        
        //#]
    }
    
    public HashMap getCenter() {
        return center;
    }
    
    public void setCenter(HashMap p_center) {
        center = p_center;
    }
    
    public HashSet getPeriphery() {
        return periphery;
    }
    
    public void setPeriphery(HashSet p_periphery) {
        periphery = p_periphery;
    }
    
    public int getRedundancy() {
        return redundancy;
    }
    
    public void setRedundancy(int p_redundancy) {
        redundancy = p_redundancy;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chemUtil\MatchedSite.java
*********************************************************************/

