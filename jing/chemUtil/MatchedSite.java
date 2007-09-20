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
    
    protected HashMap periphery = new HashMap();		//## attribute periphery 
    
    protected int redundancy = 1;		//## attribute redundancy 
    
    
    // Constructors
    
    //## operation MatchedSite(Collection) 
    /*public  MatchedSite(Collection p_gcSet) {
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
        		periphery.add(gc);
        	}
        	else {
        		throw new InvalidGraphComponentException();
        	}
        }
        //#]
    }*/
    
    public  MatchedSite() {
    }
    
    //## operation addPeriphery(Node) 
    public boolean putPeriphery(GraphComponent gc, GraphComponent superGc) {
        //#[ operation addPeriphery(Node) 
        if (gc == null || superGc == null) throw new NullPointerException("Node");
        if (contains(gc)) return false;
        else {
        	periphery.put(superGc, gc);                    
        	return true;
        }
        //#]
    }
    
    //## operation contains(Node) 
    public boolean contains(GraphComponent p_gc) {
        //#[ operation contains(Node) 
        if (p_gc == null) throw new NullPointerException("Node");
        return (containsAsCenter(p_gc) || containsAsPeriphery(p_gc));
        //#]
    }
    
    /*public boolean contains(HashSet arcs, Node mainNode){
    	
    	Iterator iter = arcs.iterator();
    	while (iter.hasNext()){
    		Arc arc = (Arc)iter.next();
    		Iterator nodeIter = arc.getNeighbor();
    		while (nodeIter.hasNext()){
        		Node node = (Node)nodeIter.next();
        		if (!contains(node) && node != mainNode) return false;
        	}
    	}
    	
    	return true;
    }*/
    
    //## operation containsAsCenter(Node) 
    public boolean containsAsCenter(GraphComponent p_gc) {
        //#[ operation containsAsCenter(Node) 
        if (p_gc == null) throw new NullPointerException("Node");
        return center.containsValue(p_gc);
        //#]
    }
    
    //## operation containsAsPeriphery(Node) 
    public boolean containsAsPeriphery(GraphComponent p_gc) {
        //#[ operation containsAsPeriphery(Node) 
        if (p_gc == null) throw new NullPointerException("Node");
        return periphery.containsValue(p_gc);
        //#]
    }
    
    public boolean isSuper(Object p_matchedSite){
    	if (this == p_matchedSite) return false;
    	
    	if (!(p_matchedSite instanceof MatchedSite)) return false;
    	MatchedSite ms = (MatchedSite)p_matchedSite;
    	
    	if (this.center.keySet().size() < ms.center.keySet().size()) return false;
    	if (this.periphery.size() < ms.periphery.size()) return false;
    	
    	Iterator centerIter = ms.getCenter().keySet().iterator();
    	
    	while (centerIter.hasNext()){
    		Integer cid = (Integer)centerIter.next();
    		GraphComponent cg1 = (GraphComponent)getCenter().get(cid);
    		GraphComponent cg2 = (GraphComponent)(ms.getCenter().get(cid));
    		if (cg1 == null) return false;
    		if (cg1 != cg2)
    			return false;
    		
    		
    	}
    	
    	Iterator peripheryIter = ms.periphery.keySet().iterator();
    	while (peripheryIter.hasNext()){
    		GraphComponent keycg = (GraphComponent)peripheryIter.next();
    		GraphComponent cg1 = (GraphComponent)getPeriphery().get(keycg);
    		GraphComponent cg2 = (GraphComponent)(ms.getPeriphery().get(keycg));
    		if (cg1 == null) return false;
    		if (cg1 != cg2)
    			return false;
    	}
    	return true;
    }
    
    //## operation equals(Object) 
    public boolean equals(Object p_matchedSite) {
        //#[ operation equals(Object) 
        if (this == p_matchedSite) return true;
        
        if (!(p_matchedSite instanceof MatchedSite)) return false;
        
        MatchedSite ms = (MatchedSite)p_matchedSite;
        
        return (center.equals(ms.center) && periphery.equals(ms.periphery));
        
    	//return false;
        
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
        if (p_ms1.isSuper(p_ms2))
        	return p_ms1;
        HashMap c1 = p_ms1.getCenter();
        HashMap c2 = p_ms2.getCenter();
        HashMap cm = merged.getCenter();
        HashMap p1 = p_ms1.getPeriphery();
        HashMap p2 = p_ms2.getPeriphery();
        HashMap pm = merged.getPeriphery();
        
        cm.putAll(c1);
        pm.putAll(p1);
        
        for (Iterator cIter = c2.keySet().iterator(); cIter.hasNext();) {
        	Integer centralID = (Integer)cIter.next();
        	GraphComponent n = (GraphComponent)c2.get(centralID);
        	if (cm.containsKey(centralID) || cm.containsValue(n)) {
        		
        		return null;
        	}
        	else if (p1.containsValue(n) || p2.containsValue(n)) {
        		return null;
        	}
        	else {
        		cm.put(centralID,n);
        	}
        }
        
        for (Iterator pIter = p2.keySet().iterator(); pIter.hasNext();) {
        	GraphComponent key = (GraphComponent)pIter.next();
        	GraphComponent n = (GraphComponent)p2.get(key);
        	if (pm.containsKey(key) || cm.containsValue(n)) {
        		return null;
        	}
        	else if (c1.containsValue(n) || c2.containsValue(n)){
        		return null;
        	}
        	else {
        		pm.put(key, n);
        	}
        }
        
        return merged;
        
        
        
        
        
        //#]
    }
    
    //## operation putCenter(Integer,Node) 
    public boolean putCenter(Integer p_centralID, GraphComponent p_gc) {
        //#[ operation putCenter(Integer,Node) 
        if (p_gc == null) throw new NullPointerException("Node");
        if (contains(p_gc)) return false;
        else {
        	center.put(p_centralID, p_gc);
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
    
    public void setCenter(LinkedHashMap p_center) {
        center = p_center;
    }
    
	public String toString(){
		return center.toString();
	}
	
    public HashMap getPeriphery() {
        return periphery;
    }
    
    /*public void setPeriphery(LinkedHashSet p_periphery) {
        periphery = p_periphery;
    }*/
    
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

