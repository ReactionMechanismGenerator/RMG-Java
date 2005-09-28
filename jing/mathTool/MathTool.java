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



package jing.mathTool;


import java.util.*;

//## package jing::mathTool 

//----------------------------------------------------------------------------
// jing\mathTool\MathTool.java                                                                  
//----------------------------------------------------------------------------

//## class MathTool 
public class MathTool {
    
    
    // Constructors
    
    //## operation MathTool() 
    private  MathTool() {
        //#[ operation MathTool() 
        //#]
    }
    
    //## operation expand(Iterator) 
    public static Collection expand(Iterator p_iter) {
        //#[ operation expand(Iterator) 
        if (!(p_iter.hasNext())) return new LinkedList();
        
        Object present = p_iter.next();
        Collection result = expand(p_iter);
        if (result.isEmpty()) {
        	if (present instanceof Collection) {
        		Iterator iter1 = ((Collection)present).iterator();
        		while (iter1.hasNext()) {
        			Object newItem = iter1.next();
        			LinkedList newSet = new LinkedList();
        			newSet.add(newItem);
        			result.add(newSet);
        		}
        	}
        	else {
        		Object newItem = present;
        		LinkedList newSet = new LinkedList();
        		newSet.add(newItem);
        		result.add(newSet);
        	}
        	return result;
        }
        
        if (present instanceof Collection) {
        	Iterator iter1 = ((Collection)present).iterator();
        	if (!(iter1.hasNext())) return result;
        	
        	LinkedList newResult = new LinkedList();
        	while (iter1.hasNext()) {
        		Object newItem = iter1.next();
        		Iterator iter2 = result.iterator();
        		while (iter2.hasNext()) {
        			LinkedList thisCollection = new LinkedList((Collection)iter2.next());
        			thisCollection.add(newItem);
        			newResult.add(thisCollection);
        		}
        	}
        	return newResult;
        }
        else {
        	Object newItem = present;
        	Iterator iter2 = result.iterator();
        	while (iter2.hasNext()) {
        		LinkedList thisCollection = new LinkedList((Collection)iter2.next());
        		thisCollection.add(newItem);
        	}
        	return result;
        }
        
        
        
        //#]
    }
    
    //## operation expandDisjointly(Iterator) 
    public static Collection expandDisjointly(Iterator p_iter) {
        //#[ operation expandDisjointly(Iterator) 
        if (!(p_iter.hasNext())) return new LinkedList();
        
        Object present = p_iter.next();
        Collection result = expandDisjointly(p_iter);
        if (result.isEmpty()) {
        	if (present instanceof Collection) {
        		Iterator iter1 = ((Collection)present).iterator();
        		while (iter1.hasNext()) {
        			Object newItem = iter1.next();
        			LinkedList newSet = new LinkedList();
        			newSet.add(newItem);
        			result.add(newSet);
        		}
        	}
        	else {
        		Object newItem = present;
        		LinkedList newSet = new LinkedList();
        		newSet.add(newItem);
        		result.add(newSet);
        	}
        	return result;
        }
        
        LinkedList newResult = new LinkedList();
        if (present instanceof Collection) {
        	Iterator iter1 = ((Collection)present).iterator();
        	if (!(iter1.hasNext())) return result;
        	
        	while (iter1.hasNext()) {
        		Object newItem = iter1.next();
        		Iterator iter2 = result.iterator();
        		while (iter2.hasNext()) {
        			LinkedList thisCollection = new LinkedList((Collection)iter2.next());
        			if (!thisCollection.contains(newItem)) {
        				thisCollection.add(newItem);
        				newResult.add(thisCollection);
        			}
        		}
        	}
        	return newResult;
        }
        else {
        	Object newItem = present;
        	Iterator iter2 = result.iterator();
        	while (iter2.hasNext()) {
        		LinkedList thisCollection = new LinkedList((Collection)iter2.next());
        		if (!thisCollection.contains(newItem)) {
        			thisCollection.add(newItem);
        			newResult.add(thisCollection);
        		}
        	}
        	return newResult;
        }
        
        
        
        //#]
    }
    
    //## operation formatDouble(double,int,int) 
    public static String formatDouble(double p_double, int p_length, int p_decimal) {
        //#[ operation formatDouble(double,int,int) 
        Double d = new Double (p_double);
        String s = d.toString();
        StringTokenizer st = new StringTokenizer(s,".");
        String result = st.nextToken();
        
        if (result.length()>p_length-p_decimal-1) throw new InvalidDoubleFormatException(s);
        result += ".";
        
        if (st.hasMoreTokens()) {
        	String dec = st.nextToken().trim();
        	if (dec.length() > p_decimal) dec = dec.substring(0,p_decimal-1);
        	result += dec;
        }
        else {
        	for (int i=0; i<p_decimal-1; i++) {
        		result += "0";
        	}
        }
        
        return result;
        //#]
    }
    
    //## operation formatInteger(int,int,String) 
    public static String formatInteger(int p_int, int p_length, String p_align) {
        //#[ operation formatInteger(int,int,String) 
        String result = Integer.toString(p_int);
        if (p_align.equals("L")) return result;
        else if (p_align.equals("R")) {
        	int space = p_length-result.length();
        	if (space < 0) throw new InvalidIntegerFormatException();
        	for (int i=0; i<space; i++) {
        		 result = " " + result;
        	}
        }
        else throw new InvalidIntegerFormatException("unknown align sign " + p_align);
        
        return result;
        	
        //#]
    }
    
    //## operation isCollectionDisjoint(Collection,Collection) 
    public static boolean isCollectionDisjoint(Collection p_c1, Collection p_c2) {
        //#[ operation isCollectionDisjoint(Collection,Collection) 
        if (p_c1 == null || p_c2 == null) throw new NullPointerException();
        
        for (Iterator iter = p_c1.iterator(); iter.hasNext(); ) {
        	Object o = iter.next();
        	if (p_c2.contains(o)) return false;
        }
        for (Iterator iter = p_c2.iterator(); iter.hasNext(); ) {
        	Object o = iter.next();
        	if (p_c1.contains(o)) return false;
        }
        
        return true;
        //#]
    }
    
    //## operation isListEquivalent(LinkedList,LinkedList) 
    public static boolean isListEquivalent(LinkedList p_list1, LinkedList p_list2) {
        //#[ operation isListEquivalent(LinkedList,LinkedList) 
        if (p_list1.size() != p_list2.size()) return false;
        
        boolean found = false;
        
        LinkedList templist = (LinkedList)(p_list2.clone());
        for (Iterator iter1 = p_list1.iterator(); iter1.hasNext();) {
        	Object o1 = iter1.next();
        	found = false;
        	for (Iterator iter2 = templist.iterator(); iter2.hasNext();) {
        		Object o2 = iter2.next();
        		if (o1.equals(o2)) { 
        			found = true;
        			iter2.remove();
        			break;
        		}
        	}
        	if (!found) return false;
        }
        
        if (templist.isEmpty()) return true;
        else return false;
        //#]
    }
    
    //## operation isSub(Object,Object) 
    public static boolean isSub(Object p_child, Object p_father) {
        //#[ operation isSub(Object,Object) 
        if (p_child == null || p_father == null) return false;
        if (p_child == p_father) return true;
        
        Object o1 = p_child;
        Object o2 = p_father;
        
        if ((o1 instanceof Collection) && (o2 instanceof Collection)) {
        	if (((Collection)o2).containsAll((Collection)o1)) {
        		return true;
        	}
        }
        // if only o2 is collection, we judge if o2 contains o1 as an element
        else if (!(o1 instanceof Collection) && (o2 instanceof Collection)) {
        	if (((Collection)o2).contains(o1)) {
        		return true;
        	}
        }
        // if o1 and o2 are both not collection, we judge if o1 isSub of o2
        else if (!(o1 instanceof Collection) && !(o2 instanceof Collection)) {
        	if (o1.equals(o2)) {
        		return true;
        	}
        }
        
        return false;
        //#]
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\mathTool\MathTool.java
*********************************************************************/

