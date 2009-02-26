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
    public static LinkedList expand(Iterator p_iter) {
        //#[ operation expand(Iterator) 
        if (!(p_iter.hasNext())) return new LinkedList();
        
        Object present = p_iter.next();
		LinkedList result = expand(p_iter);
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
    public static LinkedList expandDisjointly(Iterator p_iter) {
        //#[ operation expandDisjointly(Iterator) 
        if (!(p_iter.hasNext())) return new LinkedList();
        
        Object present = p_iter.next();
		LinkedList result = expandDisjointly(p_iter);
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
	
	 //## operation isListEquivalent(LinkedList,LinkedList) 
    public static boolean isListEqual(LinkedList p_list1, LinkedList p_list2) {
        //#[ operation isListEquivalent(LinkedList,LinkedList) 
        if (p_list1.size() != p_list2.size()) return false;
        
        boolean found = false;
        
        LinkedList templist = (LinkedList)(p_list2.clone());
        for (Iterator iter1 = p_list1.iterator(); iter1.hasNext();) {
        	Object o1 = iter1.next();
        	found = false;
        	for (Iterator iter2 = templist.iterator(); iter2.hasNext();) {
        		Object o2 = iter2.next();
        		if (o1==o2) { 
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
    
	  public static double log10(double p_data) {
	        //#[ operation log10(double) 
	        return Math.log(p_data)/Math.log(10.0);
	        //#]
	    }
}
/*********************************************************************
	File Path	: RMG\RMG\jing\mathTool\MathTool.java
*********************************************************************/

