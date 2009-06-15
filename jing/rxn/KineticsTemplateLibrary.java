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
import jing.param.Temperature;
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
    
    /**
     * Added by MRH on 11-Jun-2009
     * 	This function assigns a set of kinetic parameters to a collection of
     * 		functional groups.
     * 
     * 	The function first checks if a set of kinetic parameters is already
     * 		associated with the collection of FG nodes.  If not, the current
     * 		Kinetics will be assigned to the current collection of FGs.  If
     * 		so, RMG attempts to assign the "best" kinetic parameters to the
     * 		collection.  The "best" is defined as the set of parameters whose
     * 		temperature range includes the system temperature and has the lowest
     * 		(i.e. most confident) rank.  In the event that the ranks are equal,
     * 		the first set of kinetic parameters is used (as default).
     * 
     *	If any of the collection of FGs have multiple kinetics, RMG displays
     *		this fact to the user.
     * 
     * @param p_fgc: Collection of functional groups, whose kinetics need assigning
     * @param p_kinetics: Kinetics under consideration
     * @param p_temp: Temperature of system
     * @return
     */
    public KineticsTemplate addKinetics(LinkedHashSet p_fgc, Kinetics p_kinetics, Temperature p_temp) {
        KineticsTemplate old = getKineticsTemplate(p_fgc);
        // if there is already a number in the library, and it is not the root node, output information
        // otherwise, add the kt number into library, or replace the root number
        if (old != null && old.getKinetics().getRank() != 0) {
        	// If new kinetics data has a rank of 0, we do not want to replace the
        	//	old data with the new data
        	if (p_kinetics.getRank() == 0) return old;
//			System.out.println(" Multiple values found for: " + constructNodeEntryString(p_fgc) + ". Ignoring " + p_kinetics.getTRange() +"K data with A="+p_kinetics.getAValue());
        	// Print to screen that multiple kinetic parameters have been found for
        	//	the same set of nodes
			System.out.println(" Multiple sets of Arrhenius parameters found for: " + constructNodeEntryString(p_fgc) + ".");
        	//	Check if the reaction system temperature falls within the current Kinetics valid temperature range
        	String tRangeOld = old.getKinetics().getTRange();
        	String tRangeNew = p_kinetics.getTRange();
        	boolean isOldWithinValidTRange = determineIfTempWithinValidTempRange(p_temp,tRangeOld);
        	boolean isNewWithinValidTRange = determineIfTempWithinValidTempRange(p_temp,tRangeNew);
        	//	Determine which set of Kinetics should be stored
        	if (isOldWithinValidTRange) {
            	//	If both OLD and NEW are valid, take the one with lowest (best) rank
        		if (isNewWithinValidTRange) {
        			int oldRank = old.getKinetics().getRank();
        			int newRank = p_kinetics.getRank();
        			if (newRank < oldRank) {
        				removeKineticsTemplate(old);
        	        	KineticsTemplate kt = new KineticsTemplate(p_fgc,p_kinetics);
        	        	addKineticsTemplate(kt);
        	        	System.out.println("   Replacing Kinetics for " + constructNodeEntryString(p_fgc) +
        	        			": From " + old.getKinetics().getTRange() + "K data with A=" + old.getKinetics().getAValue() +
        	        			" to " + p_kinetics.getTRange() + "K data with A=" + p_kinetics.getAValue());
        	        	return kt;
        			} else return old;
        		//	If OLD is valid but NEW is not, return OLD
        		} else {
        			return old;
        		}
        	} else {
        		//	If NEW is valid but OLD is not, return NEW
        		if (isNewWithinValidTRange) {
    				removeKineticsTemplate(old);
    	        	KineticsTemplate kt = new KineticsTemplate(p_fgc,p_kinetics);
    	        	addKineticsTemplate(kt);
    	        	System.out.println("   Replacing Kinetics for " + constructNodeEntryString(p_fgc) +
    	        			": From " + old.getKinetics().getTRange() + "K data with A=" + old.getKinetics().getAValue() +
    	        			" to " + p_kinetics.getTRange() + "K data with A=" + p_kinetics.getAValue());
    	        	return kt;
    	        //	If both OLD and NEW are not valid, take the one with lowest (best) rank
        		} else {
        			int oldRank = old.getKinetics().getRank();
        			int newRank = p_kinetics.getRank();
        			if (newRank < oldRank) {
        				removeKineticsTemplate(old);
        	        	KineticsTemplate kt = new KineticsTemplate(p_fgc,p_kinetics);
        	        	addKineticsTemplate(kt);
                		System.out.println("  System temperature " + p_temp.getK() + 
                    			"K not within valid Temperature Range (" + p_kinetics.getTRange() +
                    			"K) for node: " + constructNodeEntryString(p_fgc));
        	        	System.out.println("   Replacing Kinetics for " + constructNodeEntryString(p_fgc) +
        	        			": From " + old.getKinetics().getTRange() + "K data with A=" + old.getKinetics().getAValue() +
        	        			" to " + p_kinetics.getTRange() + "K data with A=" + p_kinetics.getAValue());
        	        	return kt;
        			} else return old;
        		}        		
        	}
        }
        else {
        	String tRange = p_kinetics.getTRange();
        	if (!determineIfTempWithinValidTempRange(p_temp,tRange)) {
        		System.out.println("  System temperature " + p_temp.getK() + 
        			"K not within valid Temperature Range (" + p_kinetics.getTRange() +
        			"K) for node: " + constructNodeEntryString(p_fgc));
        	}
        	KineticsTemplate kt = new KineticsTemplate(p_fgc,p_kinetics);
        	addKineticsTemplate(kt);
        	return kt;
        }
        
        
        
        //#]
    }
    
    /**
     * A glorified "toString" function
     * 	The collection of functional groups are written as a string.
     * 		The nodes are written, one at a time, with a " + "
     * 		separating each node.
     * 
     * @param p_fgc: Collection of functional groups
     * @return
     */
    public String constructNodeEntryString(LinkedHashSet p_fgc) {
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
    	return s;
    }
    
    /**
     * Function determines if the Temperature p_temp falls within the
     * 	temperature range tRange.
     * 
     * @param p_temp: Reaction System Temperature
     * @param tRange: Kinetics valid temperature range
     * @return boolean
     */
    public boolean determineIfTempWithinValidTempRange(Temperature p_temp, String tRange) {
    	/*	
    	 *	The temperature range should be "tMin-tMax" or "tOnly".
    	 * 		Split the tRange variable using the regular expression "-".
    	 * 		If the length of the split is 2, we have "tMin-tMax"
    	 * 		If the length is 1, we have "tOnly" 
    	*/
    	boolean isWithinRange = false;
    	String[] minMax = tRange.split("-");
    	double sysTemp = p_temp.getK();
    	if (minMax.length == 2) {
        	double tMin = Double.parseDouble(minMax[0]);
        	double tMax = Double.parseDouble(minMax[1]);
        	if (tMin <= sysTemp && tMax >= sysTemp) {
        		isWithinRange = true;
        	}
    	} else if (minMax.length == 1) {
    		double tOnly = Double.parseDouble(minMax[0]);
    		if (tOnly == sysTemp) {
    			isWithinRange = true;
    		}
    	}
    	return isWithinRange;
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

