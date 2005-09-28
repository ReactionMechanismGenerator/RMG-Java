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
// jing\chem\ChemElement.java                                                                  
//----------------------------------------------------------------------------

/**
OVERVIEW: Properties which all atoms of the same type share.
*/
//## class ChemElement 
public class ChemElement {
    
    protected static ChemElementDictionary chemElementDictionary = ChemElementDictionary.getInstance();		//## attribute chemElementDictionary 
    
    /**
    Name of the atom
    */
    protected String name;		//## attribute name 
    
    /**
    Number of protons
    */
    protected int proton;		//## attribute proton 
    
    /**
    Valency of the atom
    */
    protected double valency;		//## attribute valency 
    
    /**
    Atomic weight of the atom
    */
    protected double weight;		//## attribute weight 
    
    
    // Constructors
    
    //## operation ChemElement() 
    private  ChemElement() {
        //#[ operation ChemElement() 
        //#]
    }
    /**
    Requires:
    Effects: privite constructor called by create() to create new object
    Modifies:
    */
    //## operation ChemElement(String,int,double,double) 
    private  ChemElement(String p_name, int p_proton, double p_valency, double p_weight) {
        //#[ operation ChemElement(String,int,double,double) 
          name = p_name;
          proton = p_proton;
          valency = p_valency;
          weight = p_weight;
        
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: call privite constructor to construct new ChemElement according to p_name.  Now, here we consider to create C(or carbon), H(or hydrogen), O(or oxygen), if pass in other names, throw UnknownSymbolException.  Notice: user can specify an chemElement list file to read in by ChemElementDictionary.  But this create method guarantees our system able to at least include C,H,O.
    Modifies: ChemElementDictionary
    */
    //## operation create(String) 
    private static ChemElement create(String p_name) throws UnknownSymbolException {
        //#[ operation create(String) 
        ChemElement ChemElement = null;
        
        if (p_name.equals("C")) {
        	ChemElement = new ChemElement("C",12,4,12);
        }
        else if (p_name.equals("H")) {
        	ChemElement = new ChemElement("H",1,1,1);
        }
        else if (p_name.equals("O")) {
        	ChemElement = new ChemElement("O",16,2,32);
        }
        else {
        	throw new UnknownSymbolException("ChemElement");
        }
        return ChemElement;
        //#]
    }
    
    //## operation getName() 
    public String getName() {
        //#[ operation getName() 
        return name;
        //#]
    }
    
    //## operation getProton() 
    public int getProton() {
        //#[ operation getProton() 
        return proton;
        //#]
    }
    
    //## operation getValency() 
    public double getValency() {
        //#[ operation getValency() 
        return valency;
        //#]
    }
    
    //## operation getWeight() 
    public double getWeight() {
        //#[ operation getWeight() 
        return weight;
        //#]
    }
    
    /**
    Requires: pass in valid string name for elements
    Effects: search chemElementDictionary by p_name to see if this chemElement is already in the system.  if it is, return the instance of existing chemElement;  if it is not, call create(p_name) to create a new chemElement.
    Modifies: 
    */
    // Argument Stringp_name : 
    /**
    the name of chemical element you want to create.
    */
    //## operation make(String) 
    public static ChemElement make(String p_name) throws UnknownSymbolException {
        //#[ operation make(String) 
        try {                                                                                     
        	String internalName = translateName(p_name);
        	ChemElement ce = chemElementDictionary.getChemElement(internalName);
        	if (ce == null) {
        		ce = create(internalName);
        		chemElementDictionary.putChemElement(ce);
        	}
        	return ce;
        }
        catch (UnknownSymbolException e) {
        	throw new UnknownSymbolException("Chemical Element: " + p_name);
        }
        
        
        
        //#]
    }
    
    //## operation translateName(String) 
    public static String translateName(String p_name) throws NullSymbolException, UnknownSymbolException {
        //#[ operation translateName(String) 
        if (p_name == null) throw new NullSymbolException("ChemElement");
        
        if ((p_name.compareToIgnoreCase("C")==0) || (p_name.compareToIgnoreCase("Carbon")==0)) {
        	return "C";
        }
        else if ((p_name.compareToIgnoreCase("H")==0) || (p_name.compareToIgnoreCase("Hydrogen")==0)) {
        	return "H";
        }
        else if ((p_name.compareToIgnoreCase("O")==0) || (p_name.compareToIgnoreCase("Oxygen")==0)) {
        	return "O";
        }
        else {
        	throw new UnknownSymbolException("ChemElement");
        }
        //#]
    }
    
    protected static ChemElementDictionary getChemElementDictionary() {
        return chemElementDictionary;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\ChemElement.java
*********************************************************************/

