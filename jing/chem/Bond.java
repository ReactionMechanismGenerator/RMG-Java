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
// jing\chem\Bond.java                                                                  
//----------------------------------------------------------------------------

//## class Bond 
public class Bond {
    
    protected static BondDictionary bondDictionary = BondDictionary.getInstance();		//## attribute bondDictionary 
    
    protected String location;		//## attribute location 
    
    protected String name;		//## attribute name 
    
    protected double order;		//## attribute order 
    
    protected int piElectrons;		//## attribute piElectrons 
    
    
    // Constructors
    
    //## operation Bond() 
    private  Bond() {
        //#[ operation Bond() 
        //#]
    }
    //## operation Bond(String,double,int) 
    private  Bond(String p_name, double p_order, int p_piElectrons) {
        //#[ operation Bond(String,double,int) 
        name = p_name;
        order = p_order;
        piElectrons = p_piElectrons;
        
        
        
        //#]
    }
    
    //## operation changeBond(int) 
    public Bond changeBond(int p_changedOrder) {
        //#[ operation changeBond(int) 
        if (isBenzene() && p_changedOrder == -1) {
          	//convertBenzene(arc,"D");
           	throw new InvalidBondMutationException("benezene bond changed");
        }
        int newOrder = (int)getOrder()+p_changedOrder;
        String newB = (new Integer(newOrder)).toString();
        Bond b = null;
        if (newB.equals("1")) 
        	b = Bond.make("S");
        else if (newB.equals("2")) 
        	b = Bond.make("D");
        else if (newB.equals("3")) 
        	b = Bond.make("T");
        else {
        	throw new InvalidBondMutationException("change bond to: " + newB);
        }
        
        return b;
        //#]
    }
    
    //## operation creat(String) 
    private static Bond creat(String p_name) throws UnknownSymbolException {
        //#[ operation creat(String) 
        Bond bond = null;
        
        if (p_name.equals("S")) {
        	bond = new Bond("S", 1, 0);
        }
        else if (p_name.equals("D")) {
        	bond = new Bond("D", 2, 2);   
        }
        else if (p_name.equals("Dcis")) {
        	bond = new Bond("Dcis", 2, 2);
        	bond.setLocation("cis");
        }
        else if (p_name.equals("Dtrans")) {
        	bond = new Bond("Dtrans", 2, 2);
        	bond.setLocation("trans");   
        }
        else if (p_name.equals("T")) {
        	bond = new Bond("T", 3, 2);
        }
        else if (p_name.equals("B")) {
        	bond = new Bond("B", 1.5, 1);
        }
        else {
        	throw new UnknownSymbolException("Bond");
        }
        
        return bond;
        
        
        
        //#]
    }
    
    //## operation getName() 
    public String getName() {
        //#[ operation getName() 
        return name;
        //#]
    }
    
    //## operation getOrder() 
    public double getOrder() {
        //#[ operation getOrder() 
        return order;
        //#]
    }
    
    //## operation getPiElectrons() 
    public int getPiElectrons() {
        //#[ operation getPiElectrons() 
        return piElectrons;
        //#]
    }
    
    //## operation isBenzene() 
    public boolean isBenzene() {
        //#[ operation isBenzene() 
        if (name != null) 
        	return name.equals("B");
        else 
        	return false;
        //#]
    }
    
    //## operation isDouble() 
    public boolean isDouble() {
        //#[ operation isDouble() 
        if (name != null)
        	return (name.equals("D") || name.equals("Dtrans") || name.equals("Dcis"));
        else 
        	return false;
        //#]
    }
    
    //## operation isSingle() 
    public boolean isSingle() {
        //#[ operation isSingle() 
        if (name != null) 
        	return name.equals("S");
        else 
        	return false;
        //#]
    }
    
    //## operation isTriple() 
    public boolean isTriple() {
        //#[ operation isTriple() 
        if (name != null) 
        	return name.equals("T");
        else 
        	return false;
        //#]
    }
    
    //## operation make(String) 
    public static Bond make(String p_name) throws UnknownSymbolException {
        //#[ operation make(String) 
        try {
        	String internalName = translateName(p_name);
        	Bond bond = bondDictionary.getBond(internalName);
        	if (bond == null) {
        		bond = creat(internalName);
        		bondDictionary.putBond(bond);
        	}
        	return bond;
        }
        catch (UnknownSymbolException e) {
        	throw new UnknownSymbolException("Bond: " + p_name);
        }
        //#]
    }
    
    //## operation translateName(String) 
    public static String translateName(String p_name) throws NullSymbolException, UnknownSymbolException {
        //#[ operation translateName(String) 
        if (p_name == null) throw new NullSymbolException("Bond");
        if (p_name.compareToIgnoreCase("S")==0 || p_name.compareToIgnoreCase("Single")==0) {
        	return "S";
        }
        else if (p_name.compareToIgnoreCase("D")==0 || p_name.compareToIgnoreCase("Double")==0) {
        	return "D";   
        }
        else if (p_name.compareToIgnoreCase("Dcis")==0 || p_name.compareToIgnoreCase("DoubleCis")==0) {
        	return "Dcis";
        }
        else if (p_name.compareToIgnoreCase("Dtrans")==0 || p_name.compareToIgnoreCase("DoubleTrans")==0) {
        	return "Dtrans";
        }
        else if (p_name.compareToIgnoreCase("T")==0 || p_name.compareToIgnoreCase("Triple")==0) {
        	return "T";
        }
        else if (p_name.compareToIgnoreCase("B")==0 || p_name.compareToIgnoreCase("Benzene")==0) {
        	return "B";
        }
        else {
        	throw new UnknownSymbolException("Bond");
        }
        
        
        //#]
    }
    
    private static BondDictionary getBondDictionary() {
        return bondDictionary;
    }
    
    public String getLocation() {
        return location;
    }
    
    public void setLocation(String p_location) {
        location = p_location;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\Bond.java
*********************************************************************/

