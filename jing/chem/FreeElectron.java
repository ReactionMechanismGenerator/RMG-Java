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
// jing\chem\FreeElectron.java                                                                  
//----------------------------------------------------------------------------

//## class FreeElectron 
public class FreeElectron {
    
    protected static FreeElectronDictionary freeElectronDictionary = FreeElectronDictionary.getInstance();		//## attribute freeElectronDictionary 
    
    protected String name;		//## attribute name 
    
    protected int order = 0;		//## attribute order 
    
    protected String spin = null;		//## attribute spin 
    
    
    // Constructors
    
    //## operation FreeElectron() 
    private  FreeElectron() {
        //#[ operation FreeElectron() 
        //#]
    }
    //## operation FreeElectron(String,int,String) 
    private  FreeElectron(String p_name, int p_order, String p_spin) {
        //#[ operation FreeElectron(String,int,String) 
        name = p_name;
        order = p_order;
        spin = p_spin;
        
        
        
        //#]
    }
    //## operation FreeElectron(String,int) 
    private  FreeElectron(String p_name, int p_order) {
        //#[ operation FreeElectron(String,int) 
        name = p_name;
        order = p_order;
        //#]
    }
    
    //## operation creat(String) 
    private static FreeElectron creat(String p_name) throws UnknownSymbolException {
        //#[ operation creat(String) 
        FreeElectron electron = null;
        
        if (p_name.equals("0")) {
        	electron = new FreeElectron("0",0);
        }
        else if (p_name.equals("1")) {
        	electron = new FreeElectron("1",1);
        }
        // if spin in not specified then by default it is a triplet
        else if (p_name.equals("2")) {
        	electron = new FreeElectron("2",2);
        }
        else if (p_name.equals("2S")) {
        	electron = new FreeElectron("2S",2,"S");
        }
        else if (p_name.equals("2T")) {
        	electron = new FreeElectron("2T",2,"T");
        }
        else if (p_name.equals("3")) {
        	electron = new FreeElectron("3",3);
        }
        else if (p_name.equals("4")) {
			electron = new FreeElectron("4",4);
        }
        else {
        	throw new UnknownSymbolException("FreeElectron");
        }
        return electron;
        
        
        
        //#]
    }
    
    //## operation getFreeElectronDictionary() 
    protected FreeElectronDictionary getFreeElectronDictionary() {
        //#[ operation getFreeElectronDictionary() 
        return freeElectronDictionary;
        //#]
    }
    
    //## operation getOrder() 
    public int getOrder() {
        //#[ operation getOrder() 
        return order;
        //#]
    }
    
    //## operation getSpin() 
    public String getSpin() {
        //#[ operation getSpin() 
        return spin;
        //#]
    }
    
    //## operation make(String) 
    public static FreeElectron make(String p_name) throws UnknownSymbolException {
        //#[ operation make(String) 
        try {
        	String internalName = translateName(p_name);
        	FreeElectron electron = freeElectronDictionary.getFreeElectron(internalName);
        	if (electron == null) {
        		electron = creat(internalName);
        		freeElectronDictionary.putFreeElectron(electron);
        	}
        	return electron;
        }
        catch (UnknownSymbolException e) {
        	throw new UnknownSymbolException("Free Electron: " + p_name);
        }
        //#]
    }
    
    //## operation translateName(String) 
    public static String translateName(String p_name) throws NullSymbolException, UnknownSymbolException {
        //#[ operation translateName(String) 
        if (p_name == null) throw new NullSymbolException("FreeElectron");
        
        if (p_name.equals("0")) {
        	return "0";
        } 
        else if (p_name.equals("1") || (p_name.compareToIgnoreCase("Mono")==0) || (p_name.compareToIgnoreCase("MonoRadical")==0)) {
        	return "1";
        }
        else if ((p_name.equals("2")) || (p_name.compareToIgnoreCase("DiRadical")==0)) {
        	return "2";
        }
        else if ((p_name.compareToIgnoreCase("2T")==0) || (p_name.compareToIgnoreCase("Triplet")==0)) {
        	return "2T";
        }
        else if ((p_name.compareToIgnoreCase("2S")==0) || (p_name.compareToIgnoreCase("Singlet")==0)) {
        	return "2S";
        }
        else if ((p_name.equals("3")) || (p_name.compareToIgnoreCase("TriRadical")==0)) {
        	return "3";
        }
        else if ((p_name.equals("4")) || (p_name.compareToIgnoreCase("TetraRadical")==0)) {
			return "4";
        }
        else {
        	throw new UnknownSymbolException("FreeElectron");
        }
        //#]
    }
    
    public String getName() {
        return name;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\FreeElectron.java
*********************************************************************/

