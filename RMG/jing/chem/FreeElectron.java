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

