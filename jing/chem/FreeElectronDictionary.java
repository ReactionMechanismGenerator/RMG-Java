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
// jing\chem\FreeElectronDictionary.java                                                                  
//----------------------------------------------------------------------------

//## class FreeElectronDictionary 
public class FreeElectronDictionary {
    
    private static FreeElectronDictionary INSTANCE = new FreeElectronDictionary();		//## attribute INSTANCE 
    
    protected HashMap dictionary;		//## attribute dictionary 
    
    
    // Constructors
    
    //## operation FreeElectronDictionary() 
    private  FreeElectronDictionary() {
        //#[ operation FreeElectronDictionary() 
        dictionary = new HashMap();
        //#]
    }
    
    //## operation getFreeElectron(String) 
    public FreeElectron getFreeElectron(String p_name) {
        //#[ operation getFreeElectron(String) 
        return (FreeElectron)(dictionary.get(p_name));
        //#]
    }
    
    //## operation getInstance() 
    protected static FreeElectronDictionary getInstance() {
        //#[ operation getInstance() 
        return INSTANCE;
        //#]
    }
    
    //## operation putFreeElectron(FreeElectron) 
    public void putFreeElectron(FreeElectron p_freeElectron) {
        //#[ operation putFreeElectron(FreeElectron) 
        dictionary.put(p_freeElectron.name, p_freeElectron);
        //#]
    }
    
    //## operation size() 
    public int size() {
        //#[ operation size() 
        return dictionary.size();
        //#]
    }
    
    public HashMap getDictionary() {
        return dictionary;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\FreeElectronDictionary.java
*********************************************************************/

