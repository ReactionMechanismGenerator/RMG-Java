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



package jing.chemUtil;


import java.util.*;

//## package jing::chemUtil 

//----------------------------------------------------------------------------
// jing\chemUtil\DummyLeaf.java                                                                  
//----------------------------------------------------------------------------

//## class DummyLeaf 
/**
 * Stores the name and depth of the leafs which are difficult to futher classify and are names as "Others" in kinetics tree.
 */
public class DummyLeaf {
    
    private int depth;		//## attribute depth 
    
    private String name;		//## attribute name 
    
    
    // Constructors
    
    /**
     * Takes in a name and depth of the dummy leaf from the thermo and kinetics heirarchy trees. 
     * More specifically all Other nodes are included as Dummy leafs.
     */
    //## operation DummyLeaf(String,int) 
    public  DummyLeaf(String p_name, int p_depth) {
        //#[ operation DummyLeaf(String,int) 
        name = p_name;
        depth = p_depth;
        
        
        
        //#]
    }
    
    public  DummyLeaf() {
    }
    
   
    
    public int getDepth() {
        return depth;
    }
    
    public String getName() {
        return name;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chemUtil\DummyLeaf.java
*********************************************************************/

