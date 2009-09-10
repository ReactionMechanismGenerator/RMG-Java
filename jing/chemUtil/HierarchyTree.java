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
import jing.chem.Matchable;

//## package jing::chemUtil 

//----------------------------------------------------------------------------
// jing\chemUtil\HierarchyTree.java                                                                  
//----------------------------------------------------------------------------

//## class HierarchyTree 
public class HierarchyTree extends Tree {
    
    
    // Constructors
    
    //## operation HierarchyTree(HierarchyTreeNode) 
    public  HierarchyTree(HierarchyTreeNode p_root) {
        //#[ operation HierarchyTree(HierarchyTreeNode) 
        super(p_root);
        p_root.setDepth(0);
        //#]
    }
    public  HierarchyTree() {
    }
    
    /**
    Requires:
    Effects: return the best-matched tree leaf for the p_element.  i.e., the returned tree node is a leaf, and the element it contains will be the best match of p_element in the whole tree.
    Modifies:
    */
    //## operation findMatchedLeaf(Matchable) 
    public HierarchyTreeNode findMatchedLeaf(Matchable p_element) {
        //#[ operation findMatchedLeaf(Matchable) 
        if (root == null) return null;
        
        return ((HierarchyTreeNode)root).findMatchedLeaf(p_element);
        
        
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: return the best-matched path when searching the best-matched tree leaf for the p_element.  The returned stack will have the best-matched leaf at the top, and the tree root at the bottom.
    Modifies:
    */
    //## operation findMatchedPath(Matchable) 
    public Stack findMatchedPath(Matchable p_element) {
       // Returns the PATH
        if (root == null) return null;
        
        Stack path = new Stack();
        ((HierarchyTreeNode)root).findMatchedPath(p_element,path);
        
        return path;
		
    }
    
    /**
    Requires:
    Effects: check
    (1) if the root is a hierarchyTreeNode
    (2) if the hierarchy relations are all correct in the tree
    Modifies:
    */
    //## operation hierarchyOk() 
    public boolean hierarchyOk() {
        //#[ operation hierarchyOk() 
        if (!(root instanceof HierarchyTreeNode)) return false;
        return ((HierarchyTreeNode)root).hierarchyOk();
        //#]
    }
    
    /**
    Requires:
    Effects: Check the repOk() for the normal tree and the hierarchyOk() for hierarhcyTree
    Modifies:
    */
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        return (super.repOk() && hierarchyOk());
        //#]
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chemUtil\HierarchyTree.java
*********************************************************************/

