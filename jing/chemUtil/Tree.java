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
// jing\chemUtil\Tree.java                                                                  
//----------------------------------------------------------------------------

//## class Tree 
/**
 * Class used to implement a tree. This class tree will be extended to form the heirarchy tree which 
 * is used to implement all the thermo and kinetics trees. 
 */
public class Tree {
    
    protected TreeNode root;
    
    // Constructors
    
    /**
     * Initializes the tree by taking the root of the tree.
     */
    //## operation Tree(TreeNode) 
    public  Tree(TreeNode p_root) {
        initRelations();
        //#[ operation Tree(TreeNode) 
        root = p_root;
        
        
        
        //#]
    }
    public  Tree() {
        initRelations();
    }
    
    /**
     * 
     * @return The root of the tree
     */
    //## operation getRoot() 
    public TreeNode getRoot() {
        //#[ operation getRoot() 
        return root;
        //#]
    }
    
    /**
    @return hight of the tree
    
    */
    //## operation height() 
    public int height() {
        //#[ operation height() 
        return root.height();
        
        
        
        //#]
    }
    
    /**
    Return true iff the tree is empty, i.e., no root node
    
    */
    //## operation isEmpty() 
    public boolean isEmpty() {
        //#[ operation isEmpty() 
        if (root == null) return true;
        else return false;
        
        
        
        //#]
    }
    
    /**
     * Calls the repOK function of the root TreeNode.
     * 
     */
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        return root.repOk();
        //#]
    }
    
    //## operation setRoot(TreeNode) 
    public void setRoot(TreeNode p_root) {
        //#[ operation setRoot(TreeNode) 
        root = p_root;
        //#]
    }
    
    /**
    Return number of nodes
    
    */
    //## operation size() 
    public int size() {
        //#[ operation size() 
        return root.size();
   
    }
    
    
    private TreeNode newRoot() {
        root = new TreeNode();
        return root;
    }
    
    public void deleteRoot() {
        root=null;
    }
    
    /**
     * Initializes the root of the tree.
     *
     */
    protected void initRelations() {
        root = newRoot();
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chemUtil\Tree.java
*********************************************************************/

