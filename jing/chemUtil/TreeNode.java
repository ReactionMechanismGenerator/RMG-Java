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
// jing\chemUtil\TreeNode.java                                                                  
//----------------------------------------------------------------------------

//## class TreeNode 
/**
 * This class is used to model the nodes of a tree. Each TreeNode can have only one father and can 
 * have many children.
 */
public class TreeNode {
    
    protected Object element;		//## attribute element 
    
    protected LinkedHashSet children;
    protected TreeNode father;
    
    // Constructors
    
    //## operation TreeNode() 
    public  TreeNode() {
        {
            children=new LinkedHashSet();
        }
        //#[ operation TreeNode() 
        //#]
    }
    /**
    Requires:
    Effects: construct a new object by setting the element and childre as the pass-ins.
    Modifies:
    */
    //## operation TreeNode(Object,HashSet) 
    public  TreeNode(Object p_element, LinkedHashSet p_children) {
        {
            children=new LinkedHashSet();
        }
        //#[ operation TreeNode(Object,HashSet) 
        element = p_element;
        children = p_children;
        //#]
    }
    /**
    Requires:
    Effects:construct a new tree node by setting the pass-in object as the element of this tree node.
    Modifies:
    */
    //## operation TreeNode(Object) 
    public  TreeNode(Object p_element) {
        {
            children=new LinkedHashSet();
        }
        //#[ operation TreeNode(Object) 
        element = p_element;
        //#]
    }
    
    //## operation __setFather(TreeNode) 
    private void __setFather(TreeNode p_TreeNode) {
        //#[ operation __setFather(TreeNode) 
        father = p_TreeNode;
        //#]
    }
    
    //## operation _setFather(TreeNode) 
    private void _setFather(TreeNode p_TreeNode) {
        //#[ operation _setFather(TreeNode) 
        if(father != null)
            father._removeChildren(this);
        __setFather(p_TreeNode);
        //#]
    }
    
    /**
    Requires:
    Effects: remove all the children of this tree node.
    Modifies: this.children
    */
    //## operation clearChildren() 
    public void clearChildren() {
        //#[ operation clearChildren() 
        // call clear() of the collection of this node's children to remove all the children
        children.clear();
        //#]
    }
    
    /**
    Requries:
    Effects: return the iterator over this node's children's collection
    Modifies:
    */
    //## operation getChildren() 
    public Iterator getChildren() {
        //#[ operation getChildren() 
        Iterator iter = children.iterator();
        return iter;
        //#]
    }
    
    /**
    Requires:
    Effects: return the number of children of this tree node
    Modifies:
    */
    //## operation getChildrenNumber() 
    public int getChildrenNumber() {
        //#[ operation getChildrenNumber() 
        return children.size();
        //#]
    }
    
    /**
    Requires:
    Effects: accessor of the element stored in this tree node.
    Modifies:
    */
    //## operation getElement() 
    public Object getElement() {
        //#[ operation getElement() 
        return element;
        //#]
    }
    
    /**
    Requries:
    Effects: return the number of levels in the tree rooted at this tree node.
    Modifies:
    */
    //## operation height() 
    public int height() {
        //#[ operation height() 
        if (this == null) return 0;
        int h = 0, htemp;
        Iterator iter = getChildren();
        while (iter.hasNext()) {
        	htemp = ((TreeNode)iter.next()).height();
        	if (h < htemp) h = htemp;
        }
        
        return 1 + h;
        
        
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: return true if this tree node has no children
    Modifies:
    */
    //## operation isLeaf() 
    public boolean isLeaf() {
        //#[ operation isLeaf() 
        if (children == null) return true;
        else return (children.size() == 0);
        
        
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: remove the pass-in tree node from the collection of this node's children
    Modifies: this.children
    */
    //## operation removeChildren(TreeNode) 
    public void removeChildren(TreeNode p_TreeNode) {
        //#[ operation removeChildren(TreeNode) 
        // call the remove() of the children collection of this tree node
        children.remove(p_TreeNode);
        //#]
    }
    
    /**
    Requires:
    Effects: return true. (no check yet)
    Modifies:
    */
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        return true;
        //#]
    }
    
    /**
    Requires: 
    Effects: mutator of the element stored in this tree node
    Modifies: this.element.
    */
    //## operation setElement(Object) 
    protected void setElement(Object p_element) {
        //#[ operation setElement(Object) 
        // set the element as the pass-in object
        element = p_element;
        //#]
    }
    
    /**
    Requires:
    Effects: calculate the number of the total tree nodes in the tree rooted at this tree node. 
    Modifies:
    */
    //## operation size() 
    public int size() {
        //#[ operation size() 
        if (this == null) return 0;
        int s = 1;
        Iterator iter = getChildren();
        while (iter.hasNext()) {
        	s += ((TreeNode)iter.next()).size();
        }
        
        return s;
        
        
        
        
        //#]
    }
    
    public void _addChildren(TreeNode p_TreeNode) {
        children.add(p_TreeNode);
    }
    
    public void addChildren(TreeNode p_TreeNode) {
        if(p_TreeNode != null)
            p_TreeNode._setFather(this);
        _addChildren(p_TreeNode);
    }
    
    public void _removeChildren(TreeNode p_TreeNode) {
        children.remove(p_TreeNode);
    }
    
    public void _clearChildren() {
        children.clear();
    }
    
    public TreeNode getFather() {
        return father;
    }
    
    public void _clearFather() {
        father = null;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chemUtil\TreeNode.java
*********************************************************************/

