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



package jing.chemUtil;


import jing.chem.*;
import java.util.*;
import jing.chem.Matchable;

//## package jing::chemUtil 

//----------------------------------------------------------------------------
// jing\chemUtil\HierarchyTreeNode.java                                                                  
//----------------------------------------------------------------------------

/**
Hierarchy tree holds a structure that each child is a substruture of the father.  So, the element of HierarchyTreeNode should implement the interface matchable.
*/
//## class HierarchyTreeNode 
public class HierarchyTreeNode extends TreeNode {
    
    /**
    depth is the distance from this node to root.  
    */
    protected int depth = -1;		//## attribute depth 
    
    /**
    This is a dummy child of this hierachy tree node.  It has to be a leaf.  It means all the match pattern other than the ones represented by all the real children in the children set.  (It is used in the case that such representation can't be easily drawn.)  Dummy child cannot have children.
    */
    protected DummyLeaf dummyChild = null;		//## attribute dummyChild 
    
    
    // Constructors
    
    //## operation HierarchyTreeNode() 
    private  HierarchyTreeNode() {
        //#[ operation HierarchyTreeNode() 
        //#]
    }
    /**
    Requires:
    Effects: construct a new HierarchyTreeNode by setting it element as the pass-in Matchable object, and setting it children as the pass-in HashSet object, if all the tree node in the collection satisfies repOk() check and isSub(this) check.  If any of the children is not an instance of HierarchyTreeNode, throw new HierarchyTreeNodeExcption; if any of the children can't pass the repOk() or isSub(this) check, throw InvalidHierarchyRelationException.
    Modifies
    */
    //## operation HierarchyTreeNode(Matchable,int,HashSet) 
    public  HierarchyTreeNode(Matchable p_element, int p_depth, LinkedHashSet p_children) {
        //#[ operation HierarchyTreeNode(Matchable,int,HashSet) 
        super(p_element);
        
        depth = p_depth;
        
        Iterator iter = children.iterator();
        while (iter.hasNext()) {
        	try {
        		HierarchyTreeNode child = (HierarchyTreeNode)iter.next();
        		if (!child.isSubAtCentralNodes(this) || !child.repOk() || child.getDepth() != getDepth()+1) throw new InvalidHierarchyRelationException();
        	}
        	catch (ClassCastException e) {
        		throw new HierarchyTreeNodeException();
        	}
        }
        
        children = p_children;
        
        
        
        
        //#]
    }
    /**
    Requires:
    Effects: use super class's constructor.
    Modifies:
    */
    //## operation HierarchyTreeNode(Matchable,int) 
    public  HierarchyTreeNode(Matchable p_element, int p_depth) {
        //#[ operation HierarchyTreeNode(Matchable,int) 
        // use a matchable element to construct
        super(p_element);
        depth = p_depth;
        //#]
    }
    
    /**
    Requires:
    Effects: if p_child is not a substructure of this HierarchyTreeNode, throw InvalidHierarchyRelationException; otherwise, add p_child as the children of this node
    Modifes: this.children
    */
    //## operation addChildren(HierarchyTreeNode) 
    public void addChildren(HierarchyTreeNode p_child) throws HierarchyTreeNodeException, InvalidHierarchyRelationException {
        //#[ operation addChildren(HierarchyTreeNode) 
			if (!(p_child.isSubAtCentralNodes(this))) {
				String s = "Father " + ((Matchable)(this.element)).getName() + " Child " + ((Matchable)(p_child.element)).getName();
				System.out.println(s);
				//throw new InvalidHierarchyRelationException(s);
			}
        
        if (p_child.getDepth() != getDepth()+1) throw new InvalidTreeNodeLevelException(); 
        
        super.addChildren(p_child);
        
        
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: calculate how many level different from this node to p_treeNode.  return = p_treeNode.depth()-this.depth()
    Modifies:
    */
    //## operation calculateDistance(HierarchyTreeNode) 
    public int calculateDistance(HierarchyTreeNode p_treeNode) {
        //#[ operation calculateDistance(HierarchyTreeNode) 
        return -(depth - p_treeNode.getDepth());
        //#]
    }
    
    /**
    Requires:
    Effects: find out if there is a leaf match between the pass-in sturcture with node or any of the node's children that is a leaf.  if there is one, return it; otherwise, return null.
    Modifies:
    */
    //## operation findMatchedLeaf(Matchable) 
    public HierarchyTreeNode findMatchedLeaf(Matchable p_element) {
        //#[ operation findMatchedLeaf(Matchable) 
        if (p_element.isSubAtCentralNodes((Matchable)element)) {
        	//if there is a match and this node is a leaf, return this node;
         	if (isLeaf()) return this;
        
        	//if there is a match and this node is not a leaf, check its children recursively
        	Iterator iter = children.iterator();
        	while (iter.hasNext()) {
        		HierarchyTreeNode node = (HierarchyTreeNode)iter.next();
        		HierarchyTreeNode match = node.findMatchedLeaf(p_element);
        		if (match != null) return match;
        	}
        	// if all other children don't match, but has a dummy child, still return this node as a leaf
        	if (dummyChild != null) return this;
        }
        
        // if not matched here, return null
        return null;
        
        
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: find out a matched path beginning with this tree node.  add the matched tree nodes to p_path.  the return value indicates if the end of matched path is a leaf.  if it is, return true, otherwise, return false.
    Modifies: p_path
    */
    // Argument Matchablep_element : 
    /**
    the element we will match with the elements stroed in the tree.
    */
    // Argument Stackp_path : 
    /**
    This stack stores the tree nodes on the found path.
    */
    //## operation findMatchedPath(Matchable,Stack) 
    public boolean findMatchedPath(Matchable p_element, Stack p_path) {
        //#[ operation findMatchedPath(Matchable,Stack) 
        if (p_element.isSubAtCentralNodes((Matchable)element)) {
        	//if there is a match and this node is a leaf, add node to the p_path;
        	p_path.push(this);
         	if (isLeaf()) return true;
        
        	//if there is a match and this node is not a leaf, check its children recursively
        	for (Iterator iter = children.iterator(); iter.hasNext(); ) {
        		HierarchyTreeNode node = (HierarchyTreeNode)iter.next();
        		boolean match = node.findMatchedPath(p_element,p_path);
        		if (match) return true;
        	}
        	
        	// if all real children don't match, but there is a dummy child, still return true;
        	if (hasDummyChild()) {
        		return true;
        	}
        }
        // if there is no match, return false to the upper level
        return false;
        
        
        
        
        //#]
    }
    
    //## operation hasDummyChild() 
    public boolean hasDummyChild() {
        //#[ operation hasDummyChild() 
        return (dummyChild != null);
        //#]
    }
    
    /**
    Requires:
    Effects: if all of these checks are true, return true; otherwise, return false;
    (1) if the element is an instance of Matchable
    (2) if the tree rooted as this node is in the right hierarchical order, i.e, for every children node x of the father node y, it should be satisfied: (x.element).isSub(y.element) return true.
    Modifies: 
    */
    //## operation hierarchyOk() 
    public boolean hierarchyOk() {
        //#[ operation hierarchyOk() 
        if (!(element instanceof Matchable)) return false;
        
        if (!isLeaf()) {
        	// if has real children, check their hierarchy
        	Iterator iter = children.iterator();
        	while (iter.hasNext()) {
        		HierarchyTreeNode child = (HierarchyTreeNode)iter.next();
        		if (!child.isSubAtCentralNodes(this)) return false;
        		else if (!child.hierarchyOk()) return false;
        	}               
        }
        return true;
        //#]
    }
    
    /**
    Requires:
    Effects: return true iff this node has no children and no dummyChild.
    Modifies:
    */
    //## operation isLeaf() 
    public boolean isLeaf() {
        //#[ operation isLeaf() 
        return (super.isLeaf() && dummyChild==null);
        //#]
    }
    
    /**
    Requires:
    Effects: return true iff this node's element is a subgraph at center of the p_father's element
    Modifies:
    */
    //## operation isSubAtCentralNodes(HierarchyTreeNode) 
    public boolean isSubAtCentralNodes(HierarchyTreeNode p_father) {
        //#[ operation isSubAtCentralNodes(HierarchyTreeNode) 
        Matchable c1 = (Matchable)this.getElement();
        Matchable c2 = (Matchable)p_father.getElement();
        
        return (c1.isSubAtCentralNodes(c2));
        //#]
    }
    
    /**
    Requires:
    Effects: return true iff:
    (1) super.repOk() is true;
    (2) hierarchyOk() is true;
    Modifies:
    */
    //## operation repOk() 
    public boolean repOk() {
        //#[ operation repOk() 
        return (super.repOk() && hierarchyOk());
        //#]
    }
    
    public int getDepth() {
        return depth;
    }
    
    public void setDepth(int p_depth) {
        depth = p_depth;
    }
    
    public DummyLeaf getDummyChild() {
        return dummyChild;
    }
    
    public void setDummyChild(DummyLeaf p_dummyChild) {
        dummyChild = p_dummyChild;
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chemUtil\HierarchyTreeNode.java
*********************************************************************/

