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
        //#[ operation findMatchedPath(Matchable) 
        if (root == null) return null;
        
        Stack path = new Stack();
        ((HierarchyTreeNode)root).findMatchedPath(p_element,path);
        
        return path;
        
        
        
        //#]
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

