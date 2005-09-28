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

//## package jing::chemUtil 

//----------------------------------------------------------------------------
// jing\chemUtil\Tree.java                                                                  
//----------------------------------------------------------------------------

//## class Tree 
public class Tree {
    
    protected TreeNode root;
    
    // Constructors
    
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
    
    //## operation getRoot() 
    public TreeNode getRoot() {
        //#[ operation getRoot() 
        return root;
        //#]
    }
    
    /**
    Requires:
    Effects: return hight of the tree
    Modifies:
    */
    //## operation height() 
    public int height() {
        //#[ operation height() 
        return root.height();
        
        
        
        //#]
    }
    
    /**
    Requires:
    Effects: return true iff the tree is empty, i.e., no root node
    Modifies:
    */
    //## operation isEmpty() 
    public boolean isEmpty() {
        //#[ operation isEmpty() 
        if (root == null) return true;
        else return false;
        
        
        
        //#]
    }
    
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
    Requires:
    Effects: return number of nodes
    Modifies:
    */
    //## operation size() 
    public int size() {
        //#[ operation size() 
        return root.size();
        
        
        
        //#]
    }
    
    public TreeNode newRoot() {
        root = new TreeNode();
        return root;
    }
    
    public void deleteRoot() {
        root=null;
    }
    
    protected void initRelations() {
        root = newRoot();
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chemUtil\Tree.java
*********************************************************************/

