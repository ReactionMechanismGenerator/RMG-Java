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


import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import jing.chemUtil.Graph;
import jing.rxnSys.CoreEdgeReactionModel;

//## package jing::chem 

//----------------------------------------------------------------------------
//jing\chem\SpeciesDictionary.java                                                                  
//----------------------------------------------------------------------------

//## class SpeciesDictionary 
public class ChemGraphDictionary {
  
  private static ChemGraphDictionary INSTANCE = new ChemGraphDictionary();		//## attribute INSTANCE 
  
  private static HashMap dictionary;		//## attribute dictionary 
  
  
  // Constructors
  
  //## operation SpeciesDictionary() 
  private  ChemGraphDictionary() {
      //#[ operation SpeciesDictionary() 
      dictionary = new HashMap();
      //#]
  }
  
  //## operation getInstance() 
  public static ChemGraphDictionary getInstance() {
      //#[ operation getInstance() 
      return INSTANCE;
      //#]
  }
  

	public static ChemGraph getChemGraphFromGraph(Graph p_g) {
		if (p_g == null) throw new NullPointerException();
		Iterator iter = dictionary.keySet().iterator();
		while(iter.hasNext()) {
			Graph g = (Graph)iter.next();
			if (g.isEquivalent(p_g)) {
				return (ChemGraph)dictionary.get(g);
			}
		}
		return null;
	}
  
  
  //## operation putSpecies(Species) 
  public void putSpecies(ChemGraph p_cg) {
      //#[ operation putSpecies(Species) 
      
      	dictionary.put(p_cg.getGraph(),p_cg);
      
      return;
      
      
      //#]
  }
  
  //## operation remove(ChemGraph) 
  public void remove(Graph p_Graph) {
      //#[ operation remove(ChemGraph) 
      if (p_Graph != null) dictionary.remove(p_Graph);
      //#]
  }
  
  //## operation size() 
  public int size() {
      //#[ operation size() 
      return dictionary.size();
      //#]
  }
  
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\SpeciesDictionary.java
*********************************************************************/

