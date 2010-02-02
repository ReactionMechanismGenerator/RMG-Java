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

