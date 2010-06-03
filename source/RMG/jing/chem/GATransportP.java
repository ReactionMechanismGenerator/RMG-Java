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
import jing.chemUtil.*;

public class GATransportP {

    public static GATransportP INSTANCE = new GATransportP();
    protected static PrimaryTransportLibrary primaryLibrary;
    protected TransportGALibrary transportLibrary;

    // Constructors
    public GATransportP() {
        initGAGroupLibrary();
        initializePrimaryTransportLibrary();
    }

    public TransportData generateTransportData(ChemGraph p_chemGraph) {
    	TransportData trans_result = primaryLibrary.getTransportData(p_chemGraph.getGraph());
    	
    	// If the chemgraph was found in a primary transport library,
    	//		return the results
    	if (trans_result != null) {
    		//p_chemGraph.fromprimarytranslibrary = true;
    		return trans_result;
    	}
    	// Else, estimate the transport properties via group additivity
    	LJData lj_result = new LJData();
    	lj_result = getGAGroup(p_chemGraph);
    	trans_result = new TransportData(lj_result,false);

        return trans_result;
    }

    public LJData getGAGroup(ChemGraph p_chemGraph) {
        LJData result = new LJData();
        result.na = p_chemGraph.getAtomNumber();
    
        Graph g = p_chemGraph.getGraph();
        HashMap oldCentralNode = (HashMap)(p_chemGraph.getCentralNode()).clone();
        
        int na=p_chemGraph.getAtomNumber();//determine the number of atoms (before saturation)
        
        // saturate radical site
        int max_radNum_molecule = ChemGraph.getMAX_RADICAL_NUM();
        int max_radNum_atom = Math.min(8,max_radNum_molecule);
        int [] idArray = new int[max_radNum_molecule];
        Atom []  atomArray = new Atom[max_radNum_molecule];
        Node [][] newnode = new Node[max_radNum_molecule][max_radNum_atom];

        int radicalSite = 0;
        Iterator iter = p_chemGraph.getNodeList();
        FreeElectron satuated = FreeElectron.make("0");
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
           	Atom atom = (Atom)node.getElement();
           	if (atom.isRadical()) {
           		radicalSite ++;
           		// save the old radical atom
           		idArray[radicalSite-1] = node.getID().intValue();
           		atomArray[radicalSite-1] = atom;
           		// new a saturated atom and replace the old one
           		Atom newAtom = new Atom(atom.getChemElement(),satuated);
           		node.setElement(newAtom);
           		node.updateFeElement();
           	}
        }

        // add H to saturate chemgraph
        Atom H = Atom.make(ChemElement.make("H"),satuated);
        Bond S = Bond.make("S");
        for (int i=0;i<radicalSite;i++) {
           	Node node = p_chemGraph.getNodeAt(idArray[i]);
           	Atom atom = atomArray[i];
           	int HNum = atom.getRadicalNumber();
           	for (int j=0;j<HNum;j++) {
           		newnode[i][j] = g.addNode(H);
           		g.addArcBetween(node,S,newnode[i][j]);
           	}
           	node.updateFgElement();
        }

        // find all the transport groups
        iter = p_chemGraph.getNodeList();
        while (iter.hasNext()) {
          	Node node = (Node)iter.next();
          	Atom atom = (Atom)node.getElement();
          	if (!(atom.getType().equals("H"))) {
           		if (!atom.isRadical()) {
           			p_chemGraph.resetThermoSite(node);
           			LJGroupData thisGAValue = new LJGroupData();
           			if (node.getInCycle()){ //depending on whether the atom is in a cycle or not, use the appropriate library
           				thisGAValue = transportLibrary.findRingGroup(p_chemGraph);
           			}
           			else {
           				thisGAValue = transportLibrary.findGroup(p_chemGraph);
           			}
                    
           			if (thisGAValue == null) {
           				System.err.println("Transport group not found: " + node.getID());
           			}
           			else {
           				result.plus(thisGAValue);
           			}
           		}
           		else {
          			System.err.println("Error: Radical detected after satuation!");
           		}
           	}
        }
        
        // recover the chem graph structure
        // recover the radical
        for (int i=0; i<radicalSite; i++) {
        	int id = idArray[i];
        	Node node = g.getNodeAt(id);
        	node.setElement(atomArray[i]);
        	node.updateFeElement();
        	int HNum = atomArray[i].getRadicalNumber();
        	//get rid of extra H
        	for (int j=0;j<HNum;j++) {
        		g.removeNode(newnode[i][j]);
        	}
        	node.updateFgElement();
        }
       
        p_chemGraph.setCentralNode(oldCentralNode);
        return result;
    }

    protected void initGAGroupLibrary() {
        transportLibrary = TransportGALibrary.getINSTANCE();
    }

    public static GATransportP getINSTANCE() {
        return INSTANCE;
    }
    
    public void initializePrimaryTransportLibrary(){
    	HashMap ptlLibrary = PrimaryTransportLibrary.library;
    	HashMap ptlDictionary = PrimaryTransportLibrary.dictionary;
    	primaryLibrary = new PrimaryTransportLibrary(ptlDictionary, ptlLibrary);
    }

}