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
import jing.param.*;
import jing.chemUtil.*;
import jing.param.*;
import jing.chemUtil.*;
import jing.mathTool.*;
import jing.param.Temperature;

//## package jing::chem

//----------------------------------------------------------------------------
// jing\chem\GATP.java
//----------------------------------------------------------------------------

/**
Our new group additivity thermo property estimator based on the tree-structured libraries.
We have basically four libraries:
(1) group library: tree structure
(2) radical library: tree structure
(3) ring correction library: list
(4) other correction: tree structure
*/
//## class GATP
public class GATP implements GeneralGAPP {

    final protected static double ENTHALPY_HYDROGEN = 52.1;		//## attribute ENTHALPY_HYDROGEN

    private static GATP INSTANCE = new GATP();		//## attribute INSTANCE

    /**
    This is a thermo library for the species whose thermal data can't be estimated by Group Additivity Method.  For example, H2, H.
    */
//    protected static HashMap library;		//## attribute library

    protected ThermoGAGroupLibrary thermoLibrary;

     protected static PrimaryThermoLibrary primaryLibrary;//svp

    // Constructors


    private  GATP() {
     //   initializeLibrary();
        initGAGroupLibrary();
        initializePrimaryThermoLibrary();//svp
    }


    public ThermoData generateThermoData(ChemGraph p_chemGraph) {
        ThermoData result = primaryLibrary.getThermoData(p_chemGraph.getGraph());
        //System.out.println(result);
        if (result != null) {
        	p_chemGraph.fromprimarythermolibrary = true;
        	return result;
        }
	  // I think this is redundant: -rwest
      //  result = getFromLibrary(p_chemGraph.getChemicalFormula());
      //  if (result != null) return result;

        result = new ThermoData();

        result.plus(getGAGroup(p_chemGraph));

        // comment out, waiting for Bill and Joanna making the right ring correction library
        result.plus(getRingCorrection(p_chemGraph));
        result.plus(getOtherCorrection(p_chemGraph));

        return result;
    }

	/*
    //## operation getFromLibrary(String)
    public ThermoData getFromLibrary(String p_chemicalFormula) {
        //#[ operation getFromLibrary(String)
        return (ThermoData)library.get(p_chemicalFormula);
        //#]
    }
	*/

    /**
    Requires: pass-in ChemGraph object repOk() == true;
    Effects: calculate thermal data:
    (1) group values
    (2) radical correction
    (3) symmetric number and optical number correction to entropy
    Modifies:
    */
    //## operation getGAGroup(ChemGraph)
    public ThermoGAValue getGAGroup(ChemGraph p_chemGraph) {
        //#[ operation getGAGroup(ChemGraph)
        ThermoData result = new ThermoData();
    
        Graph g = p_chemGraph.getGraph();
        HashMap oldCentralNode = (HashMap)(p_chemGraph.getCentralNode()).clone();

        // satuate radical site
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
           		// new a satuated atom and replace the old one
           		Atom newAtom = new Atom(atom.getChemElement(),satuated);
           		node.setElement(newAtom);
           		node.updateFeElement();
           	}
        }

        // add H to satuate chem graph
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

        // find all the thermo groups
        iter = p_chemGraph.getNodeList();
        while (iter.hasNext()) {
          	Node node = (Node)iter.next();
          	Atom atom = (Atom)node.getElement();
          	if (!(atom.getType().equals("H"))) {
           		if (!atom.isRadical()) {

           			p_chemGraph.resetThermoSite(node);
           			ThermoGAValue thisGAValue = thermoLibrary.findGAGroup(p_chemGraph);
                    //System.out.println(thisGAValue);
 
                    //2/5/09 gmagoon: for acyclic molecules, check for and include gauche and 1,5 corrections (cyclic molecules should also have these corrections, but they appear to be more complicated and depend upon the actual 3D-structure; in general, it seems there are fewer gauche corrections for cyclic than acyclic (for the same degree of carbon atom))
                                if(p_chemGraph.isAcyclic())
                                {
                                    p_chemGraph.resetThermoSite(node);
                                    ThermoGAValue thisGaucheValue = thermoLibrary.findGaucheGroup(p_chemGraph);
                                    p_chemGraph.resetThermoSite(node);
                                    ThermoGAValue thisOneFiveValue = thermoLibrary.find15Group(p_chemGraph);
                                    if(thisGaucheValue!=null)
                                        result.plus(thisGaucheValue);
                                    if(thisOneFiveValue!=null)
                                        result.plus(thisOneFiveValue);
                                }
                                if (thisGAValue == null) {
           				System.err.println("Thermo group not found: " + node.getID());
           			}
           			else {
           				//System.out.println(node.getID() + " " + thisGAValue.getName()+ "  "+thisGAValue.toString());
           				result.plus(thisGAValue);

           			}
           		}
           		else {
          			System.err.println("Error: Radical detected after satuation!");
           		}
           	}
        }

        // find the BDE for all radical groups
        for (int i=0; i<radicalSite; i++) {
          	int id = idArray[i];
           	Node node = g.getNodeAt(id);
           	Atom old = (Atom)node.getElement();
           	node.setElement(atomArray[i]);
           	node.updateFeElement();

            // get rid of the extra H at ith site
          	int HNum = atomArray[i].getRadicalNumber();
           	for (int j=0;j<HNum;j++) {
           		g.removeNode(newnode[i][j]);
           	}
           	node.updateFgElement();

           	p_chemGraph.resetThermoSite(node);
           	ThermoGAValue thisGAValue = thermoLibrary.findRadicalGroup(p_chemGraph);
           	if (thisGAValue == null) {
           		System.err.println("Radical group not found: " + node.getID());
           	}
           	else {
           		//System.out.println(node.getID() + " radical correction: " + thisGAValue.getName() + "  "+thisGAValue.toString());
           		result.plus(thisGAValue);
            }

            //recover the satuated site for next radical site calculation
          	node.setElement(old);
          	node.updateFeElement();
           	for (int j=0;j<HNum;j++) {
           		newnode[i][j] = g.addNode(H);
           		g.addArcBetween(node,S,newnode[i][j]);
           	}
           	node.updateFgElement();

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

         // subtract the enthalphy of H from the result
         int rad_number = p_chemGraph.getRadicalNumber();
         ThermoGAValue enthalpy_H = new ThermoGAValue(ENTHALPY_HYDROGEN * rad_number, 0,0,0,0,0,0,0,0,0,0,0,null);
         result.minus(enthalpy_H);

         // make the symmetric number correction to entropy

         //if (p_chemGraph.isAcyclic()){//gmagoon 7/25/09: it looks like we have cyclic symmetry number code; we may want to try commenting the cyclic check out; if it turns out that this doesn't work, however, we may wish to revert
                int sigma = p_chemGraph.getSymmetryNumber();
	         ThermoGAValue symmtryNumberCorrection = new ThermoGAValue(0,GasConstant.getCalMolK()*Math.log(sigma),0,0,0,0,0,0,0,0,0,0,null);
                 result.minus(symmtryNumberCorrection);
         //}


         p_chemGraph.setCentralNode(oldCentralNode);
         //System.out.println(result);
         return result;




        //#]
    }



    //2/5/09 gmagoon: it seems the function below is broken and not working properly
    //## operation getOtherCorrection(ChemGraph)
    public ThermoGAValue getOtherCorrection(ChemGraph p_chemGraph) {
        //#[ operation getOtherCorrection(ChemGraph)
		HashMap oldCentralNode = (HashMap)(p_chemGraph.getCentralNode()).clone();
        ThermoGAValue ga = thermoLibrary.findOtherCorrection(p_chemGraph);
        //if (ga != null) System.out.println("Other Correction: " + ga.getName());
		p_chemGraph.setCentralNode(oldCentralNode);
        return ga;
        //#]
    }

    //## operation getRingCorrection(ChemGraph)
    public ThermoGAValue getRingCorrection(ChemGraph p_chemGraph) {
        //#[ operation getRingCorrection(ChemGraph)
		 if (p_chemGraph.isAcyclic()) return null;

		 	HashMap oldCentralNode = (HashMap)(p_chemGraph.getCentralNode()).clone();
	        ChemGraph sat = p_chemGraph;
	        if (sat.isRadical()) {
				sat = ChemGraph.saturate(p_chemGraph);
	        }
	        ThermoGAValue ga = thermoLibrary.findRingCorrection(sat);
	        /*System.out.println("Ring Correction for "+ p_chemGraph.generateChemicalFormula() +" : " + ga.getName());
	        System.out.println(p_chemGraph.toStringWithoutH());*/
			p_chemGraph.setCentralNode(oldCentralNode);
	        return ga;


        //#]
    }

    //## operation initGAGroupLibrary()
    protected void initGAGroupLibrary() {
        //#[ operation initGAGroupLibrary()
        thermoLibrary = ThermoGAGroupLibrary.getINSTANCE();
        //#]
    }

/*
    public void initializeLibrary() {
        library = new HashMap();
        // put in H2
        ThermoData td_H2 = new ThermoData(0.000,31.233,6.895,6.975,6.994,7.009,7.081,7.219,7.720,0,0,0,"library value for H2");
        library.put("H2", td_H2);

        // put in H
        ThermoData td_H = new ThermoData(52.103,27.419,4.968,4.968,4.968,4.968,4.968,4.968,4.968, 0,0,0,"library value for H radical");
        library.put("H.",td_H);

    }
 */


      public void initializePrimaryThermoLibrary(){//svp

//        primaryLibrary = PrimaryThermoLibrary.getINSTANCE();
    	  HashMap ptlLibrary = PrimaryThermoLibrary.library;
    	  HashMap ptlDictionary = PrimaryThermoLibrary.dictionary;
    	  primaryLibrary = new PrimaryThermoLibrary(ptlDictionary, ptlLibrary);

      }


    protected static GATP getINSTANCE() {
        return INSTANCE;
    }

	/*
    public static HashMap getLibrary() {
        return library;
    }

    public static void setLibrary(HashMap p_library) {
        library = p_library;
    }
	 */

}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\GATP.java
*********************************************************************/

