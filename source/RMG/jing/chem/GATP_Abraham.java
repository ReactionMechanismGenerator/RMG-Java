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


public class GATP_Abraham implements GeneralAbramGAPP {
	
    private static GATP_Abraham INSTANCE = new GATP_Abraham();
    protected ThermoGAGroupLibrary thermoLibrary;
	
    // Constructors
	
    private  GATP_Abraham() {
        initGAGroupLibrary();
    }
	
	public AbramData generateAbramData(ChemGraph p_chemGraph) {
		// Generation of Abraham Solute Parameters
		AbramData result_Abraham= new AbramData();
		// get the Group-Additive parameters (SBELA)
		result_Abraham.plus(getABGroup(p_chemGraph));
		// get the McGowan Volume
		result_Abraham.V = getMcGowanVolume(p_chemGraph);
		
        return result_Abraham;
    }
	
	public double getMcGowanVolume(ChemGraph p_chemGraph) {
		// Get the McGowan's Volume
		// Returned volumes are in cm^3/mol/100 (see note below)
		// See Table 2 in Abraham & McGowan, Chromatographia Vol. 23, No. 4, p. 243. April 1987
		// doi: 10.1007/BF02311772
		// 
		// "V is scaled to have similar values to the other 
		// descriptors by division by 100 and has units of (cm3mol−1/100)."
		// the contibutions in this function are in cm3/mol, and the division by 100 is done at the very end.
		
		double Vtot = 0;
		double thisV = 0;
		Iterator iter = p_chemGraph.getNodeList();
        while (iter.hasNext()) {
			thisV = 0;
          	Node node = (Node)iter.next();
          	Atom atom = (Atom)node.getElement();
			String element = (String)atom.getType();
          	if (element.equals("C")) thisV = thisV + 16.35;
			else if (element.equals("N")) thisV = thisV + 14.39;
			else if (element.equals("O")) thisV = thisV + 12.43;
			else if (element.equals("F")) thisV = thisV + 10.48;
			else if (element.equals("H")) thisV = thisV + 8.71;
			else if (element.equals("Si")) thisV = thisV + 26.83;
			else if (element.equals("P")) thisV = thisV + 24.87;
			else if (element.equals("S")) thisV = thisV + 22.91;
			else if (element.equals("Cl")) thisV = thisV + 20.95;
			else if (element.equals("B")) thisV = thisV + 18.32;
			else if (element.equals("Ge")) thisV = thisV + 31.02;
			else if (element.equals("As")) thisV = thisV + 29.42;
			else if (element.equals("Se")) thisV = thisV + 27.81;
			else if (element.equals("Br")) thisV = thisV + 26.21;
			else if (element.equals("Sn")) thisV = thisV + 39.35;
			else if (element.equals("Te")) thisV = thisV + 36.14;
			else if (element.equals("I")) thisV = thisV + 34.53;
			// else throw an error
			Vtot = Vtot + thisV;
		}
		
		iter = p_chemGraph.getArcList();
		while (iter.hasNext()) {
			Arc arc = (Arc)iter.next();
			Vtot = Vtot - 6.56;
		}
		return Vtot/100;  // division by 100 to get units correct.
	}
	
	public AbrahamGAValue getABGroup(ChemGraph p_chemGraph) {
        
        AbramData result_abram = new AbramData();
        Graph g = p_chemGraph.getGraph();
        HashMap oldCentralNode = (HashMap)(p_chemGraph.getCentralNode()).clone();
		
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
          			AbrahamGAValue thisAbrahamValue = thermoLibrary.findAbrahamGroup(p_chemGraph);
					
                    if (thisAbrahamValue == null) {
           				System.err.println("Abraham group not found: " + node.getID());
                        System.out.println(p_chemGraph.toString());
           			}
           			else {
           				//System.out.println(node.getID() + " " + thisGAValue.getName()+ "  "+thisGAValue.toString());
                        result_abram.plus(thisAbrahamValue);
           			}
           		}
           		else {
          			System.err.println("Error: Radical detected after satuation!");
           		}
           	}
        }
		
		//        // find the BDE for all radical groups
		//        for (int i=0; i<radicalSite; i++) {
		//          	int id = idArray[i];
		//           	Node node = g.getNodeAt(id);
		//           	Atom old = (Atom)node.getElement();
		//           	node.setElement(atomArray[i]);
		//           	node.updateFeElement();
		//
		//            // get rid of the extra H at ith site
		//          	int HNum = atomArray[i].getRadicalNumber();
		//           	for (int j=0;j<HNum;j++) {
		//           		g.removeNode(newnode[i][j]);
		//           	}
		//           	node.updateFgElement();
		//
		//           	p_chemGraph.resetThermoSite(node);
		//           	ThermoGAValue thisGAValue = thermoLibrary.findRadicalGroup(p_chemGraph);
		//           	if (thisGAValue == null) {
		//           		System.err.println("Radical group not found: " + node.getID());
		//           	}
		//           	else {
		//           		//System.out.println(node.getID() + " radical correction: " + thisGAValue.getName() + "  "+thisGAValue.toString());
		//           		result.plus(thisGAValue);
		//            }
		//
		//            //recover the satuated site for next radical site calculation
		//          	node.setElement(old);
		//          	node.updateFeElement();
		//           	for (int j=0;j<HNum;j++) {
		//           		newnode[i][j] = g.addNode(H);
		//           		g.addArcBetween(node,S,newnode[i][j]);
		//           	}
		//           	node.updateFgElement();
		//
		//         }
		//
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
		//
		//         // substract the enthalphy of H from the result
		//         int rad_number = p_chemGraph.getRadicalNumber();
		//         ThermoGAValue enthalpy_H = new ThermoGAValue(ENTHALPY_HYDROGEN * rad_number, 0,0,0,0,0,0,0,0,0,0,0,null);
		//         result.minus(enthalpy_H);
		//
		//         // make the symmetric number correction to entropy
		//
		//         if (p_chemGraph.isAcyclic()){
		//			 int sigma = p_chemGraph.getSymmetryNumber();
		//	         ThermoGAValue symmtryNumberCorrection = new ThermoGAValue(0,GasConstant.getCalMolK()*Math.log(sigma),0,0,0,0,0,0,0,0,0,0,null);
		//			 result.minus(symmtryNumberCorrection);
		//         }
		
		p_chemGraph.setCentralNode(oldCentralNode);
		
		// Abraham intercepts for each descriptor
        result_abram.S=result_abram.S+0.277;
        result_abram.E=result_abram.E+0.248;
        result_abram.A=result_abram.A+0.003;
        result_abram.L=result_abram.L+0.13;
        result_abram.B=result_abram.B+0.071;
		return result_abram;
    }
	
    protected void initGAGroupLibrary() {
        thermoLibrary = ThermoGAGroupLibrary.getINSTANCE();
    }
	
	protected static GATP_Abraham getINSTANCE() {
        return INSTANCE;
    }
}
