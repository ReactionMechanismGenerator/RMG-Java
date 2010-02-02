/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jing.chem;

import java.util.*;
import jing.param.*;
import jing.chemUtil.*;
import jing.param.*;
import jing.chemUtil.*;
import jing.mathTool.*;
import jing.param.Temperature;

/**
 *
 * @author User1
 */
public class GATP_Abraham implements GeneralAbramGAPP {

    private static GATP_Abraham INSTANCE = new GATP_Abraham();		//## attribute INSTANCE

    protected static HashMap library;		//## attribute library

    protected ThermoGAGroupLibrary thermoLibrary;

    protected static PrimaryThermoLibrary primaryLibrary;//svp

    // Constructors

    //## operation GATP()
    private  GATP_Abraham() {
        //#[ operation GATP()
        initializeLibrary();
        initGAGroupLibrary();
        initializePrimaryThermoLibrary();//svp
        //#]
    }

        public AbramData generateAbramData(ChemGraph p_chemGraph) {
       
             
       // Generation of Abraham Solute Parameters
       AbramData result_Abraham= new AbramData();
       result_Abraham.plus(getABGroup(p_chemGraph));


        return result_Abraham;
        //#]
    }



        public AbrahamGAValue getABGroup(ChemGraph p_chemGraph) {
        //#[ operation getGAGroup(ChemGraph)
        
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
//         // substrate the enthalphy of H from the result
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
        //#]
    }

            //## operation initGAGroupLibrary()
    protected void initGAGroupLibrary() {
        //#[ operation initGAGroupLibrary()
        thermoLibrary = ThermoGAGroupLibrary.getINSTANCE();
        //#]
    }

    public void initializeLibrary() {
        library = new HashMap();
        // put in H2
        ThermoData td_H2 = new ThermoData(0.000,31.233,6.895,6.975,6.994,7.009,7.081,7.219,7.720,0,0,0,"library value for H2");
        library.put("H2", td_H2);

        // put in H
        ThermoData td_H = new ThermoData(52.103,27.419,4.968,4.968,4.968,4.968,4.968,4.968,4.968, 0,0,0,"library value for H radical");
        library.put("H.",td_H);
		
    }


        protected static GATP_Abraham getINSTANCE() {
        return INSTANCE;
    }

    public static HashMap getLibrary() {
        return library;
    }

    public static void setLibrary(HashMap p_library) {
        library = p_library;
    }


      public void initializePrimaryThermoLibrary(){//svp

        primaryLibrary = PrimaryThermoLibrary.getINSTANCE();

      }


}
