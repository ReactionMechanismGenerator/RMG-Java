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
public class GATP_Solvation implements GeneralSolvationGAPP {

    private static GATP_Solvation INSTANCE = new GATP_Solvation();		//## attribute INSTANCE

    protected static HashMap library;		//## attribute library

    protected ThermoGAGroupLibrary thermoLibrary;

    protected static PrimaryThermoLibrary primaryLibrary;//svp

    // Constructors


    private  GATP_Solvation() {

        //#[ operation GATP()
        //initializeLibrary();
        //initGAGroupLibrary();
        //initializePrimaryThermoLibrary();//svp
        //#]

        initializeLibrary();
        initGAGroupLibrary();
        initializePrimaryThermoLibrary();//svp

    }

        public ThermoData generateSolvThermoData(ChemGraph p_chemGraph) {


       double r_solute=p_chemGraph.getRadius();
       double r_solvent; r_solvent=3.498*Math.pow(10,-10);// 3.311;                    // Manually assigned solvent radius [=] meter Calculated using Connolly solvent excluded volume from Chem3dPro
       double r_cavity=r_solute+r_solvent;                    // Cavity radius [=] Angstrom
       double rho; rho=0.00309*Math.pow(10,30);  //0.00381;                             // number density of solvent [=] molecules/Angstrom^3   Value here is for decane using density =0.73 g/cm3
       double parameter_y=(88/21)*rho*Math.pow(r_solvent, 3);  // Parameter y from Ashcraft Thesis Refer pg no. 60
       double parameter_ymod=parameter_y/(1-parameter_y);     // parameter_ymod= y/(1-y) Defined for convenience
       double R=8.314;                                        // Gas constant units J/mol K
       double T=298;                                          // Standard state temperature

       // Definitions of K0, K1 and K2 correspond to those for K0', K1' and K2' respectively from Ashcraft's Thesis
       double K0= -R*(-Math.log(1-parameter_y)+(4.5*Math.pow(parameter_ymod, 2)));
       double K1= (R*0.5/r_solvent)*((6*parameter_ymod)+(18*Math.pow(parameter_ymod, 2)));
       double K2= -(R*0.25/Math.pow(r_solvent,2))*((12*parameter_ymod)+(18*Math.pow(parameter_ymod, 2)));

       // Basic definition of entropy change of solvation from Ashcfrat's Thesis
       double deltaS0;
       deltaS0=K0+(K1*r_cavity)+(K2*Math.pow(r_cavity,2));


       // Generation of Abraham Solute Parameters
       AbramData result_Abraham= new AbramData();
       result_Abraham= p_chemGraph.getAbramData();


       // Solute descriptors from the Abraham Model
       double S=result_Abraham.S;
       double B=result_Abraham.B;
       double E=result_Abraham.E;
       double L=result_Abraham.L;
       double A=result_Abraham.A;

        //Manually specified solvent descriptors (constants here are for decane)
        double c=0.156;         //-0.12;
        double s=0;             //0.56;
        double b=0;             //0.7;
        double e=-0.143;        //-0.2;
        double l=0.989;         //0.94;
        double a=0;             //3.56;

        double logK=c + s*S + b*B + e*E + l*L + a*A;    // Implementation of Abraham Model
        double deltaG0_octanol=-8.314*298*logK;
 //       System.out.println("The free energy of solvation in octanol at 298K w/o reference state corrections  = " + deltaG0_octanol +" J/mol for " );

        // Calculation of enthalpy change of solvation using the data obtained above
        double deltaH0=deltaG0_octanol+(T*deltaS0);
        deltaS0=deltaS0/4.18;   //unit conversion from J/mol to cal/mol
        deltaH0=deltaH0/4180;   //unit conversion from J/mol to kcal/mol

       // Generation of Gas Phase data to add to the solution phase quantities
       ThermoData solvationCorrection = new ThermoData(deltaH0, deltaS0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,"Solvation correction");
    

       // Now, solvationCorrection contains solution phase estimates of CORRECTION TO H298, S298 and all the gas phase heat capacities.
       // Assuming the solution phase heat capcities to be the same as that in the gas phase we wouls now want to pass on this
       // modified version of result to the kinetics codes. This might require reading in a keyword from the condition.txt file.
       // Exactly how this will be done is yet to be figured out.

        return solvationCorrection;
        //#]
    }


            //protected static GATP_Solvation getINSTANCE() {



        public AbrahamGAValue getABGroup(ChemGraph p_chemGraph) {
        //#[ operation getGAGroup(ChemGraph)
        
        AbramData result_abram = new AbramData();
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
          			AbrahamGAValue thisAbrahamValue = thermoLibrary.findAbrahamGroup(p_chemGraph);


                                if (thisAbrahamValue == null) {
           				System.err.println("Abraham group not found: " + node.getID());
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
//         // recover the chem graph structure
//         // recover the radical
//         for (int i=0; i<radicalSite; i++) {
//           	int id = idArray[i];
//           	Node node = g.getNodeAt(id);
//           	node.setElement(atomArray[i]);
//           	node.updateFeElement();
//           	int HNum = atomArray[i].getRadicalNumber();
//         	//get rid of extra H
//           	for (int j=0;j<HNum;j++) {
//          	g.removeNode(newnode[i][j]);
//           	}
//           	node.updateFgElement();
//         }
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


        protected static GATP_Solvation getINSTANCE() {

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
