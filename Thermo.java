import java.util.*;
import java.io.*;
import jing.chem.*;
import jing.chemParser.*;
import jing.param.*;
import jing.chemUtil.*;
//import bondGroups.*;
import jing.rxn.*;
import jing.rxnSys.*;
import jing.mathTool.*;


/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class Thermo {

  //## configuration RMG::RMG
public static void main(String[] args) {
  initializeSystemProperties();

 try {
          FileReader in = new FileReader("sandeep.txt");
          BufferedReader data = new BufferedReader(in);

          Graph g = ChemParser.readChemGraph(data);

          in.close();

         System.out.println(g);
         ChemGraph chemgraph = ChemGraph.make(g);
		 Species spe = Species.make("molecule",chemgraph);
		 System.out.println("The number of resonance isomers is " + spe.getResonanceIsomersHashSet().size());
		 System.out.println("The NASA data is \n"+ spe.getNasaThermoData());
        //int K = chemgraph.getKekule();
        int symm = chemgraph.getSymmetryNumber();
        //System.out.println("number of Kekule structures = "+K);
        System.out.println(" symmetry number = "+symm);

         Temperature T = new Temperature(298.0,"K");

         String chemicalFormula = chemgraph.getChemicalFormula();

         System.out.println(chemicalFormula + "  H=" + chemgraph.calculateH(T));


    //      Species species = Species.make(chemicalFormula, chemgraph);
          // this is equal to  System.out.println(species.toString());
    //      System.out.println(species.toStringWithoutH());
  //        species.generateResonanceIsomers();
  //        Iterator rs = species.getResonanceIsomers();
  //        while (rs.hasNext()){
  //          ChemGraph cg = (ChemGraph)rs.next();
  //          Species s = cg.getSpecies();
  //          String string = s.toStringWithoutH();
  //          System.out.print(string);
  //        }

   //       Species species = Species.make(chemicalFormula, chemgraph);
   //       Iterator iterator = species.getResonanceIsomers();
   //       System.out.println(iterator);

 }
 catch (FileNotFoundException e) {
   System.err.println("File was not found!\n");
 }
 catch(IOException e){
   System.err.println("Something wrong with ChemParser.readChemGraph");
 }
 catch(ForbiddenStructureException e){
   System.err.println("Something wrong with ChemGraph.make");
 }


System.out.println("Done!\n");

};

 public static void initializeSystemProperties() {
	 String name= "database";
	 String workingDir = System.getenv("RMG");
     System.setProperty("RMG.workingDirectory", workingDir);
 //  System.setProperty("jing.chem.ChemGraph.forbiddenStructureFile",
 //                     workingDir +
 //                     "database/forbiddenStructure/forbiddenStructure.txt");
     System.setProperty("jing.chem.ChemGraph.forbiddenStructureFile", workingDir + "/databases/"+name+"/forbiddenStructure/ForbiddenStructure.txt");
     System.setProperty("jing.chem.ThermoGAGroupLibrary.pathName", workingDir + "/databases/" + name+"/thermo");
     System.setProperty("jing.rxn.ReactionTemplateLibrary.pathName", workingDir + "/databases/" + name+"/kinetics/kinetics");

   // System.setProperty("jing.rxn.ReactionTemplateLibrary.pathName",
 //                     workingDir + "database/kinetics/kinetics");
 //  System.setProperty("jing.rxnSys.ReactionModelGenerator.conditionFile",
 //                     workingDir + "database/condition/condition.txt");

 };







}
