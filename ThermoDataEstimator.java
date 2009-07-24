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


public class ThermoDataEstimator {//gmagoon 7/24/09: based off of Thermo.java rev. 1.6; this performs original functionality described in the manual

  //## configuration RMG::RMG
  //first argument will be file to read
public static void main(String[] args) {
  initializeSystemProperties();

 try {
          FileReader in = new FileReader(args[0]);
          BufferedReader data = new BufferedReader(in);

          Graph g = ChemParser.readChemGraph(data);

          in.close();

         System.out.println(g);
         ChemGraph chemgraph = ChemGraph.make(g);
		 Species spe = Species.make("molecule",chemgraph);
		 System.out.println("The number of resonance isomers is " + spe.getResonanceIsomersHashSet().size());
		 System.out.println("The NASA data is \n"+ spe.getNasaThermoData());
		 System.out.println("ThermoData is \n" +  spe.getChemGraph().getThermoData().toString());
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
	 String name= "RMG_database";
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

