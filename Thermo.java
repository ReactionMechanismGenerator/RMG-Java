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
	LinkedHashSet speciesSet = new LinkedHashSet();
    String thermo_output = "";

 try {
          FileReader in = new FileReader("thermo_input.txt");
          BufferedReader data = new BufferedReader(in);
          
          // Read the first line of thermo_input.txt
          String line = ChemParser.readMeaningfulLine(data);
          StringTokenizer st = new StringTokenizer(line);
          // The first line should start with "Solvation", otherwise do nothing and display a message to the user
          if (st.nextToken().startsWith("Solvation")) {
        	  line = st.nextToken().toLowerCase();
        	  // The options for the "Solvation" field are "on" or "off" (as of 18May2009), otherwise do nothing and display a message to the user
        	  // Note: I use "Species.useInChI" because the "Species.useSolvation" updates were not yet committed.
        	  if (line.equals("on")) {
        		  Species.useInChI = true;
        		  thermo_output += "Solution-phase chemistry!\n\n";
        	  } else if (line.equals("off")) {
        		  Species.useInChI = false;
        		  thermo_output += "Gas-phase chemistry.\n\n";
        	  } else {
        		  System.out.println("Error in reading thermo_input.txt file:\nThe field 'Solvation' has the options 'on' or 'off'.");
        		  return;
        	  }
        	  // Read in the ChemGraphs and compute their thermo, while there are ChemGraphs to read in
        	  line = ChemParser.readMeaningfulLine(data);
        	  while (line != null) {
        		  String speciesName = line;
        		  Graph g = ChemParser.readChemGraph(data);
        		  ChemGraph cg = null;
        		  try {
        			  cg = ChemGraph.make(g);
        		  } catch (ForbiddenStructureException e) {
        			  System.out.println("Error in reading graph: Graph contains a forbidden structure.\n" + g.toString());
        			  System.exit(0);
        		  }
        		  Species species = Species.make(speciesName,cg);
        		  speciesSet.add(species);
        		  line = ChemParser.readMeaningfulLine(data);
        	  }
          } else
        	  System.out.println("Error in reading thermo_input.txt file:\nThe first line must read 'Solvation: on/off'.");

          in.close();
          
          thermo_output += "Order of entries: Name (read from thermo_input.txt) H298 S298 Cp300 Cp400 Cp500 Cp600 Cp800 Cp1000 C1500\n" +
          	"Units of H: kcal/mol\nUnits of S and all Cp: cal/mol/K\n\n";
          
          Iterator iter = speciesSet.iterator();       
          while (iter.hasNext()){
        	  Species spe = (Species)iter.next();
        	  thermo_output += spe.getName() + "\t" + spe.getChemGraph().getThermoData().toString() + "\n";
          }
          
          try {
        	  File thermoOutput = new File("thermo_output.txt");
        	  FileWriter fw = new FileWriter(thermoOutput);
        	  fw.write(thermo_output);
        	  fw.close();
          } catch (IOException e) {
        	  System.out.println("Error in writing thermo_output.txt file.");
        	  System.exit(0);
          }

//		 System.out.println("The number of resonance isomers is " + spe.getResonanceIsomersHashSet().size());
//		 System.out.println("The NASA data is \n"+ spe.getNasaThermoData());
//		 System.out.println("ThermoData is \n" +  spe.getChemGraph().getThermoData().toString());
//         System.out.println("AbramData is \n" +  spe.getChemGraph().getAbramData().toString());
//         System.out.println("UnifacData is \n" +  spe.getChemGraph().getUnifacData().toString());
//        //int K = chemgraph.getKekule();
//        int symm = chemgraph.getSymmetryNumber();
//        //System.out.println("number of Kekule structures = "+K);
//        System.out.println(" symmetry number = "+symm);

//         Temperature T = new Temperature(298.0,"K");
//
//         String chemicalFormula = chemgraph.getChemicalFormula();
//
//         System.out.println(chemicalFormula + "  H=" + chemgraph.calculateH(T));


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
// catch(ForbiddenStructureException e){
//   System.err.println("Something wrong with ChemGraph.make");
// }


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
