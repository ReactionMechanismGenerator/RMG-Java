//!*******************************************************************
//!
//!     Copyright: Jing Song, 2002, all rights reserved
//!
//!*******************************************************************

/*********************************************************************
	Rhapsody	: 4.0
	Login		: Administrator
	Component	: RMG
	Configuration 	: RMG
	Model Element	: RMG
//!	Generated Date	: Sun, 10, Nov 2002
	File Path	: RMG\RMG\MainRMG.java
*********************************************************************/



import java.util.*;
import java.io.*;
import jing.chem.*;
import jing.chemParser.*;
import jing.param.*;
import jing.chemUtil.*;
import jing.rxn.*;
import jing.rxnSys.*;
import jing.mathTool.*;

//----------------------------------------------------------------------------
// MainRMG.java
//----------------------------------------------------------------------------

public class RMG {
	
    //## configuration RMG::RMG
    public static void main(String[] args) {
	
    long begin = System.currentTimeMillis();
	//begin = getCpuTime();
	//System.out.println(getCpuTime()/1e9 + " " + (System.currentTimeMillis()-begin)/1e3);
	//System.out.println(getCpuTime()/1e9 + " " + (System.currentTimeMillis()-begin)/1e3);
	//Pressure pres = new Pressure(30,"bar");
	//System.out.println("The size of the object pressure is "+getObjectSize(pres));
	initializeSystemProperties(args[0]);
    ReactionModelGenerator rmg = new ReactionModelGenerator();
    rmg.modelGeneration();
    ReactionSystem rs = rmg.getReactionSystem();
    CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)rs.getReactionModel();

    System.out.println("\n\\\\\\\\\\\\\\\\\\\\\\\\\\    Final Reaction Model Output    \\\\\\\\\\\\\\\\\\\\\\\\\\");
    cerm.printPDepModel(rs.getPresentTemperature());
    System.out.println("Model Edge:\n");
    System.out.print(cerm.getEdge().getSpeciesNumber());
    System.out.print(" Species;\t");
    System.out.print(cerm.getEdge().getReactionNumber());
    System.out.print(" Reactions.\n");


    LinkedList speList = new LinkedList(rs.getReactionModel().getSpeciesSet());

    System.out.println("\\\\\\\\\\\\\\\\\\\\\\\\\\\\    Concentration Profile Output    \\\\\\\\\\\\\\\\\\\\\\\\\\");
    System.out.println(rs.printConcentrationProfile(speList));
    System.out.println("\\\\\\\\\\\\\\\\\\\\\\\\\\\\    Mole Fraction Profile Output    \\\\\\\\\\\\\\\\\\\\\\\\\\");
    System.out.println(rs.printMoleFractionProfile(speList));

    long end = System.currentTimeMillis();
    double min = (end-begin)/1E3/60;
    System.out.println("Running Time is: " + String.valueOf(min) + " minutes.");

    };

	public static void initializeSystemProperties(String inputfile) {
	    File f = new File(".");
	    String dir = f.getAbsolutePath();
		File chemdis = new File("chemdis");
		chemdis.mkdir();
		File therfit = new File("therfit");
		therfit.mkdir();
		File fit3p = new File("fit3p");
		fit3p.mkdir();
		File chemkin = new File("chemkin");
		chemkin.mkdir();
		File Restart = new File("Restart");
		Restart.mkdir();
		
		
		 String workingDir = System.getenv("RMG");
	     System.setProperty("RMG.workingDirectory", workingDir);
	     System.setProperty("jing.rxnSys.ReactionModelGenerator.conditionFile",inputfile);
		 try {//svp
             String initialConditionFile = System.getProperty("jing.rxnSys.ReactionModelGenerator.conditionFile");
             if (initialConditionFile == null) {
                     System.out.println("undefined system property: jing.rxnSys.ReactionModelGenerator.conditionFile");
                     System.exit(0);
             }

             FileReader in = new FileReader(initialConditionFile);
             BufferedReader reader = new BufferedReader(in);
             String line = ChemParser.readMeaningfulLine(reader);
			 	line = ChemParser.readMeaningfulLine(reader);
			 
             if (line.startsWith("Database")){
               StringTokenizer st = new StringTokenizer(line);
               String next = st.nextToken();
               String name = st.nextToken().trim();
               System.setProperty("jing.chem.ChemGraph.forbiddenStructureFile", workingDir + "/databases/"+name+"/forbiddenStructure/ForbiddenStructure.txt");
               System.setProperty("jing.chem.ThermoGAGroupLibrary.pathName", workingDir + "/databases/" + name+"/thermo");
               System.setProperty("jing.rxn.ReactionTemplateLibrary.pathName", workingDir + "/databases/" + name+"/kinetics/kinetics");
             }
  }
  catch (IOException e) {
     System.err.println("Error in read in reaction system initialization file!");
}
  };
}
/*********************************************************************
	File Path	: RMG\RMG\MainRMG.java
*********************************************************************/

