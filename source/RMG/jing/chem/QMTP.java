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
import jing.param.*;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;

//quantum mechanics thermo property estimator; analog of GATP
public class QMTP implements GeneralGAPP {

    final protected static double ENTHALPY_HYDROGEN = 52.1; //needed for HBI
    private static QMTP INSTANCE = new QMTP();		//## attribute INSTANCE
    protected static PrimaryThermoLibrary primaryLibrary;//Note: may be able to separate this out into GeneralGAPP, as this is common to both GATP and QMTP
    public static String qmfolder= "QMfiles/";
    //   protected static HashMap library;		//as above, may be able to move this and associated functions to GeneralGAPP (and possibly change from "x implements y" to "x extends y"), as it is common to both GATP and QMTP
    protected ThermoGAGroupLibrary thermoLibrary; //needed for HBI
    public static String qmprogram= "both";//the qmprogram can be "mopac", "gaussian03", "both" (MOPAC and Gaussian), or "mm4"/"mm4hr"
    public static boolean usePolar = false; //use polar keyword in MOPAC
    public static boolean useCanTherm = true; //whether to use CanTherm in MM4 cases for interpreting output via force-constant matrix; this will hopefully avoid zero frequency issues
    public static boolean useHindRot = false;//whether to use HinderedRotor scans with MM4 (requires useCanTherm=true)
    public static double R = 1.9872; //ideal gas constant in cal/mol-K (does this appear elsewhere in RMG, so I don't need to reuse it?)
    public static double Hartree_to_kcal = 627.5095; //conversion from Hartree to kcal/mol taken from Gaussian thermo white paper
    public static double Na = 6.02214179E23;//Avagadro's number; cf. http://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=physchem_in!
    public static double k = 1.3806504E-23;//Boltzmann's constant in J/K; cf. http://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=physchem_in!
    public static double h = 6.62606896E-34;//Planck's constant in J-s; cf. http://physics.nist.gov/cgi-bin/cuu/Value?h|search_for=universal_in!
    public static double c = 299792458. *100;//speed of light in vacuum in cm/s, cf. http://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=universal_in!
    public static double deltaTheta=5.0;//degree increment for rotor scans when using useHindRot
    // Constructors

    //## operation QMTP()
    private QMTP() {
       // initializeLibrary(); //gmagoon 72509: commented out in GATP, so I am mirroring the change here; other library functions below also commented out
        initializePrimaryThermoLibrary();
    }
       //## operation generateThermoData(ChemGraph)
    public ThermoData generateThermoData(ChemGraph p_chemGraph) {
        //#[ operation generateThermoData(ChemGraph)
        //first, check for thermo data in the primary thermo library and library (?); if it is there, use it
        ThermoData result = primaryLibrary.getThermoData(p_chemGraph.getGraph());
        //System.out.println(result);
        if (result != null) {
        	p_chemGraph.fromprimarythermolibrary = true;
        	return result;
        }
       
        
        //        result = getFromLibrary(p_chemGraph.getChemicalFormula());//gmagoon 72509: commented out in GATP, so I am mirroring the change here
//        if (result != null) return result;
        
        result=new ThermoData();
         
        int maxRadNumForQM = Global.maxRadNumForQM;
        if (p_chemGraph.getRadicalNumber() > maxRadNumForQM)//use HBI if the molecule has more radicals than maxRadNumForQM; this is helpful because ; also MM4 (and MM3) look like they may have issues with radicals
        {//this code is based closely off of GATP saturation (in getGAGroup()), but there are some modifications, particularly for symmetry correction
            //find the initial symmetry number
            int sigmaRadical = p_chemGraph.getSymmetryNumber();
            
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

            // add H to saturate chem graph
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

            //find the saturated symmetry number
            int sigmaSaturated = p_chemGraph.getSymmetryNumber();
            
   //         result = generateThermoData(g);//I'm not sure what GATP does, but this recursive calling will use HBIs on saturated species if it exists in PrimaryThermoLibrary
            //check the primary thermo library for the saturated graph
            result = primaryLibrary.getThermoData(p_chemGraph.getGraph());
            //System.out.println(result);
            if (result != null) {
        	p_chemGraph.fromprimarythermolibrary = true;
            }
            else{
                result=generateQMThermoData(p_chemGraph);
            }
            
            // find the BDE for all radical groups
            if(thermoLibrary == null) initGAGroupLibrary();
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

                //recover the saturated site for next radical site calculation
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
           
            
            //correct the symmetry number based on the relative radical and saturated symmetry number; this should hopefully sidestep potential complications based on the fact that certain symmetry effects could be included in HBI value itself, and the fact that the symmetry number correction for saturated molecule has already been implemented, and it is likely to be different than symmetry number considered here, since the correction for the saturated molecule will have been external symmetry number, whereas RMG's ChemGraph symmetry number estimator includes both internal and external symmetry contributions; even so, I don't know if this will handle a change from chiral to achiral (or vice versa) properly            
	    ThermoGAValue symmetryNumberCorrection = new ThermoGAValue(0,-1*GasConstant.getCalMolK()*Math.log((double)(sigmaRadical)/(double)(sigmaSaturated)),0,0,0,0,0,0,0,0,0,0,null);
            result.plus(symmetryNumberCorrection);
            
            p_chemGraph.setCentralNode(oldCentralNode);  
            
            //display corrected thermo to user
            String [] InChInames = getQMFileName(p_chemGraph);//determine the filename (InChIKey) and InChI with appended info for triplets, etc.
            String name = InChInames[0];
            String InChIaug = InChInames[1];
            System.out.println("HBI-based thermo for " + name + "("+InChIaug+"): "+ result.toString());//print result, at least for debugging purposes
        }
        else{
            result = generateQMThermoData(p_chemGraph);
        }
        
        return result;
        //#]
    }
   
    public ThermoData generateQMThermoData(ChemGraph p_chemGraph){
        //if there is no data in the libraries, calculate the result based on QM or MM calculations; the below steps will be generalized later to allow for other quantum mechanics packages, etc.
        String qmProgram = qmprogram;
	String qmMethod = "";
	if(qmProgram.equals("mm4")||qmProgram.equals("mm4hr")){
	    qmMethod = "mm4";
	    if(qmProgram.equals("mm4hr")) useHindRot=true;
	}
	else{
	    qmMethod="pm3"; //may eventually want to pass this to various functions to choose which "sub-function" to call
	}
        
        ThermoData result = new ThermoData();
	double[] dihedralMinima = null;
        
        String [] InChInames = getQMFileName(p_chemGraph);//determine the filename (InChIKey) and InChI with appended info for triplets, etc.
        String name = InChInames[0];
        String InChIaug = InChInames[1];
        String directory = qmfolder;
        File dir=new File(directory);
        directory = dir.getAbsolutePath();//this and previous three lines get the absolute path for the directory
	if(qmMethod.equals("pm3")){
	    //first, check to see if the result already exists and the job terminated successfully
	    boolean gaussianResultExists = successfulGaussianResultExistsQ(name,directory,InChIaug);
	    boolean mopacResultExists = successfulMopacResultExistsQ(name,directory,InChIaug);
	    if(!gaussianResultExists && !mopacResultExists){//if a successful result doesn't exist from previous run (or from this run), run the calculation; if a successful result exists, we will skip directly to parsing the file
		 //steps 1 and 2: create 2D and 3D mole files
		molFile p_3dfile = create3Dmolfile(name, p_chemGraph);
		 //3. create the Gaussian or MOPAC input file
		directory = qmfolder;
		dir=new File(directory);
		directory = dir.getAbsolutePath();//this and previous three lines get the absolute path for the directory
		int attemptNumber=1;//counter for attempts using different keywords
		int successFlag=0;//flag for success of Gaussian run; 0 means it failed, 1 means it succeeded
		int maxAttemptNumber=1;
		int multiplicity = p_chemGraph.getRadicalNumber()+1; //multiplicity = radical number + 1
		while(successFlag==0 && attemptNumber <= maxAttemptNumber){
		    //IF block to check which program to use
		    if (qmProgram.equals("gaussian03")){
			if(p_chemGraph.getAtomNumber() > 1){
			    maxAttemptNumber = createGaussianPM3Input(name, directory, p_3dfile, attemptNumber, InChIaug, multiplicity);
			}
			else{
			    maxAttemptNumber = createGaussianPM3Input(name, directory, p_3dfile, -1, InChIaug, multiplicity);//use -1 for attemptNumber for monoatomic case
			}
			//4. run Gaussian
			successFlag = runGaussian(name, directory);
		    }
		    else if (qmProgram.equals("mopac") || qmProgram.equals("both")){
			maxAttemptNumber = createMopacPM3Input(name, directory, p_3dfile, attemptNumber, InChIaug, multiplicity);
			successFlag = runMOPAC(name, directory);
		    }
		    else{
			System.out.println("Unsupported quantum chemistry program");
			System.exit(0);
		    }
		    //new IF block to check success
		    if(successFlag==1){
			System.out.println("Attempt #"+attemptNumber + " on species " + name + " ("+InChIaug+") succeeded.");
		    }
		    else if(successFlag==0){
			if(attemptNumber==maxAttemptNumber){//if this is the last possible attempt, and the calculation fails, exit with an error message
			    if(qmProgram.equals("both")){ //if we are running with "both" option and all keywords fail, try with Gaussian
				qmProgram = "gaussian03";
				System.out.println("*****Final MOPAC attempt (#" + maxAttemptNumber + ") on species " + name + " ("+InChIaug+") failed. Trying to use Gaussian.");
				attemptNumber=0;//this needs to be 0 so that when we increment attemptNumber below, it becomes 1 when returning to the beginning of the for loop
				maxAttemptNumber=1;
			    }
			    else{
				System.out.println("*****Final attempt (#" + maxAttemptNumber + ") on species " + name + " ("+InChIaug+") failed.");
				System.out.print(p_chemGraph.toString());
			        System.exit(0);
			//	return new ThermoData(1000,0,0,0,0,0,0,0,0,0,0,0,"failed calculation");
			    }
			}
			System.out.println("*****Attempt #"+attemptNumber + " on species " + name + " ("+InChIaug+") failed. Will attempt a new keyword.");
			attemptNumber++;//try again with new keyword
		    }
		}


	    }
	    //5. parse QM output and record as thermo data (function includes symmetry/point group calcs, etc.); if both Gaussian and MOPAC results exist, Gaussian result is used
	    if (gaussianResultExists || (qmProgram.equals("gaussian03") && !mopacResultExists)){
		result = parseGaussianPM3(name, directory, p_chemGraph);
	    }
	    else if (mopacResultExists || qmProgram.equals("mopac") || qmProgram.equals("both")){
		result = parseMopacPM3(name, directory, p_chemGraph);
	    }
	    else{
		System.out.println("Unexpected situation in QMTP thermo estimation");
		System.exit(0);
	    }
	}
	else{//mm4 case
	    //first, check to see if the result already exists and the job terminated successfully
	    boolean mm4ResultExists = successfulMM4ResultExistsQ(name,directory,InChIaug);
	    if(!mm4ResultExists){//if a successful result doesn't exist from previous run (or from this run), run the calculation; if a successful result exists, we will skip directly to parsing the file
		 //steps 1 and 2: create 2D and 3D mole files
		molFile p_3dfile = create3Dmolfile(name, p_chemGraph);
		 //3. create the MM4 input file
		directory = qmfolder;
		dir=new File(directory);
		directory = dir.getAbsolutePath();//this and previous three lines get the absolute path for the directory
		int attemptNumber=1;//counter for attempts using different keywords
		int successFlag=0;//flag for success of MM4 run; 0 means it failed, 1 means it succeeded
		int maxAttemptNumber=1;
		int multiplicity = p_chemGraph.getRadicalNumber()+1; //multiplicity = radical number + 1
		while(successFlag==0 && attemptNumber <= maxAttemptNumber){
		    maxAttemptNumber = createMM4Input(name, directory, p_3dfile, attemptNumber, InChIaug, multiplicity);
		    //4. run MM4
		    successFlag = runMM4(name, directory);
		    //new IF block to check success
		    if(successFlag==1){
			System.out.println("Attempt #"+attemptNumber + " on species " + name + " ("+InChIaug+") succeeded.");
			//run rotor calculations if necessary
			int rotors = p_chemGraph.getInternalRotor();
			if(useHindRot && rotors > 0){
			    //we should re-run scans even if pre-existing scans exist because atom numbering may change from case to case; a better solution would be to check for stored CanTherm output and use that if available
			    System.out.println("Running rotor scans on "+name+"...");
			    dihedralMinima = createMM4RotorInput(name, directory, p_chemGraph, rotors);//we don't worry about checking InChI here; if there is an InChI mismatch it should be caught
			    runMM4Rotor(name, directory, rotors);
			}
		    }
		    else if(successFlag==0){
			if(attemptNumber==maxAttemptNumber){//if this is the last possible attempt, and the calculation fails, exit with an error message
				System.out.println("*****Final attempt (#" + maxAttemptNumber + ") on species " + name + " ("+InChIaug+") failed.");
				System.out.print(p_chemGraph.toString());
			        System.exit(0);
				//return new ThermoData(1000,0,0,0,0,0,0,0,0,0,0,0,"failed calculation");
			}
			System.out.println("*****Attempt #"+attemptNumber + " on species " + name + " ("+InChIaug+") failed. Will attempt a new keyword.");
			attemptNumber++;//try again with new keyword
		    }
		}
		if(useCanTherm){
		    performCanThermCalcs(name, directory, p_chemGraph, dihedralMinima, false);
		    if (p_chemGraph.getInternalRotor()>0) performCanThermCalcs(name, directory, p_chemGraph, dihedralMinima, true);//calculate RRHO case for comparison
		}

	    }
	    //5. parse MM4 output and record as thermo data (function includes symmetry/point group calcs, etc.)
	    if(!useCanTherm) result = parseMM4(name, directory, p_chemGraph);
	    else{
		//if (qmdata==null) qmdata = getQMDataWithCClib(name, directory, p_chemGraph, true);//get qmdata if it is null (i.e. if a pre-existing successful result exists and it wasn't read in above)
		result = parseCanThermFile(name, directory, p_chemGraph);
		if (p_chemGraph.getInternalRotor()>0) parseCanThermFile(name+"_RRHO", directory, p_chemGraph);//print the RRHO result for comparison
	    }
	}
        
        return result;
    }

    protected static QMTP getINSTANCE() {
        return INSTANCE;
    }
    

    public void initializePrimaryThermoLibrary(){//svp

        primaryLibrary = PrimaryThermoLibrary.getINSTANCE();

      }
      
    //creates a 3D molFile; for monoatomic species, it just returns the 2D molFile
    public molFile create3Dmolfile(String name, ChemGraph p_chemGraph){
	//1. create a 2D file
	//use the absolute path for directory, so we can easily reference from other directories in command-line paths
	//can't use RMG.workingDirectory, since this basically holds the RMG environment variable, not the workingDirectory
	String directory = "2Dmolfiles/";
	File dir=new File(directory);
	directory = dir.getAbsolutePath();
	molFile p_2dfile = new molFile(name, directory, p_chemGraph);
	molFile p_3dfile = new molFile();//it seems this must be initialized, so we initialize to empty object
	//2. convert from 2D to 3D using RDKit if the 2D molfile is for a molecule with 2 or more atoms
	int atoms = p_chemGraph.getAtomNumber();
	if(atoms > 1){
	    int distGeomAttempts=1;
	    if(atoms > 3){//this check prevents the number of attempts from being negative
		distGeomAttempts = 5*(p_chemGraph.getAtomNumber()-3); //number of conformer attempts is just a linear scaling with molecule size, due to time considerations; in practice, it is probably more like 3^(n-3) or something like that
	    }
	    p_3dfile = embed3D(p_2dfile, distGeomAttempts);
	    return p_3dfile;
	}
	else{
	    return p_2dfile;
	}


    }
    //embed a molecule in 3D, using RDKit
    public molFile embed3D(molFile twoDmolFile, int numConfAttempts){
    //convert to 3D MOL file using RDKit script
        int flag=0;
        String directory = "3Dmolfiles/";
        File dir=new File(directory);
        directory = dir.getAbsolutePath();//this uses the absolute path for the directory
        String name = twoDmolFile.getName();
        try{   
            File runningdir=new File(directory);
	    String command="";
	    if (System.getProperty("os.name").toLowerCase().contains("windows")){//special windows case where paths can have spaces and are allowed to be surrounded by quotes
		command = "python \""+System.getProperty("RMG.workingDirectory")+"/scripts/distGeomScriptMolLowestEnergyConf.py\" ";
		String twoDmolpath=twoDmolFile.getPath();
		command=command.concat("\""+twoDmolpath+"\" ");
		command=command.concat("\""+name+".mol\" ");//this is the target file name; use the same name as the twoDmolFile (but it will be in he 3Dmolfiles folder
		command=command.concat("\""+name+".cmol\" ");//this is the target file name for crude coordinates (corresponding to the minimum energy conformation based on UFF refinement); use the same name as the twoDmolFile (but it will be in he 3Dmolfiles folder) and have suffix .cmol
		command=command.concat(numConfAttempts + " ");
		command=command.concat("\"" + System.getenv("RDBASE")+"\"");//pass the $RDBASE environment variable to the script so it can use the approprate directory when importing rdkit
	    }    
	    else{//non-Windows case
		command = "python "+System.getProperty("RMG.workingDirectory")+"/scripts/distGeomScriptMolLowestEnergyConf.py ";
		String twoDmolpath=twoDmolFile.getPath();
		command=command.concat(""+twoDmolpath+" ");
		command=command.concat(name+".mol ");//this is the target file name; use the same name as the twoDmolFile (but it will be in he 3Dmolfiles folder
		command=command.concat(name+".cmol ");//this is the target file name for crude coordinates (corresponding to the minimum energy conformation based on UFF refinement); use the same name as the twoDmolFile (but it will be in he 3Dmolfiles folder) and have suffix .cmol
		command=command.concat(numConfAttempts + " ");
		command=command.concat(System.getenv("RDBASE"));//pass the $RDBASE environment variable to the script so it can use the approprate directory when importing rdkit
	    }
	    Process pythonProc = Runtime.getRuntime().exec(command, null, runningdir);
            String killmsg= "Python process for "+twoDmolFile.getName()+" did not complete within 120 seconds, and the process was killed. File was probably not written.";//message to print if the process times out
            Thread timeoutThread = new TimeoutKill(pythonProc, killmsg, 120000L); //create a timeout thread to handle cases where the UFF optimization get's locked up (cf. Ch. 16 of "Ivor Horton's Beginning Java 2: JDK 5 Edition"); once we use the updated version of RDKit, we should be able to get rid of this
            timeoutThread.start();//start the thread
            //check for errors and display the error if there is one
            InputStream is = pythonProc.getErrorStream();
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line=null;
            while ( (line = br.readLine()) != null) {
                    line = line.trim();
                    System.err.println(line);
                    flag=1;
            }
            //if there was an error, indicate the file and InChI
            if(flag==1){
                System.out.println("RDKit received error (see above) on " + twoDmolFile.getName()+". File was probably not written.");
            }
            int exitValue = pythonProc.waitFor();
            if(timeoutThread.isAlive())//if the timeout thread is still alive (indicating that the process has completed in a timely manner), stop the timeout thread
                timeoutThread.interrupt();
        }
        catch (Exception e) {
            String err = "Error in running RDKit Python process \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }
        

        
// gmagoon 6/3/09 comment out InChI checking for now; in any case, the code will need to be updated, as it is copied from my testing code
//        //check whether the original InChI is reproduced
//        if(flag==0){
//            try{
//                File f=new File("c:/Python25/"+molfilename);
//                File newFile= new File("c:/Python25/mol3d.mol");
//                if(newFile.exists()){
//                    newFile.delete();//apparently renaming will not work unless target file does not exist (at least on Vista)
//                }
//                f.renameTo(newFile);
//                String command = "c:/Users/User1/Documents/InChI-1/cInChI-1.exe c:/Python25/mol3d.mol inchi3d.inchi /AuxNone /DoNotAddH";//DoNotAddH used to prevent adding Hs to radicals (this would be a problem for current RDKit output which doesn't use M RAD notation)
//                Process inchiProc = Runtime.getRuntime().exec(command);	
//               // int exitValue = inchiProc.waitFor();
//                Thread.sleep(200);//****update: can probably eliminate this using buffered reader
//                inchiProc.destroy();
//                
//                //read output file
//                File outputFile = new File("inchi3d.inchi");
//                FileReader fr = new FileReader(outputFile);
//                BufferedReader br = new BufferedReader(fr);
//        	String line=null;
//                String inchi3d=null;
//                while ( (line = br.readLine()) != null) {
//                        line = line.trim();
//                        if(line.startsWith("InChI="))
//                        {
//                            inchi3d=line;
//                        }
//                }
//                fr.close();
//                
//                //return file to original name:
//                File f2=new File("c:/Python25/mol3d.mol");
//                File newFile2= new File("c:/Python25/"+molfilename);
//                if(newFile2.exists()){
//                    newFile2.delete();
//                }
//                f2.renameTo(newFile2);
//                
//                //compare inchi3d with input inchi and print a message if they don't match
//                if(!inchi3d.equals(inchiString)){
//                    if(inchi3d.startsWith(inchiString)&&inchiString.length()>10){//second condition ensures 1/C does not match 1/CH4; 6 characters for InChI=, 2 characters for 1/, 2 characters for atom layer
//                        System.out.println("(probably minor) For File: "+ molfilename+" , 3D InChI (" + inchi3d+") begins with, but does not match original InChI ("+inchiString+"). SMILES string: "+ smilesString); 
//                        
//                    }
//                    else{
//                        System.out.println("For File: "+ molfilename+" , 3D InChI (" + inchi3d+") does not match original InChI ("+inchiString+"). SMILES string: "+ smilesString);
//                    }
//                }
//            }
//            catch (Exception e) {
//                String err = "Error in running InChI process \n";
//                err += e.toString();
//                e.printStackTrace();
//                System.exit(0);
//            }
//        }
        
        
        //construct molFile pointer to new file (name will be same as 2D mol file
        return new molFile(name, directory);
    }
    
    //creates Gaussian PM3 input file in directory with filename name.gjf by using OpenBabel to convert p_molfile
    //attemptNumber determines which keywords to try
    //the function returns the maximum number of keywords that can be attempted; this will be the same throughout the evaluation of the code, so it may be more appropriate to have this as a "constant" attribute of some sort
    //attemptNumber=-1 will call a special set of keywords for the monoatomic case
    public int createGaussianPM3Input(String name, String directory, molFile p_molfile, int attemptNumber, String InChIaug, int multiplicity){
        //write a file with the input keywords
        int scriptAttempts = 18;//the number of keyword permutations available; update as additional options are added
        int maxAttemptNumber=2*scriptAttempts;//we will try a second time with crude coordinates if the UFF refined coordinates do not work
        try{
            File inpKey=new File(directory+"/inputkeywords.txt");
            String inpKeyStr="%chk="+directory+"/RMGrunCHKfile.chk\n";
            inpKeyStr+="%mem=6MW\n";
            inpKeyStr+="%nproc=1\n";
            if(attemptNumber==-1) inpKeyStr+="# pm3 freq";//keywords for monoatomic case (still need to use freq keyword to get molecular mass)
            else if(attemptNumber%scriptAttempts==1) inpKeyStr+="# pm3 opt=(verytight,gdiis) freq IOP(2/16=3)";//added IOP option to avoid aborting when symmetry changes; 3 is supposed to be default according to documentation, but it seems that 0 (the default) is the only option that doesn't work from 0-4; also, it is interesting to note that all 4 options seem to work for test case with z-matrix input rather than xyz coords; cf. http://www.ccl.net/cgi-bin/ccl/message-new?2006+10+17+005 for original idea for solution
            else if(attemptNumber%scriptAttempts==2) inpKeyStr+="# pm3 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)";//use different SCF method; this addresses at least one case of failure for a C4H7J species
            else if(attemptNumber%scriptAttempts==3) inpKeyStr+="# pm3 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm";//try multiple different options (no gdiis, use calcfc, nosymm); 7/21/09: added maxcyc option to fix case of MPTBUKVAJYJXDE-UHFFFAOYAPmult3 (InChI=1/C4H10O5Si/c1-3-7-9-10(5,6)8-4-2/h4-5H,3H2,1-2H3/mult3) (file manually copied to speed things along)
            else if(attemptNumber%scriptAttempts==4) inpKeyStr+="# pm3 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm";//7/8/09: numerical frequency keyword version of keyword #3; used to address GYFVJYRUZAKGFA-UHFFFAOYALmult3 (InChI=1/C6H14O6Si/c1-3-10-13(8,11-4-2)12-6-5-9-7/h6-7H,3-5H2,1-2H3/mult3) case; (none of the existing Gaussian or MOPAC combinations worked with it)
            else if(attemptNumber%scriptAttempts==5) inpKeyStr+="# pm3 opt=(verytight,gdiis,small) freq IOP(2/16=3)";//7/10/09: somehow, this worked for problematic case of ZGAWAHRALACNPM-UHFFFAOYAF (InChI=1/C8H17O5Si/c1-3-11-14(10,12-4-2)13-8-5-7(9)6-8/h7-9H,3-6H2,1-2H3); (was otherwise giving l402 errors); even though I had a keyword that worked for this case, I manually copied the fixed log file to QMfiles folder to speed things along; note that there are a couple of very low frequencies (~5-6 cm^-1 for this case)
            else if(attemptNumber%scriptAttempts==6) inpKeyStr+="# pm3 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)";//used for troublesome C5H7J2 case (similar error to C5H7J below); calcfc is not necessary for this particular species, but it speeds convergence and probably makes it more robust for other species
            else if(attemptNumber%scriptAttempts==7) inpKeyStr+="# pm3 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)"; //use numerical frequencies; this takes a relatively long time, so should only be used as one of the last resorts; this seemed to address at least one case of failure for a C6H10JJ species; 7/15/09: maxcyc=200 added to address GVCMURUDAUQXEY-UHFFFAOYAVmult3 (InChI=1/C3H4O7Si/c1-2(9-6)10-11(7,8)3(4)5/h6-7H,1H2/mult3)...however, result was manually pasted in QMfiles folder to speed things along
            else if(attemptNumber%scriptAttempts==8) inpKeyStr+="# pm3 opt=tight freq IOP(2/16=3)";//7/10/09: this worked for problematic case of SZSSHFMXPBKYPR-UHFFFAOYAF (InChI=1/C7H15O5Si/c1-3-10-13(8,11-4-2)12-7-5-6-9-7/h7H,3-6H2,1-2H3) (otherwise, it had l402.exe errors); corrected log file was manually copied to QMfiles to speed things along; we could also add a freq=numerical version of this keyword combination for added robustness; UPDATE: see below
            else if(attemptNumber%scriptAttempts==9) inpKeyStr+="# pm3 opt=tight freq=numerical IOP(2/16=3)";//7/10/09: used for problematic case of CIKDVMUGTARZCK-UHFFFAOYAImult4 (InChI=1/C8H15O6Si/c1-4-12-15(10,13-5-2)14-7-6-11-8(7,3)9/h7H,3-6H2,1-2H3/mult4 (most other cases had l402.exe errors);  corrected log file was manually copied to QMfiles to speed things along
            else if(attemptNumber%scriptAttempts==10) inpKeyStr+="# pm3 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)";//7/8/09: similar to existing #5, but uses tight rather than verytight; used for ADMPQLGIEMRGAT-UHFFFAOYAUmult3 (InChI=1/C6H14O5Si/c1-4-9-12(8,10-5-2)11-6(3)7/h6-7H,3-5H2,1-2H3/mult3)
            else if(attemptNumber%scriptAttempts==11) inpKeyStr+="# pm3 opt freq IOP(2/16=3)"; //use default (not verytight) convergence criteria; use this as last resort
            else if(attemptNumber%scriptAttempts==12) inpKeyStr+="# pm3 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)";//to address problematic C10H14JJ case
            else if(attemptNumber%scriptAttempts==13) inpKeyStr+="# pm3 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm";// added 6/10/09 for very troublesome RRMZRNPRCUANER-UHFFFAOYAQ (InChI=1/C5H7/c1-3-5-4-2/h3H,1-2H3) case...there were troubles with negative frequencies, where I don't think they should have been; step size of numerical frequency was adjusted to give positive result; accuracy of result is questionable; it is possible that not all of these keywords are needed; note that for this and other nearly free rotor cases, I think heat capacity will be overestimated by R/2 (R vs. R/2) (but this is a separate issue)
            else if(attemptNumber%scriptAttempts==14) inpKeyStr+="# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm";// added 6/22/09 for troublesome QDERTVAGQZYPHT-UHFFFAOYAHmult3(InChI=1/C6H14O4Si/c1-4-8-11(7,9-5-2)10-6-3/h4H,5-6H2,1-3H3/mult3); key aspects appear to be tight (rather than verytight) convergence criteria, no calculation of frequencies during optimization, use of numerical frequencies, and probably also the use of opt=small
      //gmagoon 7/9/09: commented out since although this produces a "reasonable" result for the problematic case, there is a large amount of spin contamination, apparently introducing 70+ kcal/mol of instability      else if(attemptNumber==12) inpKeyStr+="# pm3 opt=(verytight,gdiis,small) freq=numerical IOP(2/16=3) IOP(4/21=200)";//7/9/09: similar to current number 9 with keyword small; this addresses case of VCSJVABXVCFDRA-UHFFFAOYAI (InChI=1/C8H19O5Si/c1-5-10-8(4)13-14(9,11-6-2)12-7-3/h8H,5-7H2,1-4H3)
            else if(attemptNumber%scriptAttempts==15) inpKeyStr+="# pm3 opt=(verytight,gdiis,calcall) IOP(2/16=3)";//used for troublesome C5H7J case; note that before fixing, I got errors like the following: "Incomplete coordinate system.  Try restarting with Geom=Check Guess=Read Opt=(ReadFC,NewRedundant) Incomplete coordinate system. Error termination via Lnk1e in l103.exe"; we could try to restart, but it is probably preferrable to have each keyword combination standalone; another keyword that may be helpful if additional problematic cases are encountered is opt=small; 6/9/09 note: originally, this had # pm3 opt=(verytight,gdiis,calcall) freq IOP(2/16=3)" (with freq keyword), but I discovered that in this case, there are two thermochemistry sections and cclib parses frequencies twice, giving twice the number of desired frequencies and hence produces incorrect thermo; this turned up on C5H6JJ isomer
   //gmagoon 7/3/09: it is probably best to retire this keyword combination in light of the similar combination below         //else if(attemptNumber==6) inpKeyStr+="# pm3 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) IOP(4/21=2)";//6/10/09: worked for OJZYSFFHCAPVGA-UHFFFAOYAK (InChI=1/C5H7/c1-3-5-4-2/h1,4H2,2H3) case; IOP(4/21) keyword was key            
            else if(attemptNumber%scriptAttempts==16) inpKeyStr+="# pm3 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm";//6/29/09: worked for troublesome ketene case: CCGKOQOJPYTBIH-UHFFFAOYAO (InChI=1/C2H2O/c1-2-3/h1H2) (could just increase number of iterations for similar keyword combination above (#6 at the time of this writing), allowing symmetry, but nosymm seemed to reduce # of iterations; I think one of nosymm or higher number of iterations would allow the similar keyword combination to converge; both are included here for robustness)
            else if(attemptNumber%scriptAttempts==17) inpKeyStr+="# pm3 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm";//7/1/09: added for case of ZWMVZWMBTVHPBS-UHFFFAOYAEmult3 (InChI=1/C4H4O2/c1-3-5-6-4-2/h1-2H2/mult3)
            else if(attemptNumber%scriptAttempts==0) inpKeyStr+="# pm3 opt=(calcall,small,maxcyc=100) IOP(2/16=3)"; //6/10/09: used to address troublesome FILUFGAZMJGNEN-UHFFFAOYAImult3 case (InChI=1/C5H6/c1-3-5-4-2/h3H,1H2,2H3/mult3)
            else throw new Exception();//this point should not be reached
           // if(multiplicity == 3) inpKeyStr+= " guess=mix"; //assumed to be triplet biradical...use guess=mix to perform unrestricted ; nevermind...I think this would only be for singlet biradicals based on http://www.gaussian.com/g_tech/g_ur/k_guess.htm
            if (usePolar) inpKeyStr += " polar";
            FileWriter fw = new FileWriter(inpKey);
            fw.write(inpKeyStr);
            fw.close();
        }
        catch(Exception e){
            String err = "Error in writing inputkeywords.txt \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }
        
        //call the OpenBabel process (note that this requires OpenBabel environment variable)
        try{ 
            File runningdir=new File(directory);
	    String command=null;
	    String molPath=null;
	    if(attemptNumber<=scriptAttempts){//use UFF-refined coordinates
		molPath = p_molfile.getPath();
	    }
	    else{//use crude coordinates
		molPath = p_molfile.getCrudePath();
	    }
	    if (System.getProperty("os.name").toLowerCase().contains("windows")){//special windows case
		command = "babel -imol \""+ molPath+ "\" -ogjf \"" + name+".gjf\" -xf inputkeywords.txt --title \""+InChIaug+"\"";
	    }
	    else{
		command = "babel -imol "+ molPath+ " -ogjf " + name+".gjf -xf inputkeywords.txt --title "+InChIaug;
	    }
	    Process babelProc = Runtime.getRuntime().exec(command, null, runningdir);
            //read in output
            InputStream is = babelProc.getInputStream();
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line=null;
            while ( (line = br.readLine()) != null) {
		//do nothing
            }
            int exitValue = babelProc.waitFor();
        }
        catch(Exception e){
            String err = "Error in running OpenBabel MOL to GJF process \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }
        return maxAttemptNumber;
    }

    //creates MM4 input file and MM4 batch file in directory with filenames name.mm4 and name.com, respectively using MoleCoor
    //attemptNumber determines which keywords to try
    //the function returns the maximum number of keywords that can be attempted; this will be the same throughout the evaluation of the code, so it may be more appropriate to have this as a "constant" attribute of some sort
    public int createMM4Input(String name, String directory, molFile p_molfile, int attemptNumber, String InChIaug, int multiplicity){
        //Step 1: write the script for MM4 batch operation
//	Example script file:
//	#! /bin/csh
//	cp testEthylene.mm4 CPD.MM4
//	cp $MM4_DATDIR/BLANK.DAT PARA.MM4
//	cp $MM4_DATDIR/CONST.MM4 .
//	$MM4_EXEDIR/mm4 <<%
//	1
//	2
//	0
//	%
//	mv TAPE4.MM4 testEthyleneBatch.out
//	mv TAPE9.MM4 testEthyleneBatch.opt
//	exit

	int scriptAttempts = 2;//the number of script permutations available; update as additional options are added
        int maxAttemptNumber=2*scriptAttempts;//we will try a second time with crude coordinates if the UFF refined coordinates do not work
        try{
	    //create batch file with executable permissions: cf. http://java.sun.com/docs/books/tutorial/essential/io/fileAttr.html#posix
            File inpKey = new File(directory+"/"+name+".com");
            String inpKeyStr="#! /bin/csh\n";
            inpKeyStr+="cp "+name+".mm4 CPD.MM4\n";
            inpKeyStr+="cp $MM4_DATDIR/BLANK.DAT PARA.MM4\n";
	    inpKeyStr+="cp $MM4_DATDIR/CONST.MM4 .\n";
	    inpKeyStr+="$MM4_EXEDIR/mm4 <<%\n";
	    inpKeyStr+="1\n";//read from first line of .mm4 file
	    if (!useCanTherm){
		if(attemptNumber%scriptAttempts==1) inpKeyStr+="2\n"; //Block-Diagonal Method then Full-Matrix Method
		else if(attemptNumber%scriptAttempts==0) inpKeyStr+="3\n"; //Full-Matrix Method only
		else throw new Exception();//this point should not be reached
	    }
	    else{//CanTherm case: write the FORCE.MAT file
		if(attemptNumber%scriptAttempts==1) inpKeyStr+="4\n"; //Block-Diagonal Method then Full-Matrix Method
		else if(attemptNumber%scriptAttempts==0) inpKeyStr+="5\n"; //Full-Matrix Method only
		else throw new Exception();//this point should not be reached
		inpKeyStr+="\n";//<RETURN> for temperature
		inpKeyStr+="4\n";//unofficial option 4 for vibrational eigenvector printout to generate Cartesian force constant matrix in FORCE.MAT file
		inpKeyStr+="0\n";//no vibrational amplitude printout
	    }
            inpKeyStr+="0\n";//terminate the job
	    inpKeyStr+="%\n";
	    inpKeyStr+="mv TAPE4.MM4 "+name+".mm4out\n";
	    inpKeyStr+="mv TAPE9.MM4 "+name+".mm4opt\n";
	    if(useCanTherm){
		inpKeyStr+="mv FORCE.MAT "+name+".fmat\n";
	    }
	    inpKeyStr+="exit\n";
	    FileWriter fw = new FileWriter(inpKey);
            fw.write(inpKeyStr);
            fw.close();
        }
        catch(Exception e){
            String err = "Error in writing MM4 script file\n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }

        //Step 2: call the MoleCoor process to create the MM4 input file from the mole file
        try{
            File runningdir=new File(directory);
	    //this will only be run on Linux so we don't have to worry about Linux vs. Windows issues
	    String command = "python "+System.getenv("RMG")+"/scripts/MM4InputFileMaker.py ";
	    //first argument: input file path; for the first attempts, we will use UFF refined coordinates; if that doesn't work, (e.g. case of cyclopropene, InChI=1/C3H4/c1-2-3-1/h1-2H,3H2 OOXWYYGXTJLWHA-UHFFFAOYAJ) we will try crude coordinates (.cmol suffix)
	    if(attemptNumber<=scriptAttempts){
		command=command.concat(p_molfile.getPath() + " ");
	    }
	    else{
		command=command.concat(p_molfile.getCrudePath() + " ");
	    }
	    //second argument: output path
	    String inpfilepath=directory+"/"+name+".mm4";
	    command=command.concat(inpfilepath+ " ");
	    //third argument: molecule name (the augmented InChI)
	    command=command.concat(InChIaug+ " ");
	    //fourth argument: PYTHONPATH
	    command=command.concat(System.getenv("RMG")+"/source/MoleCoor");//this will pass $RMG/source/MoleCoor to the script (in order to get the appropriate path for importing
	    Process molecoorProc = Runtime.getRuntime().exec(command, null, runningdir);
            //read in output
            InputStream is = molecoorProc.getInputStream();
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line=null;
            while ( (line = br.readLine()) != null) {
		//do nothing
            }
            int exitValue = molecoorProc.waitFor();
        }
        catch(Exception e){
            String err = "Error in running MoleCoor MOL to .MM4 process \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }
        return maxAttemptNumber;
    }

    //creates MM4 rotor input file and MM4 batch file in directory with filenames name.mm4roti and name.comi, respectively
    //the function returns the set of rotor dihedral angles for the minimum energy conformation
    public double[] createMM4RotorInput(String name, String directory, ChemGraph p_chemgraph, int rotors){
	//read in the optimized coordinates from the completed "normal" MM4 job; this will be used as a template for the rotor input files
	String mm4optContents = "";
	double[] dihedralMinima = new double[rotors];
	try{
	    FileReader mm4opt = new FileReader(directory+"/"+name+".mm4opt");
	    BufferedReader reader = new BufferedReader(mm4opt);
	    String line=reader.readLine();
	    while(line!=null){
		mm4optContents+=line+"\n";
		line=reader.readLine();
	    }
	    mm4opt.close();
	}
	catch(Exception e){
		String err = "Error in reading .mm4opt file\n";
		err += e.toString();
		e.printStackTrace();
		System.exit(0);
	}
	String[] lines = mm4optContents.split("[\\r?\\n]+");//split by newlines, excluding blank lines; cf. http://stackoverflow.com/questions/454908/split-java-string-by-new-line
	int indexForFirstAtom = lines.length - p_chemgraph.getAtomNumber();//assumes the last line is for the last atom
	lines[1]=lines[1].substring(0,78)+" 2";//take the first 78 characters of line 2, and append the option number for the NDRIVE option; in other words, we are replacing the NDRIVE=0 option with the desired option number
	//reconstruct mm4optContents
	mm4optContents = "";
	for(int j=0; j<lines.length;j++){
	    mm4optContents +=lines[j]+"\n";
	}
	//iterate over all the rotors in the molecule
	int i = 0;//rotor index
	LinkedHashMap rotorInfo = p_chemgraph.getInternalRotorInformation();
	Iterator iter = rotorInfo.keySet().iterator();
	while(iter.hasNext()){
	    i++;
	    int[] rotorAtoms = (int[])iter.next();
	    try{
		//write one script file for each rotor
		//Step 1: write the script for MM4 batch operation
//	Example script file:
//	#! /bin/csh
//	cp testEthylene.mm4 CPD.MM4
//	cp $MM4_DATDIR/BLANK.DAT PARA.MM4
//	cp $MM4_DATDIR/CONST.MM4 .
//	$MM4_EXEDIR/mm4 <<%
//	1
//	2
//	0
//	%
//	mv TAPE4.MM4 testEthyleneBatch.out
//	mv TAPE9.MM4 testEthyleneBatch.opt
//	exit
		//create batch file with executable permissions: cf. http://java.sun.com/docs/books/tutorial/essential/io/fileAttr.html#posix
		File inpKey = new File(directory+"/"+name+".com"+i);
		String inpKeyStr="#! /bin/csh\n";
		inpKeyStr+="cp "+name+".mm4rot"+i+" CPD.MM4\n";
		inpKeyStr+="cp $MM4_DATDIR/BLANK.DAT PARA.MM4\n";
		inpKeyStr+="cp $MM4_DATDIR/CONST.MM4 .\n";
		inpKeyStr+="$MM4_EXEDIR/mm4 <<%\n";
		inpKeyStr+="1\n";//read from first line of .mm4 file
		inpKeyStr+="3\n"; //Full-Matrix Method only
		//inpKeyStr+="2\n"; //Block-Diagonal Method then Full-Matrix Method

		inpKeyStr+="0\n";//terminate the job
		inpKeyStr+="%\n";
		inpKeyStr+="mv TAPE4.MM4 "+name+".mm4rotout"+i+"\n";
		inpKeyStr+="mv TAPE9.MM4 "+name+".mm4rotopt"+i+"\n";
		inpKeyStr+="exit\n";
		FileWriter fw = new FileWriter(inpKey);
		fw.write(inpKeyStr);
		fw.close();
	    }
	    catch(Exception e){
		String err = "Error in writing MM4 script files\n";
		err += e.toString();
		e.printStackTrace();
		System.exit(0);
	    }

	    //Step 2: create the MM4 rotor job input file for rotor i, using the output file from the "normal" MM4 job as a template
	    //Step 2a: find the dihedral angle for the minimized conformation (we need to start at the minimum energy conformation because CanTherm assumes that theta=0 is the minimum
	    //extract the lines for dihedral atoms
	    String dihedral1s = lines[indexForFirstAtom+rotorAtoms[0]-1];
	    String atom1s = lines[indexForFirstAtom+rotorAtoms[1]-1];
	    String atom2s = lines[indexForFirstAtom+rotorAtoms[2]-1];
	    String dihedral2s = lines[indexForFirstAtom+rotorAtoms[3]-1];
	    //extract the x,y,z coordinates for each atom
	    double[] dihedral1 = {Double.parseDouble(dihedral1s.substring(0,10)),Double.parseDouble(dihedral1s.substring(10,20)),Double.parseDouble(dihedral1s.substring(20,30))};
	    double[] atom1 = {Double.parseDouble(atom1s.substring(0,10)),Double.parseDouble(atom1s.substring(10,20)),Double.parseDouble(atom1s.substring(20,30))};
	    double[] atom2 = {Double.parseDouble(atom2s.substring(0,10)),Double.parseDouble(atom2s.substring(10,20)),Double.parseDouble(atom2s.substring(20,30))};
	    double[] dihedral2 = {Double.parseDouble(dihedral2s.substring(0,10)),Double.parseDouble(dihedral2s.substring(10,20)),Double.parseDouble(dihedral2s.substring(20,30))};
	    //determine the dihedral angle
	    dihedralMinima[i-1] = calculateDihedral(dihedral1,atom1,atom2,dihedral2);
	    //double dihedral = calculateDihedral(dihedral1,atom1,atom2,dihedral2);
	    //if (dihedral < 0) dihedral = dihedral + 360;//make sure dihedral is positive; this will save an extra character in the limited space for specifying starting and ending angles
	    //eventually, when problems arise due to collinear atoms (both arguments to atan2 are zero) we can iterate over other atom combinations (with atoms in each piece determined by the corresponding value in rotorInfo for atom2 and the complement of these (full set = p_chemgraph.getNodeIDs()) for atom1) until they are not collinear

	    //Step 2b: write the file for rotor i
	    try{
		FileWriter mm4roti = new FileWriter(directory+"/"+name+".mm4rot"+i);
		//mm4roti.write(mm4optContents+"\n"+String.format("  %3d  %3d  %3d  %3d     %5.1f%5.1f%5.1f", rotorAtoms[0],rotorAtoms[1],rotorAtoms[2],rotorAtoms[3], dihedral, dihedral + 360.0 - deltaTheta, deltaTheta)+"\n");//deltaTheta should be less than 100 degrees so that dihedral-deltaTheta still fits
		//mm4roti.write(mm4optContents+"\n"+String.format("  %3d  %3d  %3d  %3d     %5.1f%5.1f%5.1f", rotorAtoms[0],rotorAtoms[1],rotorAtoms[2],rotorAtoms[3], dihedral, dihedral - deltaTheta, deltaTheta)+"\n");//deltaTheta should be less than 100 degrees so that dihedral-deltaTheta still fits
		mm4roti.write(mm4optContents+String.format("  %3d  %3d  %3d  %3d     %5.1f%5.1f%5.1f", rotorAtoms[0],rotorAtoms[1],rotorAtoms[2],rotorAtoms[3], 0.0, 360.0-deltaTheta, deltaTheta)+"\n");//M1, M2, M3, M4, START, FINISH, DIFF (as described in MM4 manual); //this would require miminum to be stored and used to adjust actual angles before sending to CanTherm
		mm4roti.close();
	    }
	    catch(Exception e){
		String err = "Error in writing MM4 rotor input file\n";
		err += e.toString();
		e.printStackTrace();
		System.exit(0);
	    }
	}
        return dihedralMinima;
    }
    
    //given x,y,z (cartesian) coordinates for dihedral1 atom, (rotor) atom 1, (rotor) atom 2, and dihedral2 atom, calculates the dihedral angle (in degrees, between +180 and -180) using the atan2 formula at http://en.wikipedia.org/w/index.php?title=Dihedral_angle&oldid=373614697
    public double calculateDihedral(double[] dihedral1, double[] atom1, double[] atom2, double[] dihedral2){
	//calculate the vectors between the atoms
	double[] b1 = {atom1[0]-dihedral1[0], atom1[1]-dihedral1[1], atom1[2]-dihedral1[2]};
	double[] b2 = {atom2[0]-atom1[0], atom2[1]-atom1[1], atom2[2]-atom1[2]};
	double[] b3 = {dihedral2[0]-atom2[0], dihedral2[1]-atom2[1], dihedral2[2]-atom2[2]};
	//calculate norm of b2
	double normb2 = Math.sqrt(b2[0]*b2[0]+b2[1]*b2[1]+b2[2]*b2[2]);
	//calculate necessary cross products
	double[] b1xb2 = {b1[1]*b2[2]-b1[2]*b2[1], b2[0]*b1[2]-b1[0]*b2[2], b1[0]*b2[1]-b1[1]*b2[0]};
	double[] b2xb3 = {b2[1]*b3[2]-b2[2]*b3[1], b3[0]*b2[2]-b2[0]*b3[2], b2[0]*b3[1]-b2[1]*b3[0]};
	//compute arguments for atan2 function (includes dot products)
	double y = normb2*(b1[0]*b2xb3[0]+b1[1]*b2xb3[1]+b1[2]*b2xb3[2]);//|b2|*b1.(b2xb3)
	double x = b1xb2[0]*b2xb3[0]+b1xb2[1]*b2xb3[1]+b1xb2[2]*b2xb3[2];//(b1xb2).(b2xb3)
	double dihedral = Math.atan2(y, x);
	//return dihedral*180.0/Math.PI;//converts from radians to degrees
	return Math.toDegrees(dihedral);//converts from radians to degrees
    }

    //returns the extra Mopac keywords to use for radical species, given the spin multiplicity (radical number + 1)
    public String getMopacRadicalString(int multiplicity){
        if(multiplicity==1) return "";
        else if (multiplicity==2) return "uhf doublet";
        else if (multiplicity==3) return "uhf triplet";
        else if (multiplicity==4) return "uhf quartet";
        else if (multiplicity==5) return "uhf quintet";
        else if (multiplicity==6) return "uhf sextet";
        else if (multiplicity==7) return "uhf septet";
        else if (multiplicity==8) return "uhf octet";
        else if (multiplicity==9) return "uhf nonet";
        else{
            System.out.println("Invalid multiplicity encountered: "+multiplicity);
            System.exit(0);
        }
        
        return "this should not be returned: error associated with getMopacRadicalString()";
    }
    
    //creates MOPAC PM3 input file in directory with filename name.mop by using OpenBabel to convert p_molfile
    //attemptNumber determines which keywords to try
    //the function returns the maximum number of keywords that can be attempted; this will be the same throughout the evaluation of the code, so it may be more appropriate to have this as a "constant" attribute of some sort
    //unlike createGaussianPM3 input, this requires an additional input specifying the spin multiplicity (radical number + 1) for the species
    public int createMopacPM3Input(String name, String directory, molFile p_molfile, int attemptNumber, String InChIaug, int multiplicity){
        //write a file with the input keywords
        int scriptAttempts = 5;//the number of keyword permutations available; update as additional options are added
        int maxAttemptNumber=2*scriptAttempts;//we will try a second time with crude coordinates if the UFF refined coordinates do not work
        String inpKeyStrBoth = "";//this string will be written at both the top (for optimization) and the bottom (for thermo/force calc)
        String inpKeyStrTop = "";//this string will be written only at the top
        String inpKeyStrBottom = "";//this string will be written at the bottom
        String radicalString = getMopacRadicalString(multiplicity);
        try{
      //      File inpKey=new File(directory+"/inputkeywords.txt");
      //      String inpKeyStr="%chk="+directory+"\\RMGrunCHKfile.chk\n";
      //      inpKeyStr+="%mem=6MW\n";
      //     inpKeyStr+="%nproc=1\n";

            if(attemptNumber%scriptAttempts==1){
                inpKeyStrBoth="pm3 "+radicalString;
                inpKeyStrTop=" precise nosym";
                inpKeyStrBottom="oldgeo thermo nosym precise ";//7/10/09: based on a quick review of recent results, keyword combo #1 rarely works, and when it did (CJAINEUZFLXGFA-UHFFFAOYAUmult3 (InChI=1/C8H16O5Si/c1-4-11-14(9,12-5-2)13-8-6-10-7(8)3/h7-8H,3-6H2,1-2H3/mult3)), the grad. norm on the force step was about 1.7 (too large); I manually removed this result and re-ran...the entropy was increased by nearly 20 cal/mol-K...perhaps we should add a check for the "WARNING" that MOPAC prints out when the gradient is high; 7/22/09: for the case of FUGDBSHZYPTWLG-UHFFFAOYADmult3 (InChI=1/C5H8/c1-4-3-5(4)2/h4-5H,1-3H2/mult3), adding nosym seemed to resolve 1. large discrepancies from Gaussian and 2. negative frequencies in mass-weighted coordinates and possibly related issue in discrepancies between regular and mass-weighted coordinate frequencies
            }
            else if(attemptNumber%scriptAttempts==2){//7/9/09: used for VCSJVABXVCFDRA-UHFFFAOYAI (InChI=1/C8H19O5Si/c1-5-10-8(4)13-14(9,11-6-2)12-7-3/h8H,5-7H2,1-4H3); all existing Gaussian keywords also failed; the Gaussian result was also rectified, but the resulting molecule was over 70 kcal/mol less stable, probably due to a large amount of spin contamination (~1.75 in fixed Gaussian result vs. 0.754 for MOPAC)
                inpKeyStrBoth="pm3 "+radicalString;
                inpKeyStrTop=" precise nosym gnorm=0.0 nonr";
                inpKeyStrBottom="oldgeo thermo nosym precise ";
            }
            else if(attemptNumber%scriptAttempts==3){//7/8/09: used for ADMPQLGIEMRGAT-UHFFFAOYAUmult3 (InChI=1/C6H14O5Si/c1-4-9-12(8,10-5-2)11-6(3)7/h6-7H,3-5H2,1-2H3/mult3); all existing Gaussian keywords also failed; however, the Gaussian result was also rectified, and the resulting conformation was about 1.0 kcal/mol more stable than the one resulting from this, so fixed Gaussian result was manually copied to QMFiles folder
                inpKeyStrBoth="pm3 "+radicalString;
                inpKeyStrTop=" precise nosym gnorm=0.0";
                inpKeyStrBottom="oldgeo thermo nosym precise "; //precise appeared to be necessary for the problematic case (to avoid negative frequencies); 
            }
            else if(attemptNumber%scriptAttempts==4){//7/8/09: used for GYFVJYRUZAKGFA-UHFFFAOYALmult3 (InChI=1/C6H14O6Si/c1-3-10-13(8,11-4-2)12-6-5-9-7/h6-7H,3-5H2,1-2H3/mult3) case (negative frequency issues in MOPAC) (also, none of the existing Gaussian combinations worked with it); note that the Gaussian result appears to be a different conformation as it is about 0.85 kcal/mol more stable, so the Gaussian result was manually copied to QMFiles directory; note that the MOPAC output included a very low frequency (4-5 cm^-1)
                inpKeyStrBoth="pm3 "+radicalString;
                inpKeyStrTop=" precise nosym gnorm=0.0 bfgs";
                inpKeyStrBottom="oldgeo thermo nosym precise "; //precise appeared to be necessary for the problematic case (to avoid negative frequencies)
            }
//            else if(attemptNumber==5){
//                inpKeyStrBoth="pm3 "+radicalString;
//                inpKeyStrTop=" precise nosym gnorm=0.0 ddmin=0.0";
//                inpKeyStrBottom="oldgeo thermo nosym precise ";
//            }
//            else if(attemptNumber==6){
//                inpKeyStrBoth="pm3 "+radicalString;
//                inpKeyStrTop=" precise nosym gnorm=0.0 nonr ddmin=0.0";
//                inpKeyStrBottom="oldgeo thermo nosym precise ";
//            }
//            else if(attemptNumber==7){
//                inpKeyStrBoth="pm3 "+radicalString;
//                inpKeyStrTop=" precise nosym bfgs gnorm=0.0 ddmin=0.0";
//                inpKeyStrBottom="oldgeo thermo nosym precise ";
//            }
            else if(attemptNumber%scriptAttempts==0){//used for troublesome HGRZRPHFLAXXBT-UHFFFAOYAVmult3 (InChI=1/C3H2O4/c4-2(5)1-3(6)7/h1H2/mult3) case (negative frequency and large gradient issues)
                inpKeyStrBoth="pm3 "+radicalString;
                inpKeyStrTop=" precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000";
                inpKeyStrBottom="oldgeo thermo nosym precise ";
            }
         //   else if(attemptNumber==9){//used for troublesome CMARQPBQDRXBTN-UHFFFAOYAAmult3 (InChI=1/C3H2O4/c1-3(5)7-6-2-4/h1H2/mult3) case (negative frequency issues)
         //       inpKeyStrBoth="pm3 "+radicalString;
         //       inpKeyStrTop=" precise nosym recalc=1 dmax=0.05 gnorm=0.0 cycles=1000 t=1000";
         //       inpKeyStrBottom="oldgeo thermo nosym precise ";
         //   }
         //   else if(attemptNumber==10){//used for ATCYLHQLTOSVFK-UHFFFAOYAMmult4 (InChI=1/C4H5O5/c1-3(5)8-9-4(2)7-6/h6H,1-2H2/mult4) case (timeout issue; also, negative frequency issues); note that this is very similar to the keyword below, so we may want to consolidate
         //       inpKeyStrBoth="pm3 "+radicalString;
         //       inpKeyStrTop=" precise nosym recalc=1 dmax=0.05 gnorm=0.2 cycles=1000 t=1000";
         //       inpKeyStrBottom="oldgeo thermo nosym precise ";
         //   }
            else throw new Exception();//this point should not be reached
      //      FileWriter fw = new FileWriter(inpKey);
      //      fw.write(inpKeyStr);
      //      fw.close();
        }
        catch(Exception e){
            String err = "Error in writing inputkeywords.txt \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }
        
        String polarString = "";
        if (usePolar){
            if(multiplicity == 1) polarString = System.getProperty("line.separator") + System.getProperty("line.separator") + System.getProperty("line.separator")+ "oldgeo polar nosym precise " + inpKeyStrBoth;
            else polarString = System.getProperty("line.separator") + System.getProperty("line.separator") + System.getProperty("line.separator")+ "oldgeo static nosym precise " + inpKeyStrBoth;
        }
        
        //call the OpenBabel process (note that this requires OpenBabel environment variable)
        try{ 
            File runningdir=new File(directory);
            String inpKeyStrTopCombined = inpKeyStrBoth + inpKeyStrTop;
	    String command = null;
	    if(attemptNumber<=scriptAttempts){//use UFF-refined coordinates
		command = "babel -imol "+ p_molfile.getPath()+ " -omop " + name+".mop -xk \"" + inpKeyStrTopCombined + "\" --title \""+InChIaug+"\"";
	    }
	    else{//use the crude coordinates
		command = "babel -imol "+ p_molfile.getCrudePath()+ " -omop " + name+".mop -xk \"" + inpKeyStrTopCombined + "\" --title \""+InChIaug+"\"";
	    }
	    Process babelProc = Runtime.getRuntime().exec(command, null, runningdir);
            //read in output
            InputStream is = babelProc.getInputStream();
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line=null;
            while ( (line = br.readLine()) != null) {
                //do nothing
            }
            int exitValue = babelProc.waitFor();
            
            //append the final keywords to the end of the file just written
           // File mopacInpFile = new File(directory+"/"+name+".mop");
            FileWriter fw = new FileWriter(directory+"/"+name+".mop", true);//filewriter with append = true
            fw.write(System.getProperty("line.separator") + inpKeyStrBottom + inpKeyStrBoth + polarString);//on Windows Vista, "\n" appears correctly in WordPad view, but not Notepad view (cf. http://forums.sun.com/thread.jspa?threadID=5386822)
            fw.close();
        }
        catch(Exception e){
            String err = "Error in running OpenBabel MOL to MOP process \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }
        return maxAttemptNumber;
    }
    
    //name and directory are the name and directory for the input (and output) file;
    //input is assumed to be preexisting and have the .gjf suffix
    //returns an integer indicating success or failure of the Gaussian calculation: 1 for success, 0 for failure;
    public int runGaussian(String name, String directory){
        int flag = 0;
        int successFlag=0;
        try{ 
            String command = "g03 ";
            command=command.concat(qmfolder+"/"+name+".gjf ");//specify the input file; space is important
	    command=command.concat(qmfolder+"/"+name+".log");//specify the output file
	    Process gaussianProc = Runtime.getRuntime().exec(command);

            //check for errors and display the error if there is one
            InputStream is = gaussianProc.getErrorStream();
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line=null;
            while ( (line = br.readLine()) != null) {
                    line = line.trim();
                    System.err.println(line);
                    flag=1;
            }
            //if there was an error, indicate that an error was obtained
            if(flag==1){
                System.out.println("Gaussian process received error (see above) on " + name);
            }
            int exitValue = gaussianProc.waitFor();
        }
        catch(Exception e){
            String err = "Error in running Gaussian process \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }
        //look in the output file to check for the successful termination of the Gaussian calculation
        //failed jobs will contain the a line beginning with " Error termination" near the end of the file
        int failureFlag=0;
        String errorLine = "";//string to store the error
        try{
            FileReader in = new FileReader(directory+"/"+name+".log");
            BufferedReader reader = new BufferedReader(in);
            String line=reader.readLine();
            while(line!=null){
                if (line.startsWith(" Error termination ")){
                    failureFlag=1;
                    errorLine = line.trim();
                    System.out.println("*****Error in Gaussian log file: "+errorLine);//print the error (note that in general, I think two lines will be printed)
                }
                else if (line.startsWith(" ******")){//also look for imaginary frequencies
                    if (line.contains("imaginary frequencies")){
                        System.out.println("*****Imaginary freqencies found:");
                        failureFlag=1;
                    }
                }
                line=reader.readLine();
            }
        }
        catch(Exception e){
            String err = "Error in reading Gaussian log file \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }
        //if the failure flag is still 0, the process should have been successful
        if (failureFlag==0) successFlag=1;
        
        return successFlag;
    }

        //name and directory are the name and directory for the input (and output) file;
    //input script is assumed to be preexisting and have the .com suffix
    //returns an integer indicating success or failure of the calculation: 1 for success, 0 for failure
    public int runMM4(String name, String directory){
        int successFlag=0;
        //int flag = 0;
	try{
	    File runningDirectory = new File(qmfolder);
            String command=name+".com";
	   // File script = new File(qmfolder+command);
	   // command = "./"+command;
	   // script.setExecutable(true);
	    command = "csh "+ command;
	    Process mm4Proc = Runtime.getRuntime().exec(command, null, runningDirectory);

	    //check for errors and display the error if there is one
//	    InputStream is = mm4Proc.getErrorStream();
//	    InputStreamReader isr = new InputStreamReader(is);
//	    BufferedReader br = new BufferedReader(isr);
//	    String line=null;
//	    while ( (line = br.readLine()) != null) {
//		line = line.trim();
//		if(!line.equals("STOP   statement executed")){//string listed here seems to be typical
//		    System.err.println(line);
//		    flag=1;
//		}
//	    }
	    InputStream is = mm4Proc.getInputStream();
	    InputStreamReader isr = new InputStreamReader(is);
	    BufferedReader br = new BufferedReader(isr);
	    String line=null;
	    while ( (line = br.readLine()) != null) {
		//do nothing
	    }
	    //if there was an error, indicate that an error was obtained
//	    if(flag==1){
//		System.out.println("MM4 process received error (see above) on " + name);
//	    }


            int exitValue = mm4Proc.waitFor();
        }
        catch(Exception e){
            String err = "Error in running MM4 process \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }
        //look in the output file to check for the successful termination of the MM4 calculation (cf. successfulMM4ResultExistsQ)
        File file = new File(directory+"/"+name+".mm4out");
	int failureFlag=1;//flag (1 or 0) indicating whether the MM4 job failed
	int failureOverrideFlag=0;//flag (1 or 0) to override success as measured by failureFlag
        if(file.exists()){//if the file exists, do further checks; otherwise, we will skip to final statement and return false
            try{
                FileReader in = new FileReader(file);
                BufferedReader reader = new BufferedReader(in);
                String line=reader.readLine();
                while(line!=null){
                    String trimLine= line.trim();
                    if (trimLine.equals("STATISTICAL THERMODYNAMICS ANALYSIS")){
                       failureFlag = 0;
                    }
                    else if (trimLine.endsWith("imaginary frequencies,")){//read the number of imaginary frequencies and make sure it is zero
                        String[] split = trimLine.split("\\s+");
			if (Integer.parseInt(split[3])>0){
			    System.out.println("*****Imaginary freqencies found:");
			    failureOverrideFlag=1;
			}
                    }
		    else if (trimLine.contains("             0.0     (fir )")){
			if (useCanTherm){//zero frequencies are only acceptable when CanTherm is used
			    System.out.println("*****Warning: zero freqencies found (values lower than 7.7 cm^-1 are rounded to zero in MM4 output); CanTherm should hopefully correct this:");
			}
			else{
			    System.out.println("*****Zero freqencies found:");
			    failureOverrideFlag=1;
			}
		    }
                    line=reader.readLine();
                }
            }
            catch(Exception e){
                String err = "Error in reading MM4 output file \n";
                err += e.toString();
                e.printStackTrace();
                System.exit(0);
            }
	}
        //if the failure flag is still 0, the process should have been successful
        if(failureOverrideFlag==1) failureFlag=1; //job will be considered a failure if there are imaginary frequencies or if job terminates to to excess time/cycles
        //if the failure flag is 0 and there are no negative frequencies, the process should have been successful
        if (failureFlag==0) successFlag=1;

        return successFlag;
    }

    //name and directory are the name and directory for the input (and output) file;
    //input script is assumed to be preexisting and have the .comi suffix where i is a number between 1 and rotors
    public void runMM4Rotor(String name, String directory, int rotors){
	 for(int j=1;j<=rotors;j++){
	    try{
		File runningDirectory = new File(qmfolder);
		String command=name+".com"+j;
	       // File script = new File(qmfolder+command);
	       // command = "./"+command;
	       // script.setExecutable(true);
		command = "csh "+ command;
		Process mm4Proc = Runtime.getRuntime().exec(command, null, runningDirectory);

		//check for errors and display the error if there is one
    //	    InputStream is = mm4Proc.getErrorStream();
    //	    InputStreamReader isr = new InputStreamReader(is);
    //	    BufferedReader br = new BufferedReader(isr);
    //	    String line=null;
    //	    while ( (line = br.readLine()) != null) {
    //		line = line.trim();
    //		if(!line.equals("STOP   statement executed")){//string listed here seems to be typical
    //		    System.err.println(line);
    //		    flag=1;
    //		}
    //	    }
		InputStream is = mm4Proc.getInputStream();
		InputStreamReader isr = new InputStreamReader(is);
		BufferedReader br = new BufferedReader(isr);
		String line=null;
		while ( (line = br.readLine()) != null) {
		    //do nothing
		}
		//if there was an error, indicate that an error was obtained
    //	    if(flag==1){
    //		System.out.println("MM4 process received error (see above) on " + name);
    //	    }


		int exitValue = mm4Proc.waitFor();
	    }
	    catch(Exception e){
		String err = "Error in running MM4 rotor process \n";
		err += e.toString();
		e.printStackTrace();
		System.exit(0);
	    }
	}

        return;
    }
    
    //name and directory are the name and directory for the input (and output) file;
    //input is assumed to be preexisting and have the .mop suffix
    //returns an integer indicating success or failure of the MOPAC calculation: 1 for success, 0 for failure;
    //this function is based on the Gaussian analogue
    public int runMOPAC(String name, String directory){
        int flag = 0;
        int successFlag=0;
        try{ 
            String command = System.getenv("MOPAC_LICENSE")+"MOPAC2009.exe ";
            command=command.concat(directory+"/"+name+".mop ");//specify the input file; space is important
            command=command.concat(directory+"/"+name+".out");//specify the output file
            Process mopacProc = Runtime.getRuntime().exec(command);
            //check for errors and display the error if there is one
            InputStream is = mopacProc.getErrorStream();
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line=null;
            while ( (line = br.readLine()) != null) {
                    line = line.trim();
                    System.err.println(line);
                    flag=1;
            }
            //if there was an error, indicate that an error was obtained
            if(flag==1){
                System.out.println("MOPAC process received error (see above) on " + name);
            }
            int exitValue = mopacProc.waitFor();
        }
        catch(Exception e){
            String err = "Error in running MOPAC process \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }
        //look in the output file to check for the successful termination of the calculation (this is a trimmed down version of what appears in successfulMOPACResultExistsQ (it doesn't have the InChI check)
        File file = new File(directory+"/"+name+".out");
        int failureFlag=1;//flag (1 or 0) indicating whether the MOPAC job failed
        int failureOverrideFlag=0;//flag (1 or 0) to override success as measured by failureFlag
        if(file.exists()){//if the file exists, do further checks; otherwise, we will skip to final statement and return false
            try{
                FileReader in = new FileReader(file);
                BufferedReader reader = new BufferedReader(in);
                String line=reader.readLine();
                while(line!=null){
                    String trimLine = line.trim();
                    if (trimLine.equals("DESCRIPTION OF VIBRATIONS")){//check for this line; if it is here, check for negative frequencies
                        //if(!MopacFileContainsNegativeFreqsQ(name, directory)) failureFlag=0;
                        failureFlag=0;
                    }
                    //negative frequencies notice example:
//         NOTE: SYSTEM IS NOT A GROUND STATE, THEREFORE ZERO POINT
//         ENERGY IS NOT MEANINGFULL. ZERO POINT ENERGY PRINTED
//         DOES NOT INCLUDE THE  2 IMAGINARY FREQUENCIES
                    else if (trimLine.endsWith("IMAGINARY FREQUENCIES")){
                        System.out.println("*****Imaginary freqencies found:");
                        failureOverrideFlag=1;
                    }
                    else if (trimLine.equals("EXCESS NUMBER OF OPTIMIZATION CYCLES")){//exceeding max cycles error
                        failureOverrideFlag=1;
                    }
                    else if (trimLine.equals("NOT ENOUGH TIME FOR ANOTHER CYCLE")){//timeout error
                        failureOverrideFlag=1;
                    }
                    line=reader.readLine();
                }
            }
            catch(Exception e){
                String err = "Error in reading MOPAC output file \n";
                err += e.toString();
                e.printStackTrace();
                System.exit(0);
            }
        }
        if(failureOverrideFlag==1) failureFlag=1; //job will be considered a failure if there are imaginary frequencies or if job terminates to to excess time/cycles
        //if the failure flag is 0 and there are no negative frequencies, the process should have been successful
        if (failureFlag==0) successFlag=1;
        
        return successFlag;
    }
    
    //parse the results using cclib and return a ThermoData object; name and directory indicate the location of the Gaussian .log file
    //may want to split this into several functions
    public ThermoData parseGaussianPM3(String name, String directory, ChemGraph p_chemGraph){
//        //parse the Gaussian file using cclib
//        int natoms = 0; //number of atoms from Gaussian file; in principle, this should agree with number of chemGraph atoms
//        ArrayList atomicNumber = new ArrayList(); //vector of atomic numbers (integers) (apparently Vector is thread-safe; cf. http://answers.yahoo.com/question/index?qid=20081214065127AArZDT3; ...should I be using this instead?)
//        ArrayList x_coor = new ArrayList(); //vectors of x-, y-, and z-coordinates (doubles) (Angstroms) (in order corresponding to above atomic numbers)
//        ArrayList y_coor = new ArrayList();
//        ArrayList z_coor = new ArrayList();
//        double energy = 0; //PM3 energy (Hf298) in Hartree
//        double molmass = 0; //molecular mass in amu
//        ArrayList freqs = new ArrayList(); //list of frequencies in units of cm^-1
//        double rotCons_1 = 0;//rotational constants in (1/s)
//        double rotCons_2 = 0;
//        double rotCons_3 = 0; 
//        int gdStateDegen = p_chemGraph.getRadicalNumber()+1;//calculate ground state degeneracy from the number of radicals; this should give the same result as spin multiplicity in Gaussian input file (and output file), but we do not explicitly check this (we could use "mult" which cclib reads in if we wanted to do so); also, note that this is not always correct, as there can apparently be additional spatial degeneracy for non-symmetric linear molecules like OH radical (cf. http://cccbdb.nist.gov/thermo.asp)
//        try{   
//            File runningdir=new File(directory);
//            String command = "c:/Python25/python.exe c:/Python25/GaussianPM3ParsingScript.py ";//this should eventually be modified for added generality
//            String logfilepath=directory+"/"+name+".log";
//            command=command.concat(logfilepath);
//            Process cclibProc = Runtime.getRuntime().exec(command, null, runningdir);
//            //read the stdout of the process, which should contain the desired information in a particular format
//            InputStream is = cclibProc.getInputStream();
//            InputStreamReader isr = new InputStreamReader(is);
//            BufferedReader br = new BufferedReader(isr);
//            String line=null;
//            //example output:
////            C:\Python25>python.exe GaussianPM3ParsingScript.py TEOS.out
////            33
////            [ 6  6  8 14  8  6  6  8  6  6  8  6  6  1  1  1  1  1  1  1  1  1  1  1  1
////              1  1  1  1  1  1  1  1]
////            [[ 2.049061 -0.210375  3.133106]
////             [ 1.654646  0.321749  1.762752]
////             [ 0.359284 -0.110429  1.471465]
////             [-0.201871 -0.013365 -0.12819 ]
////             [ 0.086307  1.504918 -0.82893 ]
////             [-0.559186  2.619928 -0.284003]
////             [-0.180246  3.839463 -1.113029]
////             [ 0.523347 -1.188305 -1.112765]
////             [ 1.857584 -1.018167 -1.495088]
////             [ 2.375559 -2.344392 -2.033403]
////             [-1.870397 -0.297297 -0.075427]
////             [-2.313824 -1.571765  0.300245]
////             [-3.83427  -1.535927  0.372171]
////             [ 1.360346  0.128852  3.917699]
////             [ 2.053945 -1.307678  3.160474]
////             [ 3.055397  0.133647  3.403037]
////             [ 1.677262  1.430072  1.750899]
////             [ 2.372265 -0.029237  0.985204]
////             [-0.245956  2.754188  0.771433]
////             [-1.656897  2.472855 -0.287156]
////             [-0.664186  4.739148 -0.712606]
////             [-0.489413  3.734366 -2.161038]
////             [ 0.903055  4.016867 -1.112198]
////             [ 1.919521 -0.229395 -2.269681]
////             [ 2.474031 -0.680069 -0.629949]
////             [ 2.344478 -3.136247 -1.273862]
////             [ 1.786854 -2.695974 -2.890647]
////             [ 3.41648  -2.242409 -2.365094]
////             [-1.884889 -1.858617  1.28054 ]
////             [-1.976206 -2.322432 -0.440995]
////             [-4.284706 -1.26469  -0.591463]
////             [-4.225999 -2.520759  0.656131]
////             [-4.193468 -0.809557  1.112677]]
////            -14.1664924726
////            [    9.9615    18.102     27.0569    31.8459    39.0096    55.0091
////                66.4992    80.4552    86.4912   123.3551   141.6058   155.5448
////               159.4747   167.0013   178.5676   207.3738   237.3201   255.3487
////               264.5649   292.867    309.4248   344.6503   434.8231   470.2074
////               488.9717   749.1722   834.257    834.6594   837.7292   839.6352
////               887.9767   892.9538   899.5374   992.1851  1020.6164  1020.8671
////              1028.3897  1046.7945  1049.1768  1059.4704  1065.1505  1107.4001
////              1108.1567  1109.0466  1112.6677  1122.7785  1124.4315  1128.4163
////              1153.3438  1167.6705  1170.9627  1174.9613  1232.1826  1331.8459
////              1335.3932  1335.8677  1343.9556  1371.37    1372.8127  1375.5428
////              1396.0344  1402.4082  1402.7554  1403.2463  1403.396   1411.6946
////              1412.2456  1412.3519  1414.5982  1415.3613  1415.5698  1415.7993
////              1418.5409  2870.7446  2905.3132  2907.0361  2914.1662  2949.2646
////              2965.825   2967.7667  2971.5223  3086.3849  3086.3878  3086.6448
////              3086.687   3089.2274  3089.4105  3089.4743  3089.5841  3186.0753
////              3186.1375  3186.3511  3186.365 ]
////            [ 0.52729  0.49992  0.42466]
////note: above example has since been updated to print molecular mass; also frequency and atomic number format has been updated
//            String [] stringArray;
//            natoms = Integer.parseInt(br.readLine());//read line 1: number of atoms
//            stringArray = br.readLine().replace("[", "").replace("]","").trim().split(",\\s+");//read line 2: the atomic numbers (first removing braces)
//           // line = br.readLine().replace("[", "").replace("]","");//read line 2: the atomic numbers (first removing braces)
//           // StringTokenizer st = new StringTokenizer(line); //apprently the stringTokenizer class is deprecated, but I am having trouble getting the regular expressions to work properly
//            for(int i=0; i < natoms; i++){
//               // atomicNumber.add(i,Integer.parseInt(stringArray[i]));
//                atomicNumber.add(i,Integer.parseInt(stringArray[i]));
//            }
//            for(int i=0; i < natoms; i++){
//                stringArray = br.readLine().replace("[", "").replace("]","").trim().split("\\s+");//read line 3+: coordinates for atom i; used /s+ for split; using spaces with default limit of 0 was giving empty string
//                x_coor.add(i,Double.parseDouble(stringArray[0]));
//                y_coor.add(i,Double.parseDouble(stringArray[1]));
//                z_coor.add(i,Double.parseDouble(stringArray[2]));
//            }
//            energy = Double.parseDouble(br.readLine());//read next line: energy
//            molmass = Double.parseDouble(br.readLine());//read next line: molecular mass (in amu)
//            if (natoms>1){//read additional info for non-monoatomic species
//                stringArray = br.readLine().replace("[", "").replace("]","").trim().split(",\\s+");//read next line: frequencies
//                for(int i=0; i < stringArray.length; i++){
//                    freqs.add(i,Double.parseDouble(stringArray[i]));
//                }
//                stringArray = br.readLine().replace("[", "").replace("]","").trim().split("\\s+");//read next line rotational constants (converting from GHz to Hz in the process)
//                rotCons_1 = Double.parseDouble(stringArray[0])*1000000000;
//                rotCons_2 = Double.parseDouble(stringArray[1])*1000000000;
//                rotCons_3 = Double.parseDouble(stringArray[2])*1000000000;
//            }
//            while ( (line = br.readLine()) != null) {
//                //do nothing (there shouldn't be any more information, but this is included to get all the output)
//            }
//            int exitValue = cclibProc.waitFor();
//        }
//        catch (Exception e) {
//            String err = "Error in running ccLib Python process \n";
//            err += e.toString();
//            e.printStackTrace();
//            System.exit(0);
//        } 
//   
//        ThermoData result = calculateThermoFromPM3Calc(natoms, atomicNumber, x_coor, y_coor, z_coor, energy, molmass, freqs, rotCons_1, rotCons_2, rotCons_3, gdStateDegen);
//        System.out.println("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes
//        return result;
        
        String command = null;
        if (System.getProperty("os.name").toLowerCase().contains("windows")){//special windows case where paths can have spaces and are allowed to be surrounded by quotes
	    command = "python \""+ System.getProperty("RMG.workingDirectory")+"/scripts/GaussianPM3ParsingScript.py\" ";
	    String logfilepath="\""+directory+"/"+name+".log\"";
	    command=command.concat(logfilepath);
	    command=command.concat(" \""+ System.getenv("RMG")+"/source\"");//this will pass $RMG/source to the script (in order to get the appropriate path for importing
	}
	else{//non-Windows case
	    command = "python "+ System.getProperty("RMG.workingDirectory")+"/scripts/GaussianPM3ParsingScript.py ";
	    String logfilepath=directory+"/"+name+".log";
	    command=command.concat(logfilepath);
	    command=command.concat(" "+ System.getenv("RMG")+"/source");//this will pass $RMG/source to the script (in order to get the appropriate path for importing
	}
	ThermoData result = getPM3MM4ThermoDataUsingCCLib(name, directory, p_chemGraph, command);
        System.out.println("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes
        return result;
    }

    //parse the results using cclib and return a ThermoData object; name and directory indicate the location of the MM4 .mm4out file
    public ThermoData parseMM4(String name, String directory, ChemGraph p_chemGraph){
	String command = "python "+System.getProperty("RMG.workingDirectory")+"/scripts/MM4ParsingScript.py ";
	String logfilepath=directory+"/"+name+".mm4out";
	command=command.concat(logfilepath);
	command=command.concat(" "+ System.getenv("RMG")+"/source");//this will pass $RMG/source to the script (in order to get the appropriate path for importing
        ThermoData result = getPM3MM4ThermoDataUsingCCLib(name, directory, p_chemGraph, command);
        System.out.println("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes
        return result;
    }

    //parse the results using cclib and CanTherm and return a ThermoData object; name and directory indicate the location of the MM4 .mm4out file
    //formerly known as parseMM4withForceMat
    public QMData performCanThermCalcs(String name, String directory, ChemGraph p_chemGraph, double[] dihedralMinima, boolean forceRRHO){
	//1. parse the MM4 file with cclib to get atomic number vector and geometry
	QMData qmdata = getQMDataWithCClib(name, directory, p_chemGraph, true);
	//unpack the needed results
	double energy = qmdata.energy;
	double stericEnergy = qmdata.stericEnergy;
	ArrayList freqs = qmdata.freqs;
	//2. compute H0/E0;  note that we will compute H0 for CanTherm by H0=Hf298(harmonicMM4)-(H298-H0)harmonicMM4, where harmonicMM4 values come from cclib parsing;  298.16 K is the standard temperature used by MM4; also note that Hthermal should include the R*T contribution (R*298.16 (enthalpy vs. energy difference)) and H0=E0 (T=0)	
	double T_MM4 = 298.16;
	energy *= Hartree_to_kcal;//convert from Hartree to kcal/mol
	stericEnergy *= Hartree_to_kcal;//convert from Hartree to kcal/mol
	double Hthermal = 5./2.*R*T_MM4/1000.;//enthalpy vs. energy difference(RT)+translation(3/2*R*T) contribution to thermal energy
	//rotational contribution
	if(p_chemGraph.getAtomNumber()==1) Hthermal += 0.0;
	else if (p_chemGraph.isLinear()) Hthermal += R*T_MM4/1000.;
	else Hthermal += 3./2.*R*T_MM4/1000;
	//vibrational contribution
	if(p_chemGraph.getAtomNumber()!=1)Hthermal+=R*calcVibH(freqs, T_MM4, h, k, c)/1000.;
	energy = energy - Hthermal;
	//3. write CanTherm input file
	//determine point group using the SYMMETRY Program
	String geom = qmdata.getSYMMETRYinput();
	String pointGroup = determinePointGroupUsingSYMMETRYProgram(geom);
	double sigmaCorr = getSigmaCorr(pointGroup);
	String canInp = "Calculation: Thermo\n";
	canInp += "Trange: 300 100 13\n";//temperatures from 300 to 1500 in increments of 100
	canInp += "Scale: 1.0\n";//scale factor of 1
	canInp += "Mol 1\n";
	if(p_chemGraph.getAtomNumber()==1) canInp += "ATOM\n";
	else if (p_chemGraph.isLinear()) canInp+="LINEAR\n";
	else canInp+="NONLINEAR\n";
	canInp += "GEOM MM4File " + name+".mm4out\n";//geometry file; ***special MM4 treatment in CanTherm; another option would be to use mm4opt file, but CanTherm code would need to be modified accordingly
	canInp += "FORCEC MM4File "+name+".fmat\n";//force constant file; ***special MM4 treatment in CanTherm
	if (forceRRHO) name = name + "_RRHO"; //"_RRHO" will be appended to InChIKey for RRHO calcs (though we want to use unmodified name in getQMDataWithCClib above, and in GEOM and FORCEC sections
	canInp += "ENERGY "+ energy +" MM4\n";//***special MM4 treatment in CanTherm
	canInp+="EXTSYM "+Math.exp(-sigmaCorr)+"\n";//***modified treatment in CanTherm; traditional EXTSYM integer replaced by EXTSYM double, to allow fractional values that take chirality into account
	canInp+="NELEC 1\n";//multiplicity = 1; all cases we consider will be non-radicals
	int rotors = p_chemGraph.getInternalRotor();
	String rotInput = null;
	if(!useHindRot || rotors==0 || forceRRHO) canInp += "ROTORS 0\n";//do not consider hindered rotors
	else{
	    int rotorCount = 0;
	    canInp+="ROTORS "+rotors+ " "+name+".rotinfo\n";
	    canInp+="POTENTIAL separable mm4files_inertia\n";//***special MM4 treatment in canTherm;
	    String rotNumbersLine=""+stericEnergy;//the line that will contain V0 (kcal/mol), and all the dihedral minima (degrees)
	    rotInput = "L1: 1 2 3\n";
	    LinkedHashMap rotorInfo = p_chemGraph.getInternalRotorInformation();
	    Iterator iter = rotorInfo.keySet().iterator();
	    while(iter.hasNext()){
		rotorCount++;
		int[] rotorAtoms = (int[])iter.next();
		LinkedHashSet rotatingGroup = (LinkedHashSet)rotorInfo.get(rotorAtoms);
		Iterator iterGroup = rotatingGroup.iterator();
		rotInput += "L2: " + rotorAtoms[4] +" "+rotorAtoms[1]+ " "+ rotorAtoms[2];//note: rotorAtoms[4] is the rotorSymmetry number as estimated by calculateRotorSymmetryNumber; this will be taken into account elsewhere
		while (iterGroup.hasNext()){//print the atoms associated with rotorAtom 2
		    Integer id = (Integer)iterGroup.next();
		    if(id != rotorAtoms[2]) rotInput+= " "+id;//don't print atom2
		}
		rotInput += "\n";
		canInp+=name + ".mm4rotopt"+rotorCount+" ";// potential files will be named as name.mm4rotopti
		rotNumbersLine+=" "+dihedralMinima[rotorCount-1];
	    }
	    canInp+="\n"+rotNumbersLine+"\n";
	}
	canInp+="0\n";//no bond-additivity corrections
	try{
            File canFile=new File(directory+"/"+name+".can");
            FileWriter fw = new FileWriter(canFile);
            fw.write(canInp);
            fw.close();
	    if(rotInput !=null){//write the rotor information
		File rotFile=new File(directory+"/"+name+".rotinfo");
		FileWriter fwr = new FileWriter(rotFile);
		fwr.write(rotInput);
		fwr.close();
	    }
        }
        catch(Exception e){
            String err = "Error in writing CanTherm input \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }
	//4 call CanTherm 
	try{
	    File runningDirectory = new File(qmfolder);
            String canCommand="python " + System.getenv("RMG")+"/source/CanTherm/source/CanTherm.py "+name+".can";
	    Process canProc = Runtime.getRuntime().exec(canCommand, null, runningDirectory);
	    InputStream is = canProc.getInputStream();
	    InputStreamReader isr = new InputStreamReader(is);
	    BufferedReader br = new BufferedReader(isr);
	    String line=null;
	    while ( (line = br.readLine()) != null) {
		    //do nothing
	    }

            int exitValue = canProc.waitFor();
        }
        catch(Exception e){
            String err = "Error in running CanTherm process \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }
        return qmdata;
    }

    public ThermoData parseCanThermFile(String name, String directory, ChemGraph p_chemGraph){
	//5. read CanTherm output
	Double Hf298 = null;
	Double S298 = null;
	Double Cp300 = null;
	Double Cp400 = null;
	Double Cp500 = null;
	Double Cp600 = null;
	Double Cp800 = null;
	Double Cp1000 = null;
	Double Cp1500 = null;
	File file = new File(directory+"/"+name+".canout");
        try{
	    FileReader in = new FileReader(file);
	    BufferedReader reader = new BufferedReader(in);
	    String line=reader.readLine();
	    while(!line.startsWith("Hf298 S298 Cps:")){//get to the end of the file with the data we want
		line=reader.readLine();
	    }
	    String[] split = reader.readLine().trim().split("\\s+");//read the next line, which contains the information we want
	    Hf298 = Double.parseDouble(split[0]);
	    S298 = Double.parseDouble(split[1]);
	    Cp300 = Double.parseDouble(split[2]);
	    Cp400 = Double.parseDouble(split[3]);
	    Cp500 = Double.parseDouble(split[4]);
	    Cp600 = Double.parseDouble(split[5]);
	    Cp800 = Double.parseDouble(split[7]);
	    Cp1000 = Double.parseDouble(split[9]);
	    Cp1500 = Double.parseDouble(split[14]);
	    reader.close();
	}
	catch(Exception e){
            String err = "Error in reading CanTherm .canout file \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }

	ThermoData result = new ThermoData(Hf298,S298,Cp300,Cp400,Cp500,Cp600,Cp800,Cp1000,Cp1500,3,1,1,"MM4 calculation; includes CanTherm analysis of force-constant matrix");//this includes rough estimates of uncertainty
        System.out.println("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes

	return result;
    }

    //parse the results using cclib and return a ThermoData object; name and directory indicate the location of the MOPAC .out file
    public ThermoData parseMopacPM3(String name, String directory, ChemGraph p_chemGraph){
        String command=null;
        if (System.getProperty("os.name").toLowerCase().contains("windows")){//special windows case where paths can have spaces and are allowed to be surrounded by quotes
	    command = "python \""+System.getProperty("RMG.workingDirectory")+"/scripts/MopacPM3ParsingScript.py\" ";
	    String logfilepath="\""+directory+"/"+name+".out\"";
	    command=command.concat(logfilepath);
	    command=command.concat(" \""+ System.getenv("RMG")+"/source\"");//this will pass $RMG/source to the script (in order to get the appropriate path for importing
	}
	else{//non-Windows case
	    command = "python "+System.getProperty("RMG.workingDirectory")+"/scripts/MopacPM3ParsingScript.py ";
	    String logfilepath=directory+"/"+name+".out";
	    command=command.concat(logfilepath);
	    command=command.concat(" "+ System.getenv("RMG")+"/source");//this will pass $RMG/source to the script (in order to get the appropriate path for importing
	}
        ThermoData result = getPM3MM4ThermoDataUsingCCLib(name, directory, p_chemGraph, command);
        System.out.println("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes
        return result;
    }
    
    //separated from parseMopacPM3, since the function was originally based off of parseGaussianPM3 and was very similar (differences being command and logfilepath variables);
    public ThermoData getPM3MM4ThermoDataUsingCCLib(String name, String directory, ChemGraph p_chemGraph, String command){
        //parse the Mopac file using cclib
        int natoms = 0; //number of atoms from Mopac file; in principle, this should agree with number of chemGraph atoms
        ArrayList atomicNumber = new ArrayList(); //vector of atomic numbers (integers) (apparently Vector is thread-safe; cf. http://answers.yahoo.com/question/index?qid=20081214065127AArZDT3; ...should I be using this instead?)
        ArrayList x_coor = new ArrayList(); //vectors of x-, y-, and z-coordinates (doubles) (Angstroms) (in order corresponding to above atomic numbers)
        ArrayList y_coor = new ArrayList();
        ArrayList z_coor = new ArrayList();
        double energy = 0; //PM3 energy (Hf298) in Hartree (***note: in the case of MOPAC, the MOPAC file will contain in units of kcal/mol, but modified ccLib will return in Hartree)
        double molmass = 0; //molecular mass in amu
        ArrayList freqs = new ArrayList(); //list of frequencies in units of cm^-1
        double rotCons_1 = 0;//rotational constants in (1/s)
        double rotCons_2 = 0;
        double rotCons_3 = 0; 
        int gdStateDegen = p_chemGraph.getRadicalNumber()+1;//calculate ground state degeneracy from the number of radicals; this should give the same result as spin multiplicity in Gaussian input file (and output file), but we do not explicitly check this (we could use "mult" which cclib reads in if we wanted to do so); also, note that this is not always correct, as there can apparently be additional spatial degeneracy for non-symmetric linear molecules like OH radical (cf. http://cccbdb.nist.gov/thermo.asp)
        try{   
            Process cclibProc = Runtime.getRuntime().exec(command);
            //read the stdout of the process, which should contain the desired information in a particular format
            InputStream is = cclibProc.getInputStream();
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);

            String line=null;
            //example output:
//            C:\Python25>python.exe GaussianPM3ParsingScript.py TEOS.out
//            33
//            [ 6  6  8 14  8  6  6  8  6  6  8  6  6  1  1  1  1  1  1  1  1  1  1  1  1
//              1  1  1  1  1  1  1  1]
//            [[ 2.049061 -0.210375  3.133106]
//             [ 1.654646  0.321749  1.762752]
//             [ 0.359284 -0.110429  1.471465]
//             [-0.201871 -0.013365 -0.12819 ]
//             [ 0.086307  1.504918 -0.82893 ]
//             [-0.559186  2.619928 -0.284003]
//             [-0.180246  3.839463 -1.113029]
//             [ 0.523347 -1.188305 -1.112765]
//             [ 1.857584 -1.018167 -1.495088]
//             [ 2.375559 -2.344392 -2.033403]
//             [-1.870397 -0.297297 -0.075427]
//             [-2.313824 -1.571765  0.300245]
//             [-3.83427  -1.535927  0.372171]
//             [ 1.360346  0.128852  3.917699]
//             [ 2.053945 -1.307678  3.160474]
//             [ 3.055397  0.133647  3.403037]
//             [ 1.677262  1.430072  1.750899]
//             [ 2.372265 -0.029237  0.985204]
//             [-0.245956  2.754188  0.771433]
//             [-1.656897  2.472855 -0.287156]
//             [-0.664186  4.739148 -0.712606]
//             [-0.489413  3.734366 -2.161038]
//             [ 0.903055  4.016867 -1.112198]
//             [ 1.919521 -0.229395 -2.269681]
//             [ 2.474031 -0.680069 -0.629949]
//             [ 2.344478 -3.136247 -1.273862]
//             [ 1.786854 -2.695974 -2.890647]
//             [ 3.41648  -2.242409 -2.365094]
//             [-1.884889 -1.858617  1.28054 ]
//             [-1.976206 -2.322432 -0.440995]
//             [-4.284706 -1.26469  -0.591463]
//             [-4.225999 -2.520759  0.656131]
//             [-4.193468 -0.809557  1.112677]]
//            -14.1664924726
//            [    9.9615    18.102     27.0569    31.8459    39.0096    55.0091
//                66.4992    80.4552    86.4912   123.3551   141.6058   155.5448
//               159.4747   167.0013   178.5676   207.3738   237.3201   255.3487
//               264.5649   292.867    309.4248   344.6503   434.8231   470.2074
//               488.9717   749.1722   834.257    834.6594   837.7292   839.6352
//               887.9767   892.9538   899.5374   992.1851  1020.6164  1020.8671
//              1028.3897  1046.7945  1049.1768  1059.4704  1065.1505  1107.4001
//              1108.1567  1109.0466  1112.6677  1122.7785  1124.4315  1128.4163
//              1153.3438  1167.6705  1170.9627  1174.9613  1232.1826  1331.8459
//              1335.3932  1335.8677  1343.9556  1371.37    1372.8127  1375.5428
//              1396.0344  1402.4082  1402.7554  1403.2463  1403.396   1411.6946
//              1412.2456  1412.3519  1414.5982  1415.3613  1415.5698  1415.7993
//              1418.5409  2870.7446  2905.3132  2907.0361  2914.1662  2949.2646
//              2965.825   2967.7667  2971.5223  3086.3849  3086.3878  3086.6448
//              3086.687   3089.2274  3089.4105  3089.4743  3089.5841  3186.0753
//              3186.1375  3186.3511  3186.365 ]
//            [ 0.52729  0.49992  0.42466]
//note: above example has since been updated to print molecular mass; also frequency and atomic number format has been updated
            String [] stringArray;
            natoms = Integer.parseInt(br.readLine());//read line 1: number of atoms
            stringArray = br.readLine().replace("[", "").replace("]","").trim().split(",\\s+");//read line 2: the atomic numbers (first removing braces)
           // line = br.readLine().replace("[", "").replace("]","");//read line 2: the atomic numbers (first removing braces)
           // StringTokenizer st = new StringTokenizer(line); //apprently the stringTokenizer class is deprecated, but I am having trouble getting the regular expressions to work properly
            for(int i=0; i < natoms; i++){
               // atomicNumber.add(i,Integer.parseInt(stringArray[i]));
                atomicNumber.add(i,Integer.parseInt(stringArray[i]));
            }
            for(int i=0; i < natoms; i++){
                stringArray = br.readLine().replace("[", "").replace("]","").trim().split("\\s+");//read line 3+: coordinates for atom i; used /s+ for split; using spaces with default limit of 0 was giving empty string
                x_coor.add(i,Double.parseDouble(stringArray[0]));
                y_coor.add(i,Double.parseDouble(stringArray[1]));
                z_coor.add(i,Double.parseDouble(stringArray[2]));
            }
            energy = Double.parseDouble(br.readLine());//read next line: energy
            molmass = Double.parseDouble(br.readLine());//read next line: molecular mass (in amu)
            if (natoms>1){//read additional info for non-monoatomic species
                stringArray = br.readLine().replace("[", "").replace("]","").trim().split(",\\s+");//read next line: frequencies
                for(int i=0; i < stringArray.length; i++){
                    freqs.add(i,Double.parseDouble(stringArray[i]));
                }
                stringArray = br.readLine().replace("[", "").replace("]","").trim().split("\\s+");//read next line rotational constants (converting from GHz to Hz in the process)
                rotCons_1 = Double.parseDouble(stringArray[0])*1000000000;
                rotCons_2 = Double.parseDouble(stringArray[1])*1000000000;
                rotCons_3 = Double.parseDouble(stringArray[2])*1000000000;
            }
            while ( (line = br.readLine()) != null) {
                //do nothing (there shouldn't be any more information, but this is included to get all the output)
            }
            int exitValue = cclibProc.waitFor();
        }
        catch (Exception e) {
            String err = "Error in running ccLib Python process \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        } 
   
        ThermoData result = calculateThermoFromPM3MM4Calc(natoms, atomicNumber, x_coor, y_coor, z_coor, energy, molmass, freqs, rotCons_1, rotCons_2, rotCons_3, gdStateDegen);
        return result;
    }
    //returns a thermo result, given results from quantum PM3 calculation or MM4 calculation (originally, this was in parseGaussianPM3 function
    public ThermoData calculateThermoFromPM3MM4Calc(int natoms, ArrayList atomicNumber, ArrayList x_coor, ArrayList y_coor, ArrayList z_coor, double energy, double molmass, ArrayList freqs, double rotCons_1, double rotCons_2, double rotCons_3, int gdStateDegen){
        //determine point group using the SYMMETRY Program
        String geom = natoms + "\n";
        for(int i=0; i < natoms; i++){
            geom += atomicNumber.get(i) + " "+ x_coor.get(i) + " " + y_coor.get(i) + " " +z_coor.get(i) + "\n";
        }
       // String pointGroup = determinePointGroupUsingSYMMETRYProgram(geom, 0.01);
        String pointGroup = determinePointGroupUsingSYMMETRYProgram(geom);
        
        //calculate thermo quantities using stat. mech. equations        
        //boolean linearity = p_chemGraph.isLinear();//determine linearity (perhaps it would be more appropriate to determine this from point group?)
        boolean linearity = false;
        if (pointGroup.equals("Cinfv")||pointGroup.equals("Dinfh")) linearity=true;//determine linearity from 3D-geometry; changed to correctly consider linear ketene radical case
        //we will use number of atoms from above (alternatively, we could use the chemGraph); this is needed to test whether the species is monoatomic
        double Hf298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500;
	double sigmaCorr = getSigmaCorr(pointGroup);
        Hf298 = energy*Hartree_to_kcal;
        S298 = R*Math.log(gdStateDegen)+R*(3./2.*Math.log(2.*Math.PI*molmass/(1000.*Na*Math.pow(h,2.)))+5./2.*Math.log(k*298.15)-Math.log(100000.)+5./2.);//electronic + translation; note use of 10^5 Pa for standard pressure; also note that molecular mass needs to be divided by 1000 for kg units
        Cp300 = 5./2.*R;
        Cp400 = 5./2.*R;
        Cp500 = 5./2.*R;
        Cp600 = 5./2.*R;
        Cp800 = 5./2.*R;
        Cp1000 = 5./2.*R;
        Cp1500 = 5./2.*R;
        if(natoms>1){//include statistical correction and rotational (without symmetry number, vibrational contributions if species is polyatomic
            if(linearity){//linear case
                //determine the rotational constant (note that one of the rotcons will be zero)
                double rotCons;
                if(rotCons_1 > 0.0001) rotCons = rotCons_1;
                else rotCons = rotCons_2;
                S298 += R*sigmaCorr+R*(Math.log(k*298.15/(h*rotCons))+1)+R*calcVibS(freqs, 298.15, h, k, c);
                Cp300 += R + R*calcVibCp(freqs, 300., h, k, c); 
                Cp400 += R + R*calcVibCp(freqs, 400., h, k, c);
                Cp500 += R + R*calcVibCp(freqs, 500., h, k, c);
                Cp600 += R + R*calcVibCp(freqs, 600., h, k, c);
                Cp800 += R + R*calcVibCp(freqs, 800., h, k, c);
                Cp1000 += R + R*calcVibCp(freqs, 1000., h, k, c);
                Cp1500 += R + R*calcVibCp(freqs, 1500., h, k, c);
            }
            else{//nonlinear case
                S298 += R*sigmaCorr+R*(3./2.*Math.log(k*298.15/h)-1./2.*Math.log(rotCons_1*rotCons_2*rotCons_3/Math.PI)+3./2.)+R*calcVibS(freqs, 298.15, h, k, c);
                Cp300 += 3./2.*R + R*calcVibCp(freqs, 300., h, k, c);
                Cp400 += 3./2.*R + R*calcVibCp(freqs, 400., h, k, c);
                Cp500 += 3./2.*R + R*calcVibCp(freqs, 500., h, k, c);
                Cp600 += 3./2.*R + R*calcVibCp(freqs, 600., h, k, c);
                Cp800 += 3./2.*R + R*calcVibCp(freqs, 800., h, k, c);
                Cp1000 += 3./2.*R + R*calcVibCp(freqs, 1000., h, k, c);
                Cp1500 += 3./2.*R + R*calcVibCp(freqs, 1500., h, k, c);
            }
        }
        ThermoData result = new ThermoData(Hf298,S298,Cp300,Cp400,Cp500,Cp600,Cp800,Cp1000,Cp1500,5,1,1,"PM3 or MM4 calculation");//this includes rough estimates of uncertainty
        return result;
    }

    //gets the statistical correction for S in dimensionless units (divided by R)
    public double getSigmaCorr(String pointGroup){
        double sigmaCorr=0;
        //determine statistical correction factor for 1. external rotational symmetry (affects rotational partition function) and 2. chirality (will add R*ln2 to entropy) based on point group
        //ref: http://cccbdb.nist.gov/thermo.asp
        //assumptions below for Sn, T, Th, O, I seem to be in line with expectations based on order reported at: http://en.wikipedia.org/w/index.php?title=List_of_character_tables_for_chemically_important_3D_point_groups&oldid=287261611 (assuming order = symmetry number * 2 (/2 if chiral))...this appears to be true for all point groups I "know" to be correct
        //minor concern: does SYMMETRY appropriately calculate all Sn groups considering 2007 discovery of previous errors in character tables (cf. Wikipedia article above)
        if (pointGroup.equals("C1")) sigmaCorr=+Math.log(2.);//rot. sym. = 1, chiral
        else if (pointGroup.equals("Cs")) sigmaCorr=0; //rot. sym. = 1
        else if (pointGroup.equals("Ci")) sigmaCorr=0; //rot. sym. = 1
        else if (pointGroup.equals("C2")) sigmaCorr=0;//rot. sym. = 2, chiral (corrections cancel)
        else if (pointGroup.equals("C3")) sigmaCorr=+Math.log(2.)-Math.log(3.);//rot. sym. = 3, chiral
        else if (pointGroup.equals("C4")) sigmaCorr=+Math.log(2.)-Math.log(4.);//rot. sym. = 4, chiral
        else if (pointGroup.equals("C5")) sigmaCorr=+Math.log(2.)-Math.log(5.);//rot. sym. = 5, chiral
        else if (pointGroup.equals("C6")) sigmaCorr=+Math.log(2.)-Math.log(6.);//rot. sym. = 6, chiral
        else if (pointGroup.equals("C7")) sigmaCorr=+Math.log(2.)-Math.log(7.);//rot. sym. = 7, chiral
        else if (pointGroup.equals("C8")) sigmaCorr=+Math.log(2.)-Math.log(8.);//rot. sym. = 8, chiral
        else if (pointGroup.equals("D2")) sigmaCorr=+Math.log(2.)-Math.log(4.);//rot. sym. = 4, chiral
        else if (pointGroup.equals("D3")) sigmaCorr=+Math.log(2.)-Math.log(6.);//rot. sym. = 6, chiral
        else if (pointGroup.equals("D4")) sigmaCorr=+Math.log(2.)-Math.log(8.);//rot. sym. = 8, chiral
        else if (pointGroup.equals("D5")) sigmaCorr=+Math.log(2.)-Math.log(10.);//rot. sym. = 10, chiral
        else if (pointGroup.equals("D6")) sigmaCorr=+Math.log(2.)-Math.log(12.);//rot. sym. = 12, chiral
        else if (pointGroup.equals("D7")) sigmaCorr=+Math.log(2.)-Math.log(14.);//rot. sym. = 14, chiral
        else if (pointGroup.equals("D8")) sigmaCorr=+Math.log(2.)-Math.log(16.);//rot. sym. = 16, chiral
        else if (pointGroup.equals("C2v")) sigmaCorr=-Math.log(2.);//rot. sym. = 2
        else if (pointGroup.equals("C3v")) sigmaCorr=-Math.log(3.);//rot. sym. = 3
        else if (pointGroup.equals("C4v")) sigmaCorr=-Math.log(4.);//rot. sym. = 4
        else if (pointGroup.equals("C5v")) sigmaCorr=-Math.log(5.);//rot. sym. = 5
        else if (pointGroup.equals("C6v")) sigmaCorr=-Math.log(6.);//rot. sym. = 6
        else if (pointGroup.equals("C7v")) sigmaCorr=-Math.log(7.);//rot. sym. = 7
        else if (pointGroup.equals("C8v")) sigmaCorr=-Math.log(8.);//rot. sym. = 8
        else if (pointGroup.equals("C2h")) sigmaCorr=-Math.log(2.);//rot. sym. = 2
        else if (pointGroup.equals("C3h")) sigmaCorr=-Math.log(3.);//rot. sym. = 3
        else if (pointGroup.equals("C4h")) sigmaCorr=-Math.log(4.);//rot. sym. = 4
        else if (pointGroup.equals("C5h")) sigmaCorr=-Math.log(5.);//rot. sym. = 5
        else if (pointGroup.equals("C6h")) sigmaCorr=-Math.log(6.);//rot. sym. = 6
        else if (pointGroup.equals("C7h")) sigmaCorr=-Math.log(7.);//rot. sym. = 7
        else if (pointGroup.equals("C8h")) sigmaCorr=-Math.log(8.);//rot. sym. = 8
        else if (pointGroup.equals("D2h")) sigmaCorr=-Math.log(4.);//rot. sym. = 4
        else if (pointGroup.equals("D3h")) sigmaCorr=-Math.log(6.);//rot. sym. = 6
        else if (pointGroup.equals("D4h")) sigmaCorr=-Math.log(8.);//rot. sym. = 8
        else if (pointGroup.equals("D5h")) sigmaCorr=-Math.log(10.);//rot. sym. = 10
        else if (pointGroup.equals("D6h")) sigmaCorr=-Math.log(12.);//rot. sym. = 12
        else if (pointGroup.equals("D7h")) sigmaCorr=-Math.log(14.);//rot. sym. = 14
        else if (pointGroup.equals("D8h")) sigmaCorr=-Math.log(16.);//rot. sym. = 16
        else if (pointGroup.equals("D2d")) sigmaCorr=-Math.log(4.);//rot. sym. = 4
        else if (pointGroup.equals("D3d")) sigmaCorr=-Math.log(6.);//rot. sym. = 6
        else if (pointGroup.equals("D4d")) sigmaCorr=-Math.log(8.);//rot. sym. = 8
        else if (pointGroup.equals("D5d")) sigmaCorr=-Math.log(10.);//rot. sym. = 10
        else if (pointGroup.equals("D6d")) sigmaCorr=-Math.log(12.);//rot. sym. = 12
        else if (pointGroup.equals("D7d")) sigmaCorr=-Math.log(14.);//rot. sym. = 14
        else if (pointGroup.equals("D8d")) sigmaCorr=-Math.log(16.);//rot. sym. = 16
        else if (pointGroup.equals("S4")) sigmaCorr=-Math.log(2.);//rot. sym. = 2 ;*** assumed achiral
        else if (pointGroup.equals("S6")) sigmaCorr=-Math.log(3.);//rot. sym. = 3 ;*** assumed achiral
        else if (pointGroup.equals("S8")) sigmaCorr=-Math.log(4.);//rot. sym. = 4 ;*** assumed achiral
        else if (pointGroup.equals("T")) sigmaCorr=+Math.log(2.)-Math.log(12.);//rot. sym. = 12, *** assumed chiral
        else if (pointGroup.equals("Th")) sigmaCorr=-Math.log(12.);//***assumed rot. sym. = 12
        else if (pointGroup.equals("Td")) sigmaCorr=-Math.log(12.);//rot. sym. = 12
        else if (pointGroup.equals("O")) sigmaCorr=+Math.log(2.)-Math.log(24.);//***assumed rot. sym. = 24, chiral
        else if (pointGroup.equals("Oh")) sigmaCorr=-Math.log(24.);//rot. sym. = 24
        else if (pointGroup.equals("Cinfv")) sigmaCorr=0;//rot. sym. = 1
        else if (pointGroup.equals("Dinfh")) sigmaCorr=-Math.log(2.);//rot. sym. = 2
        else if (pointGroup.equals("I")) sigmaCorr=+Math.log(2.)-Math.log(60.);//***assumed rot. sym. = 60, chiral
        else if (pointGroup.equals("Ih")) sigmaCorr=-Math.log(60.);//rot. sym. = 60
        else if (pointGroup.equals("Kh")) sigmaCorr=0;//arbitrarily set to zero...one could argue that it is infinite; apparently this is the point group of a single atom (cf. http://www.cobalt.chem.ucalgary.ca/ps/symmetry/tests/G_Kh); this should not have a rotational partition function, and we should not use the symmetry correction in this case
        else{//this point should not be reached, based on checks performed in determinePointGroupUsingSYMMETRYProgram
            System.out.println("Unrecognized point group: "+ pointGroup);
            System.exit(0);
        }
	
	return sigmaCorr;
}
    
    
    //determine the point group using the SYMMETRY program (http://www.cobalt.chem.ucalgary.ca/ps/symmetry/)
    //required input is a line with number of atoms followed by lines for each atom including atom number and x,y,z coordinates
    //finalTol determines how loose the point group criteria are; values are comparable to those specifed in the GaussView point group interface
    //public String determinePointGroupUsingSYMMETRYProgram(String geom, double finalTol){
    public String determinePointGroupUsingSYMMETRYProgram(String geom){
        int attemptNumber = 1;
        int maxAttemptNumber = 2;
        boolean pointGroupFound=false;
        //write the input file
        try {
            File inputFile=new File(qmfolder+"symminput.txt");//SYMMETRY program directory
            FileWriter fw = new FileWriter(inputFile);
            fw.write(geom);
            fw.close();
        } catch (IOException e) {
            String err = "Error writing input file for point group calculation";
            err += e.toString();
            System.out.println(err);
            System.exit(0);
        }
        String result = "";
        String command = "";
        while (attemptNumber<=maxAttemptNumber && !pointGroupFound){
            //call the program and read the result
            result = "";
            String [] lineArray;
            try{ 
                if (System.getProperty("os.name").toLowerCase().contains("windows")){//the Windows case where the precompiled executable seems to need to be called from a batch script
		    if(attemptNumber==1) command = "\""+System.getProperty("RMG.workingDirectory")+"/scripts/symmetryDefault2.bat\" "+qmfolder+ "symminput.txt";//12/1/09 gmagoon: switched to use slightly looser criteria of 0.02 rather than 0.01 to handle methylperoxyl radical result from MOPAC
		    else if (attemptNumber==2) command = "\""+System.getProperty("RMG.workingDirectory")+"/scripts/symmetryLoose.bat\" " +qmfolder+ "symminput.txt";//looser criteria (0.1 instead of 0.01) to properly identify C2v group in VBURLMBUVWIEMQ-UHFFFAOYAVmult5 (InChI=1/C3H4O2/c1-3(2,4)5/h1-2H2/mult5) MOPAC result; C2 and sigma were identified with default, but it should be C2 and sigma*2
		    else{
			System.out.println("Invalid attemptNumber: "+ attemptNumber);
			System.exit(0);
		    }
                }
		else{//in other (non-Windows) cases, where it is compiled from scratch, we should be able to run this directly
		    if(attemptNumber==1) command = System.getProperty("RMG.workingDirectory")+"/bin/SYMMETRY.EXE -final 0.02 " +qmfolder+ "symminput.txt";//12/1/09 gmagoon: switched to use slightly looser criteria of 0.02 rather than 0.01 to handle methylperoxyl radical result from MOPAC
		    else if (attemptNumber==2) command = System.getProperty("RMG.workingDirectory")+"/bin/SYMMETRY.EXE -final 0.1 " +qmfolder+ "symminput.txt";//looser criteria (0.1 instead of 0.01) to properly identify C2v group in VBURLMBUVWIEMQ-UHFFFAOYAVmult5 (InChI=1/C3H4O2/c1-3(2,4)5/h1-2H2/mult5) MOPAC result; C2 and sigma were identified with default, but it should be C2 and sigma*2
		    else{
			System.out.println("Invalid attemptNumber: "+ attemptNumber);
			System.exit(0);
		    }
		}
                Process symmProc = Runtime.getRuntime().exec(command);
                //check for errors and display the error if there is one
                InputStream is = symmProc.getInputStream();
                InputStreamReader isr = new InputStreamReader(is);
                BufferedReader br = new BufferedReader(isr);
                String line=null;
                while ( (line = br.readLine()) != null) {
                    if(line.startsWith("It seems to be the ")){//last line, ("It seems to be the [x] point group") indicates point group
                        lineArray = line.split(" ");//split the line around spaces
                        result = lineArray[5];//point group string should be the 6th word
                    }
                }
                int exitValue = symmProc.waitFor();
            }
            catch(Exception e){
                String err = "Error in running point group calculation process using SYMMETRY \n";
                err += e.toString();
                e.printStackTrace();
                System.exit(0);
            }
            //check for a recognized point group
            if (result.equals("C1")||result.equals("Cs")||result.equals("Ci")||result.equals("C2")||result.equals("C3")||result.equals("C4")||result.equals("C5")||result.equals("C6")||result.equals("C7")||result.equals("C8")||result.equals("D2")||result.equals("D3")||result.equals("D4")||result.equals("D5")||result.equals("D6")||result.equals("D7")||result.equals("D8")||result.equals("C2v")||result.equals("C3v")||result.equals("C4v")||result.equals("C5v")||result.equals("C6v")||result.equals("C7v")||result.equals("C8v")||result.equals("C2h")||result.equals("C3h")||result.equals("C4h")||result.equals("C5h")||result.equals("C6h")||result.equals("C7h")||result.equals("C8h")||result.equals("D2h")||result.equals("D3h")||result.equals("D4h")||result.equals("D5h")||result.equals("D6h")||result.equals("D7h")||result.equals("D8h")||result.equals("D2d")||result.equals("D3d")||result.equals("D4d")||result.equals("D5d")||result.equals("D6d")||result.equals("D7d")||result.equals("D8d")||result.equals("S4")||result.equals("S6")||result.equals("S8")||result.equals("T")||result.equals("Th")||result.equals("Td")||result.equals("O")||result.equals("Oh")||result.equals("Cinfv")||result.equals("Dinfh")||result.equals("I")||result.equals("Ih")||result.equals("Kh")) pointGroupFound=true;
            else{
                if(attemptNumber < maxAttemptNumber) System.out.println("Attempt number "+attemptNumber+" did not identify a recognized point group (" +result+"). Will retry with looser point group criteria.");
                else{
                    System.out.println("Final attempt number "+attemptNumber+" did not identify a recognized point group (" +result+"). Exiting.");
                    System.exit(0);
                }
                attemptNumber++;
            }
        } 
        System.out.println("Point group: "+ result);//print result, at least for debugging purposes
        
        return result;
    }
    
    //gmagoon 6/8/09
    //calculate the vibrational contribution (divided by R, dimensionless) at temperature, T, in Kelvin to entropy
    //p_freqs in cm^-1; c in cm/s; k in J/K; h in J-s
    //ref.: http://cccbdb.nist.gov/thermo.asp
    public double calcVibS(ArrayList p_freqs, double p_T, double h, double k, double c){
        double Scontrib = 0;
        double dr;
        for(int i=0; i < p_freqs.size(); i++){
            double freq = (Double)p_freqs.get(i);
            dr = h*c*freq/(k*p_T); //frequently used dimensionless ratio
            Scontrib = Scontrib - Math.log(1.-Math.exp(-dr))+dr*Math.exp(-dr)/(1.-Math.exp(-dr));
        }
            
        return Scontrib;
    }
    
    //gmagoon 6/8/09
    //calculate the vibrational contribution (divided by R, dimensionless) at temperature, T, in Kelvin to heat capacity, Cp
    //p_freqs in cm^-1; c in cm/s; k in J/K; h in J-s
    //ref.: http://cccbdb.nist.gov/thermo.asp
    public double calcVibCp(ArrayList p_freqs, double p_T, double h, double k, double c){
        double Cpcontrib = 0;
        double dr;
        for(int i=0; i < p_freqs.size(); i++){
            double freq = (Double)p_freqs.get(i);
            dr = h*c*freq/(k*p_T); //frequently used dimensionless ratio
            Cpcontrib = Cpcontrib + Math.pow(dr, 2.)*Math.exp(-dr)/Math.pow(1.-Math.exp(-dr),2.);
        }
            
        return Cpcontrib;
    }

    //gmagoon 6/23/10
    //calculate the vibrational contribution (divided by R, units of K) at temperature, T, in Kelvin to Hthermal (ZPE not included)
    //p_freqs in cm^-1; c in cm/s; k in J/K; h in J-s
    //we need to ignore zero frequencies, as MM4 does, to avoid dividing by zero; based on equation for Ev in http://www.gaussian.com/g_whitepap/thermo.htm; however, we ignore zero point contribution to be consistent with the input that CanTherm currently takes; note, however, that this could introduce some small inaccuracies, as the frequencies may be slightly different in CanTherm vs. MM4, particularly for a frequency < 7.7 cm^-1 (treated as zero in MM4)
    public double calcVibH(ArrayList p_freqs, double p_T, double h, double k, double c){
        double Hcontrib = 0;
        double dr;
        for(int i=0; i < p_freqs.size(); i++){
            double freq = (Double)p_freqs.get(i);
	    if(freq > 0.0){//ignore zero frequencies as MM4 does
		dr = h*c*freq/(k*p_T); //frequently used dimensionless ratio
		Hcontrib = Hcontrib + dr*p_T/(Math.exp(dr) - 1.);
	    }
        }

        return Hcontrib;
    }

    
    //determine the QM filename (element 0) and augmented InChI (element 1) for a ChemGraph
    //QM filename is InChIKey appended with mult3, mult4, mult5, or mult6 for multiplicities of 3 or higher
    //augmented InChI is InChI appended with /mult3, /mult4, /mult5, or /mult6 for multiplicities of 3 or higher
    public String [] getQMFileName(ChemGraph p_chemGraph){
        String [] result = new String[2];
        result[0] = p_chemGraph.getModifiedInChIKeyAnew();//need to generate InChI and key anew because ChemGraph may have changed (in particular, adding/removing hydrogens in HBI process)
        result[1] = p_chemGraph.getModifiedInChIAnew();

        return result;
    }
    
    //returns true if a Gaussian file for the given name and directory (.log suffix) exists and indicates successful completion (same criteria as used after calculation runs); terminates if the InChI doesn't match the InChI in the file or if there is no InChI in the file; returns false otherwise
    public boolean successfulGaussianResultExistsQ(String name, String directory, String InChIaug){
        //part of the code is taken from runGaussian code above
        //look in the output file to check for the successful termination of the Gaussian calculation
        //failed jobs will contain the a line beginning with " Error termination" near the end of the file
        File file = new File(directory+"/"+name+".log");
        if(file.exists()){//if the file exists, do further checks; otherwise, we will skip to final statement and return false
            int failureFlag=0;//flag (1 or 0) indicating whether the Gaussian job failed
            int InChIMatch=0;//flag (1 or 0) indicating whether the InChI in the file matches InChIaug; this can only be 1 if InChIFound is also 1;
            int InChIFound=0;//flag (1 or 0) indicating whether an InChI was found in the log file
            int InChIPartialMatch=0;//flag (1 or 0) indicating whether the InChI in the log file is a substring of the InChI in RMG's memory
            String logFileInChI="";
            try{
                FileReader in = new FileReader(file);
                BufferedReader reader = new BufferedReader(in);
                String line=reader.readLine();
                while(line!=null){
                    if (line.startsWith(" Error termination ")) failureFlag=1;
                    else if (line.startsWith(" ******")){//also look for imaginary frequencies
                        if (line.contains("imaginary frequencies")) failureFlag=1;
                    }
                    else if(line.startsWith(" InChI=")){
                        logFileInChI = line.trim();
                        //continue reading lines until a line of dashes is found (in this way, we can read InChIs that span multiple lines)
                        line=reader.readLine();
                        while (!line.startsWith(" --------")){
                            logFileInChI += line.trim();
                            line=reader.readLine();
                        }
                        InChIFound=1;
                        if(logFileInChI.equals(InChIaug)) InChIMatch=1;
                        else if(InChIaug.startsWith(logFileInChI)) InChIPartialMatch=1;
                    }
                    line=reader.readLine();
                }
            }
            catch(Exception e){
                String err = "Error in reading preexisting Gaussian log file \n";
                err += e.toString();
                e.printStackTrace();
                System.exit(0);
            }
            //if the failure flag is still 0, the process should have been successful
            if (failureFlag==0&&InChIMatch==1){
                System.out.println("Pre-existing successful quantum result for " + name + " ("+InChIaug+") has been found. This log file will be used.");
                return true;
            }
            else if (InChIFound==1 && InChIMatch == 0){//InChIs do not match (most likely due to limited name length mirrored in log file (79 characters), but possibly due to a collision)
                if(InChIPartialMatch == 1){//case where the InChI in memory begins with the InChI in the log file; we will continue and check the input file, printing a warning if there is no match
                    File inputFile = new File(directory+"/"+name+".gjf");
                    if(inputFile.exists()){//read the Gaussian inputFile
                        String inputFileInChI="";
                        try{
                            FileReader inI = new FileReader(inputFile);
                            BufferedReader readerI = new BufferedReader(inI);
                            String lineI=readerI.readLine();
                            while(lineI!=null){
                                if(lineI.startsWith(" InChI=")){
                                    inputFileInChI = lineI.trim();
                                }
                                lineI=readerI.readLine();
                            }
                        }
                        catch(Exception e){
                            String err = "Error in reading preexisting Gaussian gjf file \n";
                            err += e.toString();
                            e.printStackTrace();
                            System.exit(0);
                        }
                        if(inputFileInChI.equals(InChIaug)){
                            if(failureFlag==0){
                                System.out.println("Pre-existing successful quantum result for " + name + " ("+InChIaug+") has been found. This log file will be used. *Note that input file was read to confirm lack of InChIKey collision (InChI probably more than 79 characters)");
                                return true;
                            }
                            else{//otherwise, failureFlag==1
                                System.out.println("Pre-existing quantum result for " + name + " ("+InChIaug+") has been found, but the result was apparently unsuccessful. The file will be overwritten with a new calculation. *Note that input file was read to confirm lack of InChIKey collision (InChI probably more than 79 characters)");
                                return false;
                            }
                        }
                        else{
                            if(inputFileInChI.equals("")){//InChI was not found in input file
                                System.out.println("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI + " . InChI could not be found in the Gaussian input file. You should manually check that the log file contains the intended species.");
                                return true;
                            }
                            else{//InChI was found but doesn't match
                                System.out.println("Congratulations! You appear to have discovered the first recorded instance of an InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Gaussian input file Augmented InChI = "+inputFileInChI);
                                System.exit(0);
                            }
                        }
                    }
                    else{
                        System.out.println("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI + " . Gaussian input file could not be found to check full InChI. You should manually check that the log file contains the intended species.");
                        return true;
                    }
                }
                else{
                    System.out.println("Congratulations! You appear to have discovered the first recorded instance of an InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI);
                    System.exit(0);
                }
            }
            else if (InChIFound==0){
                System.out.println("An InChI was not found in file: " +name+".log");
                System.exit(0);
            }
            else if (failureFlag==1){//note these should cover all possible results for this block, and if the file.exists block is entered, it should return from within the block and should not reach the return statement below
                System.out.println("Pre-existing quantum result for " + name + " ("+InChIaug+") has been found, but the result was apparently unsuccessful. The file will be overwritten with a new calculation.");
                return false;
            }
        }
        //we could print a line here for cases where the file doesn't exist, but this would probably be too verbose
        return false;
    }
    //returns true if a successful result exists (either Gaussian or MOPAC)
//    public boolean [] successfulResultExistsQ(String name, String directory, String InChIaug){
//        boolean gaussianResult=successfulGaussianResultExistsQ(name, directory, InChIaug);
//        boolean mopacResult=successfulMOPACResultExistsQ(name, directory, InChIaug);
//        return (gaussianResult || mopacResult);// returns true if either a successful Gaussian or MOPAC result exists
//    }
    
    //returns true if a MOPAC output file for the given name and directory (.out suffix) exists and indicates successful completion (same criteria as used after calculation runs); terminates if the InChI doesn't match the InChI in the file or if there is no InChI in the file; returns false otherwise
    public boolean successfulMopacResultExistsQ(String name, String directory, String InChIaug){
        //part of the code is taken from analogous code for Gaussian
        //look in the output file to check for the successful termination of the calculation (assumed to be successful if "description of vibrations appears)
        File file = new File(directory+"/"+name+".out");
        if(file.exists()){//if the file exists, do further checks; otherwise, we will skip to final statement and return false
            int failureFlag=1;//flag (1 or 0) indicating whether the MOPAC job failed
            int failureOverrideFlag=0;//flag (1 or 0) to override success as measured by failureFlag
            int InChIMatch=0;//flag (1 or 0) indicating whether the InChI in the file matches InChIaug; this can only be 1 if InChIFound is also 1;
            int InChIFound=0;//flag (1 or 0) indicating whether an InChI was found in the log file
            int InChIPartialMatch=0;//flag (1 or 0) indicating whether the InChI in the log file is a substring of the InChI in RMG's memory
            String logFileInChI="";
            try{
                FileReader in = new FileReader(file);
                BufferedReader reader = new BufferedReader(in);
                String line=reader.readLine();
                while(line!=null){
                    String trimLine= line.trim();
                    if (trimLine.equals("DESCRIPTION OF VIBRATIONS")){
                       // if(!MopacFileContainsNegativeFreqsQ(name, directory)) failureFlag=0;//check for this line; if it is here, check for negative frequencies
                       failureFlag = 0;
                    }
                                        //negative frequencies notice example:
//         NOTE: SYSTEM IS NOT A GROUND STATE, THEREFORE ZERO POINT
//         ENERGY IS NOT MEANINGFULL. ZERO POINT ENERGY PRINTED
//         DOES NOT INCLUDE THE  2 IMAGINARY FREQUENCIES
                    else if (trimLine.endsWith("IMAGINARY FREQUENCIES")){
                      //  System.out.println("*****Imaginary freqencies found:");
                        failureOverrideFlag=1;
                    }
                    else if (trimLine.equals("EXCESS NUMBER OF OPTIMIZATION CYCLES")){//exceeding max cycles error
                        failureOverrideFlag=1;
                    }
                    else if (trimLine.equals("NOT ENOUGH TIME FOR ANOTHER CYCLE")){//timeout error
                        failureOverrideFlag=1;
                    }
                    else if(line.startsWith(" InChI=")){
                        logFileInChI = line.trim();//output files should take up to 240 characters of the name in the input file
                        InChIFound=1;
                        if(logFileInChI.equals(InChIaug)) InChIMatch=1;
                        else if(InChIaug.startsWith(logFileInChI)) InChIPartialMatch=1;
                    }
                    line=reader.readLine();
                }
            }
            catch(Exception e){
                String err = "Error in reading preexisting MOPAC output file \n";
                err += e.toString();
                e.printStackTrace();
                System.exit(0);
            }
            if(failureOverrideFlag==1) failureFlag=1; //job will be considered a failure if there are imaginary frequencies or if job terminates to to excess time/cycles
            //if the failure flag is still 0, the process should have been successful
            if (failureFlag==0&&InChIMatch==1){
                System.out.println("Pre-existing successful MOPAC quantum result for " + name + " ("+InChIaug+") has been found. This log file will be used.");
                return true;
            }
            else if (InChIFound==1 && InChIMatch == 0){//InChIs do not match (most likely due to limited name length mirrored in log file (240 characters), but possibly due to a collision)
               // if(InChIPartialMatch == 1){//case where the InChI in memory begins with the InChI in the log file; we will continue and check the input file, printing a warning if there is no match
                //look in the input file if the InChI doesn't match (apparently, certain characters can be deleted in MOPAC output file for long InChIs)
                File inputFile = new File(directory+"/"+name+".mop");
                if(inputFile.exists()){//read the MOPAC inputFile
                    String inputFileInChI="";
                    try{
                        FileReader inI = new FileReader(inputFile);
                        BufferedReader readerI = new BufferedReader(inI);
                        String lineI=readerI.readLine();
                        while(lineI!=null){
                            if(lineI.startsWith("InChI=")){
                                inputFileInChI = lineI.trim();
                            }
                            lineI=readerI.readLine();
                        }
                    }
                    catch(Exception e){
                        String err = "Error in reading preexisting MOPAC input file \n";
                        err += e.toString();
                        e.printStackTrace();
                        System.exit(0);
                    }
                    if(inputFileInChI.equals(InChIaug)){
                        if(failureFlag==0){
                            System.out.println("Pre-existing successful MOPAC quantum result for " + name + " ("+InChIaug+") has been found. This log file will be used. *Note that input file was read to confirm lack of InChIKey collision (InChI probably more than 240 characters or characters probably deleted from InChI in .out file)");
                            return true;
                        }
                        else{//otherwise, failureFlag==1
                            System.out.println("Pre-existing MOPAC quantum result for " + name + " ("+InChIaug+") has been found, but the result was apparently unsuccessful. The file will be overwritten with a new calculation or Gaussian result (if available) will be used. *Note that input file was read to confirm lack of InChIKey collision (InChI probably more than 240 characters or characters probably deleted from InChI in .out file)");
                            return false;
                        }
                    }
                    else{
                        if(inputFileInChI.equals("")){//InChI was not found in input file
                            System.out.println("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI + " . InChI could not be found in the MOPAC input file. You should manually check that the output file contains the intended species.");
                            return true;
                        }
                        else{//InChI was found but doesn't match
                            System.out.println("Congratulations! You appear to have discovered the first recorded instance of an InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " MOPAC input file Augmented InChI = " + inputFileInChI + " Log file Augmented InChI = "+logFileInChI);
                            System.exit(0);
                        }
                    }
                }
                else{
                    System.out.println("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI + " . MOPAC input file could not be found to check full InChI. You should manually check that the log file contains the intended species.");
                    return true;
                }
              //  }
//                else{
//                    System.out.println("Congratulations! You appear to have discovered the first recorded instance of an InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " MOPAC output file Augmented InChI = "+logFileInChI);
//                    System.exit(0);
//                }
            }
            else if (InChIFound==0){
                System.out.println("An InChI was not found in file: " +name+".out");
                System.exit(0);
            }
            else if (failureFlag==1){//note these should cover all possible results for this block, and if the file.exists block is entered, it should return from within the block and should not reach the return statement below
                System.out.println("Pre-existing MOPAC quantum result for " + name + " ("+InChIaug+") has been found, but the result was apparently unsuccessful. The file will be overwritten with a new calculation or Gaussian result (if available) will be used.");
                return false;
            }
        }
        //we could print a line here for cases where the file doesn't exist, but this would probably be too verbose
        return false;
    }

    //returns true if an MM4 output file for the given name and directory (.mm4out suffix) exists and indicates successful completion (same criteria as used after calculation runs); terminates if the InChI doesn't match the InChI in the file or if there is no InChI in the file; returns false otherwise
    public boolean successfulMM4ResultExistsQ(String name, String directory, String InChIaug){
        //part of the code is taken from analogous code for MOPAC (first ~half) and Gaussian (second ~half)
        //look in the output file to check for the successful termination of the calculation (assumed to be successful if "description of vibrations appears)
        int failureFlag=1;//flag (1 or 0) indicating whether the MM4 job failed
	int failureOverrideFlag=0;//flag (1 or 0) to override success as measured by failureFlag
	File file = new File(directory+"/"+name+".mm4out");
	File canFile = new File(directory+"/"+name+".canout");
        int InChIMatch=0;//flag (1 or 0) indicating whether the InChI in the file matches InChIaug; this can only be 1 if InChIFound is also 1;
        int InChIFound=0;//flag (1 or 0) indicating whether an InChI was found in the log file
        int InChIPartialMatch=0;//flag (1 or 0) indicating whether the InChI in the log file is a substring of the InChI in RMG's memory
	if(useCanTherm){//if we are using CanTherm, check whether a CanTherm output file exists...if it does, we will continue on, otherwise, we will rerun calculations (including MM4 calculation) from scratch to ensure our atom numbering is consistent; note: if the .canout file exists, we still want to check for the actual MM4 file even if we are using CanTherm and reading CanTherm output because 1. it ensures we have the correct species and don't have an InChI collision 2. it is needed for getting the geometry which is needed for the symmetry number corrections applied to CanTherm output (which doesn't include symmetry number considerations)
	    if(!canFile.exists()) return false;
	}
	if(file.exists()){//if the file exists, do further checks; otherwise, we will skip to final statement and return false
            String logFileInChI="";
            try{
                FileReader in = new FileReader(file);
                BufferedReader reader = new BufferedReader(in);
                String line=reader.readLine();
                while(line!=null){
                    String trimLine= line.trim();
                    if (trimLine.equals("STATISTICAL THERMODYNAMICS ANALYSIS")){
                       failureFlag = 0;
                    }
                    else if (trimLine.endsWith("imaginary frequencies,")){//read the number of imaginary frequencies and make sure it is zero
                        String[] split = trimLine.split("\\s+");
			if (Integer.parseInt(split[3])>0){
			    System.out.println("*****Imaginary freqencies found:");
			    failureOverrideFlag=1;
			}
                    }
		    else if (trimLine.contains("             0.0     (fir )")){
			if (useCanTherm){//zero frequencies are only acceptable when CanTherm is used
			    System.out.println("*****Warning: zero freqencies found (values lower than 7.7 cm^-1 are rounded to zero in MM4 output); CanTherm should hopefully correct this:");
			}
			else{
			    System.out.println("*****Zero freqencies found:");
			    failureOverrideFlag=1;
			}
		    }
                    else if(trimLine.startsWith("InChI=")){
                        logFileInChI = line.trim();//output files should take up to about 60 (?) characters of the name in the input file
                        InChIFound=1;
                        if(logFileInChI.equals(InChIaug)) InChIMatch=1;
                        else if(InChIaug.startsWith(logFileInChI)) InChIPartialMatch=1;
                    }
                    line=reader.readLine();
                }
            }
            catch(Exception e){
                String err = "Error in reading preexisting MM4 output file \n";
                err += e.toString();
                e.printStackTrace();
                System.exit(0);
            }
            if(failureOverrideFlag==1) failureFlag=1; //job will be considered a failure if there are imaginary frequencies or if job terminates to to excess time/cycles
            //if the failure flag is still 0, the process should have been successful
            if (failureFlag==0&&InChIMatch==1){
                System.out.println("Pre-existing successful MM4 result for " + name + " ("+InChIaug+") has been found. This log file will be used.");
                return true;
            }
            else if (InChIFound==1 && InChIMatch == 0){//InChIs do not match (most likely due to limited name length mirrored in log file (79 characters), but possibly due to a collision)
                if(InChIPartialMatch == 1){//case where the InChI in memory begins with the InChI in the log file; we will continue and check the input file, printing a warning if there is no match
                    File inputFile = new File(directory+"/"+name+".mm4");
                    if(inputFile.exists()){//read the MM4 inputFile
                        String inputFileInChI="";
                        try{
                            FileReader inI = new FileReader(inputFile);
                            BufferedReader readerI = new BufferedReader(inI);
                            String lineI=readerI.readLine();
			    //InChI should be repeated after in the first line of the input file
                            inputFileInChI = lineI.trim().substring(80);//extract the string starting with character 81
                        }
                        catch(Exception e){
                            String err = "Error in reading preexisting MM4 .mm4 file \n";
                            err += e.toString();
                            e.printStackTrace();
                            System.exit(0);
                        }
                        if(inputFileInChI.equals(InChIaug)){
                            if(failureFlag==0){
                                System.out.println("Pre-existing successful MM4 result for " + name + " ("+InChIaug+") has been found. This log file will be used. *Note that input file was read to confirm lack of InChIKey collision (InChI probably more than ~60 characters)");
                                return true;
                            }
                            else{//otherwise, failureFlag==1
                                System.out.println("Pre-existing MM4 result for " + name + " ("+InChIaug+") has been found, but the result was apparently unsuccessful. The file will be overwritten with a new calculation. *Note that input file was read to confirm lack of InChIKey collision (InChI probably more than ~60 characters)");
                                return false;
                            }
                        }
                        else{
                            if(inputFileInChI.equals("")){//InChI was not found in input file
                                System.out.println("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI + " . InChI could not be found in the MM4 input file. You should manually check that the log file contains the intended species.");
                                return true;
                            }
                            else{//InChI was found but doesn't match
                                System.out.println("Congratulations! You appear to have discovered the first recorded instance of an InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " MM4 input file Augmented InChI = "+inputFileInChI);
                                System.exit(0);
                            }
                        }
                    }
                    else{
                        System.out.println("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI + " . MM4 input file could not be found to check full InChI. You should manually check that the log file contains the intended species.");
                        return true;
                    }
                }
                else{
                    System.out.println("Congratulations! You appear to have discovered the first recorded instance of an InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI);
                    System.exit(0);
                }
            }
            else if (InChIFound==0){
                System.out.println("An InChI was not found in file: " +name+".mm4out");
                System.exit(0);
            }
            else if (failureFlag==1){//note these should cover all possible results for this block, and if the file.exists block is entered, it should return from within the block and should not reach the return statement below
                System.out.println("Pre-existing MM4 result for " + name + " ("+InChIaug+") has been found, but the result was apparently unsuccessful. The file will be overwritten with a new calculation.");
                return false;
            }
        }
        //we could print a line here for cases where the file doesn't exist, but this would probably be too verbose
        return false;
    }
    
//    //checks the MOPAC file for negative frequencies
//    public boolean MopacFileContainsNegativeFreqsQ(String name, String directory){
//        boolean negativeFreq=false;
//        
//        //code below copied from parseMopacPM3()
//        String command = "c:/Python25/python.exe c:/Python25/MopacPM3ParsingScript.py ";//this should eventually be modified for added generality
//        String logfilepath=directory+"/"+name+".out";
//        command=command.concat(logfilepath);
//        
//        //much of code below is copied from calculateThermoFromPM3Calc()
//        //parse the Mopac file using cclib
//        int natoms = 0; //number of atoms from Mopac file; in principle, this should agree with number of chemGraph atoms
//        ArrayList atomicNumber = new ArrayList(); //vector of atomic numbers (integers) (apparently Vector is thread-safe; cf. http://answers.yahoo.com/question/index?qid=20081214065127AArZDT3; ...should I be using this instead?)
//        ArrayList x_coor = new ArrayList(); //vectors of x-, y-, and z-coordinates (doubles) (Angstroms) (in order corresponding to above atomic numbers)
//        ArrayList y_coor = new ArrayList();
//        ArrayList z_coor = new ArrayList();
//        double energy = 0; //PM3 energy (Hf298) in Hartree (***note: in the case of MOPAC, the MOPAC file will contain in units of kcal/mol, but modified ccLib will return in Hartree)
//        double molmass = 0; //molecular mass in amu
//        ArrayList freqs = new ArrayList(); //list of frequencies in units of cm^-1
//        double rotCons_1 = 0;//rotational constants in (1/s)
//        double rotCons_2 = 0;
//        double rotCons_3 = 0; 
//        //int gdStateDegen = p_chemGraph.getRadicalNumber()+1;//calculate ground state degeneracy from the number of radicals; this should give the same result as spin multiplicity in Gaussian input file (and output file), but we do not explicitly check this (we could use "mult" which cclib reads in if we wanted to do so); also, note that this is not always correct, as there can apparently be additional spatial degeneracy for non-symmetric linear molecules like OH radical (cf. http://cccbdb.nist.gov/thermo.asp)
//        try{   
//            File runningdir=new File(directory);
//            Process cclibProc = Runtime.getRuntime().exec(command, null, runningdir);
//            //read the stdout of the process, which should contain the desired information in a particular format
//            InputStream is = cclibProc.getInputStream();
//            InputStreamReader isr = new InputStreamReader(is);
//            BufferedReader br = new BufferedReader(isr);
//            String line=null;
//            //example output:
////            C:\Python25>python.exe GaussianPM3ParsingScript.py TEOS.out
////            33
////            [ 6  6  8 14  8  6  6  8  6  6  8  6  6  1  1  1  1  1  1  1  1  1  1  1  1
////              1  1  1  1  1  1  1  1]
////            [[ 2.049061 -0.210375  3.133106]
////             [ 1.654646  0.321749  1.762752]
////             [ 0.359284 -0.110429  1.471465]
////             [-0.201871 -0.013365 -0.12819 ]
////             [ 0.086307  1.504918 -0.82893 ]
////             [-0.559186  2.619928 -0.284003]
////             [-0.180246  3.839463 -1.113029]
////             [ 0.523347 -1.188305 -1.112765]
////             [ 1.857584 -1.018167 -1.495088]
////             [ 2.375559 -2.344392 -2.033403]
////             [-1.870397 -0.297297 -0.075427]
////             [-2.313824 -1.571765  0.300245]
////             [-3.83427  -1.535927  0.372171]
////             [ 1.360346  0.128852  3.917699]
////             [ 2.053945 -1.307678  3.160474]
////             [ 3.055397  0.133647  3.403037]
////             [ 1.677262  1.430072  1.750899]
////             [ 2.372265 -0.029237  0.985204]
////             [-0.245956  2.754188  0.771433]
////             [-1.656897  2.472855 -0.287156]
////             [-0.664186  4.739148 -0.712606]
////             [-0.489413  3.734366 -2.161038]
////             [ 0.903055  4.016867 -1.112198]
////             [ 1.919521 -0.229395 -2.269681]
////             [ 2.474031 -0.680069 -0.629949]
////             [ 2.344478 -3.136247 -1.273862]
////             [ 1.786854 -2.695974 -2.890647]
////             [ 3.41648  -2.242409 -2.365094]
////             [-1.884889 -1.858617  1.28054 ]
////             [-1.976206 -2.322432 -0.440995]
////             [-4.284706 -1.26469  -0.591463]
////             [-4.225999 -2.520759  0.656131]
////             [-4.193468 -0.809557  1.112677]]
////            -14.1664924726
////            [    9.9615    18.102     27.0569    31.8459    39.0096    55.0091
////                66.4992    80.4552    86.4912   123.3551   141.6058   155.5448
////               159.4747   167.0013   178.5676   207.3738   237.3201   255.3487
////               264.5649   292.867    309.4248   344.6503   434.8231   470.2074
////               488.9717   749.1722   834.257    834.6594   837.7292   839.6352
////               887.9767   892.9538   899.5374   992.1851  1020.6164  1020.8671
////              1028.3897  1046.7945  1049.1768  1059.4704  1065.1505  1107.4001
////              1108.1567  1109.0466  1112.6677  1122.7785  1124.4315  1128.4163
////              1153.3438  1167.6705  1170.9627  1174.9613  1232.1826  1331.8459
////              1335.3932  1335.8677  1343.9556  1371.37    1372.8127  1375.5428
////              1396.0344  1402.4082  1402.7554  1403.2463  1403.396   1411.6946
////              1412.2456  1412.3519  1414.5982  1415.3613  1415.5698  1415.7993
////              1418.5409  2870.7446  2905.3132  2907.0361  2914.1662  2949.2646
////              2965.825   2967.7667  2971.5223  3086.3849  3086.3878  3086.6448
////              3086.687   3089.2274  3089.4105  3089.4743  3089.5841  3186.0753
////              3186.1375  3186.3511  3186.365 ]
////            [ 0.52729  0.49992  0.42466]
////note: above example has since been updated to print molecular mass; also frequency and atomic number format has been updated
//            String [] stringArray;
//            natoms = Integer.parseInt(br.readLine());//read line 1: number of atoms
//            stringArray = br.readLine().replace("[", "").replace("]","").trim().split(",\\s+");//read line 2: the atomic numbers (first removing braces)
//           // line = br.readLine().replace("[", "").replace("]","");//read line 2: the atomic numbers (first removing braces)
//           // StringTokenizer st = new StringTokenizer(line); //apprently the stringTokenizer class is deprecated, but I am having trouble getting the regular expressions to work properly
//            for(int i=0; i < natoms; i++){
//               // atomicNumber.add(i,Integer.parseInt(stringArray[i]));
//                atomicNumber.add(i,Integer.parseInt(stringArray[i]));
//            }
//            for(int i=0; i < natoms; i++){
//                stringArray = br.readLine().replace("[", "").replace("]","").trim().split("\\s+");//read line 3+: coordinates for atom i; used /s+ for split; using spaces with default limit of 0 was giving empty string
//                x_coor.add(i,Double.parseDouble(stringArray[0]));
//                y_coor.add(i,Double.parseDouble(stringArray[1]));
//                z_coor.add(i,Double.parseDouble(stringArray[2]));
//            }
//            energy = Double.parseDouble(br.readLine());//read next line: energy
//            molmass = Double.parseDouble(br.readLine());//read next line: molecular mass (in amu)
//            if (natoms>1){//read additional info for non-monoatomic species
//                stringArray = br.readLine().replace("[", "").replace("]","").trim().split(",\\s+");//read next line: frequencies
//                for(int i=0; i < stringArray.length; i++){
//                    freqs.add(i,Double.parseDouble(stringArray[i]));
//                }
//                stringArray = br.readLine().replace("[", "").replace("]","").trim().split("\\s+");//read next line rotational constants (converting from GHz to Hz in the process)
//                rotCons_1 = Double.parseDouble(stringArray[0])*1000000000;
//                rotCons_2 = Double.parseDouble(stringArray[1])*1000000000;
//                rotCons_3 = Double.parseDouble(stringArray[2])*1000000000;
//            }
//            while ( (line = br.readLine()) != null) {
//                //do nothing (there shouldn't be any more information, but this is included to get all the output)
//            }
//            int exitValue = cclibProc.waitFor();
//        }
//        catch (Exception e) {
//            String err = "Error in running ccLib Python process \n";
//            err += e.toString();
//            e.printStackTrace();
//            System.exit(0);
//        } 
//       
//        //start of code "new" to this function (aside from initialization of negativeFreq)
//        if(natoms > 0){
//            for (int i=0; i<freqs.size(); i++){
//                if((Double)freqs.get(i) < 0) negativeFreq = true;
//            }
//        }
//        return negativeFreq;
//    }
        //## operation initGAGroupLibrary()
    protected void initGAGroupLibrary() {
        //#[ operation initGAGroupLibrary()
        thermoLibrary = ThermoGAGroupLibrary.getINSTANCE();
        //#]
    }

    public QMData getQMDataWithCClib(String name, String directory, ChemGraph p_chemGraph, boolean getStericEnergy){
    	String command = "python "+System.getProperty("RMG.workingDirectory")+"/scripts/MM4ParsingScript.py ";
	String logfilepath=directory+"/"+name+".mm4out";
	command=command.concat(logfilepath);
	command=command.concat(" "+ System.getenv("RMG")+"/source");//this will pass $RMG/source to the script (in order to get the appropriate path for importing
        if (getStericEnergy) command=command.concat(" 1");//option to print stericenergy before molar mass (this will only be used in useHindRot cases, but it is always read in with this function
	///////////beginning of block taken from the bulk of getPM3MM4ThermoDataUsingCCLib////////////
	//parse the file using cclib
        int natoms = 0; //number of atoms from Mopac file; in principle, this should agree with number of chemGraph atoms
        ArrayList atomicNumber = new ArrayList(); //vector of atomic numbers (integers) (apparently Vector is thread-safe; cf. http://answers.yahoo.com/question/index?qid=20081214065127AArZDT3; ...should I be using this instead?)
        ArrayList x_coor = new ArrayList(); //vectors of x-, y-, and z-coordinates (doubles) (Angstroms) (in order corresponding to above atomic numbers)
        ArrayList y_coor = new ArrayList();
        ArrayList z_coor = new ArrayList();
        double energy = 0; // energy (Hf298) in Hartree
        double molmass = 0; //molecular mass in amu
	double stericEnergy = 0;//steric energy in Hartree
        ArrayList freqs = new ArrayList(); //list of frequencies in units of cm^-1
        double rotCons_1 = 0;//rotational constants in (1/s)
        double rotCons_2 = 0;
        double rotCons_3 = 0;
        int gdStateDegen = p_chemGraph.getRadicalNumber()+1;//calculate ground state degeneracy from the number of radicals; this should give the same result as spin multiplicity in Gaussian input file (and output file), but we do not explicitly check this (we could use "mult" which cclib reads in if we wanted to do so); also, note that this is not always correct, as there can apparently be additional spatial degeneracy for non-symmetric linear molecules like OH radical (cf. http://cccbdb.nist.gov/thermo.asp)
        try{
            Process cclibProc = Runtime.getRuntime().exec(command);
            //read the stdout of the process, which should contain the desired information in a particular format
            InputStream is = cclibProc.getInputStream();
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);

            String line=null;
            //example output:
//            C:\Python25>python.exe GaussianPM3ParsingScript.py TEOS.out
//            33
//            [ 6  6  8 14  8  6  6  8  6  6  8  6  6  1  1  1  1  1  1  1  1  1  1  1  1
//              1  1  1  1  1  1  1  1]
//            [[ 2.049061 -0.210375  3.133106]
//             [ 1.654646  0.321749  1.762752]
//             [ 0.359284 -0.110429  1.471465]
//             [-0.201871 -0.013365 -0.12819 ]
//             [ 0.086307  1.504918 -0.82893 ]
//             [-0.559186  2.619928 -0.284003]
//             [-0.180246  3.839463 -1.113029]
//             [ 0.523347 -1.188305 -1.112765]
//             [ 1.857584 -1.018167 -1.495088]
//             [ 2.375559 -2.344392 -2.033403]
//             [-1.870397 -0.297297 -0.075427]
//             [-2.313824 -1.571765  0.300245]
//             [-3.83427  -1.535927  0.372171]
//             [ 1.360346  0.128852  3.917699]
//             [ 2.053945 -1.307678  3.160474]
//             [ 3.055397  0.133647  3.403037]
//             [ 1.677262  1.430072  1.750899]
//             [ 2.372265 -0.029237  0.985204]
//             [-0.245956  2.754188  0.771433]
//             [-1.656897  2.472855 -0.287156]
//             [-0.664186  4.739148 -0.712606]
//             [-0.489413  3.734366 -2.161038]
//             [ 0.903055  4.016867 -1.112198]
//             [ 1.919521 -0.229395 -2.269681]
//             [ 2.474031 -0.680069 -0.629949]
//             [ 2.344478 -3.136247 -1.273862]
//             [ 1.786854 -2.695974 -2.890647]
//             [ 3.41648  -2.242409 -2.365094]
//             [-1.884889 -1.858617  1.28054 ]
//             [-1.976206 -2.322432 -0.440995]
//             [-4.284706 -1.26469  -0.591463]
//             [-4.225999 -2.520759  0.656131]
//             [-4.193468 -0.809557  1.112677]]
//            -14.1664924726
//            [    9.9615    18.102     27.0569    31.8459    39.0096    55.0091
//                66.4992    80.4552    86.4912   123.3551   141.6058   155.5448
//               159.4747   167.0013   178.5676   207.3738   237.3201   255.3487
//               264.5649   292.867    309.4248   344.6503   434.8231   470.2074
//               488.9717   749.1722   834.257    834.6594   837.7292   839.6352
//               887.9767   892.9538   899.5374   992.1851  1020.6164  1020.8671
//              1028.3897  1046.7945  1049.1768  1059.4704  1065.1505  1107.4001
//              1108.1567  1109.0466  1112.6677  1122.7785  1124.4315  1128.4163
//              1153.3438  1167.6705  1170.9627  1174.9613  1232.1826  1331.8459
//              1335.3932  1335.8677  1343.9556  1371.37    1372.8127  1375.5428
//              1396.0344  1402.4082  1402.7554  1403.2463  1403.396   1411.6946
//              1412.2456  1412.3519  1414.5982  1415.3613  1415.5698  1415.7993
//              1418.5409  2870.7446  2905.3132  2907.0361  2914.1662  2949.2646
//              2965.825   2967.7667  2971.5223  3086.3849  3086.3878  3086.6448
//              3086.687   3089.2274  3089.4105  3089.4743  3089.5841  3186.0753
//              3186.1375  3186.3511  3186.365 ]
//            [ 0.52729  0.49992  0.42466]
//note: above example has since been updated to print molecular mass and steric energy; also frequency and atomic number format has been updated
            String [] stringArray;
            natoms = Integer.parseInt(br.readLine());//read line 1: number of atoms
            stringArray = br.readLine().replace("[", "").replace("]","").trim().split(",\\s+");//read line 2: the atomic numbers (first removing braces)
           // line = br.readLine().replace("[", "").replace("]","");//read line 2: the atomic numbers (first removing braces)
           // StringTokenizer st = new StringTokenizer(line); //apprently the stringTokenizer class is deprecated, but I am having trouble getting the regular expressions to work properly
            for(int i=0; i < natoms; i++){
               // atomicNumber.add(i,Integer.parseInt(stringArray[i]));
                atomicNumber.add(i,Integer.parseInt(stringArray[i]));
            }
            for(int i=0; i < natoms; i++){
                stringArray = br.readLine().replace("[", "").replace("]","").trim().split("\\s+");//read line 3+: coordinates for atom i; used /s+ for split; using spaces with default limit of 0 was giving empty string
                x_coor.add(i,Double.parseDouble(stringArray[0]));
                y_coor.add(i,Double.parseDouble(stringArray[1]));
                z_coor.add(i,Double.parseDouble(stringArray[2]));
            }
            energy = Double.parseDouble(br.readLine());//read next line: energy
            stericEnergy = Double.parseDouble(br.readLine());//read next line: steric energy (in Hartree)
	    molmass = Double.parseDouble(br.readLine());//read next line: molecular mass (in amu)
            if (natoms>1){//read additional info for non-monoatomic species
                stringArray = br.readLine().replace("[", "").replace("]","").trim().split(",\\s+");//read next line: frequencies
                for(int i=0; i < stringArray.length; i++){
                    freqs.add(i,Double.parseDouble(stringArray[i]));
                }
                stringArray = br.readLine().replace("[", "").replace("]","").trim().split("\\s+");//read next line rotational constants (converting from GHz to Hz in the process)
                rotCons_1 = Double.parseDouble(stringArray[0])*1000000000;
                rotCons_2 = Double.parseDouble(stringArray[1])*1000000000;
                rotCons_3 = Double.parseDouble(stringArray[2])*1000000000;
            }
            while ( (line = br.readLine()) != null) {
                //do nothing (there shouldn't be any more information, but this is included to get all the output)
            }
            int exitValue = cclibProc.waitFor();
        }
        catch (Exception e) {
            String err = "Error in running ccLib Python process \n";
            err += e.toString();
            e.printStackTrace();
            System.exit(0);
        }
	//package up the result
	QMData qmdata = new QMData(natoms, atomicNumber, x_coor, y_coor, z_coor, energy, stericEnergy, molmass, freqs, rotCons_1, rotCons_2, rotCons_3);
	return qmdata;
    }

    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\QMTP.java
*********************************************************************/
