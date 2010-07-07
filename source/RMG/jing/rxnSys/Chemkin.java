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



package jing.rxnSys;


import java.io.*;

import jing.rxn.*;
import jing.chem.*;

import java.util.*;

import jing.param.*;

import org.w3c.dom.*;
import jing.mathTool.*;
import jing.rxn.Reaction;
import jing.rxn.PDepRateConstant.Mode;

import javax.xml.parsers.*;
import javax.xml.transform.*;
import org.xml.sax.SAXException;
import javax.xml.transform.dom.*;
import javax.xml.transform.stream.*;

//## package jing::rxnSys

//----------------------------------------------------------------------------
//jing\rxnSys\Chemkin.java
//----------------------------------------------------------------------------

//## class Chemkin
public class Chemkin implements DAESolver {

  protected double atol;		//## attribute atol

  protected String reactorType;		//## attribute reactorType

  protected double rtol;		//## attribute rtol
  
  public static boolean SMILESutility = false;

  //protected String thermoHeader = "";

  // Constructors

  //## operation Chemkin(double,double,String)
  public  Chemkin(double p_rtol, double p_atol, String p_reactorType) {
      //#[ operation Chemkin(double,double,String)
      if (p_rtol < 0 || p_atol < 0) throw new InvalidChemkinParameterException("Negative rtol or atol!");
      if (p_reactorType == null) throw new NullPointerException();

      String dir = System.getProperty("RMG.workingDirectory");
      
      //create the documentTypesDefinitions
      File docFile = new File("chemkin/documentTypeDefinitions");
      docFile.mkdir();
      copyFiles(dir+"/software/reactorModel/documentTypeDefinitions/reactorInput.dtd", "chemkin/documentTypeDefinitions/reactorInput.dtd");
      copyFiles(dir+"/software/reactorModel/documentTypeDefinitions/reactorOutput.dtd", "chemkin/documentTypeDefinitions/reactorOutput.dtd");
      
      rtol = p_rtol;
      atol = p_atol;
      reactorType = p_reactorType;
      
  }
  private static void copyFiles(String string, String string2)  {
	  File src = new File(string);
	  File dest = new File(string2);
	  FileInputStream fin;
	  try {
		  fin = new FileInputStream(src);
		  FileOutputStream fout = new FileOutputStream (dest);
		  int c;
		  while ((c = fin.read()) >= 0) 
			  fout.write(c);
		  fin.close();
		  fout.close();
	  } catch (FileNotFoundException e) {
		  // TODO Auto-generated catch block
		  e.printStackTrace();
	}catch (IOException e){
		e.printStackTrace();
	}

		
  }
public  Chemkin() {
  }

	public void addConversion(double [] temp1, int temp2){
	
	}
	
	public double[] getConversion(){
		return null;
	}
	
  //## operation checkChemkinMessage()
  public void checkChemkinMessage() {
      //#[ operation checkChemkinMessage()
      try {
      	String dir = System.getProperty("RMG.workingDirectory");
      	String filename = "chemkin/chem.message";
      	FileReader fr = new FileReader(filename);
      	BufferedReader br = new BufferedReader(fr);

      	String line = br.readLine().trim();
      	if (line.startsWith("NO ERRORS FOUND ON INPUT")) {
      		return;
      	}
      	else if (line.startsWith("WARNING...THERE IS AN ERROR IN THE LINKING FILE")) {
      		System.out.println("Error in chemkin linking to reactor!");
      		System.exit(0);
      	}
      	else {
      		System.out.println("Unknown message in chem.message!");
      		System.exit(0);
      	}
       }
       catch (Exception e) {
       	System.out.println("Can't read chem.message!");
       	System.out.println(e.getMessage());
       	System.exit(0);
       }

      //#]
  }

  //## operation clean()
  public void clean() {
      //#[ operation clean()
      //#]
  }

  //## operation generateSpeciesStatus(ReactionModel,ArrayList,ArrayList,ArrayList)
  private LinkedHashMap generateSpeciesStatus(ReactionModel p_reactionModel, ArrayList p_speciesChemkinName, ArrayList p_speciesConc, ArrayList p_speciesFlux) {
      //#[ operation generateSpeciesStatus(ReactionModel,ArrayList,ArrayList,ArrayList)
      int size = p_speciesChemkinName.size();
      if (size != p_speciesConc.size() || size != p_speciesFlux.size()) throw new InvalidSpeciesStatusException();
      LinkedHashMap speStatus = new LinkedHashMap();
      for (int i=0;i<size;i++){
      	String name = (String)p_speciesChemkinName.get(i);
      	int ID = parseIDFromChemkinName(name);
      	Species spe = SpeciesDictionary.getInstance().getSpeciesFromID(ID);
      	double conc = ((Double)p_speciesConc.get(i)).doubleValue();
      	double flux = ((Double)p_speciesFlux.get(i)).doubleValue();

      	System.out.println(String.valueOf(spe.getID()) + '\t' + spe.getName() + '\t' + String.valueOf(conc) + '\t' + String.valueOf(flux));

      	if (conc < 0) {
			double aTol = ReactionModelGenerator.getAtol();
			//if (Math.abs(conc) < aTol) conc = 0;
      		//else throw new NegativeConcentrationException("species " + spe.getName() + " has negative conc: " + String.valueOf(conc));
			if (conc < -100.0 * aTol) 
				throw new NegativeConcentrationException("Species " + spe.getName() + " has negative concentration: " + String.valueOf(conc));
      	}

      	SpeciesStatus ss = new SpeciesStatus(spe, 1, conc, flux);
      	speStatus.put(spe, ss);
      }
      return speStatus;

      //#]
  }

  //## operation isPDepReaction(Reaction)
  private static boolean isPDepReaction(Reaction p_reaction) {
      //#[ operation isPDepReaction(Reaction)
      if (p_reaction instanceof PDepReaction || p_reaction instanceof ThirdBodyReaction || p_reaction instanceof TROEReaction || p_reaction instanceof LindemannReaction) return true;
      else return false;

      //#]
  }

  //## operation parseIDFromChemkinName(String)
  private int parseIDFromChemkinName(String p_name) {
      //#[ operation parseIDFromChemkinName(String)
      char [] name = p_name.toCharArray();
      int pos = -1;
      for (int i=name.length-1;i>=0; i--) {
      	if (name[i]=='(') {
      		pos = i;
      		break;
      	}
      }
      if (pos < 0) throw new InvalidSpeciesException();

      String sID = p_name.substring(pos+1,name.length-1);
      return Integer.parseInt(sID);
      //#]
  }

  //## operation readReactorOutputFile(ReactionModel)
  public SystemSnapshot readReactorOutputFile(ReactionModel p_reactionModel) {
      //#[ operation readReactorOutputFile(ReactionModel)
      try {
      	// open output file and build the DOM tree
      	String dir = System.getProperty("RMG.workingDirectory");
      	String filename = "chemkin/reactorOutput.xml";
      	File inputFile = new File(filename);

      	DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
      	factory.setValidating(true); // validate the document with the DTD
      	factory.setIgnoringElementContentWhitespace(true); // ignore whitespace
      	DocumentBuilder builder = factory.newDocumentBuilder();
      	Document doc = builder.parse(inputFile);

      	// get root element and its children
      	Element root = doc.getDocumentElement();
      	NodeList rootchildren = root.getChildNodes();

      	// header is rootchildren.item(0)

      	// get return message and check for successful run
      	Element returnmessageElement = (Element) rootchildren.item(1);
      	Text returnmessageText = (Text) returnmessageElement.getFirstChild();
      	String returnmessage = returnmessageText.toString();
	returnmessage=returnmessage.trim();
      	if (!returnmessage.contains("SUCCESSFULLY COMPLETED RUN.")) {
      		System.out.println("External reactor model failed!");
      		System.out.println("Reactor model error message: " + returnmessage);
      		System.exit(0);
      	}

      	// get outputvalues element and its children
      	Element outputvaluesElement = (Element) rootchildren.item(2);
      	NodeList children = outputvaluesElement.getChildNodes();

      	// get time
      	Element timeElement = (Element) children.item(0);
      	Text timeText = (Text) timeElement.getFirstChild();
      	double time = Double.parseDouble(timeText.getData());
      	String timeUnits = timeElement.getAttribute("units");

      	// get systemstate element and its children
      	Element systemstateElement = (Element) children.item(1);
      	NodeList states = systemstateElement.getChildNodes();

       	// get temperature and its units
      	Element temperatureElement = (Element) states.item(0);
      	String tempUnits = temperatureElement.getAttribute("units");
      	Text temperatureText = (Text) temperatureElement.getFirstChild();
      	double temp = Double.parseDouble(temperatureText.getData());
          Temperature T = new Temperature(temp, tempUnits);

      	// get pressure and its units
      	Element pressureElement = (Element) states.item(1);
      	String presUnits = pressureElement.getAttribute("units");
      	Text pressureText = (Text) pressureElement.getFirstChild();
      	double pres = Double.parseDouble(pressureText.getData());
      	Pressure P = new Pressure(pres, presUnits);

      	// get species amounts (e.g. concentrations)
      	ArrayList speciesIDs = new ArrayList();
      	ArrayList amounts = new ArrayList();
      	ArrayList fluxes = new ArrayList();
      	String amountUnits = null;
      	String fluxUnits = null;

      	// loop thru all the species
      	// begin at i=2, since T and P take already the first two position of states
      	int nSpe = (states.getLength()-2)/2;
      	int index = 0;
      	LinkedHashMap inertGas = new LinkedHashMap();
      	for (int i = 2; i < nSpe+2; i++) {
      		// get amount element and the units
      		Element amountElement = (Element) states.item(i);
       		amountUnits = amountElement.getAttribute("units");

       		Element fluxElement = (Element)states.item(i+nSpe);
       		fluxUnits = fluxElement.getAttribute("units");

         		// get speciesid and store in an array list
      		String thisSpeciesID = amountElement.getAttribute("speciesid");

      		// get amount (e.g. concentraion) and store in an array list
      		Text amountText = (Text) amountElement.getFirstChild();
      		double thisAmount = Double.parseDouble(amountText.getData());
      		if (thisAmount < 0) {
				double aTol = ReactionModelGenerator.getAtol();
				//if (Math.abs(thisAmount) < aTol) thisAmount = 0;
      			//else throw new NegativeConcentrationException("Negative concentration in reactorOutput.xml: " + thisSpeciesID);
				if (thisAmount < -100.0 * aTol)
					throw new NegativeConcentrationException("Species " + thisSpeciesID + " has negative concentration: " + String.valueOf(thisAmount));
      		}

      		// get amount (e.g. concentraion) and store in an array list
      		Text fluxText = (Text)fluxElement.getFirstChild();
      		double thisFlux = Double.parseDouble(fluxText.getData());

              if (thisSpeciesID.compareToIgnoreCase("N2")==0 || thisSpeciesID.compareToIgnoreCase("Ne")==0 || thisSpeciesID.compareToIgnoreCase("Ar")==0) {
              	inertGas.put(thisSpeciesID, new Double(thisAmount));
              }
              else {
      			speciesIDs.add(index, thisSpeciesID);
      			amounts.add(index, new Double(thisAmount));
      			fluxes.add(index, new Double(thisFlux));
      			index++;
      		}
      	}

              // print results for debugging purposes
      /**
              System.out.println(returnmessage);
              System.out.println("Temp = " + temp + " " + tempUnits);
              System.out.println("Pres = " + pres + " " + presUnits);
              for (int i = 0; i < amounts.size(); i++) {
                System.out.println(speciesIDs.get(i) + " " + amounts.get(i) + " " +
                                   amountUnits);
              }
      */
      	ReactionTime rt = new ReactionTime(time, timeUnits);
      	LinkedHashMap speStatus = generateSpeciesStatus(p_reactionModel, speciesIDs, amounts, fluxes);
      	SystemSnapshot ss = new SystemSnapshot(rt, speStatus, T, P);
      	ss.inertGas = inertGas;
      	return ss;
      }
      catch (Exception e) {
      	System.out.println("Error reading reactor model output: " + e.getMessage());
      	System.exit(0);
      	return null;

      }



      //#]
  }

  //## operation runChemkin()
  public void runChemkin() {
      //#[ operation runChemkin()
      // run chemkin
      String dir = System.getProperty("RMG.workingDirectory");

      try {
         	// system call for chemkin
         	String[] command = {dir + "/software/reactorModel/chem.exe"};
         	File runningDir = new File("chemkin");
          Process chemkin = Runtime.getRuntime().exec(command, null, runningDir);
          InputStream ips = chemkin.getInputStream();
          InputStreamReader is = new InputStreamReader(ips);
          BufferedReader br = new BufferedReader(is);
          String line=null;
          while ( (line = br.readLine()) != null) {
          	//System.out.println(line);
          }
          int exitValue = chemkin.waitFor();
      }
      catch (Exception e) {
      	System.out.println("Error in running chemkin!");
      	System.out.println(e.getMessage());
      	System.exit(0);
      }


      //#]
  }

  //## operation runReactor()
  public void runReactor() {
      //#[ operation runReactor()
      // run reactor
      String dir = System.getProperty("RMG.workingDirectory");

      try {
         	// system call for reactor
         	String[] command = {dir + "/software/reactorModel/reactor.exe"};
         	File runningDir = new File("chemkin");
          Process reactor = Runtime.getRuntime().exec(command, null, runningDir);
          InputStream ips = reactor.getInputStream();
          InputStreamReader is = new InputStreamReader(ips);
          BufferedReader br = new BufferedReader(is);
          String line=null;
          while ( (line = br.readLine()) != null) {
          	//System.out.println(line);
          }
          int exitValue = reactor.waitFor();
      }
      catch (Exception e) {
      	System.out.println("Error in running reactor!");
      	System.out.println(e.getMessage());
      	System.exit(0);
      }

      //#]
  }

  //## operation solve(boolean,ReactionModel,boolean,SystemSnapshot,ReactionTime,ReactionTime,boolean)
  public SystemSnapshot solve(boolean p_initialization, ReactionModel p_reactionModel, boolean p_reactionChanged, SystemSnapshot p_beginStatus, final ReactionTime p_beginTime, ReactionTime p_endTime, boolean p_conditionChanged) {
      //#[ operation solve(boolean,ReactionModel,boolean,SystemSnapshot,ReactionTime,ReactionTime,boolean)
      writeChemkinInputFile(p_reactionModel, p_beginStatus);
      runChemkin();
      checkChemkinMessage();

      writeReactorInputFile(p_reactionModel,p_beginTime, p_endTime, p_beginStatus);
      runReactor();
      System.out.println("After ODE: from " + p_beginTime + " to " + p_endTime);
      SystemSnapshot result = readReactorOutputFile(p_reactionModel);
      return result;
      //#]
  }
    
    
  public static String writeChemkinHeader() {
	  StringBuilder header = new StringBuilder();
	  header.append("!  This chemkin file was generated by RMG - Reaction Mechanism Generator (http://rmg.mit.edu)");
	  header.append("\n!  The java code was compiled by ant at:");
	  header.append("\n!    "+jing.param.VersionInfo.getBuildDate());
	  header.append("\n!  The git repository was on the branch:");
	  header.append("\n!    "+jing.param.VersionInfo.getBranchName());
	  header.append("\n!  And at the commit with the hash:");
	  String versionHash=jing.param.VersionInfo.getVersionHash() ;
	  header.append("\n!    "+versionHash+"");
	  header.append("\n!");
	  header.append("\n!  For details visit:");
	  if (versionHash.startsWith("*") ) // error messages should start with a *
		  header.append("\n!   http://github.com/GreenGroup/RMG-Java/");
	  else {
		  header.append("\n!   http://github.com/GreenGroup/RMG-Java/tree/"+versionHash.substring(0,6) );
		  header.append("\n!To see changes since then visit:");
		  header.append("\n!   http://github.com/GreenGroup/RMG-Java/compare/"+versionHash.substring(0,6)+"...master");
	  }
	  header.append("\n\n");
	  return header.toString();
  }

  //## operation writeChemkinElement()
  public static  String writeChemkinElement() {
      //#[ operation writeChemkinElement()
      return "ELEMENTS H C O N Ne Ar He Si S END\n";
      //#]
  }

  //## operation writeChemkinInputFile(ReactionModel,SystemSnapshot)
  public static void writeChemkinInputFile(final ReactionModel p_reactionModel, SystemSnapshot p_beginStatus) {
      StringBuilder result=new StringBuilder();
      result.append(writeChemkinHeader());
	  result.append(writeChemkinElement());
	  double start = System.currentTimeMillis();
      result.append(writeChemkinSpecies(p_reactionModel, p_beginStatus));
      result.append(writeChemkinThermo(p_reactionModel));
      Global.chemkinThermo = Global.chemkinThermo + (System.currentTimeMillis() - start)/1000/60;
	  start = System.currentTimeMillis();
      result.append(writeChemkinPdepReactions(p_reactionModel, p_beginStatus));
	  Global.chemkinReaction = Global.chemkinReaction + (System.currentTimeMillis() - start)/1000/60;

      //String dir = System.getProperty("RMG.workingDirectory");
      //if (!dir.endsWith("/")) dir += "/";
      //dir += "software/reactorModel/";
	  
	  // This defines the directory that all the chemkin files will be written to.
	  // They will also be copied to the "chemkin" directory.
	  String directory = "chemkin/"+String.valueOf(p_reactionModel.getSpeciesNumber());
	  // Make directory if it don't exist
	  File directory_object = new File(directory);
	  directory_object.mkdirs();
	  
      String chemkinFile = directory+"/chem.inp";
      try {
      	FileWriter fw = new FileWriter(chemkinFile);
      	fw.write(result.toString());
      	fw.close();
      }
      catch (Exception e) {
      	System.out.println("Error in writing chemkin input file "+directory+"/chem.inp!");
      	System.out.println(e.getMessage());
      	System.exit(0);
      }
	  // put a copy of the latest in the root chemkin folder
	  copyFiles(chemkinFile, "chemkin/chem.inp");
      
	  // write tableOfRateCoeffs.txt if running with pressure-dependence
      if (PDepRateConstant.getMode() == Mode.CHEBYSHEV ||
    		  PDepRateConstant.getMode() == Mode.PDEPARRHENIUS ||
    		  PDepRateConstant.getMode() == Mode.RATE) {
    	  StringBuilder gridOfRateCoeffs = new StringBuilder();
	      gridOfRateCoeffs.append(writeGridOfRateCoeffs(p_reactionModel));
	      String newFile = directory+"/tableOfRateCoeffs.txt";
	      try {
	    	  FileWriter fw = new FileWriter(newFile);
	    	  fw.write(gridOfRateCoeffs.toString());
	    	  fw.close();
	      }
	      catch (Exception e) {
	    	  System.out.println("Error in writing "+newFile);
	    	  System.out.println(e.getMessage());
	    	  System.exit(0);
	      }
		  // put a copy of the latest in the root chemkin folder
		  copyFiles(newFile, "chemkin/tableOfRateCoeffs.txt");
      }
      
	  // Write the transport file
	  String transportFilePath = directory+"/tran.dat";
      writeTransportFile((CoreEdgeReactionModel)p_reactionModel, transportFilePath);
	  // put a copy of the latest in the root chemkin folder
	  copyFiles(transportFilePath, "chemkin/tran.dat");
	  
  }
  
  public static void writeChemkinInputFile(ReactionSystem rs) {
	  // call the above writeChemkinInputFile method, with the appropriate parameters
	  writeChemkinInputFile(rs.reactionModel, rs.initialStatus);
   }
  

//## operation writeChemkinPdepReactions(ReactionModel, SystemSnapshot)
 public static String writeChemkinPdepReactions(ReactionModel p_reactionModel, SystemSnapshot p_beginStatus) {
	 /* 
	  Writes all reactions, not just the P-dep ones.
	  Returns the result as a string.
	  
	  First p-dep (troe, thirdbody, lindemann, and from pdep networks)
	  Then non-p-dep
	  Finally duplicates
	  */
      StringBuilder result = new StringBuilder();

      String reactionHeader = "";
      String units4Ea = ArrheniusKinetics.getEaUnits();
      if (units4Ea.equals("cal/mol")) reactionHeader = "CAL/MOL\t";
      else if (units4Ea.equals("kcal/mol")) reactionHeader = "KCAL/MOL\t";
      else if (units4Ea.equals("J/mol")) reactionHeader = "JOULES/MOL\t";
      else if (units4Ea.equals("kJ/mol")) reactionHeader = "KJOULES/MOL\t";
      else if (units4Ea.equals("Kelvins")) reactionHeader = "KELVINS\t";
      
      String units4A = ArrheniusKinetics.getAUnits();
      if (units4A.equals("moles")) reactionHeader += "MOLES\n";
      else if (units4A.equals("molecules")) reactionHeader += "MOLECULES\n";
	  
      result.append("REACTIONS\t" + reactionHeader);
      
	  LinkedList pDepList = new LinkedList();
	  LinkedList nonPDepList = new LinkedList();
	  LinkedList duplicates = new LinkedList();
	  
      CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionModel;
      //first get troe, thirdbody, and Lindemann reactions (from seed mechanism and primary kinetics) and add them to the pDepList
      for (Iterator iter = cerm.getReactionSet().iterator(); iter.hasNext(); ) {
        	Reaction r = (Reaction)iter.next();
        	/*
        	 * 1Jul2009-MRH:
        	 * 	Added extra set of parenthesis.  Before, if the rxn was reverse but an instance of
        	 * 		TROEReaction, it would also be added to the pDepList, resulting in RMG reporting
        	 * 		both rxns (forward and reverse) in the chem.inp file, w/o a DUP tag.  Furthermore,
        	 * 		both rxns were given the same set of Arrhenius parameters.  Running this in
        	 * 		Chemkin-v4.1.1 resulted in an error.
        	 */
        	if (r.isForward() && (r instanceof ThirdBodyReaction || r instanceof TROEReaction || r instanceof LindemannReaction)) {        		
        		pDepList.add(r);
        	}
        }
      
	 // then get reactions from pressure-dependent networks and add them to pDepList
      for (Iterator iter = PDepNetwork.getNetworks().iterator(); iter.hasNext(); ) {
      	PDepNetwork pdn = (PDepNetwork)iter.next();
      	for (ListIterator pdniter = pdn.getNetReactions().listIterator(); pdniter.hasNext();) {
      		PDepReaction rxn = (PDepReaction) pdniter.next();
      		if (cerm.categorizeReaction(rxn) != 1) continue;
      		
      		//check if this reaction is not already in the list and also check if this reaction has a reverse reaction
      		// which is already present in the list.
      		if (rxn.getReverseReaction() == null)
      			rxn.generateReverseReaction();
      		
      		if (!rxn.reactantEqualsProduct() && !pDepList.contains(rxn) && !pDepList.contains(rxn.getReverseReaction()) )  {
      			pDepList.add(rxn);
      		}
      	}
      }

      for (Iterator iter = p_reactionModel.getReactionSet().iterator(); iter.hasNext(); ) {
      	Reaction r = (Reaction)iter.next();
		if (r.isForward() && !(r instanceof ThirdBodyReaction) && !(r instanceof TROEReaction) && !(r instanceof LindemannReaction)) {
			nonPDepList.add(r);
      	}
      }
      // First report pressure dependent reactions
      for (Iterator iter = pDepList.iterator(); iter.hasNext();){
    	  Reaction r = (Reaction)iter.next();
			// 6Jul2009-MRH:
			//	Pass both system temperature and pressure to function toChemkinString.
    	  	//		The only PDepKineticsModel that uses the passed pressure is RATE
          result.append(r.toChemkinString(p_beginStatus.getTemperature(),p_beginStatus.getPressure())+"\n");
      }
	 // Second, report non-pressure-dependent reactions
      for (Iterator iter = nonPDepList.iterator(); iter.hasNext();){
    	  Reaction r = (Reaction)iter.next();
          result.append(r.toChemkinString(p_beginStatus.getTemperature(),p_beginStatus.getPressure())+"\n");
      }
	 // Third, report duplicate reactions
      for (Iterator iter = duplicates.iterator(); iter.hasNext();){
    	  Reaction r = (Reaction)iter.next();
          result.append(r.toChemkinString(p_beginStatus.getTemperature(),p_beginStatus.getPressure())+"\n\tDUP\n");
      }

      result.append("END\n");
      return result.toString();
  }

  //## operation writeChemkinSpecies(ReactionModel,SystemSnapshot)
  public static String writeChemkinSpecies(ReactionModel p_reactionModel, SystemSnapshot p_beginStatus) {
      //#[ operation writeChemkinSpecies(ReactionModel,SystemSnapshot)

      StringBuilder result = new StringBuilder();
	  result.append("SPECIES\n");

      CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionModel;

      // write inert gas
      for (Iterator iter = p_beginStatus.getInertGas(); iter.hasNext();) {
      	String name = (String)iter.next();
      	result.append('\t' + name + '\n');
      }

      // write species
      for (Iterator iter = cerm.getSpecies(); iter.hasNext(); ) {
      	Species spe = (Species)iter.next();
      	if (spe.getChemkinName().startsWith("SPC"))
      		result.append("\t" + spe.getChemkinName() + "\t!" + spe.getName() + "\n");
      	else
      		result.append('\t' + spe.getChemkinName() + '\n');
      }

      result.append("END\n");
      return result.toString();
  }

  //## operation writeChemkinThermo(ReactionModel)
  public static String writeChemkinThermo(ReactionModel p_reactionModel) {
      //#[ operation writeChemkinThermo(ReactionModel)
      /*
	  String thermoHeader = "! neon added by pey (20/6/04) - used thermo for Ar\n";
		thermoHeader += "Ne                120186Ne  1               G  0300.00   5000.00  1000.00      1\n";
		thermoHeader += " 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2\n";
		thermoHeader += "-0.07453750E+04 0.04366001E+02 0.02500000E+02 0.00000000E+00 0.00000000E+00    3\n";
		thermoHeader += " 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.04366001E+02                   4\n";
		thermoHeader += "N2                121286N   2               G  0300.00   5000.00  1000.00      1\n";
		thermoHeader += " 0.02926640e+02 0.01487977e-01-0.05684761e-05 0.01009704e-08-0.06753351e-13    2\n";
		thermoHeader += "-0.09227977e+04 0.05980528e+02 0.03298677e+02 0.01408240e-01-0.03963222e-04    3\n";
		thermoHeader += " 0.05641515e-07-0.02444855e-10-0.01020900e+05 0.03950372e+02                   4\n";
		thermoHeader += "Ar                120186Ar  1               G  0300.00   5000.00  1000.00      1\n";
		thermoHeader += " 0.02500000e+02 0.00000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00    2\n";
		thermoHeader += "-0.07453750e+04 0.04366001e+02 0.02500000e+02 0.00000000e+00 0.00000000e+00    3\n";
		thermoHeader += " 0.00000000e+00 0.00000000e+00-0.07453750e+04 0.04366001e+02                   4\n";
       */
    //#]
      String thermoHeader = "! The first four sets of polynomial coefficients (Ar, N2, Ne, He) are from         \n";
      thermoHeader += "! THIRD MILLENIUM IDEAL GAS AND CONDENSED PHASE THERMOCHEMICAL DATABASE FOR     \n";
      thermoHeader += "! COMBUSTION WITH UPDATES FROM ACTIVE THERMOCHENICAL TABLES                     \n";
      thermoHeader += "! Authors: Alexander Burcat and Branko Ruscic                                   \n";
      thermoHeader += "!                                                                               \n";
      thermoHeader += "! The rest of the species are estimated by RMG (http://rmg.mit.edu/)            \n";
      // thermoHeader += "! Ar HF298=0.  REF=C.E. Moore 'Atomic Energy Levels' NSRDS-NBS 35 (1971) p.211  \n";
      // thermoHeader += "! NASA Glen (former Lewis) Research Center   (1988)                             \n";
      thermoHeader += "Ar                L 6/88Ar  1               G   200.000  6000.000 1000.        1\n";
      thermoHeader += " 0.25000000E+01 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2\n";
      thermoHeader += "-0.74537500E+03 0.43796749E+01 0.25000000E+01 0.00000000E+00 0.00000000E+00    3\n";
      thermoHeader += " 0.00000000E+00 0.00000000E+00-0.74537500E+03 0.43796749E+01                   4\n";
      // thermoHeader += "! N2  HF298= 0.0 KJ  REF=TSIV  Max Lst Sq Error Cp @ 6000 K 0.29%               \n";
      thermoHeader += "N2                G 8/02N   2               G   200.000  6000.000 1000.        1\n";
      thermoHeader += " 2.95257637E+00 1.39690040E-03-4.92631603E-07 7.86010195E-11-4.60755204E-15    2\n";
      thermoHeader += "-9.23948688E+02 5.87188762E+00 3.53100528E+00-1.23660988E-04-5.02999433E-07    3\n";
      thermoHeader += " 2.43530612E-09-1.40881235E-12-1.04697628E+03 2.96747038E+00                   4\n";
      // thermoHeader += "!Ne    HF298= 0.0 KJ REF=McBride, Heimel, Ehlers & Gordon                       \n";
      // thermoHeader += "!                'Thermodynamic Properties to 6000 K...' NASA SP-3001  (1963)   \n";
      thermoHeader += "Ne                L10/90Ne  1               G    200.0   6000.00  1000.0       1\n";
      thermoHeader += " 0.25000000E+01 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2\n";
      thermoHeader += "-0.74537500E+03 0.33553227E+01 0.25000000E+01 0.00000000E+00 0.00000000E+00    3\n";
      thermoHeader += " 0.00000000E+00 0.00000000E+00-0.74537498E+03 0.33553227E+01                   4\n";
      // thermoHeader += "7440-59-7                                                                       \n";                                                                               
      // thermoHeader += "He  HF298=0.0 KJ  REF=McBride, Heimel, Ehlers & Gordon "Thermodynamic Properties\n";
      // thermoHeader += "to 6000K ..." NASA SP-3001 1963.                                                \n";                                                
      thermoHeader += "He REF ELEMENT          He  1               G   200.000  6000.000 1000.        1\n";
      thermoHeader += " 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2\n";
      thermoHeader += "-7.45375000E+02 9.28723974E-01 2.50000000E+00 0.00000000E+00 0.00000000E+00    3\n";
      thermoHeader += " 0.00000000E+00 0.00000000E+00-7.45375000E+02 9.28723974E-01 0.00000000E+00    4\n\n";
      
      StringBuilder result = new StringBuilder();
	  result.append("THERMO ALL\n");
      result.append("   300.000  1000.000  5000.000\n");
      result.append(thermoHeader);

      CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionModel;
      for (Iterator iter = cerm.getSpecies(); iter.hasNext(); ) {
      	Species spe = (Species)iter.next();

      	if (spe.getNasaThermoSource() != null) {
      		result.append("!" + spe.getNasaThermoSource() + "\n");
      	}
      	/*
      	 * MRH 2MAR2010:
      	 * Added additional line to thermochemistry portion of chemkin file
      	 * 
      	 * Seyed asked the developers to add a line to the species
      	 * thermochemistry portion of the chem.inp file.  This line holds
      	 * a blank space for each species SMILES id.
      	 */
      	if (SMILESutility)
      		result.append("! [_ SMILES=\"" + spe.getInChI() + "\" _]\n");
      	result.append(spe.getNasaThermoData() + "\n");

      }
      result.append("END\n");
      result.append("\n");

      return result.toString();
  }
  

  public static String writeGridOfRateCoeffs(ReactionModel p_reactionModel) {
	  
      StringBuilder result = new StringBuilder();
	  
	  LinkedList pDepList = new LinkedList();
      CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionModel;
      
      for (Iterator iter = PDepNetwork.getNetworks().iterator(); iter.hasNext(); ) {
    	  PDepNetwork pdn = (PDepNetwork)iter.next();
    	  for (ListIterator pdniter = pdn.getNetReactions().listIterator(); pdniter.hasNext();) {
    		  PDepReaction rxn = (PDepReaction) pdniter.next();
    		  if (cerm.categorizeReaction(rxn) != 1) continue;	
    		  //check if this reaction is not already in the list and also check if this reaction has a reverse reaction
    		  // which is already present in the list.
    		  if (rxn.getReverseReaction() == null)
    			  rxn.generateReverseReaction();      		
    		  if (!rxn.reactantEqualsProduct() && !pDepList.contains(rxn) && !pDepList.contains(rxn.getReverseReaction())) {
    			  pDepList.add(rxn);
    		  }
    	  }
      }
      
      Temperature[] tempsUsedInFame = PDepRateConstant.getTemperatures();
      int numTemps = tempsUsedInFame.length;
      Pressure[] pressUsedInFame = PDepRateConstant.getPressures();
      int numPress = pressUsedInFame.length;
      
      for (int i=0; i<numTemps; i++) {
    	  for (int j=0; j<numPress; j++) {
    		  result.append("T="+tempsUsedInFame[i].getK()+"K,P="+pressUsedInFame[j].getBar()+"bar\t");
    	  }
    	  result.append("\n");
      }
      result.append("\n");
      
      for (Iterator iter = pDepList.iterator(); iter.hasNext();){
    	  PDepReaction r = (PDepReaction)iter.next();
    	  result.append(r.toString()+"\n");
    	  double[][] rates = new double[numTemps][numPress];
    	  rates = r.getPDepRate().getRateConstants();
    	  for (int i=0; i<numTemps; i++) {
    		  for (int j=0; j<numPress; j++) {
    			  result.append(rates[i][j] + "\t");
    		  }
    		  result.append("\n");
    	  }
          result.append("\n");
      }
      return result.toString();
  }

  //## operation writeReactorInputFile(ReactionModel,ReactionTime,ReactionTime,SystemSnapshot)
  public boolean writeReactorInputFile(ReactionModel p_reactionModel, ReactionTime p_beginTime, ReactionTime p_endTime, SystemSnapshot p_beginStatus) {
      //#[ operation writeReactorInputFile(ReactionModel,ReactionTime,ReactionTime,SystemSnapshot)
      // construct "input" string
      String input = "<?xml version=\"1.0\" standalone=\"no\"?>" + "\n";

      String dir = System.getProperty("RMG.workingDirectory");
      if (!dir.endsWith("/")) dir += "/";
      String dtd = dir + "software/reactorModel/documentTypeDefinitions/reactorInput.dtd";
      input += "<!DOCTYPE reactorinput SYSTEM \"" + dtd + "\">" + "\n";

      input += "<reactorinput>" + "\n";
      input += "<header>" + "\n";
      input += "<title>Reactor Input File</title>" + "\n";
      input += "<description>RMG-generated file used to call an external reactor model</description>" + "\n";
      input += "</header>" + "\n";
      input += "<inputvalues>" + "\n";
      input += "<integrationparameters>" + "\n";
      input += "<reactortype>" + reactorType + "</reactortype>" + "\n";
      input += "<starttime units=\"" + p_beginTime.getUnit() + "\">" + MathTool.formatDouble(p_beginTime.getTime(),15,6) +  "</starttime>" + "\n";
      input += "<endtime units=\"" + p_endTime.getUnit() + "\">" + MathTool.formatDouble(p_endTime.getTime(),15,6) +  "</endtime>" + "\n";
      //      input += "<starttime units=\"" + p_beginTime.unit + "\">" + MathTool.formatDouble(p_beginTime.time,15,6) +  "</starttime>" + "\n";
      //      input += "<endtime units=\"" + p_endTime.unit + "\">" + MathTool.formatDouble(p_endTime.time,15,6) +  "</endtime>" + "\n";
      input += "<rtol>" + rtol + "</rtol>" + "\n";
      input += "<atol>" + atol + "</atol>" + "\n";
      input += "</integrationparameters>" + "\n";
      input += "<chemistry>" + "\n";
      input += "</chemistry>" + "\n";
      input += "<systemstate>" + "\n";
      input += "<temperature units=\"K\">" + MathTool.formatDouble(p_beginStatus.getTemperature().getK(),15,6) + "</temperature>" + "\n";
      input += "<pressure units=\"Pa\">" + MathTool.formatDouble(p_beginStatus.getPressure().getPa(),15,6) + "</pressure>" + "\n";
      for (Iterator iter = p_beginStatus.getSpeciesStatus(); iter.hasNext();) {
      	SpeciesStatus spcStatus = (SpeciesStatus) iter.next();
      	Species thisSpecies = spcStatus.getSpecies();
      	CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionModel;
      	if (cerm.containsAsReactedSpecies(thisSpecies)) {
      		String spcChemkinName = thisSpecies.getChemkinName();
      		double concentration = spcStatus.getConcentration();
      		input += "<amount units=\"molPerCm3\" speciesid=\"" + spcChemkinName + "\">" + concentration + "</amount>" + "\n";
      	}
      }
      for (Iterator iter = p_beginStatus.getInertGas(); iter.hasNext(); ) {
      	String name = (String)iter.next();
      	double conc = p_beginStatus.getInertGas(name);
      	if (conc != 0.0)
      		input += "<amount units=\"molPerCm3\" speciesid=\"" + name + "\">" + conc + "</amount>" + "\n";
      }
      input += "</systemstate>" + "\n";
      input += "</inputvalues>" + "\n";
      input += "</reactorinput>" + "\n";

      // write "input" string to file
      try {
      	String file = "chemkin/reactorInput.xml";
      	FileWriter fw = new FileWriter(file);
      	fw.write(input);
      	fw.close();
      	return true;
      }
      catch (Exception e) {
      	System.out.println("Error in writing reactorInput.xml!");
      	System.out.println(e.getMessage());
      	return false;
      }


      //#]
  }

  public double getAtol() {
      return atol;
  }

  public void setAtol(double p_atol) {
      atol = p_atol;
  }

  public String getReactorType() {
      return reactorType;
  }

  public void setReactorType(String p_reactorType) {
      reactorType = p_reactorType;
  }

  public double getRtol() {
      return rtol;
  }

  public void setRtol(double p_rtol) {
      rtol = p_rtol;
  }
  
  
public SystemSnapshot solve(boolean p_initialization, ReactionModel p_reactionModel, boolean p_reactionChanged, SystemSnapshot p_beginStatus, ReactionTime p_beginTime, ReactionTime p_endTime, Temperature p_temperature, Pressure p_pressure, boolean p_conditionChanged, TerminationTester tt, int iternum) {
	
	//writeChemkinInputFile(p_reactionModel, p_beginStatus);
	
	runChemkin();
	checkChemkinMessage();
	
	writeReactorInputFile(p_reactionModel, p_beginTime, p_endTime, p_beginStatus);
	runReactor();
	System.out.println("After ODE: from " + p_beginTime + " to "+ p_endTime);
	SystemSnapshot result = readReactorOutputFile(p_reactionModel);
	return result;
}

public static void setSMILES(boolean yesno) {
	SMILESutility = yesno;
}

	public static void writeTransportFile(CoreEdgeReactionModel cerm, String filepath) {
		//Write core species to filepath (eg. "chemkin/tran.dat")
		String coreSpecies ="";
		// Write the three inert gas species' transport data
		//	Data comes from CHEMKIN-v4.1.1 manual
		coreSpecies +=	"Ar                 0   136.500     3.330     0.000     0.000     0.000 !CHEMKIN-v4.1.1\n" +
						"He                 0    10.200     2.576     0.000     0.000     0.000 !CHEMKIN-v4.1.1\n" +
						"N2                 1    97.530     3.621     0.000     1.760     4.000 !CHEMKIN-v4.1.1\n";
		
		Iterator iter = cerm.getSpecies();
		
		while (iter.hasNext()){
			Species spe = (Species)iter.next();
			TransportData lj4species = spe.getChemkinTransportData();
			String whitespace = "                ";
			
			// Write the 6 transport properties
			coreSpecies += spe.getChemkinName() + whitespace.substring(spe.getChemkinName().length()) +
				"   " + lj4species.toString() + " ! " + lj4species.getSource() + "\t" + 
				lj4species.getComment() + "\n";
		}
		
		try {
			File trandat = new File(filepath);
			FileWriter fw = new FileWriter(trandat);
			fw.write(coreSpecies);
			fw.close();
		}
		catch (IOException e) {
			System.out.println("Could not write "+filepath);
			System.exit(0);
		}
	}

}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\Chemkin.java
*********************************************************************/

