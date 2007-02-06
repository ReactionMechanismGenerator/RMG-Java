//!********************************************************************************
//!
//!    RMG: Reaction Mechanism Generator
//!
//!    Copyright: Jing Song, MIT, 2002, all rights reserved
//!
//!    Author's Contact: jingsong@mit.edu
//!
//!    Restrictions:
//!    (1) RMG is only for non-commercial distribution; commercial usage
//!        must require other written permission.
//!    (2) Redistributions of RMG must retain the above copyright
//!        notice, this list of conditions and the following disclaimer.
//!    (3) The end-user documentation included with the redistribution,
//!        if any, must include the following acknowledgment:
//!        "This product includes software RMG developed by Jing Song, MIT."
//!        Alternately, this acknowledgment may appear in the software itself,
//!        if and wherever such third-party acknowledgments normally appear.
//!
//!    RMG IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED
//!    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//!    OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//!    DISCLAIMED.  IN NO EVENT SHALL JING SONG BE LIABLE FOR
//!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
//!    OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
//!    OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//!    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
//!    THE USE OF RMG, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//!
//!******************************************************************************



package jing.rxnSys;


import java.io.*;

import jing.rxn.*;
import jing.chem.*;
import java.util.*;

import jing.param.*;

import org.w3c.dom.*;
import jing.mathTool.*;
import jing.rxn.Reaction;
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

  protected static String thermoHeader = "";

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
      thermoHeader += "! neon added by pey (20/6/04) - used thermo for Ar\n";
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
      //#]
  }
  private void copyFiles(String string, String string2)  {
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
      		if (Math.abs(conc) < 1.0E-25) conc = 0;
      		else throw new NegativeConcentrationException("species " + spe.getName() + " has negative conc: " + String.valueOf(conc));
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
      if (p_reaction instanceof PDepNetReaction || p_reaction instanceof ThirdBodyReaction || p_reaction instanceof TROEReaction) return true;
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
      			if (thisAmount < atol) thisAmount = 0;
      			else throw new NegativeConcentrationException("Negative concentration in reactorOutput.xml: " + thisSpeciesID);
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

  //## operation writeChemkinElement()
  public static  String writeChemkinElement() {
      //#[ operation writeChemkinElement()
      return "ELEMENTS H C O N Ne Ar END\n";
      //#]
  }

  //## operation writeChemkinInputFile(ReactionModel,SystemSnapshot)
  public static void writeChemkinInputFile(final ReactionModel p_reactionModel, SystemSnapshot p_beginStatus) {
      //#[ operation writeChemkinInputFile(ReactionModel,SystemSnapshot)

      StringBuilder result=new StringBuilder();
	  result.append(writeChemkinElement());
	  double start = System.currentTimeMillis();
      result.append(writeChemkinSpecies(p_reactionModel, p_beginStatus));
      result.append(writeChemkinThermo(p_reactionModel));
      Global.chemkinThermo = Global.chemkinThermo + (System.currentTimeMillis() - start)/1000/60;
	  start = System.currentTimeMillis();
	  result.append(writeChemkinReactions(p_reactionModel));
	  Global.chemkinReaction = Global.chemkinReaction + (System.currentTimeMillis() - start)/1000/60;

      String dir = System.getProperty("RMG.workingDirectory");
      if (!dir.endsWith("/")) dir += "/";
      dir += "software/reactorModel/";
      String file = "chemkin/chem.inp";

      try {
      	FileWriter fw = new FileWriter(file);
      	fw.write(result.toString());
      	fw.close();
      }
      catch (Exception e) {
      	System.out.println("Error in writing chemkin input file chem.inp!");
      	System.out.println(e.getMessage());
      	System.exit(0);
      }
      
      //#]
  }
  
//## operation writeChemkinReactions(ReactionModel)
  public static String writeChemkinReactions(ReactionModel p_reactionModel) {
      //#[ operation writeChemkinReactions(ReactionModel)
      StringBuilder result = new StringBuilder();
	  result.append("REACTIONS	KCAL/MOLE\n");
      CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionModel;

      
      LinkedHashSet all = cerm.getReactedReactionSet();
	  
      HashSet hs = new HashSet();
      int numfor = 0;
      int numrev = 0;
      int numdup = 0;
      int numnorev = 0;
      for (Iterator iter = all.iterator(); iter.hasNext(); ) {
      	Reaction rxn = (Reaction)iter.next();
      	if (rxn.isForward()) {
      		result.append(" " + rxn.toChemkinString(Global.temperature) + "\n");
      		
      	}
      	
      }
     
      result.append("END\n");

      return result.toString();

      //#]
  }

  
  
  //## operation writeChemkinReactions(ReactionModel)
 public static String writeChemkinPdepReactions(ReactionModel p_reactionModel) {
      //#[ operation writeChemkinReactions(ReactionModel)

      StringBuilder result = new StringBuilder();
	  result.append("REACTIONS	KCAL/MOLE\n");

      CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionModel;

      LinkedHashSet RISet = new LinkedHashSet();
      // print normal reactions
      LinkedList rxns = cerm.generatePDepReactionSet();
      LinkedList nonPDepReactionList = (LinkedList)rxns.get(0);
      LinkedList pDepReactionList = (LinkedList)rxns.get(1);
      LinkedList all = new LinkedList();
      //all.addAll(p_reactionModel.getReactionSet());
	  all.addAll(nonPDepReactionList);
      all.addAll(pDepReactionList);

      for (Iterator iter = all.iterator(); iter.hasNext(); ) {
      	Reaction rxn = (Reaction)iter.next();
      	if (rxn.isForward()) {
      		if (rxn.hasResonanceIsomer()) RISet.add(rxn);

      		else result.append(" " + rxn.toChemkinString(Global.temperature) + "\n");

      	}
      }

      // print possible duplicated reactions
      LinkedList RIList = new LinkedList(RISet);
      while (!RIList.isEmpty()) {
      	Reaction r1 = (Reaction)RIList.removeFirst();
      	boolean found = false;
      	for (Iterator iter = RIList.iterator(); iter.hasNext();) {
      		Reaction r2 = (Reaction)iter.next();
      		if (r1.equals(r2)) {
      			if (isPDepReaction(r1)) {

      				if (!found) result .append(" " + r1.toChemkinString(Global.temperature) + '\n');

      			}
      			else if  (isPDepReaction(r2)) {

      				if (!found) result.append(" " + r2.toChemkinString(Global.temperature) + '\n');

      			}
      			else {
      				if (!found) {

      					result.append(" " + r1.toChemkinString(Global.temperature) + '\n');
      	 				result.append("\t" + "DUP\n");

      	 			}

      				result.append(" " + r2.toChemkinString(Global.temperature) + '\n');
      				result.append("\t" + "DUP\n");

      			}
      			iter.remove();
      			found = true;
      		}

		}
      	if (!found) result.append(" " + r1.toChemkinString(Global.temperature) + '\n');

      	

      }

      result.append("END\n");

      return result.toString();

      //#]
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
      	result.append('\t' + spe.getChemkinName() + '\n');
      }


      result.append("END\n");


      return result.toString();

      //#]
  }

  //## operation writeChemkinThermo(ReactionModel)
  public static String writeChemkinThermo(ReactionModel p_reactionModel) {
      //#[ operation writeChemkinThermo(ReactionModel)
      StringBuilder result = new StringBuilder();
	  result.append("THERMO ALL\n");
      result.append("   300.000  1000.000  5000.000\n");
      result.append(thermoHeader);

      CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)p_reactionModel;
      for (Iterator iter = cerm.getSpecies(); iter.hasNext(); ) {
      	Species spe = (Species)iter.next();

      	result.append(spe.getNasaThermoData() + "\n");

      }
      result.append("END\n");

      return result.toString();

      //#]
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
      		input += "<amount units=\"molPerCm3\" speciesid=\"" + spcChemkinName + "\">" + String.format("%1.6e",concentration) + "</amount>" + "\n";
      	}
      }
      for (Iterator iter = p_beginStatus.getInertGas(); iter.hasNext(); ) {
      	String name = (String)iter.next();
      	double conc = p_beginStatus.getInertGas(name);
      	if (conc != 0.0)
      		input += "<amount units=\"molPerCm3\" speciesid=\"" + name + "\">" + String.format("%1.6e",conc) + "</amount>" + "\n";
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

}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxnSys\Chemkin.java
*********************************************************************/

