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



package jing.chem;


import java.io.*;
import java.util.*;
import jing.mathTool.*;
import jing.mathTool.Queue;
import jing.chemUtil.*;
import jing.chemParser.*;
import jing.chemUtil.Node;
import jing.param.Temperature;

//## package jing::chem

//----------------------------------------------------------------------------
// jing\chem\Species.java
//----------------------------------------------------------------------------

//## class Species
public class Species {

	protected boolean GATPFitExecuted = false;
    protected int ID;		//## attribute ID

    protected static int TOTAL_NUMBER = 0;		//## attribute TOTAL_NUMBER

    protected ChemGraph chemGraph;		//## attribute chemGraph

    /**
    The initial value is the parameter for N2, now it is treated as a default value if no further detailed information provided.
    unit: cm-1
    */
    protected double deltaEDown = 461;		//## attribute deltaEDown

    protected String name = null;		//## attribute name

    protected HashSet resonanceIsomers = new HashSet();		//## attribute resonanceIsomers

    protected boolean therfitExecuted = false;		//## attribute therfitExecuted

    protected LennardJones LJ;
	protected NASAThermoData nasaThermoData;
    protected ThreeFrequencyModel threeFrequencyModel;
	//protected WilhoitThermoData wilhoitThermoData;

    // Flag to tag certain species as library only... i.e. we won't try them against RMG templates.
    // They will only react as defined in the primary reaction library.   GJB
    protected boolean IsReactive = true; 
    
    // Constructors

    //## operation Species()
    private  Species() {
        initRelations();
        //#[ operation Species()
        //#]
    }
    //## operation Species(String,ChemGraph)
    private  Species(int id, String p_name, ChemGraph p_chemGraph) {
        initRelations();
		ID = id;
        //#[ operation Species(String,ChemGraph)
        name = p_name;
        chemGraph = p_chemGraph;
        generateResonanceIsomers();
        findStablestThermoData();
        calculateLJParameters();
        selectDeltaEDown();
		generateNASAThermoDatabyGATPFit();
		//generateThreeFrequencyModel();
        //generateNASAThermoData();
        //#]
    }

    //## operation addResonanceIsomer(ChemGraph)
    public boolean addResonanceIsomer(ChemGraph p_resonanceIsomer) {
        //#[ operation addResonanceIsomer(ChemGraph)
        if (resonanceIsomers == null) resonanceIsomers = new HashSet();

        p_resonanceIsomer.setSpecies(this);
        return resonanceIsomers.add(p_resonanceIsomer);



        //#]
    }

    //## operation calculateCp(Temperature)
    public double calculateCp(Temperature p_temperature) {
        //#[ operation calculateCp(Temperature)
        return getThermoData().calculateCp(p_temperature);
        //#]
    }

    //## operation calculateG(Temperature)
    public double calculateG(Temperature p_temperature) {
        //#[ operation calculateG(Temperature)
        //return getThermoData().calculateG(p_temperature);
		return nasaThermoData.calculateFreeEnergy(p_temperature);
        //#]
    }

    //## operation calculateGLowerBound(Temperature)
    //svp
      public double calculateGLowerBound(Temperature p_temperature) {
        //#[ operation calculateGLowerBound(Temperature)
        return getThermoData().calculateGLowerBound(p_temperature);
        //#]
      }

    //## operation calculateGUpperBound(Temperature)
    //svp
      public double calculateGUpperBound(Temperature p_temperature) {
        //#[ operation calculateGUpperBound(Temperature)
        return getThermoData().calculateGUpperBound(p_temperature);
        //#]
      }


    //## operation calculateH(Temperature)
    public double calculateH(Temperature p_temperature) {
        //#[ operation calculateH(Temperature)
        //return getThermoData().calculateH(p_temperature);
        return nasaThermoData.calculateEnthalpy(p_temperature);
		
		//#]
    }

    //## operation calculateLJParameters()
    public void calculateLJParameters() {
        //#[ operation calculateLJParameters()
        int cNum = getChemGraph().getCarbonNumber();

        selectLJParametersForSpecialMolecule();
        if (cNum == 1) LJ = new LennardJones(3.758, 148.6);
        else if (cNum == 2) LJ = new LennardJones(4.443, 110.7);
        else if (cNum == 3) LJ = new LennardJones(5.118, 237.1);
        else if (cNum == 4) LJ = new LennardJones(4.687, 531.4);
        else if (cNum == 5) LJ = new LennardJones(5.784, 341.1);
        else LJ = new LennardJones(5.949, 399.3);

        return;


        //#]
    }

    //## operation calculateS(Temperature)
    public double calculateS(Temperature p_temperature) {
        //#[ operation calculateS(Temperature)
        //return getThermoData().calculateS(p_temperature);
        return nasaThermoData.calculateEntropy(p_temperature);
		//#]
    }

	 //## operation callGATPFit(String)
    private boolean callGATPFit(String p_directory) {
        //#[ operation callGATPFit(String)
        if (p_directory == null) throw new NullPointerException();

        // write GATPFit input file
		String workingDirectory = System.getProperty("RMG.workingDirectory");
        // write species name
        String ls = System.getProperty("line.separator");
        String result = "SPEC " + getChemkinName() + ls;

        // write the element
        ChemGraph cg = getChemGraph();
        int Hn = cg.getHydrogenNumber();
        int Cn = cg.getCarbonNumber();
        int On = cg.getOxygenNumber();
        /*if (Cn>0) result += "ELEM C " + MathTool.formatInteger(Cn,3,"L") + ls;
        if (Hn>0) result += "ELEM H " + MathTool.formatInteger(Hn,3,"L") + ls;
        if (On>0) result += "ELEM O " + MathTool.formatInteger(On,3,"L") + ls;*/
		
		result += "ELEM C " + MathTool.formatInteger(Cn,3,"L") + ls;
        result += "ELEM H " + MathTool.formatInteger(Hn,3,"L") + ls;
        result += "ELEM O " + MathTool.formatInteger(On,3,"L") + ls;
		
        // write H and S at 298
        ThermoData td = getThermoData();
        result += "H298 " + MathTool.formatDouble(td.getH298(), 10, 2).trim() + ls;
        result += "S298 " + MathTool.formatDouble(td.getS298(), 10, 2).trim() + ls;
		
        result += "DLTH " + MathTool.formatDouble(td.getH298(), 10, 2).trim() + ls;

        // write MW, temperature, ouput format, etc
        result += "MWEI " + MathTool.formatDouble(getMolecularWeight(), 6, 1).trim() + ls;
        result += "TEMP 1000.0" + ls;
		result += "TMIN 300.0"+ls;
		result += "TMAX 5000.0" + ls;
        result += "CHEM" + ls;
        result += "TEM2 2000.0" + ls;
		if (chemGraph.isLinear())  result += "LINEAR" + ls;
		else result += "NONLINEAR" + ls;
        result += String.valueOf(cg.getAtomNumber()) + ls;
        result += String.valueOf(getInternalRotor()) + ls;
		result += "TECP 300 " + MathTool.formatDouble(td.Cp300,10,2).trim() + ls;
		result += "TECP 400 " + MathTool.formatDouble(td.Cp400,10,2).trim() + ls;
		result += "TECP 500 " + MathTool.formatDouble(td.Cp500,10,2).trim() + ls;
		result += "TECP 600 " + MathTool.formatDouble(td.Cp600,10,2).trim() + ls;
		result += "TECP 800 " + MathTool.formatDouble(td.Cp800,10,2).trim() + ls;
		result += "TECP 1000 " + MathTool.formatDouble(td.Cp1000,10,2).trim() +ls;
		result += "TECP 1500 " + MathTool.formatDouble(td.Cp1500,10,2).trim() + ls;
		result += "END" + ls;

        // finished writing text for input file, now save result to fort.1
        String GATPFit_input_name = null;
        File GATPFit_input = null;

        GATPFit_input_name = "GATPFit/INPUT.txt";
        try {
        	GATPFit_input = new File(GATPFit_input_name);
        	FileWriter fw = new FileWriter(GATPFit_input);
        	fw.write(result);
        	fw.close();
        	}
        catch (IOException e) {
        	String err = "GATPFit input file error: " + ls;
        	err += e.toString();
        	throw new GATPFitException(err);
        }

        // call GATPFit
        boolean error = false;
        try {
        	 // system call for therfit
        	String[] command = {workingDirectory +  "/software/GATPFit/GATPFit.exe"};
			File runningDir = new File("GATPFit");
        	Process GATPFit = Runtime.getRuntime().exec(command, null, runningDir);
        	int exitValue = GATPFit.waitFor();
        	GATPFitExecuted = true;

        }
        catch (Exception e) {
        	String err = "Error in running GATPFit!" + ls;
        	err += e.toString();
        	throw new GATPFitException(err);
        }

        // return error = true, if there was a problem
        return error;


        //#]
    }


	   //## operation callTherfit(String,String)
    private boolean callTherfit(String p_directory, String p_mode) {
        //#[ operation callTherfit(String,String)
		
        if (p_directory == null || p_mode == null) throw new NullPointerException();
		String workingDirectory = System.getProperty("RMG.workingDirectory");
        // write therfit input file
        String result = p_mode.toUpperCase() + '\n';
        result += getChemkinName() + '\n';
        ThermoData td = getThermoData();

        result += MathTool.formatDouble(td.getH298(), 8, 2) + '\n';
        result += MathTool.formatDouble(td.getS298(), 8, 2) + '\n';
        result += MathTool.formatDouble(td.getCp300(), 7, 3) + '\n';
        result += MathTool.formatDouble(td.getCp400(), 7, 3) + '\n';
        result += MathTool.formatDouble(td.getCp500(), 7, 3) + '\n';
        result += MathTool.formatDouble(td.getCp600(), 7, 3) + '\n';
        result += MathTool.formatDouble(td.getCp800(), 7, 3) + '\n';
        result += MathTool.formatDouble(td.getCp1000(), 7, 3) + '\n';
        result += MathTool.formatDouble(td.getCp1500(), 7, 3) + '\n';

        ChemGraph cg = getChemGraph();
        result += "C\n";
        result += MathTool.formatInteger(cg.getCarbonNumber(),3,"R") + '\n';
        result += "H\n";
        result += MathTool.formatInteger(cg.getHydrogenNumber(),3,"R") + '\n';
        int oNum = cg.getOxygenNumber();
        if (oNum == 0) {
        	result += "0\n0\n";
        }
        else {
        	result += "O\n";
        	result += MathTool.formatInteger(oNum,3,"R") + '\n';
        }
        result += "0\n0\n";

        result += "G\n";
        result += Integer.toString(getInternalRotor()) + '\n';

        // finished writing text for input file, now save result to fort.1
        String therfit_input_name = null;
        File therfit_input = null;

		try {
        	// prepare therfit input file, "fort.1" is the input file name
        	therfit_input_name = p_directory+"/fort.1";//p_directory+"/software/therfit/fort.1";
        	therfit_input = new File(therfit_input_name);
        	FileWriter fw = new FileWriter(therfit_input);
        	fw.write(result);
        	fw.close();
        	}
        catch (IOException e) {
        	System.out.println("Wrong input file for therfit!");
        	System.exit(0);
        }

        // call therfit
        boolean error = false;
		try {
       	 // system call for therfit
			String[] command = {workingDirectory + "/software/therfit/therfit.exe"};
			File runningDir = new File(p_directory );//+ "/therfit");// "/software/therfit");
			Process therfit = Runtime.getRuntime().exec(command, null, runningDir);
			therfitExecuted = true;
			InputStream is = therfit.getErrorStream();
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			String line=null;
			while ( (line = br.readLine()) != null) {
				//System.out.println(line);
				line = line.trim();
				if (!line.startsWith("*** THRFIT Job Complete")) {
					String speName = getName();
					System.out.println("therfit error for species: " + speName+"\n"+toString());
					File newfile = new File(therfit_input_name+"."+speName);
					therfit_input.renameTo(newfile);
					error = true;
				}
			}
			int exitValue = therfit.waitFor();
			br.close();
			isr.close();
		}
		catch (Exception e) {
			System.out.println("Error in run therfit!");
			System.exit(0);
		}

        // return error = true, if there was a problem
        return error;

        //#]
    }

    //## operation doDelocalization(ChemGraph,Stack)
    private ChemGraph doDelocalization(ChemGraph p_chemGraph, Stack p_path) {
        //#[ operation doDelocalization(ChemGraph,Stack)
        if (p_path.isEmpty() || p_path.size() != 3) throw new InvalidDelocalizationPathException();

        Graph graph = Graph.copy(p_chemGraph.getGraph());

        // n1-a1-n2-a2-n3.
        Node node1 = graph.getNodeAt((Integer)p_path.pop());
        Node node2 = graph.getNodeAt((Integer)p_path.pop());
        Node node3 = graph.getNodeAt((Integer)p_path.pop());
        Arc arc1 = graph.getArcBetween(node1,node2);
        Arc arc2 = graph.getArcBetween(node2,node3);

        // deal with node1
        Atom atom1 = (Atom)node1.getElement();
        Atom newAtom1 = (Atom)atom1.changeRadical(1,null);
        node1.setElement(newAtom1);

        // deal with arc1
        Bond bond1 = (Bond)arc1.getElement();
        Bond newBond1 = bond1.changeBond(-1);
        arc1.setElement(newBond1);

        // deal with node2, actually do nothing

        // deal with arc2
        Bond bond2 = (Bond)arc2.getElement();
        Bond newBond2 = bond2.changeBond(1);
        arc2.setElement(newBond2);

        // deal with node3
        Atom atom3 = (Atom)node3.getElement();
        Atom newAtom3 = (Atom)atom3.changeRadical(-1,null);
        node3.setElement(newAtom3);

        node1.updateFgElement();
        node1.updateFeElement();
        node2.updateFgElement();
        node2.updateFeElement();
        node3.updateFgElement();
        node3.updateFeElement();

        p_path = null;

        try {
        	ChemGraph newIsomer = ChemGraph.make(graph);

        	if (addResonanceIsomer(newIsomer)) return newIsomer;
        	else return null;
        }
        catch (ForbiddenStructureException e) {
        	return null;
        }
        //#]
    }

    //## operation findAllDelocalizationPaths(Node)
    private HashSet findAllDelocalizationPaths(Node p_radical) {
        //#[ operation findAllDelocalizationPaths(Node)
        HashSet allPaths = new HashSet();
        Atom atom = (Atom)p_radical.getElement();
        if (!atom.isRadical()) return allPaths;

        Iterator iter = p_radical.getNeighbor();
        while (iter.hasNext()) {
        	Arc arc1 = (Arc)iter.next();
        	Bond bond1 = (Bond)arc1.getElement();
        	if (bond1.isSingle() || bond1.isDouble()) {
        		Node node1 = arc1.getOtherNode(p_radical);
        		Iterator iter2 = node1.getNeighbor();
        		while (iter2.hasNext()) {
        			Arc arc2 = (Arc)iter2.next();
        			if (arc2 != arc1) {
        				Bond bond2 = (Bond)arc2.getElement();
        				if (bond2.isDouble() || bond2.isTriple()) {
        					Node node2 = arc2.getOtherNode(node1);
        					Stack path = new Stack();
        					path.push(p_radical.getID());
        					path.push(node1.getID());
        					path.push(node2.getID());
        					allPaths.add(path);
        				}
        			}
        		}

        	}
        }

        return allPaths;
        //#]
    }

    //## operation findStablestThermoData()
    public void findStablestThermoData() {
        //#[ operation findStablestThermoData()
        double H = chemGraph.getThermoData().getH298();
        ChemGraph stablest = chemGraph;
        if (resonanceIsomers != null) {
        	Iterator iter = resonanceIsomers.iterator();
        	while (iter.hasNext()) {
        		ChemGraph g = (ChemGraph)iter.next();
        		double newH = g.getThermoData().getH298();
        		if (newH < H) {
        			H = newH;
        			stablest = g;
        		}
        	}
        }

        chemGraph = stablest;
        //#]
    }

	  //## operation generateNASAThermoDatabyGATPFit()
    public void generateNASAThermoDatabyGATPFit() {
        //#[ operation generateNASAThermoDatabyGATPFit()
        // get working directory
        String dir = System.getProperty("RMG.workingDirectory");

        try {
        	// prepare GATPFit input file and execute system call
        	boolean error = callGATPFit(dir);
        }
        catch (GATPFitException e) {
        	throw new NASAFittingException("Error in running GATPFit: " + e.toString());
        }

        // parse output from GATPFit, "output.txt" is the output file name
        String therfit_nasa_output = "GATPFit/OUTPUT.txt";

        try {
        	FileReader in = new FileReader(therfit_nasa_output);
        	BufferedReader data = new BufferedReader(in);

        	String line = data.readLine();
        	line = data.readLine();
        	String nasaString = "";

        	while (line != null) {
        		nasaString += line + System.getProperty("line.separator");
        		line = data.readLine();
        	}
        	nasaThermoData = new NASAThermoData(nasaString);
        }
        catch (Exception e) {
        	throw new NASAFittingException("Error in reading in GATPFit output file: " + System.getProperty("line.separator") + e.toString());
        }
        //#]
    }

	   //## operation generateNASAThermoData()
    public void generateNASAThermoData() {
        ///#[ operation generateNASAThermoData()
        // get working directory
        //String dir = System.getProperty("RMG.workingDirectory");
        String mode = "WILHOI"; // use Wilhoit polynomial to extrapolate Cp
		String dir = "therfit";
        // prepare therfit input file and execute system call
        boolean error = callTherfit(dir,mode);

        // parse output from therfit, "fort.25" is the output file name
        if (error) {
        	System.out.println("Error in generating NASA thermodata!");
        	System.exit(0);
        }

        try {
        	String therfit_nasa_output = dir+"/fort.54";

        	FileReader in = new FileReader(therfit_nasa_output);
        	BufferedReader data = new BufferedReader(in);

        	String line = data.readLine();
        	String nasaString = "";

        	while (line != null) {
        		nasaString += line + "\n";
        		line = data.readLine();
        	}
        	nasaThermoData = new NASAThermoData(nasaString);
        }
        catch (Exception e) {
        	System.out.println("Wrong output from therfit!");
        	System.exit(0);
        }


        //#]
    }


    /**
    Requires:
    Effects: generate all the possible resonance isomers for the primary chemGraph
    Modifies: this.resonanceIsomers
    */
    //## operation generateResonanceIsomers()
    public void generateResonanceIsomers() {
        //#[ operation generateResonanceIsomers()
        if (chemGraph == null) throw new NullPointerException();

        // check if the representation of chemGraph is correct
        if (!chemGraph.repOk()) {
        	resonanceIsomers = null;
        	throw new InvalidChemGraphException();
        }

        // generate RI for radical
        generateResonanceIsomersFromRadicalCenter();

        // generaate RI for O2, removed, don't allow .o-o.
        //generateResonanceIsomersForOxygen();

        if (resonanceIsomers.size() == 1) {
        	ChemGraph cg = (ChemGraph)(resonanceIsomers.iterator().next());
        	if (cg == chemGraph) resonanceIsomers.clear();
        	else addResonanceIsomer(chemGraph);
        }

        /*Graph g = Graph.copy(chemGraph.getGraph());
        // generate node-electron stucture
        int nodeNumber = g.getNodeNumber();
        int [] electronOnNode = new int[nodeNumber+1];

        Iterator nodeIter = g.getNodeList();
        while (nodeIter.hasNext()) {
        	Node node = (Node)nodeIter.next();
        	Atom atom = (Atom)node.getElement();
        	if (atom.isRadical()) {
        		int ID = node.getID().intValue();
        		electronOnNode[ID] = atom.getRadicalNumber();
        		Iterator arcIter = node.getNeighbor();
        		arcLoop:while (arcIter.hasNext()) {
        			Arc arc = (Arc)arcIter.next();
        			Bond bond = (Bond)arc.getElement();
        			if (bond.isBenzene() || bond.isTriple()) {
        				electronOnNode[ID] = 1;
        				break arcLoop;
        			}
        			else if (bond.isDouble()) {
        				electronOnNode[node.getID().intValue()-1] += 1;
        			}
        		}
        	}
        }
        */

        // combine them accordingly
        //resonanceIsomer = combineElectronBetweenNodes(g);

        return;




        //#]
    }

    //## operation generateResonanceIsomersForOxygen()
    public void generateResonanceIsomersForOxygen() {
        //#[ operation generateResonanceIsomersForOxygen()
        // form a O2 graph
        Graph g1 = new Graph();
        Node n1 = g1.addNodeAt(1,Atom.make(ChemElement.make("O"), FreeElectron.make("0")));
        Node n2 = g1.addNodeAt(2,Atom.make(ChemElement.make("O"), FreeElectron.make("0")));
        g1.addArcBetween(n1,Bond.make("D"),n2);

        // form a .o-o. graph
        Graph g2 = new Graph();
        Node n3 = g2.addNodeAt(1,Atom.make(ChemElement.make("O"), FreeElectron.make("1")));
        Node n4 = g2.addNodeAt(2,Atom.make(ChemElement.make("O"), FreeElectron.make("1")));
        g2.addArcBetween(n3,Bond.make("S"),n4);

        try {
        	if (chemGraph.getGraph().isEquivalent(g1)) {
        		ChemGraph cg = ChemGraph.make(g2);
        		addResonanceIsomer(cg);
        		}
        	else if (chemGraph.getGraph().isEquivalent(g2)) {
        		ChemGraph cg = ChemGraph.make(g1);
        		addResonanceIsomer(cg);
        	}
        	return;
        }
        catch (ForbiddenStructureException e) {
        	return;
        }



        //#]
    }

    //## operation generateResonanceIsomersFromRadicalCenter()
    private void generateResonanceIsomersFromRadicalCenter() {
        //#[ operation generateResonanceIsomersFromRadicalCenter()
        // only radical is considered here
        if (chemGraph.getRadicalNumber() <= 0) return;

        addResonanceIsomer(chemGraph);

        Queue undoChemGraph = new Queue(chemGraph.getAtomNumber());
        undoChemGraph.enqueue(chemGraph);

        HashSet processedChemGraph = new HashSet();
        while (!undoChemGraph.isEmpty()) {
        	ChemGraph cg = (ChemGraph)undoChemGraph.dequeue();
        	HashSet radicalNode = cg.getRadicalNode();
        	Iterator radicalIter = radicalNode.iterator();
        	while (radicalIter.hasNext()) {
        		Node radical = (Node)radicalIter.next();
        		int radicalNumber = ((Atom)radical.getElement()).getRadicalNumber();
        		if (radicalNumber > 0) {
        			HashSet allPath = findAllDelocalizationPaths(radical);
        			Iterator pathIter = allPath.iterator();
        			while (pathIter.hasNext()) {
        				Stack path = (Stack)pathIter.next();
        				ChemGraph newCG = doDelocalization(cg, path);
        				if (newCG!=null && !processedChemGraph.contains(newCG)) {
        					undoChemGraph.enqueue(newCG);
        				}
        			}
        		}
        	}
        	processedChemGraph.add(cg);
        }
        //#]
    }

    //## operation generateThreeFrequencyModel()
	public void generateThreeFrequencyModel() {
        //#[ operation generateThreeFrequencyModel()
        if (isTriatomicOrSmaller()) return;

        //String directory = System.getProperty("RMG.workingDirectory");
        String mode = "HOE"; // use Harmonic Oscillator Eqn to calc psuedo-freqs
		String dir = "therfit";
        // prepare therfit input file and execute system call
        boolean error = callTherfit(dir,mode);
        if (error) return;

        double [] degeneracy = new double [3];
        double [] frequency = new double [3];

        try {
        	String therfit_output = dir+"/fort.25";

        	FileReader in = new FileReader(therfit_output);
        	BufferedReader data = new BufferedReader(in);

        	String line = ChemParser.readMeaningfulLine(data);
        	if (!line.startsWith(getChemkinName())) {
        		System.out.println("Wrong output from therfit!");
        		System.exit(0);
        	}

        	int index = 0;
        	line = ChemParser.readMeaningfulLine(data);
        	while (line != null) {
        		line = line.trim();
        		if (line.startsWith("VIBRATION")) {
        			if (index > 2) break;
        			StringTokenizer token = new StringTokenizer(line);
        			String temp = token.nextToken();
        			temp = token.nextToken();
        			temp = token.nextToken();
        			temp = token.nextToken();
        			String deg = token.nextToken();
        			temp = token.nextToken();
        			temp = token.nextToken();
        	        String freq = token.nextToken();
        	        degeneracy[index] = Double.parseDouble(deg);
        	        frequency[index] = Double.parseDouble(freq);
        	        index++;
        	    }
        		line = ChemParser.readMeaningfulLine(data);
        	}

        	in.close();
        }
        catch (Exception e) {
        	System.out.println("Wrong output fort.25 from therfit!");
        	System.exit(0);
        }
        // creat threeFrequencyModel for this species
        threeFrequencyModel = new ThreeFrequencyModel(frequency, degeneracy);

        //#]
    }

	  public String getChemkinName() {
	        //#[ operation getChemkinName()
	        String chemkinName = getName() + "(" + getID() + ")";
	        if (chemkinName.length() > 16) chemkinName = "SPC(" + getID() + ")";
	        return chemkinName;

	        //#]
	    }

    //## operation getInternalRotor()
    public int getInternalRotor() {
        //#[ operation getInternalRotor()
        return getChemGraph().getInternalRotor();
        //#]
    }

    //## operation getMolecularWeight()
    public double getMolecularWeight() {
        //#[ operation getMolecularWeight()
        return getChemGraph().getMolecularWeight();
        //#]
    }

	   //## operation getNasaThermoData()
    public NASAThermoData getNasaThermoData() {
        //#[ operation getNasaThermoData()
        //if (nasaThermoData==null && !therfitExecuted) generateNASAThermoData();
		if (nasaThermoData==null) generateNASAThermoData();
        return nasaThermoData;
        //#]
    }

    //## operation getResonanceIsomers()
    public Iterator getResonanceIsomers() {
        //#[ operation getResonanceIsomers()
        return resonanceIsomers.iterator();
        //#]
    }

	public HashSet getResonanceIsomersHashSet(){
		return resonanceIsomers;
	}

    /**
    Requires:
    Effects: return the thermo data of the stablest resonance isomer.
    Modifies:
    */
    //## operation getThermoData()
    public ThermoData getThermoData() {
        //#[ operation getThermoData()
        return chemGraph.getThermoData();
        //#]
    }

    //## operation getThreeFrequencyMode()
    public ThreeFrequencyModel getThreeFrequencyMode() {
        //#[ operation getThreeFrequencyMode()
        if (threeFrequencyModel==null && !therfitExecuted) generateThreeFrequencyModel();

        return threeFrequencyModel;
        //#]
    }

    //## operation hasResonanceIsomers()
    public boolean hasResonanceIsomers() {
        //#[ operation hasResonanceIsomers()
        if (resonanceIsomers == null) return false;
        else return (resonanceIsomers.size() > 0);
        //#]
    }

    //## operation hasThreeFrequencyModel()
    public boolean hasThreeFrequencyModel() {
        //#[ operation hasThreeFrequencyModel()
        return (getThreeFrequencyModel() != null);
        //#]
    }

    //## operation isRadical()
    public boolean isRadical() {
        //#[ operation isRadical()
        return chemGraph.isRadical();
        //#]
    }

    //## operation isTriatomicOrSmaller()
    public boolean isTriatomicOrSmaller() {
        //#[ operation isTriatomicOrSmaller()
        return (getChemGraph().getAtomNumber()<=3);
        //#]
    }

    //## operation make(String,ChemGraph)
    public static Species make(String p_name, ChemGraph p_chemGraph) {
        //#[ operation make(String,ChemGraph)
        SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
        Species spe = (Species)(dictionary.getSpecies(p_chemGraph));
		
        if (spe == null) {
        	String name = p_name;
        	if (name == null || name.length()==0) {
        		name = p_chemGraph.getChemicalFormula();
        	}
			int id= ++TOTAL_NUMBER;
        	spe = new Species(id,name,p_chemGraph);
        	//spe.ID =
        	dictionary.putSpecies(spe, true);

        }
        p_chemGraph.setSpecies(spe);
        return spe;
        //#]
    }

	public static Species make(String p_name, ChemGraph p_chemGraph, int id) {
        //#[ operation make(String,ChemGraph)
        SpeciesDictionary dictionary = SpeciesDictionary.getInstance();
        Species spe = (Species)(dictionary.getSpecies(p_chemGraph));

        if (spe == null) {
        	String name = p_name;
        	if (name == null || name.length()==0) {
        		name = p_chemGraph.getChemicalFormula();
        	}
        	spe = new Species(id, name,p_chemGraph);
			if (id > TOTAL_NUMBER) TOTAL_NUMBER=id;
        	dictionary.putSpecies(spe, false);
 
        }
        p_chemGraph.setSpecies(spe);
        return spe;
        //#]
    }

    //## operation repOk()
    public boolean repOk() {
        //#[ operation repOk()
        if (name == null || name.length() == 0) return false;

        if (ID < 0 || ID > TOTAL_NUMBER) return false;

        if (chemGraph != null && !chemGraph.repOk()) return false;

        Iterator iter = resonanceIsomers.iterator();
        while (iter.hasNext()) {
        	ChemGraph cg = (ChemGraph)iter.next();
        	if (cg == null || !cg.repOk()) return false;
        }

        return true;
        //#]
    }

    //## operation selectDeltaEDown()
    public void selectDeltaEDown() {
        //#[ operation selectDeltaEDown()
        String name = getName();

        if (name.equals("AR")) deltaEDown = 374.0;
        else if (name.equals("H2")) deltaEDown = 224.0;
        else if (name.equals("O2")) deltaEDown = 517.0;
        else if (name.equals("N2")) deltaEDown = 461.0;
        else if (name.equals("He")) deltaEDown = 291.0;
        else if (name.equals("CH4")) deltaEDown = 1285.0;
        else if (name.equals("HO2")) deltaEDown = 975.0;
        else if (name.equals("H2O2")) deltaEDown = 975.0;
        else if (name.equals("CO2")) deltaEDown = 417.0;
        else if (name.equals("CO")) deltaEDown = 283.0;

        return;



        //#]
    }

    //## operation selectLJParametersForSpecialMolecule()
    public void selectLJParametersForSpecialMolecule() {
        //#[ operation selectLJParametersForSpecialMolecule()
        String name = getName();
        if (name.equals("H2")) LJ = new LennardJones(2.8327, 59.7);
        else if (name.equals("O2")) LJ = new LennardJones(3.467, 106.7);
        else if (name.equals("H2O")) LJ = new LennardJones(2.641,809.1);
        else if (name.equals("H2O2")) LJ = new LennardJones(4.196,289.3);
        else if (name.equals("CO2")) LJ = new LennardJones(3.941,195.2);
        else if (name.equals("CO")) LJ = new LennardJones(3.690,91.7);
        return;
        //#]
    }
	public String toChemkinString() {
		//#[ operation toChemkinString()
		return getChemkinName();
		//#]
	}

    //## operation toChemDisString()
    public String toChemDisString() {
        //#[ operation toChemDisString()
        return "SPC" + String.valueOf(getID());
        //#]
    }

    //## operation toString()
    public String toString() {
        //#[ operation toString()
        String s = "Species " + String.valueOf(ID) + '\t';
        s = s + "Name: " + getName() + '\n';
        ChemGraph primary = getChemGraph();
        s = s + primary.toString();

        int index=0;
        for (Iterator iter = resonanceIsomers.iterator(); iter.hasNext(); ) {
        	index++;
        	ChemGraph isomer = (ChemGraph)iter.next();
        	if (isomer!=primary) {
        		s = s + '\n';
        		s = s + "isomer" + String.valueOf(index) + ":\n";
        		s = s + isomer.toString();
        	}
        }
        return s;
    }

//		## operation toString()
	    public String toString(int i) {
	        //#[ operation toString()
			String s ="";
			ChemGraph primary = getChemGraph();
	        s = s + primary.toString(i);

	        int index=0;
	        /*for (Iterator iter = resonanceIsomers.iterator(); iter.hasNext(); ) {
	        	index++;
	        	ChemGraph isomer = (ChemGraph)iter.next();
	        	if (isomer!=primary) {
	        		s = s + '\n';
	        		s = s + "isomer" + String.valueOf(index) + ":\n";
	        		s = s + isomer.toStringWithoutH(i);
	        	}
	        }*/
	        return s;
    }

    //## operation toStringWithoutH()
    public String toStringWithoutH() {
        //#[ operation toStringWithoutH()
        String s = "Species " + String.valueOf(ID) + '\t';
        s = s + "Name: " + getName() + '\n';
        ChemGraph primary = getChemGraph();
        s = s + primary.toStringWithoutH();

        int index=0;
        for (Iterator iter = resonanceIsomers.iterator(); iter.hasNext(); ) {
        	index++;
        	ChemGraph isomer = (ChemGraph)iter.next();
        	if (isomer!=primary) {
        		s = s + '\n';
        		s = s + "isomer" + String.valueOf(index) + ":\n";
        		s = s + isomer.toStringWithoutH();
        	}
        }
        return s;




        //#]
    }

	public String toStringWithoutH(int i) {
        //#[ operation toStringWithoutH()
        String s ="";
		ChemGraph primary = getChemGraph();
        s = s + primary.toStringWithoutH(i);

        int index=0;
        /*for (Iterator iter = resonanceIsomers.iterator(); iter.hasNext(); ) {
        	index++;
        	ChemGraph isomer = (ChemGraph)iter.next();
        	if (isomer!=primary) {
        		s = s + '\n';
        		s = s + "isomer" + String.valueOf(index) + ":\n";
        		s = s + isomer.toStringWithoutH(i);
        	}
        }*/
        return s;




        //#]
    }

    public int getID() {
        return ID;
    }

    public static int getTOTAL_NUMBER() {
        return TOTAL_NUMBER;
    }

    public ChemGraph getChemGraph() {
        return chemGraph;
    }

    public void setChemGraph(ChemGraph p_chemGraph) {
        chemGraph = p_chemGraph;
    }

    public double getDeltaEDown() {
        return deltaEDown;
    }

    public String getName() {
        return name;
    }

    public void setName(String p_name) {
        name = p_name;
    }

    public boolean getTherfitExecuted() {
        return therfitExecuted;
    }

    public LennardJones getLJ() {
        return LJ;
    }

    public LennardJones newLJ() {
        LJ = new LennardJones();
        return LJ;
    }

    public void deleteLJ() {
        LJ=null;
    }

    public ThreeFrequencyModel getThreeFrequencyModel() {
        return threeFrequencyModel;
    }

    protected void initRelations() {
        LJ = newLJ();
    }

    public boolean isReactive() {
    	return IsReactive;
    }	
    
    public void setReactivity(boolean reactive) {
    	IsReactive = reactive;
    }
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\Species.java
*********************************************************************/

