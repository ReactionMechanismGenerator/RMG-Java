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


package jing.chemParser;


import java.io.*;

import jing.rxn.*;
import jing.chem.*;

import java.util.*;
import jing.mathTool.*;
import jing.chemUtil.*;
import jing.rxn.Reaction;
import jing.chemUtil.Graph;
import jing.chem.ThermoGAValue;
import jing.rxn.ArrheniusKinetics;
import jing.chemUtil.HierarchyTree;
import jing.rxn.ArrheniusEPKinetics;
import jing.rxnSys.CoreEdgeReactionModel;

//## package jing::chemParser

//----------------------------------------------------------------------------
// jing\chemParser\ChemParser.java
//----------------------------------------------------------------------------

//## class ChemParser
public class ChemParser {

    protected static String COMMENT_SIGN = "//";		//## attribute COMMENT_SIGN

    final protected static int MAX_TREE_LEVEL = 20;		//## attribute MAX_TREE_LEVEL


    // Constructors

    //## operation ChemParser()
    private  ChemParser() {
        //#[ operation ChemParser()
        //#]
    }

    //## operation extractReactionSeperator(String)
    public static final String extractReactionSeperator(String p_string) {
        //#[ operation extractReactionSeperator(String)
        if (p_string == null) throw new NullPointerException();

        StringTokenizer st = new StringTokenizer(p_string);
        while (st.hasMoreTokens()) {
        	String s = st.nextToken().trim();
        	if (s.equals("<=>")) return "<=>";
        	else if (s.equals("=")) return "=";
        	else if (s.equals("=>")) return "=>";
        	else if (s.equals("->")) return "->";
        }
        return null;
        //#]
    }

    /**
     * Recursively goes inside a directory and deltes everything inside.
     */
    public static boolean deleteDir(File dir) {
        if (dir.isDirectory()) {
            String[] children = dir.list();
            for (int i=0; i<children.length; i++) {
                boolean success = deleteDir(new File(dir, children[i]));
                if (!success) {
                    return false;
                }
            }
        }
    
        // The directory is now empty so delete it
        return dir.delete();
    }
    
    // Argument Stringp_name :
    /**
    the name of union
    */
    // Argument HashMapp_unRead :
    /**
    the map including all the union, each mapping maps the unread union group's name with their union components' names.
    */
    // Argument HashMapp_dictionary :
    /**
    The real dicationary where we can find the graph information from the name.
    */
    //## operation findUnion(String,HashMap,HashMap)
    public static void findUnion(String p_name, HashMap p_unRead, HashMap p_dictionary) {
        //#[ operation findUnion(String,HashMap,HashMap)
        FunctionalGroupCollection fgc = new FunctionalGroupCollection(p_name);

        // get the union set according to the p_name
        HashSet union = (HashSet)p_unRead.get(p_name);
        Iterator union_iter = union.iterator();
        while (union_iter.hasNext()) {
        	String fg_name = (String)union_iter.next();
        	// if the name of a union componenet is another union group, call findUnion() recusively
        	if (p_unRead.containsKey(fg_name)) {
        		findUnion(fg_name, p_unRead, p_dictionary);
        	}
        	Matchable fg = (Matchable)p_dictionary.get(fg_name);
        	if (fg == null)
        		throw new InvalidFunctionalGroupException("Unknown FunctionalGroup in findUnion(): " + fg_name);
        	fgc.addFunctionalGroups(fg);
        }
        p_dictionary.put(p_name,fgc);
        p_unRead.remove(p_name);
        return;



        //#]
    }

    //## operation parseArrheniusEPKinetics(String)
    /**
     * Input parameters updated by MRH on 8-Jun-2009
     * 	This function is called when reading a kinetics Library.txt file,
     * 	if the kinetics tag is Arrhenius_EP.  As of today, all rxn families
     * 	in the RMG database are labeled Arrhenius_EP.  Instead of passing only
     * 	the information after the node(s), i.e. Trange, A, n, alpha, E0, etc.,
     *  the entire data entry is passed.
     *  
     *  The string p_string contains the entire line of data (including the randomly
     *  assigned number and node(s)).  The integer p_keyNum contains the number of
     *  nodes associated with this data entry (e.g. 1 for Beta_Scission, 2 for 
     *  H_Abstraction, 3 for Diels_Alder_Addition, etc.).
     *  
     *  The output is of type ArrheniusEPKinetics.  The difference is that the
     *  "comments" and "source" fields are no longer null.  The "source" contains
     *  the set of nodes associated with this data, and the "comments" are
     *  the end-of-line comments, hopefully denoting where this data came from.
     **/
    public static ArrheniusEPKinetics parseArrheniusEPKinetics(String p_string, int p_keyNum) {
        //#[ operation parseArrheniusEPKinetics(String)
        StringTokenizer token = new StringTokenizer(p_string);
        //if (token.countTokens() != 10) throw new InvalidKineticsFormatException();
        /*
         * If statement commented out by MRH on 8-Jun-2009
         * 	The if statement existed because RMG was looking for 10 tokens:
         * 		(Trange, A, n, alpha, E0, dA, dn, dalpha, dE0, and Rank).
         * 		RMG now has no limit on the number of tokens, which allows
         * 		a RMG user/developer to put comments after the rank
         */
        String dummyCounter = token.nextToken();	// This should be the #. associated with the data
        // Set the source of this data as the set of nodes
        String source = "";
        for (int i=0; i<p_keyNum-1; i++) {
        	source += token.nextToken() + " ";
        }
        source += token.nextToken();
        String TRange = token.nextToken();

        double A = Double.parseDouble(token.nextToken());
        if (A<0) throw new NegativeAException("Negative A:" + String.valueOf(A));
        double n = Double.parseDouble(token.nextToken());
        double alpha = Double.parseDouble(token.nextToken());
        double E = Double.parseDouble(token.nextToken());

        String s = token.nextToken();
        double DA;
        if (s.startsWith("*")) {
        	s = s.substring(1,s.length());
        	DA = Math.abs(Double.parseDouble(s));
        	if (DA<1) throw new InvalidUncertaintyException("Multiplier Uncertainty of A (<1): " + String.valueOf(DA));
        }
        else {
        // if not multiplier uncertain, transfer it into multiplier for A
        	DA = Math.abs(Double.parseDouble(s));
        	if (DA == 0) DA = 1;
        	else if (DA > A) throw new NegativeAException("lower bound of A(<0): " + String.valueOf(A-DA));
        	else {
        		DA = A/(A-DA);
        	}
        }

        double Dn = Math.abs(Double.parseDouble(token.nextToken()));
        double Dalpha = Math.abs(Double.parseDouble(token.nextToken()));
        double DE = Math.abs(Double.parseDouble(token.nextToken()));
        int rank = Integer.parseInt(token.nextToken());

        UncertainDouble ua = new UncertainDouble(A, DA, "Multiplier");
        UncertainDouble un = new UncertainDouble(n, Dn, "Adder");
        UncertainDouble ualpha = new UncertainDouble(alpha, Dalpha, "Adder");
        UncertainDouble ue = new UncertainDouble(E, DE, "Adder");

        // Set the comments of this data to whatever (if anything) follows the rank
        String comments = "";
        while (token.hasMoreTokens())
			comments += token.nextToken() + " ";

		/*if (token.hasMoreTokens()) {
	        String beginningOfComments = token.nextToken();
	        String[] splitString = p_string.split(beginningOfComments);
	        comments = beginningOfComments + splitString[1];
        }*/
        ArrheniusEPKinetics k = new ArrheniusEPKinetics(ua, un, ualpha, ue, TRange, rank, source, comments);
//        ArrheniusEPKinetics k = new ArrheniusEPKinetics(ua, un, ualpha, ue, TRange, rank, null, null);

        return k;



        //#]
    }

    //## operation parseArrheniusKinetics(String)
    /**
     * Input parameters updated by MRH on 11-Jun-2009
     * 	This function is called when reading a kinetics Library.txt file,
     * 	if the kinetics tag is Arrhenius.  As of today, no rxn families
     * 	in the RMG database are labeled Arrhenius.  Instead of passing only
     * 	the information after the node(s), i.e. Trange, A, n, Ea, etc.,
     *  the entire data entry is passed.
     *  
     *  The string p_string contains the entire line of data (including the randomly
     *  assigned number and node(s)).  The integer p_keyNum contains the number of
     *  nodes associated with this data entry (e.g. 1 for Beta_Scission, 2 for 
     *  H_Abstraction, 3 for Diels_Alder_Addition, etc.).
     *  
     *  The output is of type ArrheniusKinetics.  The difference is that the
     *  "comments" and "source" fields are no longer null.  The "source" contains
     *  the set of nodes associated with this data, and the "comments" are
     *  the end-of-line comments, hopefully denoting where this data came from.
     **/
    public static ArrheniusKinetics parseArrheniusKinetics(String p_string, int p_keyNum) {
        //#[ operation parseArrheniusKinetics(String)
        StringTokenizer token = new StringTokenizer(p_string);
        /*
         * If statement commented out by MRH on 11-Jun-2009
         * 	The if statement existed because RMG was looking for 8 tokens:
         * 		(Trange, A, n, Ea, dA, dn, dEa, and Rank).
         * 		RMG now has no limit on the number of tokens, which allows
         * 		a RMG user/developer to put comments after the rank
         */
//        if (token.countTokens() != 8) throw new InvalidKineticsFormatException();
        String dummyCounter = token.nextToken();	// This should be the #. associated with the data
        // Set the source of this data as the set of nodes
        String source = "";
        for (int i=0; i<p_keyNum-1; i++) {
        	source += token.nextToken() + " ";
        }
        source += token.nextToken();
        String TRange = token.nextToken();

        double A = Double.parseDouble(token.nextToken());
        if (A<0) throw new NegativeAException("Negative A:" + String.valueOf(A));
        double n = Double.parseDouble(token.nextToken());
        double E = Double.parseDouble(token.nextToken());

        double DA = 1;
        String s = token.nextToken();
        if (s.startsWith("*")) {
        	s = s.substring(1,s.length());
        	DA = Math.abs(Double.parseDouble(s));
        	if (DA<1) throw new InvalidUncertaintyException("Multiplier Uncertainty of A (<1): " + String.valueOf(DA));
        }
        else {
        // if not multiplier uncertain, transfer it into multiplier for A
        DA = Math.abs(Double.parseDouble(s));
        if (DA == 0) DA = 1;
        	else if (DA > A) throw new NegativeAException("lower bound of A(<0): " + String.valueOf(A-DA));
        	else {
        		DA = A/(A-DA);
        	}
        }
        double Dn = Double.parseDouble(token.nextToken());
        double DE = Double.parseDouble(token.nextToken());
        int rank = Integer.parseInt(token.nextToken());


        UncertainDouble ua = new UncertainDouble(A, DA, "Multiplier");
        UncertainDouble un = new UncertainDouble(n, Dn, "Adder");
        UncertainDouble ue = new UncertainDouble(E, DE, "Adder");

        // Set the comments of this data to whatever (if anything) follows the rank
        String comments = "";
        if (token.hasMoreTokens()) {
	        String beginningOfComments = token.nextToken();
	        String[] splitString = p_string.split(beginningOfComments);
	        comments = beginningOfComments + splitString[1];
        }
        ArrheniusKinetics k = new ArrheniusKinetics(ua, un, ue, TRange, rank, source, comments);
//        ArrheniusKinetics k = new ArrheniusKinetics(ua, un, ue, TRange, rank, null, null);

        return k;



        //#]
    }

    /*
     * MRH 23MAR2010:
     * Commented out uncalled parseArrheniusReaction method
     */
    //## operation parseArrheniusReaction(HashMap,String,double,double)
//    public static Reaction parseCoreArrheniusReaction(SpeciesDictionary p_species, String p_reactionString, double p_AMultiplier, double p_EMultiplier) {
//        //#[ operation parseArrheniusReaction(HashMap,String,double,double)
//        //boolean isReverse = false;
//		if (p_reactionString == null) throw new NullPointerException("parseArrheniusReaction");
//
//        StringTokenizer st = new StringTokenizer(p_reactionString);
//        int size = st.countTokens();
//        if (size < 6 || size > 16) throw new InvalidReactionFormatException();
//
//        int cut = size - 3;
//        String structureString = "";
//        String arrheniusString = "";
//
//        for (int i=0; i<cut-2; i++) {
//        	structureString = structureString + st.nextToken() + " ";
//        }
//		int direction = Integer.parseInt(st.nextToken());
//		int redundancy = Integer.parseInt(st.nextToken());
//
//        for (int i=cut; i<size; i++) {
//        	arrheniusString = arrheniusString + st.nextToken() + " ";
//        }
//
//        // parse structure string
//        String sep = extractReactionSeperator(structureString);
//        boolean generateReverse;
//
//        if (sep.equals("=") || sep.equals("<=>")) generateReverse = true;
//        else if (sep.equals("=>") || sep.equals("->")) generateReverse = false;
//        else throw new InvalidStructureException("Unknown reaction seperator: " + structureString);
//
//        st = new StringTokenizer(structureString,sep);
//        if (st.countTokens() != 2) throw new InvalidStructureException("Unknown format:" + structureString);
//
//        LinkedList r = parseReactionSpecies(p_species, st.nextToken());
//        LinkedList p = parseReactionSpecies(p_species, st.nextToken());
//
//        Structure s = new Structure(r,p);
//        s.setDirection(1);
//
//        // parse kinetics
//        ArrheniusKinetics k = parseSimpleArrheniusKinetics(arrheniusString, p_AMultiplier, p_EMultiplier, r.size());
//
//        return Reaction.makeReaction(s,k,generateReverse);
//        //#]
//    }

    /*
     * MRH 23MAR2010:
     * 	Commented out uncalled parseArrheniusReaction method
     */
//	## operation parseArrheniusReaction(HashMap,String,double,double)
//    public static Reaction parseEdgeArrheniusReaction(SpeciesDictionary p_species, String p_reactionString, double p_AMultiplier, double p_EMultiplier) {
//        //#[ operation parseArrheniusReaction(HashMap,String,double,double)
//        //boolean isReverse = false;
//		if (p_reactionString == null) throw new NullPointerException("parseArrheniusReaction");
//
//        StringTokenizer st = new StringTokenizer(p_reactionString);
//        int size = st.countTokens();
//        if (size < 4 || size > 14) throw new InvalidReactionFormatException();
//
//        int cut = size - 3;
//        String structureString = "";
//        String arrheniusString = "";
//
//        for (int i=0; i<cut; i++) {
//        	structureString = structureString + st.nextToken() + " ";
//        }
//        for (int i=cut; i<size; i++) {
//        	arrheniusString = arrheniusString + st.nextToken() + " ";
//        }
//
//        // parse structure string
//        String sep = extractReactionSeperator(structureString);
//        boolean generateReverse;
//
//        if (sep.equals("=") || sep.equals("<=>")) generateReverse = true;
//        else if (sep.equals("=>") || sep.equals("->")) generateReverse = false;
//        else throw new InvalidStructureException("Unknown reaction seperator: " + structureString);
//
//        st = new StringTokenizer(structureString,sep);
//        if (st.countTokens() != 2) throw new InvalidStructureException("Unknown format:" + structureString);
//
//        LinkedList r = parseReactionSpecies(p_species, st.nextToken());
//        LinkedList p = parseReactionSpecies(p_species, st.nextToken());
//        if (generateReverse){
//			Structure s = new Structure(p,r);
//	        s.setDirection(1);
//
//	        // parse kinetics
//	        ArrheniusKinetics k = parseSimpleArrheniusKinetics(arrheniusString, p_AMultiplier, p_EMultiplier, r.size());
//
//	        Reaction forward = Reaction.makeReaction(s,k,generateReverse);
//
//			return forward.getReverseReaction();
//        }
//        else {
//			Structure s = new Structure(r,p);
//	        s.setDirection(1);
//
//	        // parse kinetics
//	        ArrheniusKinetics k = parseSimpleArrheniusKinetics(arrheniusString, p_AMultiplier, p_EMultiplier, r.size());
//
//	        //Reaction forward = Reaction.makeReaction(s,k,generateReverse);
//
//			return Reaction.makeReaction(s,k,generateReverse);
//        }
//
//        //#]
//    }

    /*
     * MRH 23MAR2010:
     * 	Commented out uncalled parseArrheniusReaction method
     */
//	## operation parseArrheniusReaction(HashMap,String,double,double)
//    public static Reaction parseArrheniusReaction(SpeciesDictionary p_species, String p_reactionString, double p_AMultiplier, double p_EMultiplier, CoreEdgeReactionModel cerm) {
//        //#[ operation parseArrheniusReaction(HashMap,String,double,double)
//        //boolean isReverse = false;
//		if (p_reactionString == null) throw new NullPointerException("parseArrheniusReaction");
//
//        StringTokenizer st = new StringTokenizer(p_reactionString);
//        int size = st.countTokens();
//        if (size < 6 || size > 16) throw new InvalidReactionFormatException();
//
//        int cut = size - 3;
//        String structureString = "";
//        String arrheniusString = "";
//
//        for (int i=0; i<cut-2; i++) {
//        	structureString = structureString + st.nextToken() + " ";
//        }
//		int direction = Integer.parseInt(st.nextToken());
//		int redundancy = Integer.parseInt(st.nextToken());
//
//        for (int i=cut; i<size; i++) {
//        	arrheniusString = arrheniusString + st.nextToken() + " ";
//        }
//
//        // parse structure string
//        String sep = extractReactionSeperator(structureString);
//        boolean generateReverse;
//
//        if (sep.equals("=") || sep.equals("<=>")) generateReverse = true;
//        else if (sep.equals("=>") || sep.equals("->")) generateReverse = false;
//        else throw new InvalidStructureException("Unknown reaction seperator: " + structureString);
//
//        st = new StringTokenizer(structureString,sep);
//        if (st.countTokens() != 2) throw new InvalidStructureException("Unknown format:" + structureString);
//
//        LinkedList r = parseReactionSpecies(p_species, st.nextToken());
//        LinkedList p = parseReactionSpecies(p_species, st.nextToken());
//		Structure s = new Structure(r,p);
//		int i = cerm.categorizeReaction(s);
//		if (i== -1){
//	        if (generateReverse){
//				s = new Structure(p,r);
//		        s.setDirection(1);
//
//		        // parse kinetics
//		        ArrheniusKinetics k = parseSimpleArrheniusKinetics(arrheniusString, p_AMultiplier, p_EMultiplier, r.size());
//				//s.setRedundancy(redundancy);;
//		        Reaction forward = Reaction.makeReaction(s,k,generateReverse);
//
//				return forward.getReverseReaction();
//	        }
//	        else {
//				s = new Structure(r,p);
//		        s.setDirection(1);
//
//		        // parse kinetics
//		        ArrheniusKinetics k = parseSimpleArrheniusKinetics(arrheniusString, p_AMultiplier, p_EMultiplier, r.size());
//				//s.setRedundancy(redundancy);
//		        //Reaction forward = Reaction.makeReaction(s,k,generateReverse);
//
//				return Reaction.makeReaction(s,k,generateReverse);
//	        }
//		}
//		else if (i==1){
//			if (direction == 1){
//				s = new Structure(r,p);
//		        s.setDirection(1);
//
//		        // parse kinetics
//		        ArrheniusKinetics k = parseSimpleArrheniusKinetics(arrheniusString, p_AMultiplier, p_EMultiplier, r.size());
//				//s.setRedundancy(redundancy);
//		        return Reaction.makeReaction(s,k,generateReverse);
//			}
//
//			if (direction == -1){
//
//				s = new Structure(p,r);
//		        s.setDirection(1);
//
//		        // parse kinetics
//		        ArrheniusKinetics k = parseSimpleArrheniusKinetics(arrheniusString, p_AMultiplier, p_EMultiplier, r.size());
//		        Reaction forwardReaction = Reaction.makeReaction(s,k,generateReverse);
//		        return forwardReaction.getReverseReaction();
//			}
//		}
//
//		return null;
//        //#]
//    }
    
    /**
     * ChemParser.parseRestartReaction(p_rxnString, coreSpcsIDs, core_edge) method
     * 		Method interprets a string (in the Reaction.toChemkinString() format) and
     * 		returns a Reaction
     * 
     * @param p_rxnString: String listing the reactants, separator ("=" or "=>"), products, Arrhenius
     * 		parameters, and comments
     * @param coreSpcsIDs: int[] listing the core species IDs
     * @param core_edge: String indicating whether the rxn of interest is a "core" or "edge" rxn
     * @return
     */
    public static Reaction parseRestartReaction(String p_rxnString, int[] coreSpcsIDs, String core_edge, String eaUnits) {
    	StringTokenizer st = new StringTokenizer(p_rxnString);
    	//	First token is the rxn structure: A+B=C+D
    	//		Note: Up to 3 reactants/products allowed
    	//			: Either "=" or "=>" will separate reactants and products
    	String structure = st.nextToken();    	
    	//	Separate the reactants from the products
    	boolean generateReverse = false;
    	String[] reactsANDprods = null;
    	reactsANDprods = structure.split("\\=>");
    	// If the length of reactsANDprods = 1, no "=>" token, meaning
    	//	reverse rxn
    	if (reactsANDprods.length == 1) {
    		reactsANDprods = null;
    		reactsANDprods = structure.split("\\=");
    		generateReverse = true;
    	}
    	
    	SpeciesDictionary sd = SpeciesDictionary.getInstance();
        LinkedList r = parseReactionSpecies(sd, reactsANDprods[0]);
        LinkedList p = parseReactionSpecies(sd, reactsANDprods[1]);

        Structure s = new Structure(r,p);
        s.setDirection(1);
        
        /*
         * All rxns in restart files are listed in the forward direction.
         * 	If the rxn is from the core, both forward and backward rxn should
         * 	be stored in the CoreEdgeReactionModel; if from the edge, only
         * 	the rxn containing core species as reactants should be stored
         * 	as edge rxns in the CoreEdgeReactionModel.
         */
        boolean storeAsForward = true;
        if (core_edge.equals("edge")) {
        	boolean foundInCore = false;
        	for (Iterator iter = r.iterator(); iter.hasNext();) {
        		Species spc = (Species)iter.next();
        		for (int i=0; i<coreSpcsIDs.length; i++) {
        			// If the species is identified as core, move on to the next species
        			if (coreSpcsIDs[i] == spc.getID()) {
        				foundInCore = true;
        				break;
        			}
        		}
        		// If the spc's ID does not match any of the core species' IDs, it is
        		//	an edge species and the reverse rxn (p-->r) should be stored as an
        		//	edge rxn in the CoreEdgeReactionModel
        		if (!foundInCore) {
        			storeAsForward = false;
        			break;
        		}
        	}
        }
    	
        //	Next three tokens are the modified Arrhenius parameters
    	double rxn_A = Double.parseDouble(st.nextToken());
    	double rxn_n = Double.parseDouble(st.nextToken());
    	double rxn_E = Double.parseDouble(st.nextToken());
    	
    	//Convert rxn_E to kcal/mol units
		if (eaUnits.equals("cal/mol"))
			rxn_E = rxn_E / 1000;
		else if (eaUnits.equals("J/mol"))
			rxn_E = rxn_E / 4.184 / 1000;
		else if (eaUnits.equals("kJ/mol"))
			rxn_E = rxn_E / 4.184;
		else if (eaUnits.equals("Kelvins"))
			rxn_E = rxn_E * 1.987;
    	
    	//	For now, restart file does not contain uncertainties
    	UncertainDouble uA = new UncertainDouble(rxn_A,0.0,"A");
    	UncertainDouble un = new UncertainDouble(rxn_n,0.0,"A");
    	UncertainDouble uE = new UncertainDouble(rxn_E,0.0,"A");    	
    	
    	//	The remaining tokens are comments
    	String comments = "";
    	if (st.hasMoreTokens()) {
    		String beginningOfComments = st.nextToken();
    		int startIndex = p_rxnString.indexOf(beginningOfComments);
    		comments = p_rxnString.substring(startIndex);
    	}
    	if (comments.startsWith("!")) comments = comments.substring(1);
//    	while (st.hasMoreTokens()) {
//    		comments += st.nextToken();
//    	}
    	
    	// Generate the kinetics (assuming a rank of 1 ... as of 7/Sept/2009, the rank
    	//	of Kinetics should not be important at this stage of the mechanism
    	//	generation)
    	ArrheniusKinetics[] k = new ArrheniusKinetics[1];
        k[0] = new ArrheniusKinetics(uA,un,uE,"",1,"",comments);
        if (storeAsForward)
        	return Reaction.makeReaction(s,k,generateReverse);
        else
        	return Reaction.makeReaction(s.generateReverseStructure(),k,generateReverse);
    }

  //## operation parseArrheniusReaction(HashMap,String,double,double)
  public static Reaction parseArrheniusReaction(HashMap p_species, String p_reactionString, double p_AMultiplier, double p_EMultiplier) {
      //#[ operation parseArrheniusReaction(HashMap,String,double,double)
      if (p_reactionString == null) throw new NullPointerException("parseArrheniusReaction");

      StringTokenizer st = new StringTokenizer(p_reactionString);
      int size = st.countTokens();
      if (size < 7 || size > 17) throw new InvalidReactionFormatException();

      int cut = size - 6;
      String structureString = "";
      String arrheniusString = "";

      for (int i=0; i<cut; i++) {
              structureString = structureString + st.nextToken() + " ";
      }
      for (int i=cut; i<size; i++) {
              arrheniusString = arrheniusString + st.nextToken() + " ";
      }

      // parse structure string
      String sep = extractReactionSeperator(structureString);
      boolean generateReverse;

      if (sep.equals("=") || sep.equals("<=>")) generateReverse = true;
      else if (sep.equals("=>") || sep.equals("->")) generateReverse = false;
      else throw new InvalidStructureException("Unknown reaction seperator: " + structureString);

      st = new StringTokenizer(structureString,sep);
      if (st.countTokens() != 2) throw new InvalidStructureException("Unknown format:" + structureString);

      LinkedList r = parseReactionSpecies(p_species, st.nextToken());
      LinkedList p = parseReactionSpecies(p_species, st.nextToken());

      Structure s = new Structure(r,p);
      s.setDirection(1);

      // parse kinetics
      //ArrheniusKinetics k = parseSimpleArrheniusKinetics(arrheniusString, p_AMultiplier, p_EMultiplier);
      ArrheniusKinetics[] k = new ArrheniusKinetics[1];
      k[0] = parsePrimaryArrheniusKinetics(arrheniusString, p_AMultiplier,p_EMultiplier, r.size());

      return Reaction.makeReaction(s,k,generateReverse);
      //#]
  }

    //## operation parseReactionSpecies(HashMap,String)
    public static LinkedList parseReactionSpecies(SpeciesDictionary p_speciesSet, String p_speciesString) {
        //#[ operation parseReactionSpecies(HashMap,String)
        if (p_speciesString == null) throw new NullPointerException();

        StringTokenizer st = new StringTokenizer(p_speciesString, "+");
        int speNum = st.countTokens();
        if (speNum > 3) throw new InvalidStructureException("too many reactants/products: " + p_speciesString);

        LinkedList reactionSpe = new LinkedList();

        for (int i = 0; i < speNum; i++) {
        	String name = st.nextToken().trim();
        	if (!name.toUpperCase().equals("M")) {
        		Species spe = (Species)p_speciesSet.getSpeciesFromNameID(name);
        		if (spe == null) throw new InvalidStructureException("unknown reactant/product: " + name);
        		reactionSpe.add(spe);
        	}
        }

        return reactionSpe;
        //#]
    }

    //## operation parsePrimaryArrheniusKinetics(String)
    //svp
    /*
     * Updated by MRH on 17Feb2010
     * Until a few days ago, RMG only allowed "A" units of "mol/cm3/s" or "mol/liter/s" when a user
     * 	supplied a primary reaction library or seed mechanism to RMG.  Very recently, MRH added the
     * 	option "molecule/cm3/s" for CFG to use.  What MRH was not aware of was that EVERY "A" factor
     * 	would then be multiplied by Avogadro's #, insteaad of just the bimolecular reactions.
     * 
     * 	This update addresses that issue: "A" factors are only multiplied by the "AMultiplier" when
     * 	the reaction is bimolecular.
     */
      public static ArrheniusKinetics parsePrimaryArrheniusKinetics(String p_string, double p_AMultiplier,double p_EMultiplier, int numReacts) {
          //#[ operation parseArrheniusKinetics(String)
          StringTokenizer token = new StringTokenizer(p_string);
          if (token.countTokens() != 6) throw new InvalidKineticsFormatException();

          double A = Double.parseDouble(token.nextToken());
          if (A<0) throw new NegativeAException("Negative A:" + String.valueOf(A));
          double n = Double.parseDouble(token.nextToken());
          double E = Double.parseDouble(token.nextToken());

          double DA = 1;
          String s = token.nextToken();
          if (s.startsWith("*")) {
                  s = s.substring(1,s.length());
                  DA = Math.abs(Double.parseDouble(s));
                  if (DA<1) throw new InvalidUncertaintyException("Multiplier Uncertainty of A (<1): " + String.valueOf(DA));
          }
          else {
          // if not multiplier uncertain, transfer it into multiplier for A
          DA = Math.abs(Double.parseDouble(s));
          if (DA == 0) DA = 1;
                  else if (DA > A) throw new NegativeAException("lower bound of A(<0): " + String.valueOf(A-DA));
                  else {
                          DA = A/(A-DA);
                  }
          }
          double Dn = Double.parseDouble(token.nextToken());
          double DE = Double.parseDouble(token.nextToken());

          UncertainDouble ua;
          if (numReacts == 1) {
        	  ua = new UncertainDouble(A, DA, "Multiplier");
          } else if (numReacts == 2) {
          	ua = new UncertainDouble(A*p_AMultiplier,1,"Multiplier");
          } else if (numReacts == 3) {
          	ua = new UncertainDouble(A*p_AMultiplier*p_AMultiplier,1,"Multiplier");
          } else {
        	  ua = null;
        	  System.err.println("In ChemParser.parseSimpleArrheniusKinetics:\nThe number of " +
        			  "reactants exceeds three.");
        	  System.exit(0);
          }
          //UncertainDouble ua = new UncertainDouble(A*p_AMultiplier, DA, "Multiplier");
          UncertainDouble un = new UncertainDouble(n, Dn, "Adder");
          UncertainDouble ue = new UncertainDouble(E*p_EMultiplier, DE, "Adder");

          ArrheniusKinetics k = new ArrheniusKinetics(ua, un, ue, "Unknown", 0, null, null);

          return k;



          //#]
      }


	 //## operation parseReactionSpecies(HashMap,String)
    public static LinkedList parseReactionSpecies(HashMap p_speciesSet, String p_speciesString) {
        //#[ operation parseReactionSpecies(HashMap,String)
        if (p_speciesString == null) throw new NullPointerException();

        StringTokenizer st = new StringTokenizer(p_speciesString, "+");
        int speNum = st.countTokens();
        if (speNum > 3) throw new InvalidStructureException("too many reactants/products: " + p_speciesString);

        LinkedList reactionSpe = new LinkedList();

        for (int i = 0; i < speNum; i++) {
        	String name = st.nextToken().trim();
        	String origname = name;
        	/*
        	 * 14Feb2010: MRH adding additional comments
        	 * 	This if statement exists to catch the "third body term" in the reaction string
        	 * 		These "M" species are not RMG species, just notation to represent a pdep reaction
        	 * 	A Lindemann/Troe reaction will have (+m) in the reaction string, thus we check
        	 * 		if the "name" is equal to "m)" after tokenizing the reaction string using the
        	 * 		plus "+" character.
        	 * A 3rdBodyReaction will have +m in the reaction string, thus we also check if the "name"
        	 * 		is equal to "m" after tokenizing the reaction string.
        	 */
        	if ((!name.equals("M"))&&(!name.equals("M)"))&&(!name.equals("m"))&&(!name.equals("m)"))) {//7/29/09 gmagoon: changed this, in consultation with MRH, to accept species that start with M...; previously, species beginning with M were not read as reactants or products
				if (name.endsWith("(")) {
					/*
					 * 14Feb2010: MRH adding additional comments
					 * If the reaction string is from a Lindemann/Troe reaction, the string (+m) will
					 * 		appear in the reacion string.  When tokenizing the string using "+", A+B(+m)
					 * 		will capture "A", "B(" and "m)".  We already handle the "m)"; see above.  We
					 * 		handle the "B(" with this if statement
					 */
        			name = name.substring(0,name.length()-1).trim();
        			origname = name;
        		}
				/*
				 * 14Feb2010: MRH adding additional comments
				 * If the species comes from a Restart file, a (#) will be tacked on at the end of
				 * 	the chemical formula.  We are getting the species using its name (see below) and
				 * 	the (#) is not part of the name.  Thus, we want to remove the (#) before searching
				 *  the speciesSet.
				 * HOWEVER, some species, e.g. CH2(S), have a () that we want to keep.
				 */
				if (name.endsWith(")")) {
					String[] tempString = name.split("\\(");
					name = tempString[0];
					for (int numSplits = 1; numSplits < tempString.length; numSplits++) {
						if (!tempString[numSplits].toLowerCase().equals(tempString[numSplits].toUpperCase())) {
							name += "(" + tempString[numSplits];
						}
					}
				}
        		Species spe = (Species)p_speciesSet.get(name);
        		if (spe == null) {
        			spe = (Species)p_speciesSet.get(origname);
        			if (spe == null) {
        				System.out.println("Error in parseReactionSpecies: RMG cannot find the following species in Dictionary: " + name);
        				throw new InvalidStructureException("unknown reactant/product: " + name);
        			}
        		}
        		reactionSpe.add(spe);
        	}
        }

        return reactionSpe;
        //#]
    }

    //## operation parseSimpleArrheniusKinetics(String,double,double)
    /*
     * MRH 17Feb2010:
     * 	BUG FIX: All pre-exponential factors, regardless of whether there were one or
     * 		two reactants, were multiplied by the AMultiplier.  The pre-exponential
     * 		factor is now altered only for a bimolecular reaction.
     */
    public static ArrheniusKinetics parseSimpleArrheniusKinetics(String p_string, double p_AMultiplier, double p_EMultiplier, int numReacts) {
        //#[ operation parseSimpleArrheniusKinetics(String,double,double)
        StringTokenizer token = new StringTokenizer(p_string);
        if (token.countTokens() != 3) throw new InvalidKineticsFormatException();

        double A = Double.parseDouble(token.nextToken());
        if (A<0) throw new NegativeAException("Negative A:" + String.valueOf(A));
        double n = Double.parseDouble(token.nextToken());
        double E = Double.parseDouble(token.nextToken());

        UncertainDouble ua;
        if (numReacts == 1) {
        	ua = new UncertainDouble(A,1,"Multiplier");
        } else if (numReacts == 2) {
        	ua = new UncertainDouble(A*p_AMultiplier,1,"Multiplier");
        } else if (numReacts == 3) {
        	ua = new UncertainDouble(A*p_AMultiplier*p_AMultiplier,1,"Multiplier");
        } else {
        	ua = null;
        	System.err.println("In ChemParser.parseSimpleArrheniusKinetics:\nThe number of " +
        			"reactants exceeds three.");
        	System.exit(0);
        }
        //UncertainDouble ua = new UncertainDouble(A*p_AMultiplier, 1, "Multiplier");
        UncertainDouble un = new UncertainDouble(n, 0, "Adder");
        UncertainDouble ue = new UncertainDouble(E*p_EMultiplier, 0, "Adder");

        ArrheniusKinetics k = new ArrheniusKinetics(ua, un, ue, "Unknown", 10000, null, null);

        return k;



        //#]
    }

    //## operation parseThermoFromLibrary(String)
    public static ThermoData parseThermoFromLibrary(String p_string) {//svp
      //#[ operation parseThermoFromLibrary(String)
      String s = p_string.trim();
      StringTokenizer token = new StringTokenizer(s);

      int data_num = token.countTokens();
      if (data_num != 12) throw new InvalidThermoFormatException();

      double H, S, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500, dH, dS, dCp;
      H = Double.parseDouble(token.nextToken());
      S = Double.parseDouble(token.nextToken());
      Cp300 = Double.parseDouble(token.nextToken());
      Cp400 = Double.parseDouble(token.nextToken());
      Cp500 = Double.parseDouble(token.nextToken());
      Cp600 = Double.parseDouble(token.nextToken());
      Cp800 = Double.parseDouble(token.nextToken());
      Cp1000 = Double.parseDouble(token.nextToken());
      Cp1500 = Double.parseDouble(token.nextToken());
      dH = Double.parseDouble(token.nextToken());
      dS = Double.parseDouble(token.nextToken());
      dCp = Double.parseDouble(token.nextToken());

      return new ThermoData(H, S, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500, dH, dS, dCp, null);
      //#]
    }


    //## operation parseThermoGAValue(String)
    public static ThermoGAValue parseThermoGAValue(String p_string) {
        //#[ operation parseThermoGAValue(String)
        String s = p_string.trim();
        StringTokenizer token = new StringTokenizer(s);

        int data_num = token.countTokens();
        if (data_num != 12) throw new InvalidThermoFormatException();

        double H, S, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500,dH,dS,dCp;//svp added dH,dS,dCp
        H = Double.parseDouble(token.nextToken());
        S = Double.parseDouble(token.nextToken());
        Cp300 = Double.parseDouble(token.nextToken());
        Cp400 = Double.parseDouble(token.nextToken());
        Cp500 = Double.parseDouble(token.nextToken());
        Cp600 = Double.parseDouble(token.nextToken());
        Cp800 = Double.parseDouble(token.nextToken());
        Cp1000 = Double.parseDouble(token.nextToken());
        Cp1500 = Double.parseDouble(token.nextToken());
        dH = Double.parseDouble(token.nextToken());
        dS = Double.parseDouble(token.nextToken());
        dCp = Double.parseDouble(token.nextToken());

        return new ThermoGAValue(H, S, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500, dH,dS,dCp,null);




        //#]
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Parser for Abraham Groups
    // Amrit Jalan, April 19, 2009

        //## operation parseAbrahamGAValue(String)
    public static AbrahamGAValue parseAbrahamGAValue(String p_string) {
        //#[ operation parseThermoGAValue(String)
        String s = p_string.trim();
        StringTokenizer token = new StringTokenizer(s);

        int data_num = token.countTokens();
        if (data_num != 5) throw new InvalidThermoFormatException();

        double S,B,E,L,A;
        S = Double.parseDouble(token.nextToken());
        B = Double.parseDouble(token.nextToken());
        E = Double.parseDouble(token.nextToken());
        L = Double.parseDouble(token.nextToken());
        A = Double.parseDouble(token.nextToken());

        return new AbrahamGAValue(S,B,E,L,A);




        //#]
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Parser for Unifac Groups
    // Amrit Jalan, April 19, 2009

        //## operation parseAbrahamGAValue(String)
    public static UnifacGAValue parseUnifacGAValue(String p_string) {
        //#[ operation parseThermoGAValue(String)
        String s = p_string.trim();
        StringTokenizer token = new StringTokenizer(s);

        int data_num = token.countTokens();
        if (data_num != 2) throw new InvalidThermoFormatException();

        double R,Q;
        R = Double.parseDouble(token.nextToken());
        Q = Double.parseDouble(token.nextToken());

        return new UnifacGAValue(R,Q);




        //#]
    }






    //## operation parseThirdBodyList(String)
    public static HashMap parseThirdBodyList(String p_string, HashMap speciesList) {
        //#[ operation parseThirdBodyList(String)
        if (p_string == null) throw new NullPointerException("read third body factor");

        HashMap thirdBodyList = new HashMap();
        StringTokenizer st = new StringTokenizer(p_string, "/");
        /*
         * MRH 23APR2010:
         * Allowing RMG to handle general third-body reactions,
         * 	e.g. A+M=B+C+M (where no specific third-body colliders are given)
         * 
         * The p_string should either contain a list of third-bodies and their
         * 	collision efficiencies (e.g. H2/2.0/ H2O/15.4/ CO/0.75/ etc.)
         * OR
         * the next reaction line (e.g. O+O+M=O2+M)
         * 
         * If the former, at least two tokens will exist; if the latter, only
         * one.  For the latter, we want to skip over the while loop and return
         * an empty thirdBodyList.
         */
        if (st.countTokens() == 1) return thirdBodyList;
        while (st.hasMoreTokens()) {
        	String name = st.nextToken().trim();
        	/*
        	 * When reading in third-bodies, RMG needs to recognize that these species
        	 * 	names may have changed, i.e. CH4 is now CH4(1)
        	 * This is straightforward, except for the "inert gas" species - 
        	 * 	N2, He, Ne, and Ar - because the chemgraphs should be supplied
        	 * 	in the dictionary.txt file
        	 */
        	if (!name.toLowerCase().equals("n2") && !name.toLowerCase().equals("ar") && !name.toLowerCase().equals("he") && !name.toLowerCase().equals("ne")) {
        		Species species = (Species)speciesList.get(name);
        		if (species == null) {
        			System.out.println("Error in reading third-body colliders: ");
        			System.out.println("The species '" + name + "' is not defined in the species.txt file.");
        			System.out.println("Terminating RMG run");
        			System.exit(0);
        		}
        		name = species.getChemkinName();
        	}
        	Double factor = Double.valueOf(st.nextToken().trim());
            thirdBodyList.put(name, factor);
        }

        return thirdBodyList;
        //#]
    }

    //## operation readActions(BufferedReader)
    public static LinkedList readActions(BufferedReader p_reader) {
        //#[ operation readActions(BufferedReader)
        LinkedList action = new LinkedList();

        String line = ChemParser.readUncommentLine(p_reader);
        while (line != null) {
        	StringTokenizer token = new StringTokenizer(line);
        	String index = token.nextToken();
        	String type = token.nextToken();
        	String act = token.nextToken();

        	LinkedList site = new LinkedList();
        	Object element;
        	act = removeBrace(act);
        	token = new StringTokenizer(act,",");

        	if (type.compareToIgnoreCase("CHANGE_BOND")==0) {
        		String site1 = token.nextToken();
        		String order = token.nextToken();
        		String site2 = token.nextToken();

        		Integer s1 = new Integer(site1.substring(1,site1.length()));
        		element = new Integer(order);
        		Integer s2 = new Integer(site2.substring(1,site2.length()));

        		site.add(s1);
        		site.add(s2);
        	}
        	else if (type.compareToIgnoreCase("BREAK_BOND")==0) {
        		String site1 = token.nextToken();
        		String bond = token.nextToken();
        		String site2 = token.nextToken();

        		Integer s1 = new Integer(site1.substring(1,site1.length()));
        		element = Bond.make(bond);
        		Integer s2 = new Integer(site2.substring(1,site2.length()));

        		site.add(s1);
        		site.add(s2);
        	}
        	else if (type.compareToIgnoreCase("FORM_BOND")==0) {
        		String site1 = token.nextToken();
        		String bond = token.nextToken();
        		String site2 = token.nextToken();

        		Integer s1 = new Integer(site1.substring(1,site1.length()));
        		element = Bond.make(bond);
        		Integer s2 = new Integer(site2.substring(1,site2.length()));

        		site.add(s1);
        		site.add(s2);
        	}
        	else if (type.compareToIgnoreCase("GAIN_RADICAL")==0) {
        		String site1 = token.nextToken();
        		String order = token.nextToken();

        		Integer s1 = new Integer(site1.substring(1,site1.length()));
        		element = new Integer(order);

        		site.add(s1);
        	}
        	else if (type.compareToIgnoreCase("LOSE_RADICAL")==0) {
        		String site1 = token.nextToken();
        		String order = token.nextToken();

        		Integer s1 = new Integer(site1.substring(1,site1.length()));
        		element = new Integer(order);

        		site.add(s1);
        	}
        	else {
        		throw new InvalidActionException("Unknown Action type: " + type);
        	}

        	Action thisAction = new Action(type,site,element);
        	action.add(thisAction);

        	line = ChemParser.readUncommentLine(p_reader);
        }
        return action;



        //#]
    }

    //## operation readBond(String)
    public static Object readBond(String p_name) throws InvalidGraphFormatException {
        //#[ operation readBond(String)
        String bondType = removeBrace(p_name);
        StringTokenizer bondToken = new StringTokenizer(bondType, ",");
        int bondNum = bondToken.countTokens();
        try {
        	if (bondNum <= 0) {
        		throw new InvalidGraphFormatException("bond: " + bondType);
        	}
        	else if (bondNum == 1) {
        		String b = bondToken.nextToken();
        		Bond bond = Bond.make(b);
        		return bond;
        	}
        	else {
        		HashSet bondList = new HashSet();
        		while (bondToken.hasMoreTokens()) {
        			String b = bondToken.nextToken();
        			Bond bond = Bond.make(b);
        			bondList.add(bond);
        		}
        		return bondList;
        	}
        }
        catch (UnknownSymbolException e) {
        	throw new InvalidGraphFormatException("bond: " + bondType);
        }




        //#]
    }

    //## operation readChemGraph(BufferedReader)
    public static Graph readChemGraph(BufferedReader p_reader) throws IOException {
        //#[ operation readChemGraph(BufferedReader)
        try {
        	Graph g = new Graph();

        	int centralID = -1;

        	/*
        	 * 14MAR2010: Switched from readUncommentedLine() to readMeaningfulLine()
        	 */
        	//String line = readUncommentLine(p_reader);
        	String line = readMeaningfulLine(p_reader);
			while (line != null) {
         		StringTokenizer token = new StringTokenizer(line);
        		// read in ID
        		String index = token.nextToken();
        		if (index.endsWith(".")) {index = index.substring(0,index.length()-1);}
        		Integer ID = new Integer(index);

                // read in central ID and/or name
                int thisCentralID = -1;
                String centralIndex = token.nextToken();
                String name;
                if (!centralIndex.equals("*")) {
                	name = centralIndex;
                }
                else {
                	if (centralID == -1) centralID = 1;
                	else centralID++;
                	thisCentralID = centralID;
                	name = token.nextToken();
                }

        		// read in the radical number
        		String radical = token.nextToken();
           		Atom atom = (Atom)ChemParser.readChemNodeElement(name,radical);
                Node presentNode = null;

           		if (thisCentralID>0) presentNode = g.addNodeAt(ID.intValue(),atom,thisCentralID);
           		else presentNode = g.addNodeAt(ID.intValue(),atom);

        		// read in the bonds connected to present node site
        		while (token.hasMoreTokens()) {
        			String bondPair = token.nextToken();

        			bondPair = removeBrace(bondPair);

        			StringTokenizer bondPairToken = new StringTokenizer(bondPair,",");

        			Integer nodeID = new Integer(bondPairToken.nextToken());
        			String bondType = bondPairToken.nextToken();

        			Node otherNode = g.getNodeAt(nodeID.intValue());
        			// (1) if the other node is already in the graph, make/add the bond between those two nodes, otherwise, do nothing
        			//		wait until the other node added, and then add the arc.
        			// (2) this requires the adj list should be symetric, say, each bond appears twice. (each node linked by that bond
        			//      should has the bond recorded in the adjlist of the node.
        			if (otherNode != null) {
        				g.addArcBetween(presentNode,Bond.make(bondType),otherNode);
        			}
        		}
        		line = readUncommentLine(p_reader);
        	}
			
        	if (g.isEmpty()) g = null;
        	else g.identifyFgElement();
        	return g;
        }
        catch (Exception e) {
			if (e.getMessage() == null)
				throw new IOException("Couldn't read ChemGraph: " + e.getClass().getName());
			else
				throw new IOException("Couldn't read ChemGraph: " + e.getMessage());
        }





        //#]
    }
    
    public static Graph readAdjList(String adjlist){
    	//Initialize Graph
    	Graph g = new Graph();
    	int centralID = -1;

    	// Separate adjacency list into the individual nodes (lines)
    	String[] eachLine = adjlist.split("[\r]");
    	
    	for (int i=0; i<eachLine.length; i++) {
    		StringTokenizer token = new StringTokenizer(eachLine[i]);
    		// read in ID
    		String index = token.nextToken();
    		if (index.endsWith("."))
    			index = index.substring(0,index.length()-1);
    		Integer ID = new Integer(index);
    		
    		// read in central ID and/or name
    		int thisCentralID = -1;
    		String centralIndex = token.nextToken();
    		String name;
    		if (!centralIndex.equals("*"))
    			name = centralIndex;
    		else {
    			if (centralID == -1) centralID = 1;
    			else centralID++;
    			thisCentralID = centralID;
    			name = token.nextToken();
    		}
    		
    		// read in the radical number
    		String radical = token.nextToken();
    		Atom atom = (Atom)ChemParser.readChemNodeElement(name,radical);
    		Node presentNode = null;
    		
    		if (thisCentralID>0) presentNode = g.addNodeAt(ID.intValue(),atom,thisCentralID);
    		else presentNode = g.addNodeAt(ID.intValue(),atom);

    		// read in the bonds connected to present node site
    		while (token.hasMoreTokens()) {
    			String bondPair = token.nextToken();
    			bondPair = removeBrace(bondPair);
    			
    			StringTokenizer bondPairToken = new StringTokenizer(bondPair,",");
    			
    			Integer nodeID = new Integer(bondPairToken.nextToken());
    			String bondType = bondPairToken.nextToken();
    			
    			Node otherNode = g.getNodeAt(nodeID.intValue());
    			// (1) if the other node is already in the graph, make/add the bond between those two nodes, otherwise, do nothing
    			//		wait until the other node added, and then add the arc.
    			// (2) this requires the adj list should be symetric, say, each bond appears twice. (each node linked by that bond
    			//      should has the bond recorded in the adjlist of the node.
    			if (otherNode != null) {
    				g.addArcBetween(presentNode,Bond.make(bondType),otherNode);
    			}
    		}
    	}
    	
    	if (g.isEmpty()) g = null;
    	else g.identifyFgElement();
    	return g;
    }

    //## operation readChemNodeElement(String,String)
    public static Object readChemNodeElement(String p_name, String p_radical) throws InvalidGraphFormatException, NullSymbolException {
        //#[ operation readChemNodeElement(String,String)
        if (p_name == null || p_radical == null) throw new NullSymbolException();

        // read in radical
        FreeElectron fe = null;
        p_radical = removeBrace(p_radical);
        try {
        	fe = FreeElectron.make(p_radical);
        }
        catch (UnknownSymbolException e) {
        	throw new InvalidGraphFormatException("free electron: " + p_radical);
        }

        // read in atom
        String aList = removeBrace(p_name);

        StringTokenizer aListToken = new StringTokenizer(aList,",");
        int atomNum = aListToken.countTokens();
        try {
        	if (atomNum <= 0) {
        		throw new InvalidGraphFormatException("atom: " + aList);
        	}
        	else if (atomNum == 1) {
        		String nextToken = aListToken.nextToken();
        		ChemNodeElement atom = null;
        		try {
        			ChemElement ce = ChemElement.make(nextToken);
        			atom = Atom.make(ce,fe);
        		}
        		catch (UnknownSymbolException e) {
					String Internalname = FGElement.translateName(nextToken);
        			FGElement fge = FGElement.make(Internalname);
        			atom = FGAtom.make(fge,fe);
        		}
        		return atom;
        	}
        	else {
        		HashSet atomList = new HashSet();
        		ChemNodeElement atom = null;
        		while (aListToken.hasMoreTokens()) {
        			String nextToken = aListToken.nextToken();
        			try {
        				ChemElement ce = ChemElement.make(nextToken);
        				atom = Atom.make(ce,fe);
        			}
        			catch (UnknownSymbolException e) {
						String Internalname = FGElement.translateName(nextToken);
        				FGElement fge = FGElement.make(Internalname);
        				atom = FGAtom.make(fge,fe);
        			}
        			atomList.add(atom);
        		}
        		return atomList;
        	}
        }
        catch (UnknownSymbolException e) {
        	throw new InvalidGraphFormatException("ChemNodeElement: " + aList);
        }



        //#]
    }

    /**
    Requires:
    Effects: read and return one graph from the p_reader
    Modifies: p_reader
    */
    //## operation readFGGraph(BufferedReader)
    public static Graph readFGGraph(BufferedReader p_reader) throws InvalidGraphFormatException, IOException {
        //#[ operation readFGGraph(BufferedReader)
        Graph g = new Graph();

        String line=readUncommentLine(p_reader);

        try {
        	int centralID = -1;
        	while (line != null) {
        		StringTokenizer token = new StringTokenizer(line);
        		// read in ID
        		String index = token.nextToken();
				if (index.equalsIgnoreCase("END")) {
					throw new InvalidGraphFormatException("Please leave an empty line after each chemgraph, including before the 'END' string");
				}
        		if (index.endsWith(".")) {index = index.substring(0,index.length()-1);}
        		Integer ID = new Integer(index);

        		// read in central ID and/or name
        	    int thisCentralID = -1;
        	    String centralIndex = token.nextToken();
        	    String name;
        	    if (!centralIndex.startsWith("*")) {
        	      	name = centralIndex;
        	    }
        	    else {
        	    	if (centralIndex.equals("*")) {
        		    	if (centralID == -1) centralID = 1;
        		      	else centralID++;
        		       	thisCentralID = centralID;
        			}
        	 		else {
        	 			Integer central = new Integer(centralIndex.substring(1,centralIndex.length()));
        	 			thisCentralID = central.intValue();
        	 		}
        	 		name = token.nextToken();
        		}

        		// read in the radical number, note, the number can't be more than one, like: {#1, #2, ...}
        		// throw exception if unknown symbols are caught
        		String radical = token.nextToken();

        	 	// form atom nodes with the information of name, central id, and radical number
        		Object cne = readChemNodeElement(name,radical);
        		Node presentNode = g.addNodeAt(ID.intValue(),cne,thisCentralID);

        		// read in the bonds connected to present node site
        		while (token.hasMoreTokens()) {
        			String bondPair = token.nextToken();
        			bondPair = removeBrace(bondPair);
        			StringTokenizer bondPairToken = new StringTokenizer(bondPair,",");
        	  		Integer nodeID = new Integer(bondPairToken.nextToken());
        			String bondType = bondPairToken.nextToken();
        	        while (bondPairToken.hasMoreTokens()) {
        	        	bondType = bondType + "," + bondPairToken.nextToken();
        	        }

        			Node otherNode = g.getNodeAt(nodeID.intValue());
        			// (1) if the other node is already in the graph, make/add the bond between those two nodes, otherwise, do nothing
        			//		wait until the other node added, and then add the arc.
        			// (2) this requires the adj list should be symetric, say, each bond appears twice. (each node linked by that bond
        			//      should has the bond recorded in the adjlist of the node.
        			if (otherNode != null) {
        				// if there are more than one possible bond, add a set of bonds
        				Object bond = readBond(bondType);
        				g.addArcBetween(presentNode,bond,otherNode);
        			}
        		}
        		line=readUncommentLine(p_reader);
        	}


        	if (g.isEmpty()) g = null;
        	else g.identifyFgElement();
        	return g;
        }
        catch (InvalidGraphFormatException e) {
			System.out.println("Error on line: "+line);
			throw e;
        }




        //#]
    }

    //## operation readFunctionalGroup(String)
    public static Collection readFunctionalGroup(String p_fileName) {
        //#[ operation readFunctionalGroup(String)
        try {
        	FileReader in = new FileReader(p_fileName);

        	BufferedReader reader = new BufferedReader(in);
        	HashSet fgList = new HashSet();

        	String line = readMeaningfulLine(reader);

        	while (line != null) {
        		StringTokenizer st = new StringTokenizer(line);
        		String name = st.nextToken();
        		Graph g = readFGGraph(reader);
        		FunctionalGroup fg = FunctionalGroup.make(name, g);
        		fgList.add(fg);
        	}
        	if (fgList.isEmpty()) fgList = null;
        	return fgList;
        }
        catch (IOException e) {
        	System.err.println("Can't read species!");
        	return null;
        }
        //#]
    }

    /**
    Requires:
    Effects: read hierarchy tree structure
    Modifies:
    */
    //## operation readHierarchyTree(BufferedReader,HashMap,int)
    public static HierarchyTree readHierarchyTree(BufferedReader p_reader, HashMap p_dictionary, int p_level) throws IOException, NotInDictionaryException {
        //#[ operation readHierarchyTree(BufferedReader,HashMap,int)
        try {
        	HierarchyTreeNode root = (HierarchyTreeNode)readHierarchyTreeNode(p_reader,p_level,p_dictionary);
        	if (root == null) return null;
        	else return new HierarchyTree(root);
        }
        catch (IOException e) {
        	throw new IOException(e.getMessage());
        }
        catch (NotInDictionaryException e) {
        	throw new NotInDictionaryException("During Reading Tree: " + e.getMessage());
        }





        //#]
    }

    //## operation readHierarchyTreeNode(BufferedReader,int,HashMap)
    public static Object readHierarchyTreeNode(BufferedReader p_reader, int p_level, HashMap p_dictionary) throws IOException {
        //#[ operation readHierarchyTreeNode(BufferedReader,int,HashMap)
        try {
        	p_reader.mark(1000);
        	String line = p_reader.readLine();
        	read: while (line != null) {
        		
        		line = line.trim();
        		while (line.startsWith("/") || (line.length() == 0)) {
        			p_reader.mark(1000);
        			line = p_reader.readLine();
        			if (line == null) break read;
        			line = line.trim();
        		}
        		
        		StringTokenizer token = new StringTokenizer(line, ":");
        		int data_num = token.countTokens();
        		if (data_num != 2) {
        			throw new IOException("Invalid tree line");
        		}
        		// read in level
        		String l = token.nextToken().trim();
        		String level = l.substring(1,l.length());
        		int present_level = Integer.parseInt(level);
        		if (p_level != present_level) {
        			p_reader.reset();
        			return null;
        		}
        		String other = token.nextToken();
        		token = new StringTokenizer(other);
        		data_num = token.countTokens();
        		if (data_num < 1) {
        			throw new IOException("Invalid tree line");
        		}

        		String name = token.nextToken().trim();
        		Matchable node_element = (Matchable)p_dictionary.get(name);
        		if (node_element == null) {
        			if (name.startsWith("Other")) {
        				DummyLeaf dl = new DummyLeaf(name,p_level);
        				p_dictionary.put(name,name);
        				return dl;
        			}
        			else throw new NotInDictionaryException(name + " is not found in functional group dictionary!");
        		}
        		/*if (name.equals("Cyclopropane"))
        			System.out.println(name);*/
        		HierarchyTreeNode node = new HierarchyTreeNode(node_element,p_level);

        		while (true) {
        			// read the next level nodes (children)
        			
        			Object child = readHierarchyTreeNode(p_reader,present_level + 1, p_dictionary);
        			if (child == null) break;  // can't find children
        		    if (child instanceof HierarchyTreeNode) {
        		    	node.addChildren((HierarchyTreeNode)child);
        		    }
        		    else if (child instanceof DummyLeaf) {
        		    	node.setDummyChild((DummyLeaf)child);
        		    }
        		    else throw new InvalidChildException();
        		}

        		return node;
        	}
        	return null;
        }
        catch (IOException e) {
        	throw new IOException();
        }




        //#]
    }

    /**
    Requires:
    Effects: return the next uncommented and non-empty line.  (skip comment line and empty line).  return null, if it is the end of file
    Modifies: p_reader
    */
    //## operation readMeaningfulLine(BufferedReader)
    public static String readMeaningfulLine(BufferedReader p_reader) {
        //#[ operation readMeaningfulLine(BufferedReader)
        if (p_reader == null) return null;

        String line = null;

        try {
        	do {
        		line = p_reader.readLine();
        		if (line == null) return null;
        		line = line.trim();
        	} while (line.startsWith("//") || line.length() == 0);

        	return line;
        }
        catch (IOException e) {
        	return null;
        }
        //#]
    }

    /**
    Requires:
    Effects: read in species name and adjlist from a species file one by one, return the collection of all the species.  If nothing is read in or if there is any excetption caught, return null.
    Modifies: the species dictionary
    */
    //## operation readSpecies(String)
    public static Collection readSpecies(String p_fileName) {
        //#[ operation readSpecies(String)
        try {
        	FileReader in = new FileReader(p_fileName);

        	BufferedReader reader = new BufferedReader(in);
        	HashSet speciesList = new HashSet();

        	String line = readMeaningfulLine(reader);
        	while (line != null) {
        		StringTokenizer st = new StringTokenizer(line);
        		String name = st.nextToken();
        		Graph g = readChemGraph(reader);
        		if (g == null) throw new NullGraphException();
        		ChemGraph cg = ChemGraph.make(g);
        		Species spe = Species.make(name,cg);
        		speciesList.add(spe);
        		line = readMeaningfulLine(reader);
        	}
        	if (speciesList.isEmpty()) speciesList = null;
        	return speciesList;
        }
        catch (IOException e) {
        	System.err.println("Can't read species!");
        	System.err.println(e.getMessage());
        	return null;
        }
        catch (ForbiddenStructureException e) {
        	System.err.println("Can't read species!");
        	System.err.println("Forbidden Structure:\n" + e.getMessage());
        	return null;
        }
        //#]
    }

    /**
    Requires:
    Effects: return the next uncommented line in the p_reader.  if it is the end of file or an empty line, return null.
    Modifies:
    */
    //## operation readUncommentLine(BufferedReader)
    public static String readUncommentLine(BufferedReader p_reader) {
        //#[ operation readUncommentLine(BufferedReader)
        if (p_reader == null) return null;

        String line = null;

        try {
        	do {
        		line = p_reader.readLine();
        		if (line == null) return null;
        		line = line.trim();
        		if (line.length() == 0) return null;
        	} while (line.startsWith("//"));

        	return line;
        }
        catch (IOException e) {
        	return null;
        }
        //#]
    }

	
    public static HashSet readUnion(String p_string) throws InvalidUnionFormatException {
        try {
        	HashSet result = new HashSet();
        	if (p_string == null) return result;
        	p_string = p_string.trim();
        	String prefix = p_string.substring(0,5);
        	if (p_string.toLowerCase().startsWith("union")) {
        		p_string = p_string.substring(5,p_string.length());
			} 
			else if (p_string.startsWith("OR")) {
				p_string = p_string.substring(2,p_string.length());
			}
			else return result; // which is currently empty
			
        	p_string = ChemParser.removeBrace(p_string);
        	StringTokenizer token = new StringTokenizer(p_string,",");
        	while (token.hasMoreTokens()) {
        		String name = token.nextToken();
        		name = name.trim();
        		result.add(name);
        	}
        	
        	return result;
        }
        catch (Exception e) {
        	throw new InvalidUnionFormatException();
        }
    }

    //## operation removeBrace(String)
    public static String removeBrace(String p_name) {
        //#[ operation removeBrace(String)
        p_name = p_name.trim();
        if (p_name.startsWith("{") || p_name.startsWith("(") || p_name.startsWith("[")) {
        	int trim = 0;
        	if (p_name.endsWith("}") || p_name.endsWith(")") || p_name.endsWith("]")) {
        		trim = 1;
        	}
        	else if (p_name.endsWith("},") || p_name.endsWith("),") || p_name.endsWith("],")) {
        		trim = 2;
        	}
        	p_name = p_name.substring(1,p_name.length()-trim);
        }

        return p_name;



        //#]
    }

    //## operation writeBond(Object)
    public static String writeBond(Object p_bond) {
        //#[ operation writeBond(Object)
        if (p_bond instanceof Bond) {
        	return ((Bond)p_bond).getName();
        }
        else if (p_bond instanceof Collection) {
        	String s = "{";
        	Iterator iter = ((Collection)p_bond).iterator();
        	while (iter.hasNext()) {
        		Bond bond = (Bond)iter.next();
        		s = s + bond.getName() + ",";
        	}
        	s = s.substring(0,s.length()-1);
        	s = s + "}";
        	return s;
        }
        else throw new InvalidBondException();




        //#]
    }

    //## operation writeChemNodeElement(Object)
    public static String writeChemNodeElement(Object p_chemNodeElement) {
        //#[ operation writeChemNodeElement(Object)
        Object o = p_chemNodeElement;

        if (o instanceof Atom) {
        	return ((Atom)o).getChemElement().getName();
        }
        else if (o instanceof FGAtom) {
        	return ((FGAtom)o).getFgElement().getName();
        }
        else if (o instanceof Collection) {
        	String s = "{";
        	Iterator iter = ((Collection)o).iterator();
        	while (iter.hasNext()) {
        		Object thisObject = iter.next();
        		if (thisObject instanceof Atom) {
        			s = s + ((Atom)thisObject).getChemElement().getName() + ",";
        		}
        		else if (thisObject instanceof FGAtom) {
        			s = s + ((FGAtom)thisObject).getFgElement().getName() + ",";
        		}
        	}
        	s = s.substring(0,s.length()-1);
        	s = s + "}";
        	return s;
        }
        else throw new InvalidChemNodeElementException();





        //#]
    }

}
/*********************************************************************
	File Path	: RMG\RMG\jing\chemParser\ChemParser.java
*********************************************************************/

