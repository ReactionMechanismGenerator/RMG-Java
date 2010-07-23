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



package jing.rxn;


import jing.chem.*;
import java.util.*;
import jing.param.*;
import jing.mathTool.*;
import jing.chem.Species;
import jing.chem.ChemGraph;
import jing.chemUtil.Graph;
import jing.param.Temperature;

//## package jing::rxn

//----------------------------------------------------------------------------
// jing\rxn\Structure.java
//----------------------------------------------------------------------------

/**
the basic structure of a reaction:
A + B -> C + D
A, B, C, and D are species.
*/
//## class Structure
public class Structure {

    protected static int MAX_PRODUCT_NUMBER = 3;		//## attribute MAX_PRODUCT_NUMBER

    protected static int MAX_REACTANT_NUMBER = 3;		//## attribute MAX_REACTANT_NUMBER

    protected int direction;		//## attribute direction

    protected int redundancy = 1;		//## attribute redundancy

    protected LinkedList products;
    protected LinkedList reactants;

    // Constructors

    //## operation Structure(LinkedList,LinkedList)
    public  Structure(LinkedList p_reactants, LinkedList p_products) {
        {
            reactants=new LinkedList();
        }
        {
            products=new LinkedList();
        }
        //#[ operation Structure(LinkedList,LinkedList)
        reactants = p_reactants;
        products = p_products;



        //#]
    }

//	## operation Structure(LinkedList,LinkedList)
    public  Structure(LinkedList p_reactants, LinkedList p_products, int p_direction) {
        {
            reactants=new LinkedList();
        }
        {
            products=new LinkedList();
        }
        //#[ operation Structure(LinkedList,LinkedList)
        reactants = p_reactants;
        products = p_products;
        direction = p_direction;


        //#]
    }

    //## operation Structure(LinkedList)
    public  Structure(LinkedList p_reactants) {
        {
            reactants=new LinkedList();
        }
        {
            products=new LinkedList();
        }
        //#[ operation Structure(LinkedList)
        reactants = p_reactants;
        products = null;
        //#]
    }
    public  Structure() {
        {
            products=new LinkedList();
        }
        {
            reactants=new LinkedList();
        }
    }

    //## operation calculateHrxn(Temperature)
    public double calculateHrxn(Temperature p_temperature) {
        //#[ operation calculateHrxn(Temperature)
        double DH = 0;

        for (Iterator iter = getReactants(); iter.hasNext(); ) {
        	Object o = iter.next();
			Species spe = null;
			if (o instanceof Species)
				spe = (Species)o;
			else {
				ChemGraph cg = (ChemGraph)o;
				spe = cg.getSpecies();
			}
				
        	DH -= spe.calculateH(p_temperature);
        }

        for (Iterator iter = getProducts(); iter.hasNext(); ) {
        	Object o = iter.next();
			Species spe = null;
			if (o instanceof Species)
				spe = (Species)o;
			else {
				ChemGraph cg = (ChemGraph)o;
	        	spe = cg.getSpecies();
			}
        	DH += spe.calculateH(p_temperature);
        }

        return DH;
        //#]
    }

    //## operation calculateKeq(Temperature)
    public double calculateKeq(Temperature p_temperature) {
        //#[ operation calculateKeq(Temperature)
        double DG = 0;
        for (Iterator iter = getReactants(); iter.hasNext(); ) {
			Object o = iter.next();
			Species spe = null;
			if (o instanceof Species)
				spe = (Species)o;
			else {
				ChemGraph cg = (ChemGraph)o;
				spe = cg.getSpecies();
			}
        	DG -= spe.calculateG(p_temperature);
        }
        for (Iterator iter = getProducts(); iter.hasNext(); ) {
			Object o = iter.next();
			Species spe = null;
			if (o instanceof Species)
				spe = (Species)o;
			else {
				ChemGraph cg = (ChemGraph)o;
				spe = cg.getSpecies();
			}
        	DG += spe.calculateG(p_temperature);
        }
        double Keq;
        double T = p_temperature.getStandard();
        double R = GasConstant.getKcalMolK();
        Keq = Math.exp(-DG/R/T);
        int deltaN = getDeltaN();
        Keq *= Math.pow(GasConstant.getCCAtmMolK()*p_temperature.getK(),-deltaN);
        return Keq;
        //#]
    }

    //## operation calculateKeqLowerBound(Temperature)
      public double calculateKeqLowerBound(Temperature p_temperature) {//svp
        //#[ operation calculateKeqLowerBound(Temperature)
        double DG = 0;
        for (Iterator iter = getReactants();iter.hasNext();) {
			Object o = iter.next();
			Species spe = null;
			if (o instanceof Species)
				spe = (Species)o;
			else {
				ChemGraph cg = (ChemGraph)o;
				spe = cg.getSpecies();
			}
          DG -= spe.calculateGUpperBound(p_temperature);
        }
        for (Iterator iter = getProducts(); iter.hasNext();) {
			Object o = iter.next();
			Species spe = null;
			if (o instanceof Species)
				spe = (Species)o;
			else {
				ChemGraph cg = (ChemGraph)o;
				spe = cg.getSpecies();
			}
          DG += spe.calculateGLowerBound(p_temperature);
        }
        double Keq;
        double T = p_temperature.getStandard();
        double R = GasConstant.getKcalMolK();
        Keq = Math.exp(-DG/R/T);
        int deltaN = getDeltaN();
        Keq *= Math.pow(GasConstant.getCCAtmMolK()*p_temperature.getK(),-deltaN);
        return Keq;
        //#]
      }

      //## operation calculateKeqUpperBound(Temperature)
        public double calculateKeqUpperBound(Temperature p_temperature) {//svp
          //#[ operation calculateKeqUpperBound(Temperature)
          double DG = 0;
          for (Iterator iter = getReactants();iter.hasNext();) {
			  Object o = iter.next();
				Species spe = null;
				if (o instanceof Species)
					spe = (Species)o;
				else {
					ChemGraph cg = (ChemGraph)o;
					spe = cg.getSpecies();
				}
            DG -= spe.calculateGUpperBound(p_temperature);
          }
          for (Iterator iter = getProducts(); iter.hasNext();) {
			  Object o = iter.next();
				Species spe = null;
				if (o instanceof Species)
					spe = (Species)o;
				else {
					ChemGraph cg = (ChemGraph)o;
					spe = cg.getSpecies();
				}
            DG += spe.calculateGLowerBound(p_temperature);
          }
          double Keq;
          double T = p_temperature.getStandard();
          double R = GasConstant.getKcalMolK();
          Keq = Math.exp(-DG/R/T);
          int deltaN = getDeltaN();
          Keq *= Math.pow(GasConstant.getCCAtmMolK()*p_temperature.getK(),-deltaN);
          return Keq;
          //#]
        }


    //## operation calculateSrxn(Temperature)
    public double calculateSrxn(Temperature p_temperature) {
        //#[ operation calculateSrxn(Temperature)
        double DS = 0;

        for (Iterator iter = getReactants(); iter.hasNext(); ) {
			Object o = iter.next();
			Species spe = null;
			if (o instanceof Species)
				spe = (Species)o;
			else {
				ChemGraph cg = (ChemGraph)o;
				spe = cg.getSpecies();
			}
        	DS -= spe.calculateS(p_temperature);
        }

        for (Iterator iter = getProducts(); iter.hasNext(); ) {
			Object o = iter.next();
			Species spe = null;
			if (o instanceof Species)
				spe = (Species)o;
			else {
				ChemGraph cg = (ChemGraph)o;
				spe = cg.getSpecies();
			}
        	DS += spe.calculateS(p_temperature);
        }

        return DS;
        //#]
    }

    //## operation clearChemGraph()
    /*public void clearChemGraph() {
        //#[ operation clearChemGraph()
        for (Iterator iter = getReactants(); iter.hasNext(); ) {
        	ChemGraph cg = (ChemGraph)iter.next();
        	Species s = cg.getSpecies();
        	ChemGraph newCG = null;
        	if (s.hasResonanceIsomers()) {
        		for (Iterator riIter = s.getResonanceIsomers(); riIter.hasNext(); ) {
        			ChemGraph ri = (ChemGraph)riIter.next();
        			if (cg.equals(ri)) {
        				newCG = ri;
        				break;
        			}
        		}
        	}
        	else {
        		newCG = s.getChemGraph();
        	}
        	if (newCG == null) {
        		String error = "ChemGraph's Species doesn't include this ChemGraph structure!\n";
        		error = error + "Species info:\n " + s.toString();
        		error = error + "ChemGraph info:\n" + cg.toString();
        		throw new InvalidChemGraphException(error);
        	}
        	if (cg != newCG) {
        		reactants.set(reactants.indexOf(cg),newCG);
        		cg = null;
        	}
        }

        for (Iterator iter = getProducts(); iter.hasNext(); ) {
        	ChemGraph cg = (ChemGraph)iter.next();
        	Species s = cg.getSpecies();
        	ChemGraph newCG = null;
        	if (s.hasResonanceIsomers()) {
        		for (Iterator riIter = s.getResonanceIsomers(); riIter.hasNext(); ) {
        			ChemGraph ri = (ChemGraph)riIter.next();
        			if (cg.equals(ri)) {
        				newCG = ri;
        				break;
        			}
        		}
        	}
        	else {
        		newCG = s.getChemGraph();
        	}
        	if (newCG == null) {
        		String error = "ChemGraph's Species doesn't include this ChemGraph structure!\n";
        		error = error + "Species info:\n " + s.toString();
        		error = error + "ChemGraph info:\n" + cg.toString();
        		throw new InvalidChemGraphException(error);
        	}
        	if (newCG != cg) {
        		products.set(products.indexOf(cg),newCG);
        		cg = null;
        	}
        }

        /*Iterator iter = getReactants();
        while (iter.hasNext()) {
        	ChemGraph cg = (ChemGraph)iter.next();
        	ChemGraph newCg = ChemGraph.getChemGraphFromDictionary(cg);
        	if (newCg == null) ChemGraph.putChemGraphInDictionary(cg);
        	else {
        		reactants.set(reactants.indexOf(cg),newCg);
        		cg = null;
        	}
        }

        iter = getProducts();
        while (iter.hasNext()) {
        	ChemGraph cg = (ChemGraph)iter.next();
        	ChemGraph newCg = ChemGraph.getChemGraphFromDictionary(cg);
        	if (newCg == null) ChemGraph.putChemGraphInDictionary(cg);
        	else {
        		products.set(products.indexOf(cg),newCg);
        		cg = null;
        	}
        }
        return;


        //#]
    }*/

    /**
    Requires:
    Effects: check if this sturcture contains the pass-in species.  if it contains the species, either as reactants, or as products, return true; otherwise, return false.
    Modifies:
    */
    //## operation contains(Species)
    public boolean contains(Species p_species) {
        //#[ operation contains(Species)
        if (containsAsReactants(p_species) || containsAsProducts(p_species)) return true;
        return false;



        //#]
    }

    /**
    Requires:
    Effects: check if this structure contains the species as products.  if it does, return true; otherwise, return false
    Modifies:
    */
    //## operation containsAsProducts(Species)
    public boolean containsAsProducts(Species p_species) {
        //#[ operation containsAsProducts(Species)
        if (products.contains(p_species)) return true;
        else return false;
        //#]
    }

    /**
    Requires:
    Effects: check if this structure contains the species as reactants.  if it does, return true; otherwise, return false
    Modifies:
    */
    //## operation containsAsReactants(Species)
    public boolean containsAsReactants(Species p_species) {
        //#[ operation containsAsReactants(Species)
        if (reactants.contains(p_species)) return true;
        else return false;



        //#]
    }

    /**
    Requires:
    Effects: compare if two structures are the same.
    the order doesn't matter, i.e., A+B->C+D is equal to B+A->D+C, etc.
    Modifies:
    */
    //## operation equals(Object)
	public boolean equalsAsChemGraph(Object p_structure) {
        //check if p_structure is a structure
		if (!(p_structure instanceof Structure)) return false;
		
		//check if p_structure is equal to my structure
		if (p_structure == this) return true;
		
		Structure structure = (Structure)p_structure;
		
		//check if each of the reactant of this is present in the reactant of structure.
		Iterator iter = reactants.iterator();
		LinkedList p_structureReactants = (LinkedList)structure.reactants.clone();
		
		while (iter.hasNext()){
			ChemGraph cg = (ChemGraph)iter.next();
			if (cgPresentInList(cg,p_structureReactants))
				continue;
			else
				return false;
		}
		
		
		//check if each of the productt of this is present in the reactant of structure.
		iter = products.iterator();
		LinkedList p_structureProducts = (LinkedList)structure.products.clone();
		while (iter.hasNext()){
			ChemGraph cg = (ChemGraph)iter.next();
			if (cgPresentInList(cg,p_structureProducts))
				continue;
			else
				return false;
		}
		
		return true;
        //#]
    }
	
	public static boolean cgPresentInList( ChemGraph cg, LinkedList l){
		if (!cg.getSpecies().hasResonanceIsomers()) {
			boolean present = false;
			Iterator iter = l.iterator();
			while (iter.hasNext()){
				ChemGraph cg2 = (ChemGraph)iter.next();
				if (cg2.getSpecies().equals(cg.getSpecies())){
					l.remove(cg2);
					return true;
				}
			}
			return false;
		}
		else {
			boolean present = false;
			Iterator iter = l.iterator();
			while (iter.hasNext()){
				ChemGraph cg2 = (ChemGraph)iter.next();
				if (cg2.getSpecies().hasResonanceIsomers()){
					if (cg2.equals(cg)){
						l.remove(cg2);
						return true;
					}
				}
			}
			return false;
		}
	}
	
    public boolean equals(Object p_structure) {
        //#[ operation equals(Object)
        if (this == p_structure) return true;

        if (!(p_structure instanceof Structure)) return false;

        Structure structure = (Structure)p_structure;

		 boolean requal = isCGListEquivalentAsSpecies(reactants, structure.reactants);
	        
			if (requal){
				return isCGListEquivalentAsSpecies(products, structure.products);
			}
			else 
				return false;
		
        // compare rule: if the two sets have the same size, and one of them contains all the elements in the other set,
        // these two sets are the same.
        /*if (equalsAsSpecies(structure)){
			boolean requal = MathTool.isListEquivalent(reactants, structure.reactants);
	        boolean pequal = MathTool.isListEquivalent(products, structure.products);

			//boolean requal = MathTool.isListEqual(reactants, structure.reactants);
	        //boolean pequal = MathTool.isListEqual(products, structure.products);

	        return requal&&pequal;
        }
        else 
			return false;
		*/
        //#]
    }

	
    //## operation equals(Object)
    /*public boolean isDuplicate(Object p_structure) {
        //#[ operation equals(Object)
        if (this == p_structure) return true;

        if (!(p_structure instanceof Structure)) return false;

        Structure structure = (Structure)p_structure;

        // compare rule: if the two sets have the same size, and one of them contains all the elements in the other set,
        // these two sets are the same.
		 //if (equalsAsSpecies((Structure)p_structure)){
				//boolean requal = MathTool.isListEquivalent(reactants, structure.reactants);
		        //boolean pequal = MathTool.isListEquivalent(products, structure.products);

				boolean requal = MathTool.isListEqual(reactants, structure.reactants);
		        boolean pequal = MathTool.isListEqual(products, structure.products);

		        return requal&&pequal;
        //#]
    }*/
	
	   //## operation equalsAsSpecies(Structure)
    /*public boolean equalsAsSpecies(Structure p_structure) {
        //#[ operation equalsAsSpecies(Structure)
        if (this == p_structure) return true;

        // compare rule: if the two sets have the same size, and one of them contains all the elements in the other set,
        // these two sets are the same.
        boolean requal = isCGListEquivalentAsSpecies(reactants, p_structure.reactants);
        
		if (requal){
			return isCGListEquivalentAsSpecies(products, p_structure.products);
		}
		else 
			return false;
        //return requal&&pequal;
        //#]
    }*/

    //## operation generateReverseStructure()
    public Structure generateReverseStructure() {
        //#[ operation generateReverseStructure()
        Structure newS = new Structure(getProductList(), getReactantList());
        newS.setDirection(-getDirection());

        return newS;


        //#]
    }

	 //## operation generateSpeciesStructure()
    /*public Structure generateSpeciesStructure() {
        //#[ operation generateSpeciesStructure()
        LinkedList r = new LinkedList();
        for (Iterator iter = getReactants(); iter.hasNext(); ) {
        	ChemGraph reactant = (ChemGraph)iter.next();
        	Species spe = reactant.getSpecies();
        	r.add(spe);
        }

        LinkedList p = new LinkedList();
        for (Iterator iter = getProducts(); iter.hasNext(); ) {
        	ChemGraph product = (ChemGraph)iter.next();
        	Species spe = product.getSpecies();
        	p.add(spe);
        }

        Structure s = new Structure(r,p);
        s.setDirection(getDirection());
        s.setRedundancy(getRedundancy());

        return s;
        //#]
    }*/

    //## operation getDeltaN()
    public int getDeltaN() {
        //#[ operation getDeltaN()
        return getProductNumber()-getReactantNumber();
        //#]
    }

    //## operation getProductList()
    public LinkedList getProductList() {
        //#[ operation getProductList()
        return products;
        //#]
    }

    /**
    Requires:
    Effects: return the number of species in the itsProducts linkedlist.
    Modifies:
    */
    //## operation getProductNumber()
    public int getProductNumber() {
        //#[ operation getProductNumber()
        return products.size();
        //#]
    }

    //## operation getReactantList()
    public LinkedList getReactantList() {
        //#[ operation getReactantList()
        return reactants;
        //#]
    }

    /**
    Requires:
    Effects: return the number of speices in the itsReactants linkedlist.
    Modifies:
    */
    //## operation getReactantNumber()
    public int getReactantNumber() {
        //#[ operation getReactantNumber()
        return reactants.size();
        //#]
    }

    //## operation hashCode()
    public int hashCode() {
//		#[ operation hashCode()
        int hash = 0;
        Iterator r = getReactants();
        while (r.hasNext()) {
        	Object reactant = r.next();
        	hash += reactant.hashCode();
        }
        Iterator p = getProducts();
        while (p.hasNext()) {
        	Object product = p.next();
        	hash += product.hashCode();
        }
        return hash;
        //#]
    }

    //## operation increaseRedundancy(int)
    public void increaseRedundancy(int p_number) {
        //#[ operation increaseRedundancy(int)
        redundancy += p_number;


        //#]
    }

	   //## operation isCGListEquivalentAsSpecies(LinkedList,LinkedList)
    public static boolean isCGListEquivalentAsSpecies(LinkedList p_list1, LinkedList p_list2) {
        //#[ operation isCGListEquivalentAsSpecies(LinkedList,LinkedList)
        if (p_list1.size() != p_list2.size()) return false;

        boolean found = false;

        LinkedList templist = (LinkedList)(p_list2.clone());
        for (Iterator iter1 = p_list1.iterator(); iter1.hasNext();) {
			Species s1 = (Species)iter1.next();
        	found = false;
        	for (Iterator iter2 = templist.iterator(); iter2.hasNext();) {
				Species s2 = (Species)iter2.next();
        		if (s1==s2) {
        			found = true;
        			iter2.remove();
        			break;
        		}
        	}
        	if (!found) return false;
        }

        if (templist.isEmpty()) return true;
        else return false;
        //#]
    }
    
    /**
     *	isSpeciesListEquivalentToChemGraphListAsChemGraphs.java
     *		Compares a list of species to a list of chemgraphs and determines
     *		if the lists are equivalent.
     *
     *		The first LinkedList passed to the function is a list of species;
     *		the second is a list of chemgraphs
     *
     *		This function is primarily used to determine if a just formed
     *		rxn is equivalent to any of the rxns stored in the primary
     *		kinetic libraries.
     **/
    public static boolean isSpeciesListEquivalentToChemGraphListAsChemGraphs(LinkedList p_list1, LinkedList p_list2) {
        if (p_list1.size() != p_list2.size()) return false;

        LinkedList templist = (LinkedList)(p_list2.clone());
        
        // For each element in the FIRST list
        for (Iterator iter1 = p_list1.iterator(); iter1.hasNext();) {
            boolean matchFound = false;
        	// Grab the element
			Species sp1 = (Species)iter1.next();	// This list contains species
			
			/*
			 * If the species has resonance isomers, iterate over all
			 * 	chemgraphs and check if one of them matches the just-formed
			 * 	chemgraph (this chemgraph will be equal to the "mostStable"
			 * 	chemgraph, as computed by RMG.
			 */
			if (sp1.hasResonanceIsomers()) {
				for (Iterator iter2 = templist.iterator(); iter2.hasNext();) {
					ChemGraph cg2 = (ChemGraph)iter2.next();	// This list contains chemgraphs
					Iterator resonanceIsomerIter = sp1.getResonanceIsomers();
					while (resonanceIsomerIter.hasNext()) {
						ChemGraph cg1 = (ChemGraph)resonanceIsomerIter.next();
							if (cg1.equals(cg2)) {
								// If the two chemgraphs are equivalent, remove the current
								//	chemgraph from the Iterator iter2 so we no longer
								//	try to match that chemgraph
								matchFound = true;
								iter2.remove();
								break;	// while loop
							}
						}
					if (matchFound) break; // for loop
				}
				// If we could not match the species' chemgraph against
				//	any of the chemgraphs' graph, the list is not equivalent
				if (!matchFound) return false;
			} 
			
			else {
				ChemGraph cg1 = sp1.getChemGraph();
	        	for (Iterator iter2 = templist.iterator(); iter2.hasNext();) {
					ChemGraph cg2 = (ChemGraph)iter2.next();	// This list contains chemgraphs
	        		if (cg1.equals(cg2)) {
	        			// If the two chemgraphs are equivalent, remove the current
	        			//	chemgraph from the Iterator iter2 so we no longer
	        			//	try to match that chemgraph
	        			matchFound = true;
	        			iter2.remove();
	        			break;
	        		}
	        	}
	        	// If we could not match the species' chemgraph against
	        	//	any of the chemgraphs' graph, the list is not equivalent
	        	if (!matchFound) return false;
			}
        }

        // If we found a match for all entries in the list, the
        //	templist will be empty.
        if (templist.isEmpty()) return true;
        else return false;
        //#]
    }

    //## operation isBackward()
    public boolean isBackward() {
        //#[ operation isBackward()
        return (direction == -1);
        //#]
    }

    //## operation isForward()
    public boolean isForward() {
        //#[ operation isForward()
        return (direction == 1);
        //#]
    }

    //## operation reactantEqualsProduct()
    public boolean reactantEqualsProduct() {
        //#[ operation reactantEqualsProduct()
        return MathTool.isListEquivalent(reactants, products);
        //#]
    }

    /**
    Requires:
    Effects: check if this structure is valid, which means if the structure is chemically balanced
    (1) atom balance
    (2) valance banlance
    Modifies:
    */
    //## operation repOk()
    public boolean repOk() {
        //#[ operation repOk()
        // check the number of reactants and products
        int reactant_num = getReactantNumber();
        int product_num = getProductNumber();
        if (reactant_num > MAX_REACTANT_NUMBER || reactant_num < 0) {
        	throw new InvalidReactantNumberException();
        }
        if (product_num > MAX_REACTANT_NUMBER || product_num < 0) {
        	throw new InvalidProductNumberException();
        }

        // add other checkings
        /*
         * MRH 22JUL2010: Checking whether reaction balances
         */
        int numC=0; int numH=0; int numO=0; int numS=0; int numSi=0;
        LinkedList reactants = getReactantList();
        for (int i=0; i<reactants.size(); i++) {
        	ChemGraph cg = ((Species)reactants.get(i)).getChemGraph();
        	numC += cg.getCarbonNumber();
        	numH += cg.getHydrogenNumber();
        	numO += cg.getOxygenNumber();
        	numS += cg.getSulfurNumber();
        	numSi += cg.getSiliconNumber();
        }
        LinkedList products = getProductList();
        for (int j=0; j<products.size(); j++) {
        	ChemGraph cg = ((Species)products.get(j)).getChemGraph();
        	numC -= cg.getCarbonNumber();
        	numH -= cg.getHydrogenNumber();
        	numO -= cg.getOxygenNumber();
        	numS -= cg.getSulfurNumber();
        	numSi -= cg.getSiliconNumber();
        }
        if (numC!=0 || numH!=0 || numO!=0 || numS!=0 || numSi!=0) {
        	System.out.println("Reaction is not balanced: " + toString());
        	return false;
        }
        
        return true;

        //#]
    }

    //## operation toChemkinString()
	 public StringBuilder toChemkinString(boolean p_includeReverse) {
	        //#[ operation toChemkinString(boolean)
	        StringBuilder s = new StringBuilder();
	        int index = 0;
	        int r_num = getReactantNumber();
	        int p_num = getProductNumber();

	        Iterator iter1 = getReactants();
	        while (iter1.hasNext()) {
	        	index++;
				Object o = iter1.next();
				Species spe ;
				if (o instanceof Species){
					spe = (Species)o;
				}
				else {
					spe = ((ChemGraph)o).getSpecies();
				}
	        	if (index < r_num) s = s.append(spe.toChemkinString() + "+");
	        	else s.append(spe.toChemkinString());
	        }
	        if (p_includeReverse)
	        	s.append("=");
	        else
	        	s.append("=>");

	        index = 0;
	        Iterator iter2 = getProducts();
	        while (iter2.hasNext()) {
	        	index++;
	        	Object o = iter2.next();
				Species spe ;
				if (o instanceof Species){
					spe = (Species)o;
				}
				else {
					spe = ((ChemGraph)o).getSpecies();
				}
	        	if (index < p_num) s.append(spe.toChemkinString() + "+");
	        	else s.append(spe.toChemkinString());
	        }

	        return s;
	        //#]
	    }

//	## operation toChemkinString()
//	## operation toChemkinString()
	 public StringBuilder toRestartString(boolean p_includeReverse) {
	        //#[ operation toChemkinString(boolean)
	        String s="";
	        int index = 0;
	        int r_num = getReactantNumber();
	        int p_num = getProductNumber();

	        Iterator iter1 = getReactants();
	        while (iter1.hasNext()) {
	        	index++;
				Species spe = (Species)iter1.next();
	        	if (index < r_num) s = s + spe.getName()+"("+spe.getID()+")" + "+";
	        	else s += spe.getName()+"("+spe.getID()+")";
	        }
	        if (p_includeReverse)
	        	s += "=";
	        else
	        	s += "=>";

	        index = 0;
	        Iterator iter2 = getProducts();
	        while (iter2.hasNext()) {
	        	index++;
				Species spe = (Species)iter2.next();
	        	if (index < p_num) s = s + spe.getName()+"("+spe.getID()+")" + "+";
	        	else s += spe.getName()+"("+spe.getID()+")";
	        }

	        StringBuilder sb = new StringBuilder(s);
	        
	        return sb;
	        //#]
	    }


    //## operation toString()
    /*public String toString() {
        //#[ operation toString()
        String s="";
        int index = 0;
        int r_num = getReactantNumber();
        int p_num = getProductNumber();

        Iterator iter1 = getReactants();
        while (iter1.hasNext()) {
        	index++;
        	ChemGraph cg = (ChemGraph)iter1.next();
        	int id = cg.getSpecies().getID();
        	if (index < r_num) {
        		s = s + cg.getSpecies().getName() + "(" + String.valueOf(id) + ")" + " + ";
        	}
        	else {
        		s = s + cg.getSpecies().getName() + "(" + String.valueOf(id) + ") -> ";
        	}
        }

        index = 0;
        Iterator iter2 = getProducts();
        while (iter2.hasNext()) {
        	index++;
        	ChemGraph cg = (ChemGraph)iter2.next();
        	int id = cg.getSpecies().getID();
        	if (index < p_num) {
        		s = s + cg.getSpecies().getName() + "(" + String.valueOf(id) + ")" + " + ";
        	}
        	else s = s + cg.getSpecies().getName() + "(" + String.valueOf(id) + ")" + '\t';
        }

        s = s + "Direction = " + String.valueOf(direction) + '\t';
        s = s + "Redundancy = " + String.valueOf(redundancy);

        return s;
        //#]
    }*/

    public static int getMAX_PRODUCT_NUMBER() {
        return MAX_PRODUCT_NUMBER;
    }

    public static int getMAX_REACTANT_NUMBER() {
        return MAX_REACTANT_NUMBER;
    }

    public int getDirection() {
        return direction;
    }

    public void setDirection(int p_direction) {
        direction = p_direction;
    }

    public int getRedundancy() {
        return redundancy;
    }

    public void setRedundancy(int p_redundancy) {
        redundancy = p_redundancy;
    }

    public ListIterator getProducts() {
        ListIterator iter=products.listIterator(0);
        return iter;
    }

    public void addProducts(ChemGraph p_ChemGraph) {
        products.add(p_ChemGraph);
    }

    public void removeProducts(ChemGraph p_ChemGraph) {
        products.remove(p_ChemGraph);
    }

    public void clearProducts() {
        products.clear();
    }
	
	public String toString() {
		return toChemkinString(true).toString();
	}

    public ListIterator getReactants() {
        ListIterator iter=reactants.listIterator(0);
        return iter;
    }

    public void addReactants(ChemGraph p_ChemGraph) {
        reactants.add(p_ChemGraph);
    }

    public void removeReactants(ChemGraph p_ChemGraph) {
        reactants.remove(p_ChemGraph);
    }

    public void clearReactants() {
        reactants.clear();
    }

}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\Structure.java
*********************************************************************/

