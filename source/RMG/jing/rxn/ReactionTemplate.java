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


import java.io.*;
import jing.chem.*;
import java.util.*;

import jing.param.*;
import jing.rxnSys.PrimaryKineticLibrary;
import jing.rxnSys.ReactionModelGenerator;
import jing.mathTool.*;
import jing.chemUtil.*;
import jing.chemParser.*;
import jing.chem.Species;
import jing.chem.ChemGraph;

//## package jing::rxn 

//----------------------------------------------------------------------------
//jing\rxn\ReactionTemplate.java                                                                  
//----------------------------------------------------------------------------

/**
This is the reaction template that generates a family of reactions. 
ItsStructureTemplates could be a bunch of single Structure templates, although for many cases, it will be only one for each reaction template.
reaction templates are mutable.
*/
//## class ReactionTemplate 
public class ReactionTemplate {
  
  /**
  Reaction direction:
  1: forward reaction
  2: backward reaction
  0: undertemined 
  */
  protected int direction = 0;		//## attribute direction 
  
  protected LinkedHashMap fgDictionary = new LinkedHashMap();		//## attribute fgDictionary 
  
  protected String name = "";		//## attribute name 
  
  //protected HashMap reactionDictionaryByReactant = new HashMap();		//## attribute reactionDictionaryByReactant 
  
  protected LinkedHashMap reactionDictionaryByStructure = new LinkedHashMap();		//## attribute reactionDictionaryByStructure 
  
  protected KineticsTemplateLibrary kineticsTemplateLibrary;
  protected ReactionAdjList reactionAdjList;
  protected ReactionTemplate reverseReactionTemplate;
  protected StructureTemplate structureTemplate;
  protected LinkedHashMap forbiddenStructures = new LinkedHashMap();
  
  // Constructors
  
  /**
  default constructor
  */
  //## operation ReactionTemplate() 
  public  ReactionTemplate() {
      //#[ operation ReactionTemplate() 
      //#]
  }
  
	/*public int getNumberOfReactions(){
		return reactionDictionaryByStructure.size();
	}*/
  //## operation ableToGeneratePDepWellPaths() 
  public boolean ableToGeneratePDepWellPaths() {
      //#[ operation ableToGeneratePDepWellPaths() 
      return hasOneReactant();
      //#]
  }
  
  //## operation addReaction(TemplateReaction) 
  public void addReaction(TemplateReaction p_templateReaction) {
      //#[ operation addReaction(TemplateReaction) 
      if (p_templateReaction == null) return;
      
      if (p_templateReaction.getReactionTemplate() != this) return;
      
      Structure s = p_templateReaction.getStructure();
      reactionDictionaryByStructure.put(s,p_templateReaction);
      
      return;     
      //#]
  }

  public void removeFromReactionDictionaryByStructure(Structure s) {
      if(reactionDictionaryByStructure.containsKey(s)) reactionDictionaryByStructure.remove(s);
      return;
  }
  
  //## operation calculateDepth(HashSet) 
  private static final int calculateDepth(HashSet p_treeNodeSet) {
      //#[ operation calculateDepth(HashSet) 
      if (p_treeNodeSet==null) throw new NullPointerException();
      
      int dep = 0;
      for (Iterator iter = p_treeNodeSet.iterator(); iter.hasNext();) {
      	Object n = iter.next();
      	if (n instanceof DummyLeaf) {
      		dep += Math.abs(((DummyLeaf)n).getDepth());
      	}
      	else if (n instanceof HierarchyTreeNode) {
      		dep += Math.abs(((HierarchyTreeNode)n).getDepth());
      	}
      	else throw new InvalidKineticsKeyException();
      }
      return dep;
      //#]
  }
  
  //## operation calculateDistance(HashSet,HashSet) 
  private static final int calculateDistance(HashSet p_treeNodeSet1, HashSet p_treeNodeSet2) {
      //#[ operation calculateDistance(HashSet,HashSet) 
      return Math.abs(calculateDepth(p_treeNodeSet1)-calculateDepth(p_treeNodeSet2));
      //#]
  }
  
  //## operation extractKineticsTemplateKey(HashSet) 
  private static final LinkedHashSet extractKineticsTemplateKey(LinkedHashSet p_treeNode) {
      //#[ operation extractKineticsTemplateKey(HashSet) 
      LinkedHashSet key = new LinkedHashSet();
      for (Iterator iter = p_treeNode.iterator(); iter.hasNext(); ) {
      	Object n = iter.next();
      	if (n instanceof DummyLeaf) {
      		String name = ((DummyLeaf)n).getName();
      		key.add(name);
      	}
      	else if (n instanceof HierarchyTreeNode) {
      		Matchable fg = (Matchable)((HierarchyTreeNode)n).getElement();
      		key.add(fg);
      	}
      	else {
      		throw new InvalidKineticsKeyException();
      	}
      }
      
      return key;
      //#]
  }
  
  //## operation fillKineticsBottomToTop() 
  public void fillKineticsBottomToTop() {
      //#[ operation fillKineticsBottomToTop() 
      LinkedHashSet rootSet = new LinkedHashSet();//6/26/09 gmagoon: changed from HashSet to LinkedHashSet to make behavior deterministic
      
      Iterator tree_iter = getReactantTree();
      while (tree_iter.hasNext()) {
      	HierarchyTree tree = (HierarchyTree)tree_iter.next();
      	rootSet.add(tree.getRoot());
      }
      fillKineticsByAverage(rootSet);
      System.gc();
      //#]
  }
  
  /**
  Requires:
  Effects:
  Modifies:
  */
  //## operation fillKineticsByAverage(HashSet) 
  private Kinetics fillKineticsByAverage(LinkedHashSet p_treeNodeSet) {
      //#[ operation fillKineticsByAverage(HashSet) 
      // check if k is already in the library
      // if it is, don't do any average;
      // if it isn't, average all the k below this tree node set to get an estimated k for this tree node set;
      
      LinkedHashSet key = extractKineticsTemplateKey(p_treeNodeSet);
      Kinetics k = kineticsTemplateLibrary.getKinetics(key);
      
      if (k==null || (k!=null && k.getRank()==0)) {
      	boolean allLeaf=true;
      	LinkedList fgc = new LinkedList();
      	Iterator iter = p_treeNodeSet.iterator();
      	while (iter.hasNext()) {
      		Object node = iter.next();
      		Stack path = new Stack();
              if (node instanceof DummyLeaf) {
              	path.push(node);
              }
              else if (node instanceof HierarchyTreeNode) {
      		  	HierarchyTreeNode n = (HierarchyTreeNode)node;
      			if (n.isLeaf()) {
               		path.push(n);
      			}
      			else {
      				allLeaf = false;
      				if (n.hasDummyChild()) path.push(n.getDummyChild());
      				for (Iterator child_iter = n.getChildren(); child_iter.hasNext(); ) {
      					HierarchyTreeNode child = (HierarchyTreeNode)child_iter.next();
      					path.push(child);
      				}
      			}
              }
              else {
              	throw new InvalidKineticsKeyException();
              }
      
      		fgc.add(path);
      	}
      	if (allLeaf) return null;
      
      	Collection allPossibleTreeNodeKeySet = MathTool.expand(fgc.iterator());
      
      	// this is the set that we could find for all the combinatorial generated from the p_fgc
      	LinkedHashSet kSet = new LinkedHashSet();
      
      	for (Iterator key_iter = allPossibleTreeNodeKeySet.iterator(); key_iter.hasNext(); ) {
      		LinkedHashSet keySet = new LinkedHashSet((Collection)key_iter.next());
      		Kinetics thisK = fillKineticsByAverage(keySet);
      		if (thisK!=null) kSet.add(thisK);
      	}
      
          k = ArrheniusKinetics.average(kSet);
      	if (k==null) return null;
      	kineticsTemplateLibrary.addKinetics(key,k);
      }
      
      return k;
      //#]
  }
  
  //## operation findClosestRateConstant(Collection) 
  public Kinetics[] findClosestRateConstant(LinkedList p_matchedPathSet) {
      //#[ operation findClosestRateConstant(Collection) 
      LinkedHashSet exactTreeNode = new LinkedHashSet();
      LinkedHashSet exactKey = new LinkedHashSet();
      // deal with the top tree node(leaf)
      // if it has a dummy child, at it into the matched path
      // put the tree node into the exact tree node  set
      for (Iterator pathIter = p_matchedPathSet.iterator(); pathIter.hasNext();) {
      	Stack s = (Stack)pathIter.next();
      	HierarchyTreeNode n = (HierarchyTreeNode)(s.peek());
      	if (n.hasDummyChild()) {
      		DummyLeaf dl = n.getDummyChild();
       		s.push(dl);
        		exactTreeNode.add(dl);
        		exactKey.add(dl.getName());
        	}
        	else {
        		exactTreeNode.add(n);
       		exactKey.add(n.getElement());
        	}
      }
      // find if exact kt exists.  if it does, return it.
      KineticsTemplate kt = kineticsTemplateLibrary.getKineticsTemplate(exactKey);
      int rNum = getReactantNumber();
      if (kt!=null) {
    	  Kinetics[] k_closest = new Kinetics[1];
    	  k_closest[0] = kt.kinetics;
    	  return k_closest;
//    	  return kt.kinetics;
      }
      
      Collection allPossibleTreeNodeSet = MathTool.expand(p_matchedPathSet.iterator());
      
      HashSet bestKineticsTemplateSet = new HashSet();
      LinkedHashSet bestKineticsSet = new LinkedHashSet();
      // find the closest distance
      int closest = Integer.MAX_VALUE;
      for (Iterator iter = allPossibleTreeNodeSet.iterator(); iter.hasNext(); ) {
         	LinkedHashSet treeNode = new LinkedHashSet((Collection)iter.next());
         	LinkedHashSet key = extractKineticsTemplateKey(treeNode);
        	kt = kineticsTemplateLibrary.getKineticsTemplate(key);
         	if (kt != null) {
      	  	Kinetics k = kt.getKinetics();
         		int distance = calculateDistance(exactTreeNode,treeNode);
         		if (distance < closest) {
         			closest = distance;
         			bestKineticsTemplateSet.clear();
         			bestKineticsTemplateSet.add(kt);
         			bestKineticsSet.clear();
         			bestKineticsSet.add(k);
         		}
         		else if (distance == closest) {
         			bestKineticsTemplateSet.add(kt);
         			bestKineticsSet.add(k);
         		}
         	}
      }
 /*     if (bestKineticsSet.size() == 0) {
    	 	System.out.println("problem with rate constant"); 
    	 	System.out.println(reactionAdjList.productNumber);
    	 	System.out.println(reactionAdjList.reactantNumber);
    	  
      }*/
      if (bestKineticsSet.size() == 0) throw new RateConstantNotFoundException();
      
      // get averaged k with the closest distance
      Kinetics newK = ArrheniusKinetics.average(bestKineticsSet);
      //KineticsTemplate newKT = kineticsTemplateLibrary.addKinetics(exactKey, newK);
      //RateConstant rc = new RateConstant(newKT, bestKineticsTemplateSet, closest);
      Kinetics[] k_closest = new Kinetics[1];
      k_closest[0] = newK;
      return k_closest;
//      return newK; 
      
      
      
      //#]
  }
  
  //## operation findExactRateConstant(Collection) 
  public Kinetics[] findExactRateConstant(Collection p_matchedPathSet) {
      //#[ operation findExactRateConstant(Collection) \
		
      LinkedHashSet fgc = new LinkedHashSet();
      for (Iterator iter = p_matchedPathSet.iterator(); iter.hasNext();) {
      	Stack s = (Stack)iter.next();
      	HierarchyTreeNode node = (HierarchyTreeNode)s.peek();
      	if (node.hasDummyChild()) fgc.add(node.getDummyChild().getName());
      	else {
      		Matchable fg = (Matchable)node.getElement();
      		fgc.add(fg);
				
      	}
      }
      KineticsTemplate kt = kineticsTemplateLibrary.getKineticsTemplate(fgc);
      
      if (kt==null) return null;
      else {
			//kt.kinetics.setSource(kt.kinetics.toChemkinString());
    	  Kinetics[] k_exact = new Kinetics[1];
    	  k_exact[0] = kt.kinetics;
    	  return k_exact;
//			return kt.kinetics;
      }
      
      
      
      
      //#]
  }
  
  /**
  Requires:
  Effects: call itsKineticsTemplateLibrary.findKinetics() to find out the kinetics for the structure, return the found kinetics or null if nothing is found.
  Modifies:
  */
  //## operation findRateConstant(Structure) 
  public Kinetics[] findRateConstant(Structure p_structure) {
		double pT = System.currentTimeMillis();
		
  	  /* 
  	   *If a primary reaction library exists, check the
  	   *	current reaction against that list before attempting
  	   *	to estimate via searching the tree
  	   */
		Kinetics[] k = null;
  	  if (doesPrimaryKineticLibraryExist()) {
  		  k = getPrimaryKineticRate(p_structure);
  	  }
  	  if (k != null) {
  		  setRateCoefficientSource(k);
  		  p_structure.setDirection(getPrimaryKineticDirection(p_structure));
  		  return k;
  	  }
		
      //#[ operation findRateConstant(Structure) 
      // look for kinetics in kinetics template libarry
      LinkedList reactants = null;
      if (isForward()) {
      	p_structure.setDirection(1);
      	reactants = p_structure.reactants;
      }
      else if (isBackward()) {
      	p_structure.setDirection(-1);
      }
      else if (isNeutral()) {
      	boolean thermoConsistence = true;
      	// in H abstraction, we allow biradical abstract H from a molecule, but the reverse is now allowed
      	// therefore, such H abs reactions will be all set as forward reaction
      	if (name.equals("H_Abstraction")) {
      		Iterator iter = p_structure.reactants.iterator();
      		while (iter.hasNext()) {
      			ChemGraph cg = (ChemGraph)iter.next();
      			int rNum = cg.getRadicalNumber();
      			if (rNum >= 2) {
      				thermoConsistence = false;
      				reactants = p_structure.reactants;
      				p_structure.setDirection(1);
      			}
      		}
      	}
      
          if (thermoConsistence) {
      		Temperature T = new Temperature(298, "K");
      		// to avoid calculation error's effect, lower the threshold for forward reaction
      		//if (p_structure.calculateKeq(T)>0.999)  {
      		/*if (p_structure.calculateKeq(T)>0.999 && p_structure.calculateKeq(T) <1.001) {
      			System.out.println(p_structure.toChemkinString(true));
      		}*/
      		if (p_structure.calculateKeq(T)>0.999) {
      			p_structure.setDirection(1);
      	 		reactants = p_structure.reactants;
      		}
      		else {
      			p_structure.setDirection(-1);
      		}
      		// for intra h migration, set the ROO. as the forward
        		if (name.equals("intra_H_migration")) {
              ChemGraph rcg = (ChemGraph)(p_structure.getReactants().next());
          	HashSet rrad = rcg.getRadicalNode();
           	Atom rra = (Atom)( (Node) ( (rrad.iterator()).next())).getElement();
              ChemGraph pcg = (ChemGraph)(p_structure.getProducts().next());
            	HashSet prad = pcg.getRadicalNode();
             	Atom pra = (Atom)( (Node) ( (prad.iterator()).next())).getElement();
              if (rra.isOxygen() && pra.isCarbon()) {
              	p_structure.setDirection(1);
               	reactants = p_structure.reactants;
              }
              else if (pra.isOxygen() && rra.isCarbon())
                  p_structure.setDirection(-1);
              }
      
      	}
      }
      else {
      	throw new InvalidReactionTemplateDirectionException();
      }
      
      if (p_structure.isForward()) {
      	LinkedList fg = structureTemplate.getMatchedFunctionalGroup(reactants);
			if (fg == null) {
				Global.RT_findRateConstant += (System.currentTimeMillis()-pT)/1000/60;
				return null;
			}
			String comments = getKineticsComments(fg);
      	k = findExactRateConstant(fg);
      	if (k==null) {
				k = findClosestRateConstant(fg);
				k[0].setSource(name + " estimate: (" + k[0].getSource() + ")");
      	}
      	else k[0].setSource(name  + " exact: ");
			k[0].setComments(comments);
			Global.RT_findRateConstant += (System.currentTimeMillis()-pT)/1000/60;
      	return k;
      }
      else {
			Global.RT_findRateConstant += (System.currentTimeMillis()-pT)/1000/60;
			return null;
      }
  }
  
  
  /**
  Requires:
  Effects: call itsKineticsTemplateLibrary.findKinetics() to find out the kinetics for the structure, return the found kinetics or null if nothing is found.
  Modifies:
  */
  //## operation findRateConstant(Structure) 
  public Kinetics[] findReverseRateConstant(Structure p_structure) {
		double pT = System.currentTimeMillis();
      //#[ operation findRateConstant(Structure)
		
	  	  /* 
	  	   *If a primary kinetic library exists, check the
	  	   *	current reaction against that list before attempting
	  	   *	to estimate via searching the tree
	  	   */
			Kinetics[] k = null;
	  	  if (doesPrimaryKineticLibraryExist()) {
	  		  k = getPrimaryKineticRate(p_structure);
	  	  }
	  	  if (k != null) {
	  		  setRateCoefficientSource(k);
	  		  p_structure.setDirection(1);
	  		  return k;
	  	  }
		
      // look for kinetics in kinetics template libarry
      LinkedList reactants = null;
      if (isForward()) {
      	p_structure.setDirection(1);
      	reactants = p_structure.reactants;
      }
      else if (isBackward()) {
      	p_structure.setDirection(-1);
      }
      else if (isNeutral()) {
      	boolean thermoConsistence = true;
      	// in H abstraction, we allow biradical abstract H from a molecule, but the reverse is now allowed
      	// therefore, such H abs reactions will be all set as forward reaction
      	if (name.equals("H_Abstraction")) {
      		Iterator iter = p_structure.reactants.iterator();
      		while (iter.hasNext()) {
      			ChemGraph cg = (ChemGraph)iter.next();
      			int rNum = cg.getRadicalNumber();
      			if (rNum >= 2) {
      				thermoConsistence = false;
      				reactants = p_structure.reactants;
      				p_structure.setDirection(1);
      			}
      		}
      	}
      
          if (thermoConsistence) {
      		Temperature T = new Temperature(298, "K");
      		// to avoid calculation error's effect, lower the threshold for forward reaction
      		//if (p_structure.calculateKeq(T)>0.999)  {
      		/*if (p_structure.calculateKeq(T)>0.999 && p_structure.calculateKeq(T) <1.001) {
      			System.out.println(p_structure.toChemkinString(true));
      		}*/
      		if (p_structure.calculateKeq(T)>0.999) {
      			p_structure.setDirection(1);
      	 		reactants = p_structure.reactants;
      		}
      		else {
      			p_structure.setDirection(-1);
      		}
      		// for intra h migration, set the ROO. as the forward
        		if (name.equals("intra_H_migration")) {
              ChemGraph rcg = (ChemGraph)(p_structure.getReactants().next());
          	HashSet rrad = rcg.getRadicalNode();
           	Atom rra = (Atom)( (Node) ( (rrad.iterator()).next())).getElement();
              ChemGraph pcg = (ChemGraph)(p_structure.getProducts().next());
            	HashSet prad = pcg.getRadicalNode();
             	Atom pra = (Atom)( (Node) ( (prad.iterator()).next())).getElement();
              if (rra.isOxygen() && pra.isCarbon()) {
              	p_structure.setDirection(1);
               	reactants = p_structure.reactants;
              }
              else if (pra.isOxygen() && rra.isCarbon())
                  p_structure.setDirection(-1);
              }
      
      	}
      }
      else {
      	throw new InvalidReactionTemplateDirectionException();
      }
      
      if (p_structure.isForward()) {
      	LinkedList fg = structureTemplate.getMatchedFunctionalGroup(reactants);
			if (fg == null) {
				Global.RT_findRateConstant += (System.currentTimeMillis()-pT)/1000/60;
				return null;
			}
			String comments = getKineticsComments(fg);
      	k = findExactRateConstant(fg);
      	if (k==null) {
      		try{
				k = findClosestRateConstant(fg);
				k[0].setSource(name + " estimate: (" + k[0].getSource() + ")");
      		}
      		catch (RateConstantNotFoundException e) {
              	return k;
              }

      	}
      	else k[0].setSource(name  + " exact: ");
			k[0].setComments(comments);
			Global.RT_findRateConstant += (System.currentTimeMillis()-pT)/1000/60;
      	return k;
      }
      else {
			Global.RT_findRateConstant += (System.currentTimeMillis()-pT)/1000/60;
			return null;
      }
  }
  
  public void setRateCoefficientSource(Kinetics[] k) {
	  for (int numK=0; numK<k.length; numK++) {
  		  k[numK].setFromPrimaryKineticLibrary(true);
  		  String currentSource = k[numK].getSource();
  		  if (this.direction == -1) {
  			  if (!currentSource.contains(this.reverseReactionTemplate.name))
  				  k[numK].setSource(this.reverseReactionTemplate.name+" "+k[numK].getSource());
  		  }
  		  else	{
  			  if (!currentSource.contains(this.name))
  				  k[numK].setSource(this.name+" "+k[numK].getSource());
  		  }
	  }
  }
  
  /**
   * Given a structure, function searches the primary reaction library
   * 	to determine if user supplied a set of kinetic parameters for
   * 	this reaction.  If so, RMG uses this set of parameters.  If not,
   * 	RMG estimates the kinetic parameters by traversing the tree.
   * 
   * @param p_structure: Structure of the reaction
   * @return
   */
  private Kinetics[] getPrimaryKineticRate(Structure p_structure) {
	  boolean equivReactants = false;
	  LinkedHashSet primaryKinLibrary = getPrimaryKinLibrary();
	  Iterator iter = primaryKinLibrary.iterator();
	  LinkedList reactants_asdefinedbyRMG = p_structure.getReactantList();
	  LinkedList products_asdefinedbyRMG = p_structure.getProductList();
	  while (iter.hasNext()) {
		  Reaction rxn = (Reaction)iter.next();
		  LinkedList reactants_asdefinedinPRL = rxn.getReactantList();
		  /*
		   * The LinkedList "reactants_asdefinedbyRMG" contains species generated by RMG.
		   * 	The LinkedList contains chemgraphs (and species, indirectly).
		   * The LinkedList "reactants_asdefinedinPRL" contains the species read in by RMG,
		   * 	from the user-defined primary reaction library.  When reading in the user-
		   * 	specified adjacency lists, RMG makes a Species (meaning all resonance 
		   * 	isomers have been explored / stored).
		   * Thus, both reactants_asdefinedbyRMG and reactants_asdefinedinPRL have all
		   * 	of the resonance isomers defined ... if you know where to look
		   */
		  equivReactants = Structure.isSpeciesListEquivalentToChemGraphListAsChemGraphs(reactants_asdefinedinPRL,reactants_asdefinedbyRMG);
		  if (equivReactants) {
			  LinkedList products_asdefinedinPRL = rxn.getProductList();
			  if (Structure.isSpeciesListEquivalentToChemGraphListAsChemGraphs(products_asdefinedinPRL,products_asdefinedbyRMG)) {
				  if (rxn instanceof ThirdBodyReaction || rxn instanceof TROEReaction || rxn instanceof LindemannReaction)
					  System.out.println("RMG is only utilizing the high-pressure limit parameters for PKL reaction: " + rxn.toString());
				  return rxn.getKinetics();
			  }
		  }
		
   }
	  return null;
  }
  
  private int getPrimaryKineticDirection(Structure p_structure) {
	  boolean equivReactants = false;
	  LinkedHashSet primaryKinLibrary = getPrimaryKinLibrary();
	  Iterator iter = primaryKinLibrary.iterator();
	  LinkedList p_reactants = p_structure.getReactantList();
	  LinkedList p_products = p_structure.getProductList();
	  while (iter.hasNext()) {
		  Reaction rxn = (Reaction)iter.next();
		  LinkedList reactants = rxn.getReactantList();
		  equivReactants = Structure.isSpeciesListEquivalentToChemGraphListAsChemGraphs(reactants, p_reactants);
		  if (equivReactants) {
			  LinkedList products = rxn.getProductList();
			  if (Structure.isSpeciesListEquivalentToChemGraphListAsChemGraphs(products,p_products)) {
				  if (rxn instanceof ThirdBodyReaction || rxn instanceof TROEReaction || rxn instanceof LindemannReaction)
					  System.out.println("RMG is only utilizing the high-pressure limit parameters for PKL reaction: " + rxn.toString());
				  return rxn.getStructure().direction;
			  }
		  }
	  }
	  return -1;
  }

  private String getKineticsComments(Collection p_matchedPathSet) {
		StringBuilder comment = new StringBuilder();
		for (Iterator iter = p_matchedPathSet.iterator(); iter.hasNext();) {
      	Stack s = (Stack)iter.next();
      	HierarchyTreeNode node = (HierarchyTreeNode)s.peek();
      	if (node.hasDummyChild()) comment.append((node.getDummyChild().getName()+"    "));
      	else {
      		if (node.getElement() instanceof FunctionalGroup)
      			comment.append(((FunctionalGroup)node.getElement()).getName() + "     ");
      		else
      			comment.append(((FunctionalGroupCollection)node.getElement()).getName() + "     ");
      	}
      }
		return comment.toString();
	}

	//## operation findUnion(String,HashMap) 
  private void findUnion(String p_name, HashMap p_unRead) {
      //#[ operation findUnion(String,HashMap) 
      HashSet union = (HashSet)p_unRead.get(p_name);
      FunctionalGroupCollection fgc = new FunctionalGroupCollection(p_name);
      Iterator union_iter = union.iterator();
      while (union_iter.hasNext()) {
      	String fg_name = (String)union_iter.next();
      	if (p_unRead.containsKey(fg_name)) {
      		findUnion(fg_name, p_unRead);
      	}
      	Matchable fg = (Matchable)fgDictionary.get(fg_name);
      	if (fg == null)
      		throw new InvalidFunctionalGroupException("Unknown FunctionalGroup in findUnion(): " + fg_name);
      	fgc.addFunctionalGroups(fg);
      }
      fgDictionary.put(p_name,fgc);
      p_unRead.remove(p_name);
      return;

      //#]
  }
  
  /**
  Requires:
  Effects: Return the reversed reaction template of this reaction template if this reaction template's direction is forward.  Reverse reaction template in such aspects:
  (1) structure template
  (2) kinetics template
  (3) direction
  Note: if the reversereaction template is generated successfully, set this.itsReverseReactionTemplate to this new generated reaction template. 
  Modifies: this.itsReverseReactionTemplate.
  */
  // Argument Stringp_name : 
  /**
  pass in the name of the reverse reaction template
  */
  //## operation generateReverseReactionTemplate(String) 
  public ReactionTemplate generateReverseReactionTemplate(String p_name) {
      //#[ operation generateReverseReactionTemplate(String) 
      // guarantee validity of the template and its direction
      if (direction != 1 || !repOk()) {
      	return null;
      }
      if (name.equals("Diels-Alder reaction"))
    	  System.out.println("stop");
      ReactionTemplate reversetemplate = new ReactionTemplate();
      
      // set direction
      reversetemplate.direction = -1;
      
      // set name
      reversetemplate.name = p_name;
      
      // set the structure template
      reversetemplate.structureTemplate = structureTemplate.generateReverse(reactionAdjList);
      
      // set the reaction adjList
      reversetemplate.reactionAdjList = reactionAdjList.generateReverse();
      
      // set the kinetics template library to the forward kinetics template library to keep thermo consistency
      reversetemplate.kineticsTemplateLibrary = kineticsTemplateLibrary;
      
      // set itsReverseReactionTemplate to this generated reversetemplate
      reversetemplate.reverseReactionTemplate = this;
      this.reverseReactionTemplate = reversetemplate;
      
      return reversetemplate;
      
      
      
      
      
      //#]
  }
  
  //## operation getProductNumber() 
  public int getProductNumber() {
      //#[ operation getProductNumber() 
      return reactionAdjList.getProductNumber();
      //#]
  }
  
  //## operation getReactantNumber() 
  public int getReactantNumber() {
      //#[ operation getReactantNumber() 
      return structureTemplate.getReactantNumber();
      //#]
  }
  
  //## operation getReactantTree() 
  public Iterator getReactantTree() {
      //#[ operation getReactantTree() 
      return getStructureTemplate().getReactantTree();
      //#]
  }
  
  //## operation getReactantTreeNumber() 
  public int getReactantTreeNumber() {
      //#[ operation getReactantTreeNumber() 
      return structureTemplate.getReactantTreeNumber();
      //#]
  }
  
  //## operation getReactionFromStructure(Structure) 
  public TemplateReaction getReactionFromStructure(Structure p_structure) {
      //#[ operation getReactionFromStructure(Structure) 
      return (TemplateReaction)reactionDictionaryByStructure.get(p_structure);
      //#]
  }
  
  
  //## operation hasOneReactant() 
  public boolean hasOneReactant() {
      //#[ operation hasOneReactant() 
      return (getReactantNumber() == 1);
      //#]
  }
  
  //## operation hasTwoReactants() 
  public boolean hasTwoReactants() {
      //#[ operation hasTwoReactants() 
      return (getReactantNumber() == 2);
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
  
  //## operation isNeutral() 
  public boolean isNeutral() {
      //#[ operation isNeutral() 
      return (direction == 0);
      //#]
  }
  
  //## operation isReverse(ReactionTemplate) 
  public boolean isReverse(ReactionTemplate p_reactionTemplate) {
      //#[ operation isReverse(ReactionTemplate) 
      if (direction == 0 || reverseReactionTemplate == null) return false;
      return reverseReactionTemplate.equals(p_reactionTemplate);
      //#]
  }
  
  /**
  Requires:
  Effects: react p_reactant according to this reaction template to form a set of reactions.  (generate product according to the reaction adjlist, find out kinetics sccording to the kinticstemplatelibrary).  return the set of reactions. 
  Modifies:
  */
  //## operation reactOneReactant(Species) 
  public LinkedHashSet reactOneReactant(Species p_reactant) {
      //#[ operation reactOneReactant(Species) 
      LinkedList r = new LinkedList();
      r.add(p_reactant);
		
      /*HashSet reaction = getReactionSetFromReactant(r);
      if (reaction != null)
      	return reaction;*/
      
      ChemGraph rg;
      LinkedHashSet reaction = new LinkedHashSet();
      
      if (!p_reactant.hasResonanceIsomers()) {
      	// react the main structre
      	rg = p_reactant.getChemGraph();
      	reaction = reactOneReactant(rg);
      }
      else {
      	// react the resonance isomers
      	Iterator iter = p_reactant.getResonanceIsomers();
      	while (iter.hasNext()) {
      		rg = (ChemGraph)iter.next();
      		LinkedHashSet more = reactOneReactant(rg);
      		reaction.addAll(more);
      	}
      }
      
      // put in reactionDictionaryByReactant
      //reactionDictionaryByReactant.put(r,reaction);
      
      // return the reaction set
      return reaction;
      
      
      
      
      
      //#]
  }
  
	   //## operation reactOneReactant(ChemGraph) 
  protected LinkedHashSet reactOneReactant(ChemGraph p_chemGraph) {
      //#[ operation reactOneReactant(ChemGraph) 
	  LinkedHashSet reaction_set = new LinkedHashSet();
	//  if (name.equals("Intra_R_Add_Endocyclic") && !p_chemGraph.isAcyclic()) return reaction_set; //commented out by gmagoon 7/29/09: according to Sandeep this may be an unneeded artifact of an old bug that may have been addressed with the use of ForbiddenGroups for Intra_R_Add_Endocyclic and Intra_R_Add_Exocyclic; however, it should be kept in mind that when this is commented out, it is probably more likely that one will obtain reactions for which the group identification (and hence kinetics estimates) are not reproducible due to identification of different ring paths in the same molecule
      

	  
	
      LinkedHashSet allReactionSites = structureTemplate.identifyReactedSites(p_chemGraph,1);
	  //System.out.println("Species: "+p_chemGraph.toString());
      
      if (allReactionSites.isEmpty()) return reaction_set;
      
      // add present chemgraph into the reactant linked list
      LinkedList reactant = new LinkedList();
	  LinkedList reactantSp = new LinkedList();
      reactant.add(p_chemGraph);
	  reactantSp.add(p_chemGraph.getSpecies());
      
      LinkedHashMap structureMap = new LinkedHashMap();
      LinkedHashMap rateMap = new LinkedHashMap();
		LinkedHashMap reactionMap = new LinkedHashMap();
      
      for (Iterator iter = allReactionSites.iterator(); iter.hasNext(); ) {
		  
          MatchedSite ms = (MatchedSite)iter.next();
          HashMap site = ms.getCenter();
       	  int redundancy = ms.getRedundancy();
		  //System.out.println(ms.toString());
          // reset the reacted site for rg in reactant linkedlist
          p_chemGraph.resetReactedSite(site);
          
          boolean forbidden = false;
          Iterator forbiddenIter = forbiddenStructures.values().iterator();
          while (forbiddenIter.hasNext()){
        	  Matchable fg = (Matchable)forbiddenIter.next();
        	  if (p_chemGraph.isSubAtCentralNodes(fg)){
        		  forbidden = true;
        		  break;
        	  }
          }
          if (forbidden) break;
          
          // react reactant to form a new structure
          try {
          	LinkedList product = reactionAdjList.reactChemGraph(reactant);
			LinkedList productSp = new LinkedList();
				SpeciesDictionary sd = SpeciesDictionary.getInstance();
				for (int i=0; i< product.size(); i++){
					String name = null;
					if (((ChemGraph)product.get(i)).getSpecies() == null){
						Species sp = Species.make(name, ((ChemGraph)product.get(i)));
						productSp.add(sp);
					}
				}
				double pt = System.currentTimeMillis();
  	        boolean rpsame = MathTool.isListEquivalent(reactantSp, productSp);
				Global.checkReactionReverse = Global.checkReactionReverse + (System.currentTimeMillis()-pt)/1000/60;
  	        if (!rpsame) {
      			Structure structure = new Structure(reactant,product);
				Kinetics[] k = findRateConstant(structure);
				Structure structureSp = new Structure(reactantSp, productSp);
				
				
				
				structureSp.direction = structure.direction;
      			structure.setRedundancy(redundancy);
					Reaction old_reaction = (Reaction)reactionMap.get(structureSp);
					if (old_reaction == null){
						
						TemplateReaction r = TemplateReaction.makeTemplateReaction(structureSp,k, this,structure);
						
						structure = null;
						if (r != null)
							reactionMap.put(structureSp,r);
					}
					else {
						if (k == null)
							old_reaction.addAdditionalKinetics(null,redundancy);
						else {
							for (int i=0; i<k.length; i++) {
								old_reaction.addAdditionalKinetics(k[i], redundancy);
							}
						}
						//old_reaction.getStructure().increaseRedundancy(redundancy);
						structureSp = null;
						structure = null;
					}
      		}
          }
          catch (ForbiddenStructureException e) {
          	// do the next reaction site
          }
      	catch (InvalidProductNumberException e) {
      		// do the next reaction site
      	}
      }
      
      for (Iterator mapIter = reactionMap.values().iterator(); mapIter.hasNext(); ) {
      	
      	Reaction reaction = (Reaction)mapIter.next();
      	//System.out.println(reaction.toString());
      	if (!reaction.repOk()) throw new InvalidTemplateReactionException(reaction.toString());
			//reaction.getStructure().clearChemGraph();
		reaction.setFinalized(true);
      	reaction_set.add(reaction);
      }
      
      return reaction_set;
      //#]
  }
	
//	## operation reactTwoReactants(ChemGraph,ChemGraph) 
  /*protected TemplateReaction calculateForwardRateConstant(ChemGraph chemGraph1, ChemGraph chemGraph2, Structure p_structure) {
      //#[ operation reactTwoReactants(ChemGraph,ChemGraph) 
      
      ChemGraph r1 = chemGraph1;
      ChemGraph r2 = chemGraph2;
      if (r1 == r2) {
      	try{
      		r2 = ChemGraph.copy(r2);
      	}
      	catch (ForbiddenStructureException e) {
      		return null;
      	}
      }
      
		TemplateReaction reverseReaction = null;
      /*TemplateReaction reverseReaction = getReactionFromStructure(p_structure);
		if (reverseReaction != null){
			if (!reverseReaction.isForward()) throw new InvalidReactionDirectionException();
			return reverseReaction;
		}
      HashSet allReactionSites1 = structureTemplate.identifyReactedSites(r1,1);
      HashSet allReactionSites2 = structureTemplate.identifyReactedSites(r2,2);
      
      LinkedList reactant = new LinkedList();
      reactant.add(r1);
      reactant.add(r2);
      
		HashSet rateSet = new HashSet();
      
      for (Iterator iter1 = allReactionSites1.iterator(); iter1.hasNext(); ) {
      	MatchedSite ms1 = (MatchedSite)iter1.next();
      	HashMap site1 = ms1.getCenter();
          r1.resetReactedSite(site1);
      
          int redundancy1 = ms1.getRedundancy();
      
      	for (Iterator iter2 = allReactionSites2.iterator(); iter2.hasNext(); ) {
      		MatchedSite ms2 = (MatchedSite)iter2.next();
      		HashMap site2 = (HashMap)ms2.getCenter();
      	    r2.resetReactedSite(site2);
      
      	    int redundancy2 = ms2.getRedundancy();
      		int redundancy = redundancy1*redundancy2;
      
           	try {
      	     	LinkedList product = reactionAdjList.reactChemGraph(reactant);
					SpeciesDictionary sd = SpeciesDictionary.getInstance();
					for (int i=0; i< product.size(); i++){
						String name = null;
						if (((ChemGraph)product.get(i)).getSpecies() == null){
							Species sp = Species.make(name, ((ChemGraph)product.get(i)));
						}
					}
      	        boolean rpsame = MathTool.isListEquivalent(reactant, product);
      	        if (!rpsame) {
          			Structure structure = new Structure(reactant,product);
						if (structure.equals(p_structure)){
							Kinetics k = findRateConstant(p_structure);
							if (reverseReaction == null){
								reverseReaction = TemplateReaction.makeTemplateReaction(p_structure,k,this);
							}
							else {
								//p_structure.increaseRedundancy(redundancy);
								reverseReaction.addAdditionalKinetics(k,redundancy);
								//structure = null;
							}
								
	        			}
						/*else {
							SpeciesDictionary sd = SpeciesDictionary.getInstance();
							while (!product.isEmpty())
								sd.remove(((ChemGraph)product.remove()));
							
						}
          			
      	        }
      			
      		}
      		catch (ForbiddenStructureException e) {
      			// do the next reaction site
      		}
      		catch (InvalidProductNumberException e) {
      			// do the next reaction site
      		}
				catch (InvalidChemGraphException e){
					
				}
      	}
      }
      //reverseReaction.getStructure().clearChemGraph();
      return reverseReaction;
      //#]
  }*/
  
//	## operation reactTwoReactants(ChemGraph,ChemGraph) 
  protected TemplateReaction calculateForwardRateConstant(ChemGraph p_chemGraph, Structure p_structure) {
      //#[ operation reactTwoReactants(ChemGraph,ChemGraph) 
      
		
		TemplateReaction reverseReaction = null;
		
		LinkedHashSet rateSet = new LinkedHashSet();
		LinkedHashSet allReactionSites = structureTemplate.identifyReactedSites(p_chemGraph,1);
      
      if (allReactionSites.isEmpty()) return reverseReaction;
      
      // add present chemgraph into the reactant linked list
      LinkedList reactant = new LinkedList();
	  LinkedList reactantSp = new LinkedList();
	  reactantSp.add(p_chemGraph.getSpecies());
      reactant.add(p_chemGraph);
   
      
      for (Iterator iter = allReactionSites.iterator(); iter.hasNext(); ) {
          MatchedSite ms = (MatchedSite)iter.next();
          HashMap site = ms.getCenter();
       	int redundancy = ms.getRedundancy();
      
          // reset the reacted site for rg in reactant linkedlist
          p_chemGraph.resetReactedSite(site);
          // react reactant to form a new structure
          try {
          	LinkedList product = reactionAdjList.reactChemGraph(reactant);
			LinkedList productSp = new LinkedList();
				SpeciesDictionary sd = SpeciesDictionary.getInstance();
				for (int i=0; i< product.size(); i++){
					String name = null;
					if (((ChemGraph)product.get(i)).getSpecies() == null){
						Species sp = Species.make(name, ((ChemGraph)product.get(i)));
						productSp.add(sp);
					}
				}
      	    boolean rpsame = MathTool.isListEquivalent(reactantSp, productSp);
      	    if (!rpsame) {
      			Structure structure = new Structure(reactant,product);
				Structure structureSp = new Structure(reactantSp,productSp);
      			if (structure.equalsAsChemGraph(p_structure)){
						Kinetics[] k = findRateConstant(p_structure);
						structureSp.direction = p_structure.direction;
						if (reverseReaction == null){
							reverseReaction = TemplateReaction.makeTemplateReaction(structureSp,k,this,p_structure);
						}
						else {
							//p_structure.increaseRedundancy(redundancy);
							for (int i=0; i<k.length; i++) {
								reverseReaction.addAdditionalKinetics(k[i],redundancy);
							}
							//structure = null;
						}
							
      			}
					/*else {
						SpeciesDictionary sd = SpeciesDictionary.getInstance();
						while (!product.isEmpty())
							sd.remove(((ChemGraph)product.remove()));
						
					}*/
      		}
          }
          catch (ForbiddenStructureException e) {
          	// do the next reaction site
          }
      	catch (InvalidProductNumberException e) {
      		// do the next reaction site
      	}
      }
		//reverseReaction.getStructure().clearChemGraph();
	  reverseReaction.setFinalized(true);
      return reverseReaction;
      //#]
  }
 
  /**
   * reactTwoReactants
   * 	17-Jun-2009 MRH
   * Function "reacts" two reactants (in the form of ChemGraphs).
   * 
   * @param cg1: ChemGraph of species1
   * @param rs1: Collection of reactive sites for cg1
   * @param cg2: ChemGraph of species2
   * @param rs2: Collection of reactive sties for cg2
   * @return
   */
  protected LinkedHashSet reactTwoReactants(ChemGraph cg1, LinkedHashSet rs1, ChemGraph cg2, LinkedHashSet rs2) {

	  LinkedHashSet reaction_set = new LinkedHashSet();
      if (rs1.isEmpty() || rs2.isEmpty()) return reaction_set;
      
      ChemGraph r1 = cg1;
      ChemGraph r2;
      if (cg1 == cg2) {
      	try{
      		r2 = ChemGraph.copy(cg2);
      	}
      	catch (ForbiddenStructureException e) {
      		return reaction_set;
      	}
      }
      else {
      	r2 = cg2;
      }     
      
      LinkedList reactant = new LinkedList();
	  LinkedList reactantSp = new LinkedList();
      reactant.add(r1);
      reactant.add(r2);
	  reactantSp.add(r1.getSpecies());
	  reactantSp.add(r2.getSpecies());
      
	  LinkedHashMap structureMap = new LinkedHashMap();
	  LinkedHashMap rateMap = new LinkedHashMap();
	  LinkedHashMap reactionMap = new LinkedHashMap();
      
      for (Iterator iter1 = rs1.iterator(); iter1.hasNext(); ) {
      	MatchedSite ms1 = (MatchedSite)iter1.next();
      	HashMap site1 = ms1.getCenter();
          r1.resetReactedSite(site1);
      
          boolean forbidden1 = false;
          Iterator forbiddenIter = forbiddenStructures.values().iterator();
          while (forbiddenIter.hasNext()){
        	  Matchable fg = (Matchable)forbiddenIter.next();
        	  if (r1.isSubAtCentralNodes(fg)){
        		  forbidden1 = true;
        		  break;
        	  }
          }
          if (forbidden1) continue;
          
          int redundancy1 = ms1.getRedundancy();
      
      	for (Iterator iter2 = rs2.iterator(); iter2.hasNext(); ) {
      		MatchedSite ms2 = (MatchedSite)iter2.next();
      		HashMap site2 = (HashMap)ms2.getCenter();
      	    r2.resetReactedSite(site2);
      
      	    boolean forbidden2 = false;
      	    forbiddenIter = forbiddenStructures.values().iterator();
      	    while (forbiddenIter.hasNext()){
      	    	Matchable fg = (Matchable)forbiddenIter.next();
      	    	if (r2.isSubAtCentralNodes(fg)){
      	    		forbidden2 = true;
      	    		break;
      	    	}
 	         }
 	         if (forbidden2) continue;
          
      	    int redundancy2 = ms2.getRedundancy();
      		int redundancy = redundancy1*redundancy2;
      
           	try {
      	     	LinkedList product = reactionAdjList.reactChemGraph(reactant);
				LinkedList productSp = new LinkedList();
//				SpeciesDictionary sd = SpeciesDictionary.getInstance();
				for (int i=0; i< product.size(); i++){
					String name = null;
					if (((ChemGraph)product.get(i)).getSpecies() == null){
						Species sp = Species.make(name, ((ChemGraph)product.get(i)));
						productSp.add(sp);
					}
				}
				double pt = System.currentTimeMillis();
      	        boolean rpsame = MathTool.isListEquivalent(reactantSp, productSp);
				Global.checkReactionReverse = Global.checkReactionReverse + (System.currentTimeMillis()-pt)/1000/60;
      	        if (!rpsame) {
          			Structure structure = new Structure(reactant,product);
					Structure structureSp = new Structure(reactantSp,productSp);
					Kinetics[] k = findRateConstant(structure);
					structureSp.direction = structure.direction;
          			structure.setRedundancy(redundancy);
  					Reaction old_reaction = (Reaction)reactionMap.get(structureSp);
  					if (old_reaction == null){		
						 
						TemplateReaction r= TemplateReaction.makeTemplateReaction(structureSp,k, this,structure);
						if (r != null)
							reactionMap.put(structureSp,r);
						structure = null;
  					}
  					else {
  						if (k == null)
  							old_reaction.addAdditionalKinetics(null,redundancy);
  						else {
	  						for (int i=0; i<k.length; i++) {
	  							old_reaction.addAdditionalKinetics(k[i], redundancy);
	  						}
  						}
						//old_reaction.getStructure().increaseRedundancy(redundancy);
						structure = null;
  					}
      	        }
      			
      		}
      		catch (ForbiddenStructureException e) {
      			// do the next reaction site
      		}
      		catch (InvalidProductNumberException e) {
      			// do the next reaction site
      		}
			catch (InvalidChemGraphException e){
				
			}
      	}
      }
      
      for (Iterator mapIter = reactionMap.values().iterator(); mapIter.hasNext(); ) {
      	
      	Reaction reaction = (Reaction)mapIter.next();
      	//System.out.println(reaction.toString());
      	if (!reaction.repOk()) throw new InvalidTemplateReactionException(reaction.toString());
      	//reaction.getStructure().clearChemGraph();
		reaction.setFinalized(true);
      	reaction_set.add(reaction);
      }
      
      return reaction_set;
  }
  
  /**
  Requires:
  Effects: read in the reaction templated defined in the pass-in directory
  Modifies:
  */
  //## operation read(String,String) 
  public void read(String p_reactionTemplateName, String p_directoryName) {
      //#[ operation read(String,String) 
      String directoryName;
      
      if (!(p_directoryName.endsWith("/"))) {
      	directoryName = p_directoryName + "/";
      }
      else {
      	directoryName = p_directoryName;
      }
      
      setName(p_reactionTemplateName);
      System.out.println("Reading forwards template:   "+p_reactionTemplateName);
      
      String ReactionAdjListName = directoryName + "reactionAdjList.txt";
      String DictionaryName = directoryName + "dictionary.txt";
      String TreeName = directoryName + "tree.txt";
      String LibraryName = directoryName + "rateLibrary.txt";
      String ForbiddenName = directoryName + "forbiddenGroups.txt";
      
      try {
      	readFGDictionary(DictionaryName);
      	readForbiddenStructures(ForbiddenName);
      	String reverseRTName = readReactionAdjList(ReactionAdjListName);
      	readTree(TreeName);
      	readLibrary(LibraryName);
      	fillKineticsBottomToTop();
      	if (reverseRTName != null && reverseRTName.compareToIgnoreCase("none")!=0) {
      		System.out.println("Generating reverse template: "+reverseRTName);
      		reverseReactionTemplate = generateReverseReactionTemplate(reverseRTName);
      		reverseReactionTemplate.forbiddenStructures = forbiddenStructures;
      	}
      }
      catch (Exception e) {
      	System.err.println("Error in read in reaction template: " + name);
      	System.err.println(e.getMessage());
      	System.exit(0);
      }
      
      return;
      //#]
  }
  
  public void readForbiddenStructures(String p_fileName) throws  IOException {
      //#[ operation readFGDictionary(String) 
      try {
				if (!(new File(p_fileName)).exists()) {
					// System.out.println("forbiddenStructures file does not exist");
					return;
				}
      	FileReader in = new FileReader(p_fileName);
      	BufferedReader data = new BufferedReader(in);
      	HashMap unRead = new HashMap();
      	String fgname = null;
      
      	// step 1: read in structure
      	String line = ChemParser.readMeaningfulLine(data);
      	read: while (line != null) {
      		StringTokenizer token = new StringTokenizer(line);
      		fgname = token.nextToken();
      		data.mark(10000);
      		line = ChemParser.readMeaningfulLine(data);
      		if (line == null) break read;
      		line = line.trim();
      		String prefix = line.substring(0,5);
      		if (prefix.compareToIgnoreCase("union") == 0) {
      			HashSet union = ChemParser.readUnion(line);
       			unRead.put(fgname,union);
      		}
      		else {
      			data.reset();
      			Graph fgGraph = null;
      			try {
      				fgGraph = ChemParser.readFGGraph(data);
      			}
      			catch (InvalidGraphFormatException e) {
      				throw new InvalidFunctionalGroupException(fgname + ": " + e.getMessage());
      			}
      			if (fgGraph == null)
                      throw new InvalidFunctionalGroupException(fgname);
      
      			FunctionalGroup fg = FunctionalGroup.make(fgname, fgGraph);
      			Object old = forbiddenStructures.get(fgname);
      			if (old == null) {
      				forbiddenStructures.put(fgname,fg);
      			}
      			else {
      				FunctionalGroup oldFG = (FunctionalGroup)old;
      				if (!oldFG.equals(fg)) throw new ReplaceFunctionalGroupException(fgname);
      			}
              }
      		line = ChemParser.readMeaningfulLine(data);
      	}
      
      	while (!unRead.isEmpty()) {
      		fgname = (String)(unRead.keySet().iterator().next());
      		ChemParser.findUnion(fgname,unRead,forbiddenStructures);
      	}
      
          in.close();
      	return;
      }
      catch (Exception e) {
    	  System.out.println("Failed to read forbiddenStructures file");
      	//throw new IOException(e.getMessage());
      }
      
      
      
      
      
      
      //#]
  }
  
  
  //## operation readFGDictionary(String) 
  public void readFGDictionary(String p_fileName) throws  IOException {
      //#[ operation readFGDictionary(String) 
      try {
      	FileReader in = new FileReader(p_fileName);
      	BufferedReader data = new BufferedReader(in);
      	HashMap unRead = new HashMap();
      	String fgname = null;
      
      	// step 1: read in structure
      	String line = ChemParser.readMeaningfulLine(data);
      	read: while (line != null) {
      		StringTokenizer token = new StringTokenizer(line);
      		fgname = token.nextToken();
			
			if (fgname.toLowerCase().startsWith("others")) {
				System.out.println("Skipping dictionary definition of group "+fgname+" because its begins with 'others' and that has special meaning.");
				// gobble up the rest
				while (line!=null) {
					line=ChemParser.readUncommentLine(data);
				}
				// now get an unblank line
				line = ChemParser.readMeaningfulLine(data);
				continue read;
			}
			
      		data.mark(10000);
      		line = ChemParser.readMeaningfulLine(data);
      		if (line == null) break read;
      		line = line.trim();
      		if (line.toLowerCase().startsWith("union") || line.startsWith("OR") ) {
      			HashSet union = ChemParser.readUnion(line);
       			unRead.put(fgname,union);
      		}
      		else {
      			data.reset();
      			Graph fgGraph = null;
      			try {
      				fgGraph = ChemParser.readFGGraph(data);
      			}
      			catch (InvalidGraphFormatException e) {
      				throw new InvalidFunctionalGroupException(fgname + ": " + e.getMessage());
      			}
      			if (fgGraph == null)
                      throw new InvalidFunctionalGroupException(fgname);
      
      			FunctionalGroup fg = FunctionalGroup.make(fgname, fgGraph);
      			Object old = fgDictionary.get(fgname);
      			if (old == null) {
      				fgDictionary.put(fgname,fg);
      			}
      			else {
      				FunctionalGroup oldFG = (FunctionalGroup)old;
      				if (!oldFG.equals(fg)) throw new ReplaceFunctionalGroupException(fgname);
      			}
              }
      		line = ChemParser.readMeaningfulLine(data);
      	}
      
      	while (!unRead.isEmpty()) {
      		fgname = (String)(unRead.keySet().iterator().next());
      		ChemParser.findUnion(fgname,unRead,fgDictionary);
      	}
      
          in.close();
      	return;
      }
      catch (Exception e) {
      	throw new IOException(e.getMessage());
      }
  }
  
  //## operation readLibrary(String) 
  public void readLibrary(String p_fileName) throws IOException, InvalidKineticsFormatException, InvalidFunctionalGroupException {
      //#[ operation readLibrary(String) 
      try {
      	FileReader in = new FileReader(p_fileName);
      	BufferedReader data = new BufferedReader(in);
      
      	// step 1: read in kinetics format
      	String line = ChemParser.readMeaningfulLine(data);
      	if (line == null) throw new InvalidKineticsFormatException();
      	String format;
      	if (line.compareToIgnoreCase("Arrhenius") == 0) format = "Arrhenius";
      	else if (line.compareToIgnoreCase("Arrhenius_EP") == 0) format = "Arrhenius_EP";
      	else throw new InvalidKineticsFormatException("unknown rate constant type: " + line);
      
      	// step 2 read in content
      	kineticsTemplateLibrary = new KineticsTemplateLibrary();
      
      	int treeNum = getReactantTreeNumber();
      	int reactantNum = getReactantNumber();
      	int keyNum = Math.max(treeNum,reactantNum);
      	line = ChemParser.readMeaningfulLine(data);
      	while (line != null) {
              StringTokenizer token = new StringTokenizer(line);
              String ID = token.nextToken();
      
              // read in the names of functional group defining this kinetics
         		LinkedHashSet fgc = new LinkedHashSet();
      		for (int i=0; i<keyNum; i++) {
      			String r = token.nextToken().trim();
      	        Object fg = fgDictionary.get(r);
              	if (fg == null) {
              		throw new InvalidKineticsKeyException("unknown fg name: " + r);
              	}
                  else fgc.add(fg);
      		}
      		// read in the
      		/*
      		 * Commented out by MRH 15-Jun-2009
      		 * 	We used to pass only the TRange, modified Arrhenius parameters,
      		 * 		their uncertainties, and the data's rank to parseArrheniusKinetics.
      		 * 		Now, pass the entire line of data, so that any comments
      		 * 		written after the rank may be read in 
      		 */
//              String kinetics = token.nextToken();
//              while (token.hasMoreTokens()) {
//              	kinetics = kinetics + " " + token.nextToken();
//              }
      		Kinetics k;
//      		if (format.equals("Arrhenius")) k = ChemParser.parseArrheniusKinetics(kinetics);
//      		else if (format.equals("Arrhenius_EP")) k = ChemParser.parseArrheniusEPKinetics(kinetics);
//      		else throw new InvalidKineticsFormatException("Invalid rate constant format: " + kinetics);
      		if (format.equals("Arrhenius")) k = ChemParser.parseArrheniusKinetics(line,keyNum);
      		else if (format.equals("Arrhenius_EP")) k = ChemParser.parseArrheniusEPKinetics(line,keyNum);
      		else throw new InvalidKineticsFormatException("Invalid rate constant format: " + line);
      
          	/*
          	 * Added by MRH on 24-Jun-2009
          	 * Restoring feature for RMG to choose the best kinetics based on rank
          	 * 	and valid temperature range.  We 'get' the static Temperature variable
          	 * 	temp4BestKinetics from ReactionModelGenerator and pass it to a new
          	 * 	addKinetics function, that accepts the Temperature as an input 
          	 */
          	Temperature firstTempEncountered = ReactionModelGenerator.getTemp4BestKinetics();
      		kineticsTemplateLibrary.addKinetics(fgc,k,firstTempEncountered);
      		line = ChemParser.readMeaningfulLine(data);
       	}
      
          in.close();
          return;
      
      }
      catch (Exception e) {
      	throw new IOException("Error reading rate library. The error message is:" + '\n' + e.getMessage());
      }
      //#]
  }
  
  //## operation readReactionAdjList(String) 
  public String readReactionAdjList(String p_fileName) throws InvalidReactionAdjListFormatException, IOException, NotInDictionaryException {
      //#[ operation readReactionAdjList(String) 
      try {
      	FileReader in = new FileReader(p_fileName);
      	BufferedReader data = new BufferedReader(in);
      
      	// step 1: read in structure
      	String line = ChemParser.readMeaningfulLine(data);
      	if (line == null) throw new InvalidReactionAdjListFormatException();
      
      	StringTokenizer token = new StringTokenizer(line);
      	String fgName1 = token.nextToken();
      	int reactantNum = 1;
          String fgName2 = null;
          String seperator = token.nextToken();
          if (seperator.equals("+")) {
      	   	fgName2 = token.nextToken();
             	seperator = token.nextToken();
             	reactantNum++;
          }
      
          if (!seperator.equals("->")) throw new InvalidReactionAdjListFormatException();
          String pName1 = token.nextToken();
          int productNum = 1;
//          if (token.hasMoreTokens()) {
//          	seperator = token.nextToken();
//          	if (seperator.equals("+")) productNum++;
//          }
		/* Updated by MRH on 1-Jun-2009
			ReactionTemplate.java only looked for one or two products.
			A. Jalan found a rxn to place into the database with the following recipe:
				R1-H + R2-OOH -> R1 + R2-O + H2O
			RMG now looks for as many products as the user specifies in the recipe.
		*/
          while (token.hasMoreTokens()) {
        	  seperator = token.nextToken();
        	  if (seperator.equals("+")) productNum++;
          }
          // step 2 get reactant structure from dictionary to form structureTemplate
      	Matchable r1 = (Matchable)fgDictionary.get(fgName1);
      	if (r1 == null)	throw new NotInDictionaryException("During reading reactionAdjList: " + fgName1 + " is not found in functional group dictionary!");
      	Matchable r2 = null;
      	if (fgName2 != null) {
      		r2 = (Matchable)fgDictionary.get(fgName2);
      		if (r2 == null)	throw new NotInDictionaryException("During reading reactionAdjList: " + fgName2 + " is not found in functional group dictionary!");
      	}
      	structureTemplate = new StructureTemplate(r1,r2);
      
          // step 3 read in direction and corresponding infor about reverse and thermo consisitence
          line = ChemParser.readMeaningfulLine(data);
          line = line.trim();
          int direction = 0;
          if (line.compareToIgnoreCase("forward")==0)	direction = 1;
          else if (line.compareToIgnoreCase("thermo_consistence")==0) direction = 0;
          setDirection(direction);
          String reverseRT = null;
      
          if (direction == 1) {
          	// read in reverse reaction family's name
          	line = ChemParser.readMeaningfulLine(data).trim();
          	StringTokenizer dst = new StringTokenizer(line,":");
          	String sign = dst.nextToken().trim();
          	if (!sign.startsWith("reverse")) throw new InvalidReactionAdjListFormatException("Unknown reverse reaction family name!");
          	reverseRT = dst.nextToken().trim();
          }
      
      	// step 4: read in actions in reaction template
      	reactionAdjList = new ReactionAdjList(reactantNum,productNum);
      	while (line != null) {
      		while (!line.startsWith("Action")) line = ChemParser.readMeaningfulLine(data);
      		LinkedList actions = ChemParser.readActions(data);
      		reactionAdjList.setActions(actions);
      		line = ChemParser.readMeaningfulLine(data);
      	}
      	in.close();
      	return reverseRT;
      }
      catch (IOException e) {
      	throw new IOException("Reaction AdjList: " + e.getMessage());
      }
      catch (InvalidActionException e) {
      	throw new IOException("Reaction AdjList: " + e.getMessage());
      }
      
      
      
      
      //#]
  }
  
  //## operation readTree(String) 
  public void readTree(String p_fileName) throws IOException, NotInDictionaryException {
      //#[ operation readTree(String) 
      try {
      	FileReader in = new FileReader(p_fileName);
      	BufferedReader data = new BufferedReader(in);
      
          LinkedHashSet treeSet = new LinkedHashSet();
          HierarchyTree tree = ChemParser.readHierarchyTree(data, fgDictionary, 1);
          while (tree != null) {
          	treeSet.add(tree);
      	    tree = ChemParser.readHierarchyTree(data, fgDictionary, 1);
          }
          setReactantTree(treeSet);
          in.close();
          return;
      }
      catch (IOException e) {
      	throw new IOException("kinetics tree");
      }
      catch (NotInDictionaryException e) {
      	throw new NotInDictionaryException("kinetics tree reading: "+name+'\n'+ e.getMessage());
      }
      //#]
  }
  
  /**
  Check if the present reaction template is valid.  things to check:
  (1) if the structure template of this reaction template is valid
  (2) if the kinetics template of this reaction template is valid
  (3) if all the structure templates in kinetics are the subgraph of the structure template of this reaction template 
  (4) if the sum set of all the structure templates in kinetics is equal to the sturcture template of this reaction template
   
  */
  //## operation repOk() 
  public boolean repOk() {
      //#[ operation repOk() 
      return true;
      //#]
  }
  
  //## operation resetReactedSitesBondDissociation(LinkedList) 
  private void resetReactedSitesBondDissociation(LinkedList p_reactants) {
      //#[ operation resetReactedSitesBondDissociation(LinkedList) 
      if (p_reactants.size() != 2) return;
      Iterator iter = p_reactants.iterator();
      while (iter.hasNext()) {
      	Graph g = ((ChemGraph)iter.next()).getGraph();
      	Node n = g.getCentralNodeAt(2);
      	if (n != null) {
      		g.clearCentralNode();
      		g.setCentralNode(1,n);
      	}
      }
      return;
      
      
      
      
      //#]
  }
  
  //## operation setReactantTree(HashSet) 
  public void setReactantTree(LinkedHashSet p_treeSet) {
      //#[ operation setReactantTree(HashSet) 
      structureTemplate.setReactantTree(p_treeSet);
      
      
      
      
      
      //#]
  }
  
  //## operation toString() 
  public String toString() {
      //#[ operation toString() 
      String s = "Reaction Template Name: " + name + '\n';
      s = s + "Direction: ";
      if (direction == 1) s = s+"forward" +'\n';
      else if (direction == -1) s = s+"backward" +'\n';
      else s = s+"unknown" +'\n';
      
      return s;
      //#]
  }
  
  public int getDirection() {
      return direction;
  }
  
  public void setDirection(int p_direction) {
      direction = p_direction;
  }
  
  public HashMap getFgDictionary() {
      return fgDictionary;
  }
  
  public void setFgDictionary(LinkedHashMap p_fgDictionary) {
      fgDictionary = p_fgDictionary;
  }
  
  public String getName() {
      return name;
  }
  
  public void setName(String p_name) {
      name = p_name;
  }
  
  public KineticsTemplateLibrary getKineticsTemplateLibrary() {
      return kineticsTemplateLibrary;
  }
  
  public void setKineticsTemplateLibrary(KineticsTemplateLibrary p_KineticsTemplateLibrary) {
      kineticsTemplateLibrary = p_KineticsTemplateLibrary;
  }
  
  public ReactionAdjList getReactionAdjList() {
      return reactionAdjList;
  }
  
  public void setReactionAdjList(ReactionAdjList p_ReactionAdjList) {
      reactionAdjList = p_ReactionAdjList;
  }
  
  public ReactionTemplate getReverseReactionTemplate() {
      return reverseReactionTemplate;
  }
  
  public StructureTemplate getStructureTemplate() {
      return structureTemplate;
  }
  
  public void setStructureTemplate(StructureTemplate p_StructureTemplate) {
      structureTemplate = p_StructureTemplate;
  }
  
  public LinkedHashSet getPrimaryKinLibrary(){
      return PrimaryKineticLibrary.getReactionSet();
  }
  
  public boolean doesPrimaryKineticLibraryExist() {
	  if (PrimaryKineticLibrary.size() != 0) return true;
	  else return false;
  }
  
}
/*********************************************************************
	File Path	: RMG\RMG\jing\rxn\ReactionTemplate.java
*********************************************************************/

