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

//## package jing::chem

//----------------------------------------------------------------------------
//jing\chem\ChemElement.java
//----------------------------------------------------------------------------

/**
OVERVIEW: Properties which all atoms of the same type share.
*/
//## class ChemElement
public class ChemElement {

  protected static ChemElementDictionary chemElementDictionary = ChemElementDictionary.getInstance();		//## attribute chemElementDictionary

  /**
  Name of the atom
  */
  protected String name;		//## attribute name

  /**
  Number of protons
  */
  protected int proton;		//## attribute proton

  /**
  Valency of the atom
  */
  protected double valency;		//## attribute valency

  /**
  Atomic weight of the atom
  */
  protected double weight;		//## attribute weight


  // Constructors

  //## operation ChemElement()
  private  ChemElement() {
      //#[ operation ChemElement()
      //#]
  }
  /**
  Requires:
  Effects: privite constructor called by create() to create new object
  Modifies:
  */
  //## operation ChemElement(String,int,double,double)
  private  ChemElement(String p_name, int p_proton, double p_valency, double p_weight) {
      //#[ operation ChemElement(String,int,double,double)
        name = p_name;
        proton = p_proton;
        valency = p_valency;
        weight = p_weight;



      //#]
  }

  /**
  Requires:
  Effects: call privite constructor to construct new ChemElement according to p_name.  Now, here we consider to create C(or carbon), H(or hydrogen), O(or oxygen), if pass in other names, throw UnknownSymbolException.  Notice: user can specify an chemElement list file to read in by ChemElementDictionary.  But this create method guarantees our system able to at least include C,H,O.
  Modifies: ChemElementDictionary
  */
  //## operation create(String)
  private static ChemElement create(String p_name) throws UnknownSymbolException {
      //#[ operation create(String)
      ChemElement ChemElement = null;

      if (p_name.equals("C")) {
      	ChemElement = new ChemElement("C",6,4,12);
      }
      else if (p_name.equals("Cl")){//svp
        ChemElement = new ChemElement("Cl",17,1,35.5);
      }
      else if (p_name.equals("H")) {
      	ChemElement = new ChemElement("H",1,1,1);
      }
      else if (p_name.equals("O")) {
      	ChemElement = new ChemElement("O",8,2,16);
      }
      // Added by MRH on 18-Jun-2009
      //	Hardcoding Si and S into RMG-java
      else if (p_name.equals("Si")) {
    	  ChemElement = new ChemElement("Si",14,4,28.086);
      }
      else if (p_name.equals("S")) {
    	  ChemElement = new ChemElement("S",16,2,32.064);
      }
      else {
      	throw new UnknownSymbolException("ChemElement");
      }
      return ChemElement;
      //#]
  }

  //## operation getName()
  public String getName() {
      //#[ operation getName()
      return name;
      //#]
  }

  //## operation getProton()
  public int getProton() {
      //#[ operation getProton()
      return proton;
      //#]
  }

  //## operation getValency()
  public double getValency() {
      //#[ operation getValency()
      return valency;
      //#]
  }

  //## operation getWeight()
  public double getWeight() {
      //#[ operation getWeight()
      return weight;
      //#]
  }

  /**
  Requires: pass in valid string name for elements
  Effects: search chemElementDictionary by p_name to see if this chemElement is already in the system.  if it is, return the instance of existing chemElement;  if it is not, call create(p_name) to create a new chemElement.
  Modifies:
  */
  // Argument String p_name :
  /**
  the name of chemical element you want to create.
  */
  //## operation make(String)
  public static ChemElement make(String p_name) throws UnknownSymbolException {
      //#[ operation make(String)
      try {
      	String internalName = translateName(p_name);
      	ChemElement ce = chemElementDictionary.getChemElement(internalName);
      	if (ce == null) {
      		ce = create(internalName);
      		chemElementDictionary.putChemElement(ce);
      	}
      	return ce;
      }
      catch (UnknownSymbolException e) {
      	throw new UnknownSymbolException("Chemical Element: " + p_name);
      }



      //#]
  }

  //## operation translateName(String)
  public static String translateName(String p_name) throws NullSymbolException, UnknownSymbolException {
      //#[ operation translateName(String)
      if (p_name == null) throw new NullSymbolException("ChemElement");

      if ((p_name.compareToIgnoreCase("C")==0) || (p_name.compareToIgnoreCase("Carbon")==0)) {
      	return "C";
      }
      else if ((p_name.compareToIgnoreCase("Cl")==0) || (p_name.compareToIgnoreCase("Chlorine")==0)){//svp
        return "Cl";
      }
      else if ((p_name.compareToIgnoreCase("H")==0) || (p_name.compareToIgnoreCase("Hydrogen")==0)) {
      	return "H";
      }
      else if ((p_name.compareToIgnoreCase("O")==0) || (p_name.compareToIgnoreCase("Oxygen")==0)) {
      	return "O";
      }
      // Added by MRH on 18-Jun-2009
      //	Hardcoding Si and S into RMG-java
      else if ((p_name.compareToIgnoreCase("Si")==0) || (p_name.compareToIgnoreCase("Silicon")==0)) {
    	  return "Si";
      }
      else if ((p_name.compareToIgnoreCase("S")==0) || (p_name.compareToIgnoreCase("Sulfur")==0) || (p_name.compareToIgnoreCase("Sulphur")==0)) {
    	  return "S";
      }
      else {
      	throw new UnknownSymbolException("ChemElement");
      }
      //#]
  }

  protected static ChemElementDictionary getChemElementDictionary() {
      return chemElementDictionary;
  }

}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\ChemElement.java
*********************************************************************/

