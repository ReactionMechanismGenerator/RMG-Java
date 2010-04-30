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
//jing\chem\Atom.java
//----------------------------------------------------------------------------

//## class Atom
public class Atom implements ChemNodeElement {

  protected static HashMap dictionary = new HashMap();		//## attribute dictionary

  protected String name;		//## attribute name

  protected ChemElement chemElement;
  protected FreeElectron freeElectron;

  // Constructors

  //## operation Atom()
  protected  Atom() {
      //#[ operation Atom()
      //#]
  }
  //## operation Atom(ChemElement,FreeElectron)
  protected  Atom(ChemElement p_chemElement, FreeElectron p_electron) throws InvalidFreeElectronException {
      //#[ operation Atom(ChemElement,FreeElectron)
      chemElement = p_chemElement;
      freeElectron = p_electron;

      name = p_chemElement.getName();

      if (p_electron != null) {
      	int order = p_electron.getOrder();
      	if (order == 0) {
      		name = name;
      	}
      	else if (order == 1) {
      		name = name + " mono rad";
      	}
      	else if (order == 2) {
      		name = name + " bi-rad " + p_electron.spin;
      	}
      	else if (order == 3) {
      		name = name + " tri-rad";
      	}
      	else if (order == 4) {
			name = name + "tetra-rad";
      	}
      	else {
      		throw new InvalidFreeElectronException();
      	}
      }
      else throw new InvalidFreeElectronException();




      //#]
  }

  //## operation changeRadical(int,String)
  public ChemNodeElement changeRadical(int p_radical, String p_spin) throws InvalidFreeElectronException {
      //#[ operation changeRadical(int,String)
      Atom newAtom = null;
      FreeElectron fe = getFreeElectron();
      if (fe == null) {
      	throw new InvalidFreeElectronException();
      }
      else {
      	int order = fe.getOrder()+p_radical;
      	if (order < 0) throw new InvalidFreeElectronException();
      	String name = String.valueOf(order);
      	if (order==2) {
      		// here, set a default value for 1-centered biradical to be 2T.
      		String spin = "T";
      		if (p_spin != null) spin = p_spin;
      		name = name + spin;
      	}
      	FreeElectron newfe = FreeElectron.make(name);
      	newAtom = Atom.make(getChemElement(),newfe);
      }
      return newAtom;


      //#]
  }

  //## operation geChemElementName()
  public String geChemElementName() {
      //#[ operation geChemElementName()
      return chemElement.getName();
      //#]
  }

  //## operation getRadicalNumber()
  public int getRadicalNumber() {
      //#[ operation getRadicalNumber()
      if (freeElectron == null) throw new InvalidFreeElectronException();
      else return freeElectron.getOrder();
      //#]
  }

  //gmagoon 4/30/10: modified version of getRadicalNumber; same as getRadicalNumber, except it will not count 2S as radical
  public int getUnpairedRadicalNumber() {
      //#[ operation getRadicalNumber()
      if (freeElectron == null) throw new InvalidFreeElectronException();
      else if (freeElectron.getName().equals("2S")) return 0;
      else return freeElectron.getOrder();
      //#]
  }

  //## operation getType()
  public String getType() {
      //#[ operation getType()
      return chemElement.getName();
      //#]
  }

  //## operation getValency()
  public double getValency() throws InvalidChemNodeElementException {
      //#[ operation getValency()
      double val = -1;
      if (chemElement != null) {
      	val = chemElement.getValency();
      	if (freeElectron !=null) {
      		val -= freeElectron.getOrder();
      	}
      	else throw new InvalidFreeElectronException();
      }
      else {
      	throw new InvalidChemNodeElementException();
      }

      return val;
      //#]
  }

  //## operation getWeight()
  public double getWeight() {
      //#[ operation getWeight()
      return chemElement.getWeight();
      //#]
  }

  //## operation isAny()
  public boolean isAny() {
      //#[ operation isAny()
      return false;
      //#]
  }

  //## operation isBiradical()
  public boolean isBiradical() {
      //#[ operation isBiradical()
      if (freeElectron == null) throw new InvalidFreeElectronException();
      if (freeElectron.getOrder() == 2) return true;
      else return false;


      //#]
  }

  //## operation isCarbon()
  public boolean isCarbon() {
      //#[ operation isCarbon()
      return (getType().equals("C"));
      //#]
  }

  //## operation isChlorine()
 public boolean isChlorine() {//svp
     //#[ operation isChlorine()
     return (getType().equals("Cl"));
     //#]
 }


  //## operation isHydrogen()
  public boolean isHydrogen() {
      //#[ operation isHydrogen()
      return (getType().equals("H"));
      //#]
  }

  //## operation isNonH()
  public boolean isNonH() {
      //#[ operation isNonH()
      return false;
      //#]
  }

  //## operation isOxygen()
  public boolean isOxygen() {
      //#[ operation isOxygen()
      return (getType().equals("O"));
      //#]
  }

  //## operation isRadical()
  public boolean isRadical() {
      //#[ operation isRadical()
      if (freeElectron == null) throw new InvalidFreeElectronException();
      if (freeElectron.getOrder() == 0) return false;
      else return true;



      //#]
  }
  
  // Added by MRH on 18-Jun-2009
  //	Hardcoding Si and S into RMG-java
  public boolean isSilicon() {
	  return (getType().equals("Si"));
  }
  
  public boolean isSulfur() {
	  return (getType().equals("S"));
  }

  //## operation make(ChemElement,FreeElectron)
  public static Atom make(ChemElement p_chemElement, FreeElectron p_freeElectron) {
      //#[ operation make(ChemElement,FreeElectron)
      Atom atom = new Atom(p_chemElement, p_freeElectron);
      Atom old = (Atom)(dictionary.get(atom.getName()));
      if (old == null) {
      	dictionary.put(atom.getName(), atom);
      	return atom;
      }
      else {
      	atom = null;
      	return old;
      }


      //#]
  }

  public static HashMap getDictionary() {
      return dictionary;
  }

  public static void setDictionary(HashMap p_dictionary) {
      dictionary = p_dictionary;
  }

  public String getName() {
      return name;
  }

  public void setName(String p_name) {
      name = p_name;
  }

  public ChemElement getChemElement() {
      return chemElement;
  }

  public FreeElectron getFreeElectron() {
      return freeElectron;
  }

}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\Atom.java
*********************************************************************/

