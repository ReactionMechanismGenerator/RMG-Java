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
//jing\chem\FGElement.java
//----------------------------------------------------------------------------

//## class FGElement
public class FGElement {

  protected static FGElementDictionary fGElementDictionary = FGElementDictionary.getInstance();		//## attribute fGElementDictionary

  protected String name;		//## attribute name

  protected String type;		//## attribute type


  // Constructors

  //## operation FGElement()
  protected  FGElement() {
      //#[ operation FGElement()
      //#]
  }
  //## operation FGElement(String)
  protected  FGElement(String p_name) throws UnknownSymbolException {
      //#[ operation FGElement(String)
      name = p_name;

      if (p_name.startsWith("H") || p_name.startsWith("h")) {
      	type = "H";
      }
      else if (p_name.startsWith("Cl") || p_name.startsWith("cl") || p_name.startsWith("Chl") || p_name.startsWith("chl")){//svp
        type = "Cl";
      }
      else if (p_name.startsWith("C") || p_name.startsWith("c")) {
      	type = "C";
      }
      else if (p_name.startsWith("O") || p_name.startsWith("o")) {
      	type = "O";
      }
      else if (p_name.compareToIgnoreCase("R") == 0) {
      	type = "R";
      }
      else if (p_name.compareToIgnoreCase("R!H") == 0) {
      	type = "R!H";
      }
      /*
       * WARNING: The elseif statement involving Si must precede the elseif
       * 	statement involving S.  If not, a "Si" would be recognized as a
       * 	"S", at least with the current syntax.  This is also true for
       * 	the "Cl"/"C" elseif statements above
       */
      else if (p_name.startsWith("Si") || p_name.startsWith("si")) {
    	  type = "Si";
      }
      else if (p_name.startsWith("S") || p_name.startsWith("s")) {
    	  type = "S";
      }
      else {
      	throw new UnknownSymbolException("FGElement type: " + p_name);
      }




      //#]
  }

  //## operation create(String)
  public static FGElement create(String p_name) {
      //#[ operation create(String)
      FGElement fge = null;
      if (p_name.equals("H")) {
      	fge = new FGElement("H");
      }
      else if (p_name.equals("Cs")) {
      	fge = new FGElement("Cs");
      }
      else if (p_name.equals("Cd")) {
      	fge = new FGElement("Cd");
      }
      else if (p_name.equals("Cdd")) {
      	fge = new FGElement("Cdd");
      }
      else if (p_name.equals("Ct")) {
      	fge = new FGElement("Ct");
      }
      else if (p_name.equals("Cb")) {
      	fge = new FGElement("Cb");
      }
      else if (p_name.equals("Cbf")) {
      	fge = new FGElement("Cbf");
      }
      else if (p_name.equals("CO")) {
      	fge = new FGElement("CO");
      }
      else if (p_name.equals("Cl")){//svp
        fge = new FGElement("Cl");
      }
      else if (p_name.equals("Os")) {
      	fge = new FGElement("Os");
      }
      else if (p_name.equals("Oa")) {
      	fge = new FGElement("Oa");
      }
      else if (p_name.equals("Od")) {
      	fge = new FGElement("Od");
      }
      else if (p_name.equals("R")) {
      	fge = new FGElement("R");
      }
      else if (p_name.equals("R!H")) {
      	fge = new FGElement("R!H");
      }
      else if (p_name.equals("Sis")) {
    	  fge = new FGElement("Sis");
      }
      else if (p_name.equals("Sid")) {
    	  fge = new FGElement("Sid");
      }
      else if (p_name.equals("Sidd")) {
    	  fge = new FGElement("Sidd");
      }
      else if (p_name.equals("Sit")) {
    	  fge = new FGElement("Sit");
      }
      else if (p_name.equals("Ss")) {
    	  fge = new FGElement("Ss");
      }
      else if (p_name.equals("Sa")) {
    	  fge = new FGElement("Sa");
      }
      else if (p_name.equals("Sd")) {
    	  fge = new FGElement("Sd");
      }
      else throw new UnknownSymbolException("FGElement: " + p_name);

      return fge;
      //#]
  }

  //## operation isAny()
  public boolean isAny() {
      //#[ operation isAny()
      return (getType().equals("R"));
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
      return (getType().equals("R!H"));
      //#]
  }

  //## operation isOxygen()
  public boolean isOxygen() {
      //#[ operation isOxygen()
      return (getType().equals("O"));
      //#]
  }
  
  public boolean isSilicon() {
      return (getType().equals("Si"));
  }
  
  public boolean isSulfur() {
      return (getType().equals("S"));
  }

  //## operation make(String)
  public static FGElement make(String p_name) {
      //#[ operation make(String)
      try {
      	//String internalName = translateName(p_name);
      	//FGElement fge = fGElementDictionary.getFGElement(internalName);
		FGElement fge = fGElementDictionary.getFGElement(p_name);
      	
		if (fge == null) {
      		//fge = FGElement.create(internalName);
			fge = FGElement.create(p_name);
			fGElementDictionary.putFGElement(fge);
      	}
      	return fge;
      }
      catch (UnknownSymbolException e) {
      	throw new UnknownSymbolException("FG Element: " + p_name);
      }


      //#]
  }

  //## operation translateName(String)
  public static String translateName(String p_name) {
      //#[ operation translateName(String)
      if (p_name == null) throw new NullSymbolException("FGElement");

      if ((p_name.compareToIgnoreCase("H")==0) || (p_name.compareToIgnoreCase("Hydrogen")==0)) {
      	return "H";
      }
      else if ((p_name.compareToIgnoreCase("Cs")==0)) {
      	return "Cs";
      }
      else if ((p_name.compareToIgnoreCase("Cd")==0)) {
      	return "Cd";
      }
      else if ((p_name.compareToIgnoreCase("Cdd")==0)) {
      	return "Cdd";
      }
      else if ((p_name.compareToIgnoreCase("Ct")==0)) {
      	return "Ct";
      }
      else if ((p_name.compareToIgnoreCase("Cb")==0)) {
      	return "Cb";
      }
      else if ((p_name.compareToIgnoreCase("Cbf")==0)) {
      	return "Cbf";
      }
      else if ((p_name.compareToIgnoreCase("CO")==0)) {
      	return "CO";
      }
      else if ((p_name.compareToIgnoreCase("Cl")==0)){//svp
        return "Cl";
      }
      else if ((p_name.compareToIgnoreCase("Os")==0)) {
      	return "Os";
      }
      else if ((p_name.compareToIgnoreCase("Oa")==0)) {
      	return "Oa";
      }
      else if ((p_name.compareToIgnoreCase("Od")==0)) {
      	return "Od";
      }
      else if (p_name.compareToIgnoreCase("R")==0) {
      	return "R";
      }
      else if (p_name.compareToIgnoreCase("R!H")==0 || p_name.compareToIgnoreCase("R|H")==0) {
      	return "R!H";
      }
      else if ((p_name.compareToIgnoreCase("Sis")==0)) {
    	  return "Sis";
      }
      else if ((p_name.compareToIgnoreCase("Sid")==0)) {
    	  return "Sid";
      }
      else if ((p_name.compareToIgnoreCase("Sidd")==0)) {
    	  return "Sidd";
      }
      else if ((p_name.compareToIgnoreCase("Sit")==0)) {
    	  return "Sit";
      }
      else if ((p_name.compareToIgnoreCase("Ss")==0)) {
    	  return "Ss";
      }
      else if ((p_name.compareToIgnoreCase("Sa")==0)) {
        	return "Sa";
      }
      else if ((p_name.compareToIgnoreCase("Sd")==0)) {
    	  return "Sd";
      }
      else {
      	throw new UnknownSymbolException("FGElement");
      }
      //#]
  }

  public static FGElementDictionary getFGElementDictionary() {
      return fGElementDictionary;
  }

  public static void setFGElementDictionary(FGElementDictionary p_fGElementDictionary) {
      fGElementDictionary = p_fGElementDictionary;
  }

  public String getName() {
      return name;
  }

  public void setName(String p_name) {
      name = p_name;
  }

  public String getType() {
      return type;
  }

  public void setType(String p_type) {
      type = p_type;
  }

}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\FGElement.java
*********************************************************************/

