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

  //## operation make(String)
  public static FGElement make(String p_name) {
      //#[ operation make(String)
      try {
      	String internalName = translateName(p_name);
      	FGElement fge = fGElementDictionary.getFGElement(internalName);
      	if (fge == null) {
      		fge = FGElement.create(internalName);
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

