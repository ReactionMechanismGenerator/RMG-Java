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
import jing.rxn.*;
import jing.rxnSys.*;
import java.util.*;
import jing.chemUtil.*;
import jing.chemParser.*;

//## package jing::chem

//----------------------------------------------------------------------------
// jing\chem\PrimaryThermoLibrary.java
//----------------------------------------------------------------------------


//svp, 8/26/05
//## class PrimaryThermoLibrary
class PrimaryThermoLibrary {

  protected HashMap library;

  protected HashMap dictionary;

  private static PrimaryThermoLibrary INSTANCE = new PrimaryThermoLibrary();		//## attribute INSTANCE

  //## operation PrimaryThermoLibrary()
  private PrimaryThermoLibrary(){
    //#[ operation PrimaryThermoLibrary()
    library = new HashMap();
    dictionary = new HashMap();

    String directory = System.getProperty("jing.chem.PrimaryThermoLibrary.pathName");
    String separator = System.getProperty("file.separator");
    String dictionaryFile = directory + "/Dictionary.txt";
    String libraryFile = directory + "/Library.txt";
    //System.out.println(dictionaryFile);
    //System.out.println(libraryFile);
    try{
      read(dictionaryFile, libraryFile);
    }
    catch (IOException e){
      System.out.println("Can't read primary thermo library files!");
    }
    //#]
 }


  //## operation getThermoData(ChemGraph)
  public ThermoData getThermoData(Graph p_graph){
    //#[ operation getThermoData(ChemGraph)
    ThermoData td = (ThermoData)library.get(p_graph);
    if (td != null){
      return td;
    }
    Iterator iter = library.keySet().iterator();
    while (iter.hasNext()){
      Graph g = (Graph)iter.next();
      if (g.isEquivalent(p_graph)){
        td = (ThermoData)library.get(g);
        return td;
      }
    }
    return null;
    //#]
  }

  //## operation read(String, String)
  public void read(String p_dictionary, String p_library) throws IOException, FileNotFoundException {
    //#[ operation read(String)
      dictionary = readDictionary(p_dictionary);
      library = readLibrary(p_library, dictionary);
      //#]
  }

  //## operation readDictionary(String)
  //svp
  public HashMap readDictionary(String p_fileName) throws FileNotFoundException, IOException{
    //#[ operation readDictionary(String)
    try{
      FileReader in = new FileReader(p_fileName);
      BufferedReader data = new BufferedReader(in);

      String line = ChemParser.readMeaningfulLine(data);

      read: while(line != null){
		  //System.out.println(line);
        StringTokenizer st = new StringTokenizer(line);
        String name = st.nextToken();

        data.mark(10000);
        line = ChemParser.readMeaningfulLine(data);
        if (line == null) break read;
        line = line.trim();
        data.reset();
        Graph graph = null;

          graph = ChemParser.readChemGraph(data);

        Object old = dictionary.get(name);
        if (old == null){
          dictionary.put(name, graph);
        }
        else{
          Graph oldGraph = (Graph)old;
          if (!oldGraph.equals(graph)) {
            System.out.println("Can't replace graph in primary thermo library!");
            System.exit(0);
          }

        }
        line = ChemParser.readMeaningfulLine(data);
      }
      in.close();
      return dictionary;
    }
    catch (FileNotFoundException e){
        throw new FileNotFoundException(p_fileName);
    }
    catch (IOException e){
      throw new IOException(p_fileName + ":" + e.getMessage());
    }
    //#]
  }


  //## operation readLibrary(String)
//svp
public HashMap readLibrary(String p_thermoFileName, HashMap p_dictionary) throws IOException {
  //#[ operation readLibrary(String)
  try{
    FileReader in = new FileReader(p_thermoFileName);
    BufferedReader data = new BufferedReader(in);
    HashMap library = new HashMap();
    String line = ChemParser.readMeaningfulLine(data);
    while (line != null){
      StringTokenizer token = new StringTokenizer(line);
      String name = token.nextToken();
      String thermo = token.nextToken();
      try{
        double H = Double.parseDouble(thermo);
        thermo = thermo.concat(" ");
        for (int i = 0; i < 11; i++){
          thermo = thermo.concat(token.nextToken());
          thermo = thermo.concat(" ");
        }
          ThermoData thermoData = ChemParser.parseThermoFromLibrary(thermo);
          String comments = "";
          while(token.hasMoreTokens()){
            comments = comments + " " + token.nextToken();
          }
          ThermoData newThermoData = new ThermoData(name, thermoData, comments);
          Graph g = (Graph)dictionary.get(name);
          Object old = library.get(g);
          if (old == null){
            library.put(g, newThermoData);
          }
          else{
            ThermoData oldThermoData = (ThermoData)old;
            if (!oldThermoData.equals(newThermoData)) {
              System.out.println("Can't replace thermoData in primary thermo library!");
              System.out.println("old thermodata:"+oldThermoData.getName());
              System.out.println("new thermodata:"+newThermoData.getName());
              System.exit(0);
            }

          }

        }

      catch (NumberFormatException e){
        Object o = p_dictionary.get(thermo);
        if (o == null){
          System.out.println(name + ": "+thermo);
        }
      }
      line = ChemParser.readMeaningfulLine(data);
    }
    in.close();
      return library;
  }
  catch (Exception e){
     throw new IOException("Can't read thermo in primary thermo library!");
  }
  //#]
}

protected static PrimaryThermoLibrary getINSTANCE() {
    return INSTANCE;
}


}
/*********************************************************************
        File Path	: RMG\RMG\jing\chem\PrimaryThermoLibrary.java
*********************************************************************/

