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


import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import jing.chemUtil.Graph;
import jing.rxnSys.CoreEdgeReactionModel;

//## package jing::chem 

//----------------------------------------------------------------------------
// jing\chem\SpeciesDictionary.java                                                                  
//----------------------------------------------------------------------------

//## class SpeciesDictionary 
public class SpeciesDictionary {
    
    private static SpeciesDictionary INSTANCE = new SpeciesDictionary();		//## attribute INSTANCE 
    
    private static LinkedHashMap dictionary;		//## attribute dictionary 
    
    
    // Constructors
    
    //## operation SpeciesDictionary() 
    private  SpeciesDictionary() {
        //#[ operation SpeciesDictionary() 
        dictionary = new LinkedHashMap();
        //#]
    }
    
    //## operation getInstance() 
    public static SpeciesDictionary getInstance() {
        //#[ operation getInstance() 
        return INSTANCE;
        //#]
    }
    
    //## operation getSpecies(ChemGraph) 
    public static Species getSpecies(ChemGraph p_chemGraph) {
        //#[ operation getSpecies(ChemGraph) 
        return (Species)dictionary.get(p_chemGraph);
        
        
        
        //#]
    }
    
    //## operation getSpeciesFromID(int) 
    public static Species getSpeciesFromID(int p_id) {
        //#[ operation getSpeciesFromID(int) 
        if (p_id < 0) return null;
        
        Iterator iter = dictionary.values().iterator();
        while (iter.hasNext()) {
        	Species spe = (Species)iter.next();
        	if (spe.getID() == p_id) {
        		return spe;
        	}
        }
        
        return null;
        
        
        //#]
    }
    
    //## operation getSpeciesFromName(String) 
    public static Species getSpeciesFromName(String p_name) {
        //#[ operation getSpeciesFromName(String) 
        if (p_name == null) throw new NullPointerException();
        
        Iterator iter = dictionary.values().iterator();
        while (iter.hasNext()) {
        	Species spe = (Species)iter.next();
        	if (spe.getName().compareToIgnoreCase(p_name) == 0) {
        		return spe;
        	}
        }
        
        return null;
        
        
        
        //#]
    }
	
	public static Species getSpeciesFromGraph(Graph g) {
		if (g == null) throw new NullPointerException();
		Iterator iter = dictionary.keySet().iterator();
		while(iter.hasNext()) {
			ChemGraph cg = (ChemGraph)iter.next();
			if (cg.graph.isEquivalent(g)) {
				return (Species)dictionary.get(cg);
			}
		}
		return null;
	}
    
	 public static Species getSpeciesFromChemkinName(String p_name) {
	        //#[ operation getSpeciesFromName(String) 
	        if (p_name == null) throw new NullPointerException();
	        
	        Iterator iter = dictionary.values().iterator();
	        while (iter.hasNext()) {
	        	Species spe = (Species)iter.next();
	        	if (spe.getChemkinName().compareToIgnoreCase(p_name) == 0) {
	        		return spe;
	        	}
	        }
	        
	        return null;
	        
       
	        //#]
	    }
	
	    //## operation getSpeciesSetFromName(String) 
	    public HashSet getSpeciesSetFromName(String p_name) {
	        //#[ operation getSpeciesSetFromName(String) 
	        if (p_name == null) throw new NullPointerException();
	        
	        HashSet speSet = new HashSet();
	        Iterator iter = dictionary.values().iterator();
	        while (iter.hasNext()) {
	        	Species spe = (Species)iter.next();
	        	if (spe.getName().compareToIgnoreCase(p_name) == 0) {
	        		speSet.add(spe);
	        	}
	        }
	        
	        return speSet;
	        
	        //#]
	    }
		
    //## operation getSpeciesSet() 
    public HashSet getSpeciesSet() {
        //#[ operation getSpeciesSet() 
        return new HashSet(dictionary.values());
        //#]
    }
    
    //## operation putSpecies(Species) 
    public void putSpecies(Species p_species, boolean write) {
        //#[ operation putSpecies(Species) 
		String restartFileContent="";
		if (write){
			try{
				File allSpecies = new File ("Restart/allSpecies.txt");
				FileWriter fw = new FileWriter(allSpecies, true);
				//Species species = (Species) iter.next();
				restartFileContent = restartFileContent + p_species.getChemkinName() + " \n ";// + 0 + " (mol/cm3) \n";
				restartFileContent = restartFileContent + p_species.toString(1) + "\n\n";
				
				//restartFileContent += "\nEND";
				fw.write(restartFileContent);
				fw.close();
			}
			catch (IOException e){
				System.out.println("Could not write the restart edgespecies file");
	        	System.exit(0);
			}
		}
        if (p_species.hasResonanceIsomers()) {
        	Iterator iter = p_species.getResonanceIsomers();
        	while (iter.hasNext()) {
        		Object key = iter.next();
        		dictionary.put(key,p_species);
        	}
        }
        else {
        	dictionary.put(p_species.getChemGraph(),p_species);
        }
        return;
        
        
        //#]
    }
    
    //## operation remove(ChemGraph) 
    public void remove(ChemGraph p_chemGraph) {
        //#[ operation remove(ChemGraph) 
        if (p_chemGraph != null) dictionary.remove(p_chemGraph);
        //#]
    }
    
    //## operation size() 
    public int size() {
        //#[ operation size() 
        return dictionary.size();
        //#]
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\SpeciesDictionary.java
*********************************************************************/

