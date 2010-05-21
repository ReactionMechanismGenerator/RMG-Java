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
    
    private static LinkedList cache = new LinkedList();  // a shortlist of recently requested species used as a cache
  
    
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
        
        /* In its simplest form this routine is just: return (Species)dictionary.get(p_chemGraph);
           This implementation uses a cache consisting of the two previously gotten species
           this short list of two is checked before the full list (of potentially thousands) of species
           because very often the species being checked is the same as the penultimate species to be checked.
           This speeds things up considerably.  rwest 2009/05/08
           NB. putSpecies() has also been modified to update the cache when new species are added.
         */
        
        Species species = null;
        
        // check the cache
        Iterator iter = cache.iterator();
        while (iter.hasNext() && species==null) {
        	Species spe = (Species)iter.next();
            
            // see if 'spe' matches our p_chemGraph
            if (spe.hasResonanceIsomers()) {
                Iterator iter2 = spe.getResonanceIsomers();
                while (iter2.hasNext() && species==null) {
                    ChemGraph resonantisomer = (ChemGraph)iter2.next();
                    if (resonantisomer.equals(p_chemGraph)) species=spe;
                }
            }
        	else if (spe.getChemGraph().equals(p_chemGraph) ) {
                species=spe;
        	}
        }   
        // if species still = null then it wasn't in the cache
        
        // if we found it in the cache, then we want to move it to the end of the list so it appears as the most recent item
        if (species != null) {
            iter.remove(); // remove it from the cache (wherever it was located)
            cache.add(species); // and put it straight back again
            return species;
        }
        // if we are still here, we haven't found the species yet
       
        // now check the dictionary
        species = (Species)dictionary.get(p_chemGraph);
        
        
        // for diagnostic purposes, print stuff to a log file if we didn't get a cache hit
        /*
        try{
            File consideredSpecies = new File ("Restart/consideredSpecies.txt");
            FileWriter fw = new FileWriter(consideredSpecies, true);
            if (species==null) fw.write( "not in dictionary or cache:");
            else fw.write( "not in cache: ");
            fw.close();
        }
        catch (IOException e){
            System.out.println("Could not write the restart consideredSpecies file");
            System.exit(0);
        }
         */
        
        
        // update the cache
        if (species != null) {
            cache.add(species); // add to the end of the list
            if (cache.size() > 2)  
                cache.remove(); // remove from the front of the list
        }
        
        return species;
        
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
    
    public static Species getSpeciesFromNameID(String p_name) {
    	if (p_name == null) throw new NullPointerException();
        
        Iterator iter = dictionary.values().iterator();
        while (iter.hasNext()) {
        	Species spe = (Species)iter.next();
        	if ((spe.getName()+"("+spe.getID()+")").compareToIgnoreCase(p_name) == 0) {
        		return spe;
        	}
        }
        
        return null;
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
        
        /*  // we can't read the restart file, so no point writing it (which currently takes a LOT of time). rwest
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
        */
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
        
        // update the cache
        cache.add(p_species); // add to the end of the list
        if (cache.size() > 2) 
            cache.remove(); // remove from the front of the list

    }
    
    //## operation remove(ChemGraph) 
    public void remove(ChemGraph p_chemGraph) {
        if (p_chemGraph != null) dictionary.remove(p_chemGraph);
		// why don't we waint to throw an exception if we have a null pointer?
    }
	
    //remove all the mappings from ChemGraphs to a particular species
    public void remove(Species p_spe) {
		if (p_spe.hasResonanceIsomers()) {
			Iterator iter = p_spe.getResonanceIsomers();
			while(iter.hasNext()){
			    ChemGraph cg = (ChemGraph)iter.next();
			    cg.setSpecies(null);
			    remove(cg);
			}
		}
		else {
		    p_spe.getChemGraph().setSpecies(null);
		    remove( p_spe.getChemGraph() );
		}
		
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

