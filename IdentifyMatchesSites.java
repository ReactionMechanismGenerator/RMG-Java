import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;

import jing.chem.ChemGraph;
import jing.chem.ForbiddenStructureException;
import jing.chem.FunctionalGroup;
import jing.chem.InvalidChemGraphException;
import jing.chem.InvalidFunctionalGroupException;
import jing.chemParser.ChemParser;
import jing.chemParser.InvalidGraphFormatException;
import jing.chemUtil.Graph;
import jing.chemUtil.MatchedSite;


public class IdentifyMatchesSites {
	 public static void main(String[] args) {
	    	File funcGroups = new File ("IdentifySites.txt");
	    	try{
	    		BufferedReader reader = new BufferedReader(new FileReader(funcGroups));
	    		String line = ChemParser.readMeaningfulLine(reader);
	    		String fgname1 = line;
	    		Graph fgGraph1 = null;
	  			try {
	  				fgGraph1 = ChemParser.readFGGraph(reader);
	  			}
	  			catch (InvalidGraphFormatException e) {
	  				throw new InvalidFunctionalGroupException(fgname1 + ": " + e.getMessage());
	  			}
	  			FunctionalGroup fg1 = FunctionalGroup.make(fgname1, fgGraph1);
	  			
	  			line = ChemParser.readMeaningfulLine(reader);
	    		
	  			String fgname2 = line;
	    		Graph Graph2 = null;
	  			try {
	  				Graph2 = ChemParser.readChemGraph(reader);
	  			}
	  			catch (InvalidGraphFormatException e) {
	  				throw new InvalidFunctionalGroupException(fgname2 + ": " + e.getMessage());
	  			}
	  			ChemGraph cg2 = null;
				try {
					cg2 = ChemGraph.make(Graph2);
				} catch (InvalidChemGraphException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (ForbiddenStructureException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
	  			
	  			HashSet matchedSites = cg2.identifyReactionMatchedSite(fg1);
	  			//boolean isSub = cg2.isSubAtCentralNodes(fg1);
	  			//System.out.println("fg2 is a sub of fg1: "+ isSub);
	  			Iterator iter = matchedSites.iterator();
	  			while (iter.hasNext()) {
	  				MatchedSite o = (MatchedSite) iter.next();
	  				System.out.println("fg2 is a sub of fg1: "+ o.toString());
	  				
	  			}
	  			//System.out.println("fg2 is not sub on fg1");
	  			
	    	}
	    	catch (IOException e){
	    		System.out.println(e.toString());
	    	}
	    	
	    }
}
