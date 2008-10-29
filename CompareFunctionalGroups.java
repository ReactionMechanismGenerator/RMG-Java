import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import jing.chem.ChemGraph;
import jing.chem.FunctionalGroup;
import jing.chem.InvalidFunctionalGroupException;
import jing.chemParser.ChemParser;
import jing.chemParser.InvalidGraphFormatException;
import jing.chemUtil.Graph;


public class CompareFunctionalGroups {

	 //## configuration RMG::RMG
    public static void main(String[] args) {
    	File funcGroups = new File ("FunctionalGroups.txt");
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
    		Graph fgGraph2 = null;
  			try {
  				fgGraph2 = ChemParser.readFGGraph(reader);
  			}
  			catch (InvalidGraphFormatException e) {
  				throw new InvalidFunctionalGroupException(fgname2 + ": " + e.getMessage());
  			}
  			FunctionalGroup fg2 = FunctionalGroup.make(fgname2, fgGraph2);
  			
  			boolean isSub = fgGraph2.isSub(fg1.getGraph());
  			//boolean isSub = fg2.isSub(fg1);
  			System.out.println("fg2 is a sub of fg1: "+isSub);
  			
    	}
    	catch (IOException e){
    		System.out.println(e.toString());
    	}
    	
    }
	
}
