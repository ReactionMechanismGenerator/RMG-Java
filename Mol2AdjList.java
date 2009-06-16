import jing.chem.Species;

public class Mol2AdjList {

	/**
	 * Mol2ChemGraph converts a .mol file to an RMG adjacency list.
	 * 	This program was originally created because the InChI executable was
	 * 	not stable on Linux systems.  Rather than go from the InChI string to
	 * 	the adjacency list in "one step" (via the poorly named InChI2ChemGraph
	 * 	function), the user now supplies the .mol file; creating a .mol file
	 * 	from the InChI on Linux is normally stable.
	 */
	public static void main(String[] args) {
		String cgString = Species.mol2AdjList(args[0], args[0]);
		System.out.println(cgString);
	}

}
