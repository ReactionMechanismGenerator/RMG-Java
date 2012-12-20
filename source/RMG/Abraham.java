// //////////////////////////////////////////////////////////////////////////////
//
// RMG - Reaction Mechanism Generator
//
// Copyright (c) 2002-2011 Prof. William H. Green (whgreen@mit.edu) and the
// RMG Team (rmg_dev@mit.edu)
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// //////////////////////////////////////////////////////////////////////////////
import java.util.*;
import java.io.*;
import jing.chem.*;
import jing.chemParser.*;
import jing.param.*;
import jing.chemUtil.*;
// import bondGroups.*;
import jing.rxn.*;
import jing.rxnSys.*;
import jing.mathTool.*;

// Amrit Jalan
// June 18, 2009
// Returns the Abraham Parameters for the list of chemgraphs in Abraham_input.txt for the solvent parameters defined at
// the head of Abraham_input.txt
public class Abraham {
    // ## configuration RMG::RMG
    public static void main(String[] args) {
// initializeSystemProperties();
        RMG.globalInitializeSystemProperties();
        LinkedHashSet speciesSet = new LinkedHashSet();
        String abraham_output = "";
        double c = 0;
        double a = 0;
        double b = 0;
        double l = 0;
        double s = 0;
        double es = 0;
        String solvent = "";
        LinkedHashMap speciesFromInputFile = new LinkedHashMap();
        try {
            FileReader in = new FileReader("Abraham_input.txt");
            BufferedReader data = new BufferedReader(in);
            // Read the first line of Abraham_input.txt
            String line = ChemParser.readMeaningfulLine(data, true);
            StringTokenizer st = new StringTokenizer(line);
            // The first line should start with "SolventParameters", otherwise do nothing and display a message to the
// user
            if (st.nextToken().startsWith("SolventParameters")) {
                c = Double.parseDouble(st.nextToken());
                es = Double.parseDouble(st.nextToken());
                s = Double.parseDouble(st.nextToken());
                a = Double.parseDouble(st.nextToken());
                b = Double.parseDouble(st.nextToken());
                l = Double.parseDouble(st.nextToken());
                solvent = st.nextToken();
            } else
                System.out
                        .println("Error in reading Abraham_input.txt file:\nThe first line must read 'SolventParameters:'.");
            // Read in the ChemGraphs and compute their thermo, while there are ChemGraphs to read in
            line = ChemParser.readMeaningfulLine(data, true);
            while (line != null) {
                String speciesName = line;
                Graph g = ChemParser.readChemGraph(data);
                ChemGraph cg = null;
                try {
                    cg = ChemGraph.make(g);
                } catch (ForbiddenStructureException e) {
                    System.out
                            .println("Error in reading graph: Graph contains a forbidden structure.\n"
                                    + g.toString());
                    System.exit(0);
                }
                ReactionModelGenerator
                        .addChemGraphToListIfNotPresent_ElseTerminate(
                                speciesFromInputFile, cg, speciesName);
                Species species = Species.make(speciesName, cg);
                speciesSet.add(species);
                line = ChemParser.readMeaningfulLine(data, true);
            }
            in.close();
            abraham_output += "Solvent:" + "\t" + solvent + "\n\n";
            // abraham_output += "Name" + "\t" +"S"+"\t"+"B"+"\t"+"E"+"\t"+"L"+"\t"+"A" + "\n";
            Iterator iter = speciesSet.iterator();
            while (iter.hasNext()) {
                Species spe = (Species) iter.next();
                ChemGraph cg = spe.getChemGraph();
                // Generation of Abraham Solute Parameters
                AbramData result_Abraham = new AbramData();
                result_Abraham = cg.getAbramData();
                // Solute descriptors from the Abraham Model
                double S = result_Abraham.S;
                double B = result_Abraham.B;
                double E = result_Abraham.E;
                double L = result_Abraham.L;
                double A = result_Abraham.A;
                double logK = c + s * S + b * B + es * E + l * L + a * A; // Implementation of Abraham Model
                // abraham_output += spe.getName() + "\t" +S+"\t"+B+"\t"+E+"\t"+L+"\t"+A + "\n";
                abraham_output += spe.getName() + "\t" + logK + "\n";
            }
            try {
                File abrahamOutput = new File("Abraham_output.txt");
                FileWriter fw = new FileWriter(abrahamOutput);
                fw.write(abraham_output);
                fw.close();
            } catch (IOException e) {
                System.out.println("Error in writing Abraham_output.txt file.");
                System.exit(0);
            }
        } catch (FileNotFoundException e) {
            System.err.println("File was not found!\n");
        } catch (IOException e) {
            System.err.println("Something wrong with ChemParser.readChemGraph");
        }
        System.out.println("Done!\n");
    };
}
