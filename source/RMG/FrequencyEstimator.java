import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.ListIterator;
import jing.chem.ChemGraph;
import jing.chem.ForbiddenStructureException;
import jing.chem.FrequencyGroups;
import jing.chem.InvalidChemGraphException;
import jing.chem.Species;
import jing.chem.SpectroscopicData;
import jing.chemParser.ChemParser;
import jing.chemUtil.Graph;
import jing.rxnSys.ReactionModelGenerator;

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
/**
 * This class can be used to estimate the frequencies for a set of ChemGraphs using Frankie. To use this class, pass on
 * the command line the path to a file containing a set of adjacency lists.
 * 
 * @author jwallen
 */
public class FrequencyEstimator {
    public static void main(String[] args) {
// initializeSystemProperties();
        RMG.globalInitializeSystemProperties();
        File GATPFit = new File(System.getProperty("RMG.GATPFitDir"));
        GATPFit.mkdir();
        File frankie = new File(System.getProperty("RMG.frankieOutputDir"));
        frankie.mkdir();
        LinkedList<ChemGraph> graphList = new LinkedList<ChemGraph>();
        LinkedList<String> nameList = new LinkedList<String>();
        LinkedHashMap speciesFromInputFile = new LinkedHashMap();
        File file = new File(args[0]);
        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));
            String line = ChemParser.readMeaningfulLine(reader, true);
            if (line.toLowerCase().startsWith("database")) {
                RMG.extractAndSetDatabasePath(line);
            } else {
                System.err.println("FrequencyEstimator: Could not"
                        + " locate the Database field");
                System.exit(0);
            }
            // Read adjacency lists from file until an exception is thrown
            line = ChemParser.readMeaningfulLine(reader, true);
            while (line != null) {
                nameList.add(line);
                Graph g = ChemParser.readChemGraph(reader);
                ChemGraph cg = ChemGraph.make(g);
                ReactionModelGenerator
                        .addChemGraphToListIfNotPresent_ElseTerminate(
                                speciesFromInputFile, cg, "");
                graphList.add(cg);
                line = ChemParser.readMeaningfulLine(reader, true);
            }
        } catch (InvalidChemGraphException e) {
            e.printStackTrace();
        } catch (ForbiddenStructureException e) {
            e.printStackTrace();
        } catch (IOException e) {
            System.out.println(e.toString());
        }
        FrequencyGroups freqGroups = FrequencyGroups.getINSTANCE();
        // // the following line won't compile, with error:
        // // generateFreqData(jing.chem.ChemGraph,jing.chem.ThermoData) in jing.chem.FrequencyGroups cannot be applied
// to (jing.chem.ChemGraph)
        // // rwest commenting out from the CVS version as it's blocking my automated build process
        // freqGroups.generateFreqData(cg2);
        int counter = 0;
        for (ListIterator<ChemGraph> iter = graphList.listIterator(); iter
                .hasNext();) {
            ChemGraph cg = iter.next();
            Species spec = Species.make(nameList.get(counter), cg);
            SpectroscopicData data = freqGroups.generateFreqData(spec);
            System.out.println("");
            System.out.println(nameList.get(counter));
            System.out.println(cg.toString());
            System.out.println(cg.getFreqComments());
            System.out.println("Data:");
            System.out.print("Vibrations (cm^-1): ");
            for (int i = 0; i < data.getVibrationCount(); i++)
                System.out.print(data.getVibration(i) + "\t");
            System.out.println();
            System.out.print("Rotations (cm^-1): ");
            for (int i = 0; i < data.getRotationCount(); i++)
                System.out.print(data.getRotation(i) + "\t");
            System.out.println();
            System.out.print("Hindered rotor frequencies (cm^-1): ");
            for (int i = 0; i < data.getHinderedCount(); i++)
                System.out.print(data.getHinderedFrequency(i) + "\t");
            System.out.println();
            System.out.print("Hindered rotor barriers (cm^-1): ");
            for (int i = 0; i < data.getHinderedCount(); i++)
                System.out.print(data.getHinderedBarrier(i) + "\t");
            System.out.println();
            ++counter;
        }
    }
}
