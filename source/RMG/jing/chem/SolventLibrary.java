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
/*
 * Created by Amrit Jalan on December 22, 2010
 */
package jing.chem;

import java.io.*;
import java.util.*;
import jing.chemUtil.*;
import jing.chemParser.*;
import jing.rxnSys.Logger;

/**
 * @author amrit
 */
public class SolventLibrary {
    protected static LinkedHashMap library;

    public SolventLibrary(String name, String location) {
        library = new LinkedHashMap();
        readSolventLibrary(name, location);
    }

// Reading in the Solvent library
    public void readSolventLibrary(String name, String path) {
        String solventLibraryDirectory = System
                .getProperty("jing.chem.SolventLibrary.pathName");
        String libraryFile = solventLibraryDirectory + "/" + path
                + "/Library.txt";
        try {
            read(libraryFile, name);
        } catch (IOException e) {
            String message = String.format(
                    "RMG cannot read Solvent Library %s: %s", name,
                    e.getMessage());
            throw new RuntimeException(message);
        }
    }

    public void read(String p_library, String p_name) throws IOException,
            FileNotFoundException {
        String source = "Solvent Library: " + p_name;
        Logger.info("Reading " + source);
        library = readLibrary(p_library, source);
    }

    public LinkedHashMap readLibrary(String p_transportFileName, String source)
            throws IOException {
        try {
            FileReader in = new FileReader(p_transportFileName);
            BufferedReader data = new BufferedReader(in);
            LinkedHashMap tempLibrary = new LinkedHashMap();
            String line = ChemParser.readMeaningfulLine(data, true);
            while (line != null) {
                StringTokenizer token = new StringTokenizer(line);
                String name = token.nextToken();
                String abram = token.nextToken();
                try {
                    SolventData solventData = new SolventData(
                            Double.parseDouble(abram), Double.parseDouble(token
                                    .nextToken()), Double.parseDouble(token
                                    .nextToken()), Double.parseDouble(token
                                    .nextToken()), Double.parseDouble(token
                                    .nextToken()), Double.parseDouble(token
                                    .nextToken()), Double.parseDouble(token
                                    .nextToken()), Double.parseDouble(token
                                    .nextToken()), Double.parseDouble(token
                                    .nextToken()), Double.parseDouble(token
                                    .nextToken()), Double.parseDouble(token
                                    .nextToken()), Double.parseDouble(token
                                    .nextToken()), Double.parseDouble(token
                                    .nextToken()));
                    String comments = "";
                    while (token.hasMoreTokens()) {
                        comments = comments + " " + token.nextToken();
                    }
                    SolventData newSolventData = new SolventData(name,
                            solventData, comments);
                    tempLibrary.put(name, newSolventData);
                    library.put(name, newSolventData);
                } catch (NumberFormatException e) {
                }
                line = ChemParser.readMeaningfulLine(data, true);
            }
            in.close();
            return library;
        } catch (Exception e) {
            Logger.logStackTrace(e);
            throw new IOException(
                    "Can't read solvent descriptors in solvent library!");
        }
    }

    public static SolventData getSolventData(String solv_name) {
        if (library == null)
            return null;
        SolventData solv = (SolventData) library.get(solv_name);
        if (solv != null) {
            return solv;
        }
        Iterator iter = library.keySet().iterator();
        while (iter.hasNext()) {
            String name = (String) iter.next();
            if (name.equals(solv_name)) {
                solv = (SolventData) library.get(name);
                return solv;
            }
        }
        if (solv == null) {
            Logger.error("Solvent does not exist in RMG solvent library");
            Logger.error("Stopping.");
            throw new Error();
        }
        return null;
    }
}
