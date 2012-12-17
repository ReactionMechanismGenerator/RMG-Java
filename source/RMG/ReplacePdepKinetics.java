/**
 * ReplacePdepKinetics.class Input: [0] chem.inp file with CHEB p-dep nomenclature [1] options ARH for Arrhenius and
 * PLOG for PLOG mechanism [2] number (the pressure of interest) required for ARH option, for PLOG option the 4
 * pressures (0.1,1,10,100 atm) are hard coded Output: chem.inp file with Arrhenius/PLOG nomenclature This program takes
 * a CHEMKIN chem.inp file as input. All lines until the REACTIONS line are read-in, saved, and immediately skipped
 * over. Once the REACTIONS line is read-in, the program looks for the string "1.0E0 0.0 0.0". NOTE: MRH chose this
 * string as the indication that a CHEB p-dep rate coefficient is present. As of ~ May 2011, RMG would produce this
 * string for all CHEB p-dep kinetics. If this has changed, please update this code accordingly. Upon finding this
 * string, the program reads in the relevant parameters (including the TMAX, TMIN, PMAX, PMIN, nT, nP, and the CHEB
 * coefficients). These are sent to a function that computes the k(T,P) for a given pressure (the pressure of interest
 * supplied by the user); the temperatures used are spaced linearly on an inverse temperature space, spanning 1/TMAX to
 * 1/TMIN, with a total of 8 temperatures. NOTE: The "8" is presently hard-coded. The 8 k(T,P_of_interest) are then fit
 * via least-squares to a modified Arrhenius format. Finally, the entire CHEB formatted rate coefficient in the chem.inp
 * file is replaced with the computed modified Arrhenius format, including adding a comment that states these parameters
 * were fit from a CHEB at a certain pressure. This program takes ~10 minutes to re-cast a 3.5MB chem.inp file with
 * >7,000 CHEB rate coefficients into Arrhenius form.
 */
import Jama.Matrix;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;
import jing.chemParser.ChemParser;
import jing.param.GasConstant;

public class ReplacePdepKinetics {
    public static void main(String args[]) {
        String plogChemInpFile = "";
        String currentLine = "";
        // Read in the current chem.inp file
        try {
            FileReader fr_chemInpFile = new FileReader(args[0]);
            BufferedReader br_chemInpFile = new BufferedReader(fr_chemInpFile);
            String line = readMeaningfulAndEmptyLine(br_chemInpFile, false);
            plogChemInpFile += line + "\n";
            // Read in chem.inp until REACTIONS is found
            while (!line.startsWith("REACTIONS")) {
                line = readMeaningfulAndEmptyLine(br_chemInpFile, false);
                plogChemInpFile += line + "\n";
            }
            // Skip over the REACTIONS line
            line = readMeaningfulAndEmptyLine(br_chemInpFile, false);
            plogChemInpFile += line + "\n";
            // Find all instances of Chebyshev polynomials
            while (line != null && !line.equals("END")) {
                // Debugging Line to see which lines have been read in
                System.out.println(line);
                if (!line.startsWith("!")
                        && !line.toLowerCase().contains("low")
                        && !line.toLowerCase().contains("troe")
                        && !line.toLowerCase().contains("dup")
                        && !line.toLowerCase().contains("plog")
                        && !line.contains("CHEB")) {
                    // Find all Chebyshev rate coefficients
                    if (line.contains("1.0E0 0.0 0.0")) {
                        // Split the reaction string on (+m), this is done to separate the reactions and products in
// chebychev reaction
                        String[] wholeLine = line.split("\\(\\+m\\)");
                        if (args[1].equals("PLOG")) {
                            plogChemInpFile += wholeLine[0] + wholeLine[1]
                                    + "\t" + wholeLine[2] + "\n";
                        } else if (args[1].equals("ARH")) {
                            currentLine = wholeLine[0] + wholeLine[1] + "\t";
                        }
                        line = ChemParser.readUncommentLine(br_chemInpFile);
                        // Some pdep kinetics comments extend to multiple lines
                        while (line.startsWith("!")) {
                            line = ChemParser.readUncommentLine(br_chemInpFile);
                        }
                        // Read in the TCHEB and PCHEB information
                        StringTokenizer st_cheb = new StringTokenizer(line, "/");
                        String currentToken = st_cheb.nextToken(); // TCHEB
                        StringTokenizer temperatureInfo = new StringTokenizer(
                                st_cheb.nextToken());
                        double TMIN = Double.parseDouble(temperatureInfo
                                .nextToken());
                        double TMAX = Double.parseDouble(temperatureInfo
                                .nextToken());
                        currentToken = st_cheb.nextToken(); // PCHEB
                        StringTokenizer pressureInfo = new StringTokenizer(
                                st_cheb.nextToken());
                        double PMIN = Double.parseDouble(pressureInfo
                                .nextToken());
                        double PMAX = Double.parseDouble(pressureInfo
                                .nextToken());
                        // Read the line containing the nT and nP (basis set) info
                        line = ChemParser.readUncommentLine(br_chemInpFile);
                        st_cheb = new StringTokenizer(line, "/");
                        currentToken = st_cheb.nextToken();
                        StringTokenizer st_values = new StringTokenizer(
                                st_cheb.nextToken());
                        int nT = Integer.parseInt(st_values.nextToken());
                        int nP = Integer.parseInt(st_values.nextToken());
                        // Extract the coefficients
                        double[] coeffs = new double[nT * nP];
                        int coeffCounter = 0;
                        while (coeffCounter < nT * nP) {
                            line = ChemParser.readUncommentLine(br_chemInpFile);
                            String[] splitSlashes = line.split("/");
                            StringTokenizer st_coeffs = new StringTokenizer(
                                    splitSlashes[1]);
                            while (st_coeffs.hasMoreTokens()) {
                                coeffs[coeffCounter] = Double
                                        .parseDouble(st_coeffs.nextToken()
                                                .trim());
                                ++coeffCounter;
                            }
                        }
                        // Replace these kinetics with pdep arrhenius kinetics
                        if (args[1].equals("PLOG")) {
                            double[] P = new double[4];
                            P[0] = 0.1;
                            P[1] = 1;
                            P[2] = 10;
                            P[3] = 100;
                            for (int i = 0; i < P.length; i++) {
                                plogChemInpFile += convertCHEB2PLOG(coeffs, nT,
                                        nP, TMIN, TMAX, PMIN, PMAX, args[1],
                                        P[i])
                                        + "\n";
                                // System.out.println(convertCHEB2PLOG(coeffs,nT,nP,TMIN,TMAX,PMIN,PMAX,args[1],P[i]));
                            }
                        } else if (args[1].equals("ARH")) {
                            currentLine += convertCHEB2PLOG(coeffs, nT, nP,
                                    TMIN, TMAX, PMIN, PMAX, args[1],
                                    Double.parseDouble(args[2]));
                            plogChemInpFile += currentLine
                                    + "\t!Arrhenius fit to " + args[2] + " atm";
                        }
                    } else {
                        plogChemInpFile += line + "\n";
                    }
                } else {
                    plogChemInpFile += line + "\n";
                }
                line = readMeaningfulAndEmptyLine(br_chemInpFile, false);
            }
        } catch (FileNotFoundException e) {
            System.err.println("File was not found: " + args[0] + "\n");
        }
        // Write the output files
        try {
            String filename = null;
            if (args[1].equals("ARH")) {
                filename = new String("ARH" + args[2] + "_");
            } else if (args[1].equals("PLOG")) {
                filename = new String("PLOG" + "_");
            }
            File rxns = new File(filename + args[0]);
            FileWriter fw_rxns = new FileWriter(rxns);
            fw_rxns.write(plogChemInpFile);
            fw_rxns.close();
        } catch (IOException e) {
            System.out.println("Could not write new (CHEB-free) chem.inp file");
            System.exit(0);
        }
        System.out.println("Complete.");
    }

    public static String readMeaningfulAndEmptyLine(BufferedReader p_reader,
            boolean toTrimOrNotToTrim) {
        if (p_reader == null)
            return null;
        String line = null;
        try {
            do {
                line = p_reader.readLine();
                if (line == null)
                    return null;
                if (toTrimOrNotToTrim)
                    line = line.trim();
            } while (line.startsWith("//"));
            return line;
        } catch (IOException e) {
            return null;
        }
    }

    public static String convertCHEB2PLOG(double[] a, int nT, int nP,
            double Tlow, double Thigh, double Plow, double Phigh,
            String option, double P_of_interest) {
        String string2return = "";
        double[] T = new double[8];
        double[] Tbar = new double[T.length];
        double[] P = new double[1];
        double[] Pbar = new double[P.length];
        // Getting the normalized cheb T points
        for (int i = 0; i < T.length; i++) {
            // As cheb do not work well for T > 2000, we will evaluate in range of 500 to 2000
            double Thigh_plog = 2000;
            double Tlow_plog = 500;
            T[i] = 1 / ((1 / Thigh_plog) + i * (1 / Tlow_plog - 1 / Thigh_plog)
                    / (T.length - 1));
            Tbar[i] = (2 / T[i] - 1 / Tlow - 1 / Thigh)
                    / (1 / Thigh - 1 / Tlow);
        }
        // Getting normalized cheb P points
        P[0] = P_of_interest;
        Pbar[0] = (2 * Math.log10(P[0]) - Math.log10(Plow) - Math.log10(Phigh))
                / (Math.log10(Phigh) - Math.log10(Plow));
        // Evaluating k_cheb(Tcheb,Pcheb) to get k(T,P) values
        double[][] logk = new double[T.length][P.length];
        double[][] k_rxn = new double[T.length][P.length];
        for (int j = 0; j < P.length; j++) {
            for (int i = 0; i < T.length; i++) {
                for (int k = 0; k < nT; k++) {
                    for (int l = 0; l < nP; l++) {
                        // This gets us the log10_k values
                        logk[i][j] += a[k * nP + l] * phi(Tbar[i], k)
                                * phi(Pbar[j], l);
                    }
                }
                // This will give us the actual k values
                k_rxn[i][j] = Math.pow(10, logk[i][j]);
            }
            // Construct matrix X and vector y
            double[][] y = new double[T.length][1];
            double[][] X = new double[T.length][3];
            for (int m = 0; m < T.length; m++) {
                y[m][0] = Math.log(k_rxn[m][j]);
                X[m][0] = 1;
                X[m][1] = Math.log(T[m]);
                X[m][2] = 1 / T[m];
            }
            // Solve least-squares problem using inv(XT*X)*(XT*y)
            Matrix X_matrix = new Matrix(X);
            Matrix y_matrix = new Matrix(y);
            Matrix b_matrix = X_matrix.solve(y_matrix);
// string2return += "PLOG / " + P[j] + "\t" + Math.exp(b_matrix.get(0,0)) +
// "\t" + b_matrix.get(1,0) + "\t" + -GasConstant.getCalMolK()*b_matrix.get(2,0) +
// "\t/\n";
            // Computing log fitted rate constants
            Matrix logk_fit = X_matrix.times(b_matrix);
            double[] k_fit = new double[T.length];
            double[] res = new double[T.length];
            double rss = 0.0;
            double rmse = 0.0;
            // Computing the residual sum of squares (rss)
            for (int m = 0; m < T.length; m++) {
                k_fit[m] = logk_fit.get(m, 0);
                // Divide by ln(10), to get back to base 10 for RMSE
                res[m] = (y[m][0] - k_fit[m]) / Math.log(10);
                rss += res[m] * res[m];
            }
            // Calculating the root mean square error
            rmse = Math.sqrt(rss / T.length);
            System.out.println("Root mean square error " + rmse);
            if (option.equals("ARH")) {
                string2return += String.format(
                        "%4.3e\t%5.3f\t%8.2f\t!RMSE of fit %2.2f",
                        Math.exp(b_matrix.get(0, 0)), b_matrix.get(1, 0),
                        -GasConstant.getCalMolK() * b_matrix.get(2, 0), rmse);
            } else if (option.equals("PLOG")) {
                string2return += String
                        .format("PLOG / %6.2f %4.3e\t%5.3f\t%8.2f /\t!RMSE of fit %2.2f",
                                P[0], Math.exp(b_matrix.get(0, 0)),
                                b_matrix.get(1, 0), -GasConstant.getCalMolK()
                                        * b_matrix.get(2, 0), rmse);
            }
        }
        return string2return;
    }

    public static double phi(double argument, int order) {
        if (argument > 1.0)
            argument = 1.0;
        return Math.cos((order) * Math.acos(argument));
    }
}
