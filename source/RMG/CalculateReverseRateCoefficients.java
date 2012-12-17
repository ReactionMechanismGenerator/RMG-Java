/**
 * CalculateReverseRateCoefficients.class Input: [0] chem.inp file with CHEB p-dep nomenclature [1] Pressure at which
 * you want to calculate the reverse rate coefficient Output: reverseRateCoefficients.txt file with reverse rate
 * coefficients in Arrhenius form Note: This will calculate reverse rate coefficients only for PDEP reactions in PLOG
 * and CHEB format Note: For PLOGS it will reverse all the pressures as given in the forward reaction and not use the
 * user defined pressure Note: The reverse rate constants are valid in range of 300 - 2000 K (beyond this range the
 * k_CHEB have large errors) The k values are calculated at 20 temperature points between 2000 and 300 K
 */
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;
import Jama.Matrix;
import jing.chemParser.ChemParser;
import jing.param.GasConstant;

public class CalculateReverseRateCoefficients {
    // Assumes less than 401 species are present in the mechanism
    static String[] names = new String[400];
    static double[][] coeffs = new double[400][14];
    static int numR = 0;
    static int numP = 0;
    static StringBuilder reverseRateCoefficients = new StringBuilder();

    public static void main(String args[]) throws IOException {
        // Temperature is assumed to have units in Kelvin
        String temperatureString = "";
        // Number of temperature points, hard coded to 20
        double[] T = new double[20];
        for (int i = 0; i < T.length; i++) {
            // As cheb do not work well for T > 2000, we will evaluate in range of 300 to 2000
            double Thigh_plog = 2000;
            double Tlow_plog = 300;
            // Discretizing Temperature points, in inverse space
            T[i] = 1 / ((1 / Thigh_plog) + i * (1 / Tlow_plog - 1 / Thigh_plog)
                    / (T.length - 1));
        }
        // Read in argument 1 as Pressure to evaluate the CHEB polynomial at
        double Pressure = Double.parseDouble(args[1]);
        // Try to Read in the chem.inp file or throw file not found error
        try {
            FileReader fr_chemfile = new FileReader(args[0]);
            BufferedReader br_chemfile = new BufferedReader(fr_chemfile);
            String line = ChemParser.readMeaningfulLine(br_chemfile, true);
            // Continue reading in the file until "THERMO" is read in
            boolean foundThermo = false;
            while (!foundThermo) {
                line = ChemParser.readMeaningfulLine(br_chemfile, true);
                if (line.startsWith("THERMO")) {
                    foundThermo = true;
                    line = ChemParser.readMeaningfulLine(br_chemfile, true); // This line contains the global Tmin,
// Tmax, Tmid
                    line = ChemParser.readMeaningfulLine(br_chemfile, true); // This line should have thermo (or
// comments)
                }
            }
            // Species counter
            int counter = 0;
            // Read in all of thermo species
            while (line != null && !line.equals("END")) {
                if (!line.startsWith("!")) {
                    // Thermo data for each species stored in 4 consecutive lines
                    // The species name are stored in first 16 characters of the line in a chemkin file
                    String speciesName = line.substring(0, 16);
                    StringTokenizer st_chemkinname = new StringTokenizer(
                            speciesName);
                    // Store the name of the current species
                    names[counter] = st_chemkinname.nextToken().trim();
                    // Read in the thermo coefficients for the species
                    // The for loop reads in line 1 to 3 of NASA polynomial
                    for (int numLines = 0; numLines < 2; ++numLines) {
                        line = ChemParser
                                .readMeaningfulLine(br_chemfile, false);
                        for (int numcoeffs = 0; numcoeffs < 5; ++numcoeffs) {
                            /*
                             * Each coefficient in the NASA polynomial is 15 characters long hence substring(0,15) is
                             * 1st coefficient, substring(15,30) is second coefficient and so on
                             */
                            coeffs[counter][5 * numLines + numcoeffs] = Double
                                    .parseDouble(line.substring(15 * numcoeffs,
                                            15 * (numcoeffs + 1)).trim());
                        }
                    }
                    // The for loop reads in line 4 of NASA polynomial which contains only 4 coefficients
                    line = ChemParser.readMeaningfulLine(br_chemfile, false);
                    for (int numcoeffs = 0; numcoeffs < 4; ++numcoeffs) {
                        coeffs[counter][10 + numcoeffs] = Double
                                .parseDouble(line.substring(15 * numcoeffs,
                                        15 * (numcoeffs + 1)).trim());
                    }
                    // Finished reading in thermo for species, increment counter
                    ++counter;
                }
                line = ChemParser.readMeaningfulLine(br_chemfile, true);
            }
            // Continue reading in lines till REACTIONS is found
            while (!line.startsWith("REACTIONS")) {
                // Find "REACTIONS" line
                line = ChemParser.readMeaningfulLine(br_chemfile, true);
            }
            // Determine what units Ea is in
            StringTokenizer st = new StringTokenizer(line);
            // If no units are explicitly mentioned then it is assumed to be in cal/mol
            double R = 1.987;
            while (st.hasMoreTokens()) {
                String nextToken = st.nextToken().toLowerCase();
                if (nextToken.equals("kcal/mol"))
                    R = 1.987e-3;
                else if (nextToken.equals("kj/mol"))
                    R = 8.314e-3;
                else if (nextToken.equals("j/mol"))
                    R = 8.314;
            }
            // next line should have kinetic info
            line = ChemParser.readMeaningfulLine(br_chemfile, true);
            // Read in the lines till END is reached
            while (line != null && !line.equals("END")) {
                // By default we will set boolean reversible as true, this will be only set false only for irreversible
// reactions
                boolean reversible = true;
                /*
                 * If line is a comment or is a line giving additional kinetics skip over that line, we want to search
                 * for lines containing reaction
                 */
                if (!line.startsWith("!")
                        && !line.toLowerCase().contains("low")
                        && !line.toLowerCase().contains("troe")
                        && !line.toLowerCase().contains("dup")
                        && !line.toLowerCase().contains("plog")
                        && !line.contains("CHEB")) {
                    // Print reaction line to console
                    System.out.println(line);
                    // Initializing variables
                    String rxnString = "";
                    String fullRxnString = line;
                    String shortRxnString = "";
                    // Max react index and product index is set to 3
                    int[] reactsIndex = new int[3];
                    int[] prodsIndex = new int[3];
                    double A = 0.0;
                    double n = 0.0;
                    double E = 0.0;
                    double[] logk = new double[T.length];
                    boolean chebyshevRate = false;
                    boolean rmgRate = false;
                    boolean plogRate = false;
                    // Find all Chebyshev/plog rate coefficients
                    // 1.0E0 0.0 0.0 is identifier for cheb/plog polynomials
                    if (line.contains("1.0E0 0.0 0.0")
                            || line.contains("1.0	0.0	0.0")) {
                        st = new StringTokenizer(line);
                        rxnString = st.nextToken();
                        // Gets the reaction string without the kinetics and the comments
                        shortRxnString = rxnString;
                        // Split string to get reaction and products
                        String[] reactsANDprods = rxnString.split("=");
                        // If chebychev coefficients
                        if (shortRxnString.contains("(+m)")) {
                            chebyshevRate = true;
                            // Determine the reactants
                            reactsIndex = determineSpeciesIndex(reactsANDprods[0]
                                    .substring(0,
                                            reactsANDprods[0].length() - 4));
                            numR = determineNumberOfSpecies(reactsANDprods[0]
                                    .substring(0,
                                            reactsANDprods[0].length() - 4));
                            // Determine the products
                            prodsIndex = determineSpeciesIndex(reactsANDprods[1]
                                    .substring(0,
                                            reactsANDprods[1].length() - 4));
                            numP = determineNumberOfSpecies(reactsANDprods[1]
                                    .substring(0,
                                            reactsANDprods[1].length() - 4));
                            line = ChemParser.readUncommentLine(br_chemfile); // TCHEB & PCHEB info
                            // Skip over the comments of PDEP network
                            while (line.startsWith("!")) {
                                line = ChemParser
                                        .readUncommentLine(br_chemfile);
                            }
                            // Get Tmin and Tmax of CHEB poly
                            StringTokenizer st_cheb = new StringTokenizer(line,
                                    "/");
                            String currentToken = st_cheb.nextToken(); // TCHEB
                            StringTokenizer st_values = new StringTokenizer(
                                    st_cheb.nextToken());
                            double Tmin = Double.parseDouble(st_values
                                    .nextToken());
                            double Tmax = Double.parseDouble(st_values
                                    .nextToken());
                            // Get Pmin and Pmax of CHEB poly
                            currentToken = st_cheb.nextToken(); // PCHEB
                            st_values = new StringTokenizer(st_cheb.nextToken());
                            double Pmin = Double.parseDouble(st_values
                                    .nextToken());
                            double Pmax = Double.parseDouble(st_values
                                    .nextToken());
                            // Get the number of coeff in T and P
                            line = ChemParser.readUncommentLine(br_chemfile); // # of basis set info
                            st_cheb = new StringTokenizer(line, "/");
                            currentToken = st_cheb.nextToken();
                            st_values = new StringTokenizer(st_cheb.nextToken());
                            int nT = Integer.parseInt(st_values.nextToken());
                            int nP = Integer.parseInt(st_values.nextToken());
                            // Extract the coefficients
                            double[] coeffs = new double[nT * nP];
                            int coeffCounter = 0;
                            while (coeffCounter < nT * nP) {
                                line = ChemParser
                                        .readUncommentLine(br_chemfile);
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
                            /**
                             * CALCULATE FORWARD RATE COFFICIENT
                             */
                            // NOTE: The values returned are ln(kf)
                            double[] logkf = calculateCHEBrate(R, coeffs, nT,
                                    nP, Tmin, Tmax, Pmin, Pmax, Pressure, T);
                            /**
                             * ARHENIUS PARAMETERS OF REVERSE RATE COEFFICIENT
                             */
                            String ReverseRateCoeff = "";
                            ReverseRateCoeff = calculate_reverserate_coeff(T,
                                    reactsIndex, prodsIndex, reversible, logkf,
                                    fullRxnString, shortRxnString);
                            // Append it to string for writing to file
                            reverseRateCoefficients.append(shortRxnString
                                    + "\t" + "REV/ " + ReverseRateCoeff + "\n");
                        } else {
                            plogRate = true; // Determine the reactants
                            reverseRateCoefficients.append(shortRxnString
                                    + "\n");
                            reactsIndex = determineSpeciesIndex(reactsANDprods[0]);
                            numR = determineNumberOfSpecies(reactsANDprods[0]);
                            // Determine the products
                            prodsIndex = determineSpeciesIndex(reactsANDprods[1]);
                            numP = determineNumberOfSpecies(reactsANDprods[1]);
                            line = ChemParser.readUncommentLine(br_chemfile); // Begin reading PLOG line
                            while (line != null && line.startsWith("PLOG")) {
                                StringTokenizer st_plog = new StringTokenizer(
                                        line, "/");
                                String temporaryString = st_plog.nextToken();
                                StringTokenizer st_plog_info = new StringTokenizer(
                                        st_plog.nextToken());
                                // Get the pressure for PLOG
                                String plog_pressure = st_plog_info.nextToken();
                                // Extract the Arrhenius parameters
                                double plog_A = Double.parseDouble(st_plog_info
                                        .nextToken());
                                double plog_n = Double.parseDouble(st_plog_info
                                        .nextToken());
                                double plog_E = Double.parseDouble(st_plog_info
                                        .nextToken());
                                /**
                                 * CALCULATE FORWARD RATE CONSTANTS *
                                 */
                                // NOTE: The values returned are ln(kf)
                                double[] logkf = calculateArrheniusrate(R,
                                        plog_A, plog_n, plog_E, T);
                                /**
                                 * CALCULATE REVERSE RATE COEFFICIENTS
                                 */
                                String ReverseRateCoeff = "";
                                ReverseRateCoeff = calculate_reverserate_coeff(
                                        T, reactsIndex, prodsIndex, reversible,
                                        logkf, fullRxnString, shortRxnString);
                                reverseRateCoefficients.append("PLOG / "
                                        + plog_pressure + " "
                                        + ReverseRateCoeff + "\n");
                                // Read in Next line
                                line = ChemParser
                                        .readUncommentLine(br_chemfile);
                            }
                        }
                    } else if (line.contains("=") && !line.contains("(+M)")
                            && !line.contains("+m") && !line.contains("+M=")) {
                        rmgRate = true;
                        shortRxnString = line;
                        st = new StringTokenizer(line);
                        shortRxnString = st.nextToken();
                        if (shortRxnString.contains("=>"))
                            reversible = false;
                        A = Double.parseDouble(st.nextToken());
                        n = Double.parseDouble(st.nextToken());
                        E = Double.parseDouble(st.nextToken());
                        /**
                         * CALCULATE FORWARD RATE CONSTANTS *
                         */
                        // NOTE: The values returned are ln(kf)
                        double[] logkf = calculateArrheniusrate(R, A, n, E, T);
                        /**
                         * CALCULATE REVERSE RATE COEFFICIENTS
                         */
                        String[] reactsANDprods = shortRxnString.split("=");
                        // Determine the reactants
                        reactsIndex = determineSpeciesIndex(reactsANDprods[0]);
                        numR = determineNumberOfSpecies(reactsANDprods[0]);
                        // Determine the products
                        prodsIndex = determineSpeciesIndex(reactsANDprods[1]);
                        numP = determineNumberOfSpecies(reactsANDprods[1]);
                        String ReverseRateCoeff = "";
                        ReverseRateCoeff = calculate_reverserate_coeff(T,
                                reactsIndex, prodsIndex, reversible, logkf,
                                fullRxnString, shortRxnString);
                        reverseRateCoefficients.append(shortRxnString + "\t"
                                + "REV/ " + ReverseRateCoeff + "\n");
                    }
                }
                // Read in next line (outermost while loop)
                line = ChemParser.readMeaningfulLine(br_chemfile, true);
            }
        } catch (FileNotFoundException e) {
            System.err.println("File was not found: " + args[0] + "\n");
        }
        // Try writing the output file
        try {
            File reversek = new File("reverseRateCoefficients.txt");
            FileWriter fw_rxns = new FileWriter(reversek);
            fw_rxns.write(reverseRateCoefficients.toString());
            fw_rxns.close();
        } catch (IOException e) {
            System.out
                    .println("Could not write reverseRateCoefficients.txt files");
            System.exit(0);
        }
    }

    private static double[] calculateArrheniusrate(double R, double A,
            double n, double E, double[] T) {
        // To calculate the rate for reaction with Arrhenius parameters
        double logk[] = new double[T.length];
        for (int k = 0; k < T.length; k++) {
            logk[k] = Math.log(A * Math.pow(T[k], n) * Math.exp(-E / R / T[k]));
        }
        // ****** NOTE:The k values returned are ln(k) *******
        return logk;
    }

    private static String calculate_reverserate_coeff(double[] T,
            int[] reactsIndex, int[] prodsIndex, boolean reversible,
            double[] logkf, String fullRxnString, String shortRxnString) {
        // To calculate the reverse rate coefficient, and its arrhenius parameters
        double[] logKeq = new double[T.length];
        double[] logkr = new double[T.length];
        for (int iii = 0; iii < T.length; iii++) {
            double H_RT = 0;
            double S_R = 0;
            int coeffsCounter = 0;
            double Temperature = T[iii];
            // Check temperature to see which set of NASA polynomials to use
            if (Temperature < 1000.0)
                coeffsCounter = 0;
            else
                coeffsCounter = -7;
            // H/RT(species)= a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T
            // S/R(species) = a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5 *T^4/4 + a7
            // H/RT[reaction] = - sum(H/RT(reactants)) + sum(H/RT(products))
            // S/R[reaction] = - sum(S/R(reactants)) + sum(S/R(products))
            for (int numReacts = 0; numReacts < numR; numReacts++) {
                H_RT -= coeffs[reactsIndex[numReacts]][coeffsCounter + 7]
                        + coeffs[reactsIndex[numReacts]][coeffsCounter + 8]
                        * Temperature / 2
                        + coeffs[reactsIndex[numReacts]][coeffsCounter + 9]
                        * Temperature * Temperature / 3
                        + coeffs[reactsIndex[numReacts]][coeffsCounter + 10]
                        * Temperature * Temperature * Temperature / 4
                        + coeffs[reactsIndex[numReacts]][coeffsCounter + 11]
                        * Temperature * Temperature * Temperature * Temperature
                        / 5
                        + coeffs[reactsIndex[numReacts]][coeffsCounter + 12]
                        / Temperature;
                S_R -= coeffs[reactsIndex[numReacts]][coeffsCounter + 7]
                        * Math.log(Temperature)
                        + coeffs[reactsIndex[numReacts]][coeffsCounter + 8]
                        * Temperature
                        + coeffs[reactsIndex[numReacts]][coeffsCounter + 9]
                        * Temperature * Temperature / 2
                        + coeffs[reactsIndex[numReacts]][coeffsCounter + 10]
                        * Temperature * Temperature * Temperature / 3
                        + coeffs[reactsIndex[numReacts]][coeffsCounter + 11]
                        * Temperature * Temperature * Temperature * Temperature
                        / 4
                        + coeffs[reactsIndex[numReacts]][coeffsCounter + 13];
            }
            for (int numProds = 0; numProds < numP; numProds++) {
                H_RT += coeffs[prodsIndex[numProds]][coeffsCounter + 7]
                        + coeffs[prodsIndex[numProds]][coeffsCounter + 8]
                        * Temperature / 2
                        + coeffs[prodsIndex[numProds]][coeffsCounter + 9]
                        * Temperature * Temperature / 3
                        + coeffs[prodsIndex[numProds]][coeffsCounter + 10]
                        * Temperature * Temperature * Temperature / 4
                        + coeffs[prodsIndex[numProds]][coeffsCounter + 11]
                        * Temperature * Temperature * Temperature * Temperature
                        / 5 + coeffs[prodsIndex[numProds]][coeffsCounter + 12]
                        / Temperature;
                S_R += coeffs[prodsIndex[numProds]][coeffsCounter + 7]
                        * Math.log(Temperature)
                        + coeffs[prodsIndex[numProds]][coeffsCounter + 8]
                        * Temperature
                        + coeffs[prodsIndex[numProds]][coeffsCounter + 9]
                        * Temperature * Temperature / 2
                        + coeffs[prodsIndex[numProds]][coeffsCounter + 10]
                        * Temperature * Temperature * Temperature / 3
                        + coeffs[prodsIndex[numProds]][coeffsCounter + 11]
                        * Temperature * Temperature * Temperature * Temperature
                        / 4 + coeffs[prodsIndex[numProds]][coeffsCounter + 13];
            }
            // G/RT[reaction] = H/RT[reaction] - S/R[reaction]
            // exp(-G/RT[reaction] = exp(-H/RT[reaction] + S/R[reaction])
            // Keq = (RT/P)^-(numP-numR)*exp(-G/RT[reaction])
            // ln(Keq) = -H/RT[reaction] + S/R[reaction] + ln((RT/P)^(numR-numP))
            // RT/P = (8.314 J/molK* T K)/(10^5 N/m2)*10^6 (cm3/m3)
            logKeq[iii] = (-H_RT + S_R)
                    + Math.log(Math.pow(
                            ((8.314 * Temperature / Math.pow(10, 5)) * Math
                                    .pow(10, 6)), (numR - numP)));
            // Check if reaction is reversible or not
            if (reversible) {
                logkr[iii] = (logkf[iii] - logKeq[iii]);
            }
        }
        // Construct matrix X and vector y
        double[][] y = new double[T.length][1];
        double[][] X = new double[T.length][3];
        for (int m = 0; m < T.length; m++) {
            y[m][0] = logkr[m];
            X[m][0] = 1;
            X[m][1] = Math.log(T[m]);
            X[m][2] = 1 / T[m];
        }
        // Solve least-squares problem using inv(XT*X)*(XT*y)
        Matrix X_matrix = new Matrix(X);
        Matrix y_matrix = new Matrix(y);
        Matrix b_matrix = X_matrix.solve(y_matrix);
        // Computing log fitted rate constants
        Matrix logk_fit = X_matrix.times(b_matrix);
        // for (int i=0;i<T.length;i++){
        // System.out.println(T[i]);
        // System.out.println(logkr[i]);
        // System.out.println(logk_fit.get(i,0));
        // }
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
        String string2return = "";
        string2return = String.format(
                "%4.3e\t%5.3f\t%8.2f\t!RMSE of fit %2.2e",
                Math.exp(b_matrix.get(0, 0)), b_matrix.get(1, 0),
                -GasConstant.getCalMolK() * b_matrix.get(2, 0), rmse);
        return string2return;
    }

    private static double[] calculateCHEBrate(double R, double[] a, int nT,
            int nP, double Tlow, double Thigh, double Plow, double Phigh,
            double pressure, double[] T) {
        // TO calculate the rate values for CHEB at given Pressure as given by the user
        double[] Tbar = new double[T.length];
        double[] P = new double[1];
        double[] Pbar = new double[P.length];
        // Getting the normalized cheb T points
        for (int i = 0; i < T.length; i++) {
            Tbar[i] = (2 / T[i] - 1 / Tlow - 1 / Thigh)
                    / (1 / Thigh - 1 / Tlow);
        }
        // Getting normalized cheb P points
        P[0] = pressure;
        Pbar[0] = (2 * Math.log10(P[0]) - Math.log10(Plow) - Math.log10(Phigh))
                / (Math.log10(Phigh) - Math.log10(Plow));
        // Evaluating k_cheb(Tcheb,Pcheb) to get k(T,P) values
        double[] logk = new double[T.length];
        double[] logk_rxn = new double[T.length];
        for (int i = 0; i < T.length; i++) {
            for (int k = 0; k < nT; k++) {
                for (int l = 0; l < nP; l++) {
                    // This gets us the log10_k values
                    logk[i] += a[k * nP + l] * phi(Tbar[i], k)
                            * phi(Pbar[0], l);
                }
            }
            // This will give us the ln(k) values
            logk_rxn[i] = Math.log(Math.pow(10, logk[i]));
        }
        // ****** NOTE:The k values returned are ln(k) *******
        return logk_rxn;
    }

    public static double phi(double argument, int order) {
        if (argument > 1.0)
            argument = 1.0;
        return Math.cos((order) * Math.acos(argument));
    }

    public static int determineNumberOfSpecies(String reactsORprods) {
        StringTokenizer st_reacts = new StringTokenizer(reactsORprods, "+");
        int numSpecies = 0;
        while (st_reacts.hasMoreTokens()) {
            ++numSpecies;
            String tempString = st_reacts.nextToken();
            if (tempString.startsWith("2"))
                ++numSpecies;
        }
        return numSpecies;
    }

    public static int[] determineSpeciesIndex(String reactsORprods) {
        if (reactsORprods.startsWith(">"))
            reactsORprods = reactsORprods.substring(1, reactsORprods.length());
        int[] speciesIndex = new int[3];
        int speciesCounter = 0;
        StringTokenizer st_reacts = new StringTokenizer(reactsORprods, "+");
        while (st_reacts.hasMoreTokens()) {
            String reactantString = st_reacts.nextToken();
            boolean groupedSpecies = false;
            if (reactantString.startsWith("2")) {
                reactantString = reactantString.substring(1,
                        reactantString.length());
                groupedSpecies = true;
            }
            boolean found = false;
            for (int numSpecies = 0; numSpecies < names.length; numSpecies++) {
                if (reactantString.equals(names[numSpecies])) {
                    speciesIndex[speciesCounter] = numSpecies;
                    ++speciesCounter;
                    if (groupedSpecies) {
                        speciesIndex[speciesCounter] = numSpecies;
                        ++speciesCounter;
                    }
                    found = true;
                    break;
                }
            }
            if (!found) {
                System.err.println("Could not find thermo for species: "
                        + reactantString);
// System.exit(0);
            }
        }
        return speciesIndex;
    }
}
