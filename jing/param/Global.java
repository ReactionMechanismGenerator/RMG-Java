package jing.param;

public class Global {
	public static StringBuilder diagnosticInfo = new StringBuilder();
	public static StringBuilder enlargerInfo = new StringBuilder("Species \t singleReaction \t  doubleReaction \t longestTime \t longestTemplate \t H_Abstractiontimes \n");
	public static double RT_identifyReactedSites = 0;
	public static double RT_reactChemGraph = 0;
	public static double RT_findRateConstant = 0;
	public static double tAtInitialization;
	public static double makeSpecies = 0;
	public static double checkReactionReverse = 0;
	public static double makeTR = 0;
        //10/25/07 gmagoon: commenting out global temp/pressure parameters, which should not be used if code is applied to systems with multiple temperatures or pressures
	//public static Temperature temperature = new Temperature();
	//public static Pressure pressure = new Pressure();
        //10/29/07 gmagoon: adding new global temperature variables for the high and low temperatures of the input range
        public static Temperature lowTemperature;
        public static Temperature highTemperature;
        
	public static double solvertime= 0; //the time taken by just daspk.
	public static double writeSolverFile = 0;
	public static double readSolverFile = 0;
	public static double solverPrepossesor = 0;
	public static double transferReaction = 0;
	public static double speciesStatusGenerator = 0;
	public static int solverIterations = 0;
	
	public static double moveUnreactedToReacted = 0;
	
	public static double getReacFromStruc = 0;
	public static double generateReverse = 0;
	
	public static double chemkinThermo = 0;
	public static double chemkinReaction = 0;

}
