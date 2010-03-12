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
	public static Pressure lowPressure;
	public static Pressure highPressure;	
        
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
        public static int maxRadNumForQM;
        //5/13/08 gmagoon: added variables temporarily for automatic time-stepping timing; 6/25/08: commented out (along with timing in JDASSL)
     //   public static double JDASSLtime= 0;
     //   public static double fortrantime= 0;
     //   public static double timesteppingtime= 0;
	
	// 18-Jun-2009 MRH
	//	Used these for profiling purposes for the functions
	//		reactTwoReactants and identifyReactiveSites
//	public static int identifyReactiveSitesCount;
//	public static double RT_reactTwoReactants;
}
