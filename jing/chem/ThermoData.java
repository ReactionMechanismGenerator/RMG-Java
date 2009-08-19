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




package jing.chem;


import java.util.*;
import jing.param.*;
import jing.param.Temperature;

//## package jing::chem

//----------------------------------------------------------------------------
// jing\chem\ThermoData.java
//----------------------------------------------------------------------------

/**
Mutable objects.  Becareful when use this.plus() or this.minus(), this.sum(), since this will be changed after such operations.
*/
//## class ThermoData
public class ThermoData extends ThermoGAValue {


    // Constructors

    //## operation ThermoData()
    public  ThermoData() {
        //#[ operation ThermoData()
        super();
        //#]
    }
    //## operation ThermoData(double,double,double,double,double,double,double,double,double,String)
    public  ThermoData(double p_H298, double p_S298, double p_Cp300, double p_Cp400, double p_Cp500, double p_Cp600, double p_Cp800, double p_Cp1000, double p_Cp1500, double p_dH,double p_dS,double p_dCp,String p_comments) {
        //#[ operation ThermoData(double,double,double,double,double,double,double,double,double,String)
        super(p_H298,p_S298,p_Cp300,p_Cp400,p_Cp500,p_Cp600,p_Cp800,p_Cp1000,p_Cp1500,p_dH,p_dS,p_dCp,p_comments);



        //#]
    }
    //## operation ThermoData(ThermoGAValue)
    public  ThermoData(ThermoGAValue p_ga) {
        //#[ operation ThermoData(ThermoGAValue)
        super(p_ga);
        //#]
    }

    //## operation ThermoData(String, ThermoData, String)
    //svp
    public ThermoData(String p_name, ThermoData p_td, String p_comments){
      //#[ operation ThermoData(String, ThermoData, String)
      super(p_name, p_td, p_comments);
      //#]
    }
    
    public ThermoData(String p_name, ThermoData p_td, String p_comments, String p_source) {
    	super(p_name, p_td, p_comments, p_source);
    }


    //## operation calculateCp(double)
    public double calculateCp(double p_T) {
        //#[ operation calculateCp(double)
        double T = p_T;

        if (300<=T && T<=400) {
                return (Cp300+(Cp400-Cp300)*(T-300.0)/100.0);
        }
        else if (T<=500) {
                return (Cp400+(Cp500-Cp400)*(T-400.0)/100.0);
        }
        else if (T<=600) {
                return (Cp500+(Cp600-Cp500)*(T-500.0)/100.0);
        }
        else if (T<=800) {
                return (Cp600+(Cp800-Cp600)*(T-600.0)/200.0);
        }
        else if (T<=1000) {
                return (Cp800+(Cp1000-Cp800)*(T-800.0)/200.0);
        }
        else if (T<=1500) {
                return (Cp1000+(Cp1500-Cp1000)*(T-1000.0)/500.0);
        }
        else {
                throw new TemperatureOutOfRangeException("Thermo: 300K~1500K;");
        }
        //#]
    }

    //## operation calculateCp(Temperature)
    public double calculateCp(Temperature p_temperature) {
        //#[ operation calculateCp(Temperature)
        double T = p_temperature.getK();

        return calculateCp(T);
        //#]
    }

    //## operation calculateG(Temperature)
    public double calculateG(Temperature p_temperature) {
        //#[ operation calculateG(Temperature)
        double T = p_temperature.getK();

        return (calculateH(p_temperature)*1000 - T*calculateS(p_temperature))/1000;
        //#]
    }

    //## operation calculateH(Temperature)
    public double calculateH(Temperature p_temperature) {
        //#[ operation calculateH(Temperature)
        double T = p_temperature.getK();
        double H = H298*1000.0;

        double a350 = getCpSlope(Cp300,Cp400,300,400);
        double b350 = getCpIntecept(Cp300,Cp400,300,400);
        double a450 = getCpSlope(Cp400,Cp500,400,500);
        double b450 = getCpIntecept(Cp400,Cp500,400,500);
        double a550 = getCpSlope(Cp500,Cp600,500,600);
        double b550 = getCpIntecept(Cp500,Cp600,500,600);
        double a700 = getCpSlope(Cp600,Cp800,600,800);
        double b700 = getCpIntecept(Cp600,Cp800,600,800);
        double a900 = getCpSlope(Cp800,Cp1000,800,1000);
        double b900 = getCpIntecept(Cp800,Cp1000,800,1000);
        double a1250 = getCpSlope(Cp1000,Cp1500,1000,1500);
        double b1250 = getCpIntecept(Cp1000,Cp1500,1000,1500);

        if (T < 298) throw new TemperatureOutOfRangeException();
        if (T > 300) {
                if (T < 400)
                        H += a350*(T*T-300.0*300.0)/2.0+b350*(T-300.0);
                else
                        H += a350*(400*400-300.0*300.0)/2.0+b350*(400.0-300.0);
        }
        if (T > 400) {
                if (T < 500)
                        H += a450*(T*T-400.0*400.0)/2.0+b450*(T-400.0);
                else
                        H += a450*(500.0*500.0-400.0*400.0)/2.0+b450*(500.0-400.0);
        }
        if (T > 500) {
                if (T < 600)
                        H += a550*(T*T-500.0*500.0)/2.0+b550*(T-500.0);
                else
                        H += a550*(600.0*600.0-500.0*500.0)/2.0+b550*(600.0-500.0);
        }
        if (T > 600) {
                if (T < 800)
                        H += a700*(T*T-600.0*600.0)/2.0+b700*(T-600.0);
                else
                        H += a700*(800.0*800.0-600.0*600.0)/2.0+b700*(800.0-600.0);
        }
        if (T > 800) {
                if (T < 1000)
                        H += a900*(T*T-800.0*800.0)/2.0+b900*(T-800.0);
                else
                        H += a900*(1000.0*1000.0-800.0*800.0)/2.0+b900*(1000.0-800.0);
        }
        if (T > 1000) {
                if (T > 1500) throw new TemperatureOutOfRangeException();
                H += a1250*(T*T-1000.0*1000.0)/2.0+b1250*(T-1000.0);
        }

        return H/1000.0;



        //#]
    }

    //## operation calculateGUpperBound(Temperature)
    //svp
    public double calculateGUpperBound(Temperature p_temperature) {
      //#[ operation calculateGUpperBound(Temperature)
      double T = p_temperature.getK();
      return (calculateHUpperBound(p_temperature)*1000 - T*calculateSLowerBound(p_temperature))/1000;
      //#]
    }

  //## operation calculateGLowerBound(Temperature)
  //svp
    public double calculateGLowerBound(Temperature p_temperature) {
      //#[ operation calculateGLowerBound(Temperature)
      double T = p_temperature.getK();
      return (calculateHLowerBound(p_temperature)*1000 - T*calculateSUpperBound(p_temperature))/1000;
      //#]
    }


    //## operation calculateHUpperBound(Temperature)
    //svp
      public double calculateHUpperBound(Temperature p_temperature) {
        //#[ operation calculateHUpperBound(Temperature)
        return calculateH(p_temperature)+Math.pow(dH,0.5);
        //#]
      }

    //## operation calculateHLowerBound(Temperature)
    //svp
      public double calculateHLowerBound(Temperature p_temperature) {
        //#[ operation calculateHLowerBound(Temperature)
        return calculateH(p_temperature)-Math.pow(dH,0.5);
        //#]
      }


    //## operation calculateS(Temperature)
    public double calculateS(Temperature p_temperature) {
        //#[ operation calculateS(Temperature)
        double T = p_temperature.getK();
        double S = S298;

        double a350 = getCpSlope(Cp300,Cp400,300,400);
        double b350 = getCpIntecept(Cp300,Cp400,300,400);
        double a450 = getCpSlope(Cp400,Cp500,400,500);
        double b450 = getCpIntecept(Cp400,Cp500,400,500);
        double a550 = getCpSlope(Cp500,Cp600,500,600);
        double b550 = getCpIntecept(Cp500,Cp600,500,600);
        double a700 = getCpSlope(Cp600,Cp800,600,800);
        double b700 = getCpIntecept(Cp600,Cp800,600,800);
        double a900 = getCpSlope(Cp800,Cp1000,800,1000);
        double b900 = getCpIntecept(Cp800,Cp1000,800,1000);
        double a1250 = getCpSlope(Cp1000,Cp1500,1000,1500);
        double b1250 = getCpIntecept(Cp1000,Cp1500,1000,1500);

        if (T < 298) throw new TemperatureOutOfRangeException();
        if (T > 300) {
                if (T < 400)
                        S += a350*(T-300.0)+b350*(Math.log(T/300.0));
                else
                        S += a350*(400.0-300.0)+b350*(Math.log(400.0/300.0));
        }
        if (T > 400) {
                if (T < 500)
                        S += a450*(T-400.0)+b450*(Math.log(T/400.0));
                else
                        S += a450*(500.0-400.0)+b450*(Math.log(500.0/400.0));
        }
        if (T > 500) {
                if (T < 600)
                        S += a550*(T-500.0)+b550*(Math.log(T/500.0));
                else
                        S += a550*(600.0-500.0)+b550*(Math.log(600.0/500.0));
        }
        if (T > 600) {
                if (T < 800)
                        S += a700*(T-600.0)+b700*(Math.log(T/600.0));
                else
                        S += a700*(800.0-600.0)+b700*(Math.log(800.0/600.0));
        }
        if (T > 800) {
                if (T < 1000)
                        S += a900*(T-800.0)+b900*(Math.log(T/800.0));
                else
                        S += a900*(1000.0-800.0)+b900*(Math.log(1000.0/800.0));
        }
        if (T > 1000) {
                if (T > 1500) throw new TemperatureOutOfRangeException();
                S += a1250*(T-1000.0)+b1250*(Math.log(T/1000.0));
        }

        return S;



        //#]
    }

    //## operation calculateSUpperBound(Temperature)
    //svp
      public double calculateSUpperBound(Temperature p_temperature) {
        //#[ operation calculateSUpperBound(Temperature)
        return calculateS(p_temperature)+Math.pow(dS,0.5);
        //#]
      }

    //## operation calculateSLowerBound(Temperature)
    //svp
      public double calculateSLowerBound(Temperature p_temperature) {
        //#[ operation calculateSLowerBound(Temperature)
        return calculateS(p_temperature)-Math.pow(dS,0.5);
        //#]
      }


    //## operation copy()
    public ThermoData copy() {
        //#[ operation copy()
        if (this == null) return null;

        return new ThermoData(H298,S298,Cp300,Cp400,Cp500,Cp600,Cp800,Cp1000,Cp1500,dH,dS,dCp,comments);



        //#]
    }

    //## operation getCpIntecept(double,double,double,double)
    public static final double getCpIntecept(double p_cp1, double p_cp2, double p_t1, double p_t2) {
        //#[ operation getCpIntecept(double,double,double,double)
        return (p_cp1*p_t2-p_cp2*p_t1)/(p_t2-p_t1);
        //#]
    }

    //## operation getCpSlope(double,double,double,double)
    private static final double getCpSlope(double p_cp1, double p_cp2, double p_t1, double p_t2) {
        //#[ operation getCpSlope(double,double,double,double)
        return (p_cp2-p_cp1)/(p_t2-p_t1);
        //#]
    }

    //## operation minus(ThermoGAValue)
    public void minus(ThermoGAValue p_thermoData) {
        //#[ operation minus(ThermoGAValue)
        H298 -= p_thermoData.H298;
        S298 -= p_thermoData.S298;
        Cp300 -= p_thermoData.Cp300;
        Cp400 -= p_thermoData.Cp400;
        Cp500 -= p_thermoData.Cp500;
        Cp600 -= p_thermoData.Cp600;
        Cp800 -= p_thermoData.Cp800;
        Cp1000 -= p_thermoData.Cp1000;
        Cp1500 -= p_thermoData.Cp1500;
        dH -= Math.pow(p_thermoData.dH,2);//svp
        dS -= Math.pow(p_thermoData.dS,2);//svp
        dCp -= Math.pow(p_thermoData.dCp,2);//svp




        //#]
    }

    //## operation multiply(int)
    public void multiply(int p_multiplier) {
        //#[ operation multiply(int)
        if (p_multiplier == 1) return;

        H298 *= p_multiplier;
        S298 *= p_multiplier;
        Cp300 *= p_multiplier;
        Cp400 *= p_multiplier;
        Cp500 *= p_multiplier;
        Cp600 *= p_multiplier;
        Cp800 *= p_multiplier;
        Cp1000 *= p_multiplier;
        Cp1500 *= p_multiplier;
        dH *= p_multiplier;//svp
        dS *= p_multiplier;//svp
        dCp *= p_multiplier;//svp


        name = name + "redundancy = " + p_multiplier + '\n';
        //#]
    }

    //## operation plus(ThermoGAValue)
    public void plus(ThermoGAValue p_thermoData) {
        //#[ operation plus(ThermoGAValue)
        if (p_thermoData == null) return;

        H298 += p_thermoData.H298;
        S298 += p_thermoData.S298;
        Cp300 += p_thermoData.Cp300;
        Cp400 += p_thermoData.Cp400;
        Cp500 += p_thermoData.Cp500;
        Cp600 += p_thermoData.Cp600;
        Cp800 += p_thermoData.Cp800;
        Cp1000 += p_thermoData.Cp1000;
        Cp1500 += p_thermoData.Cp1500;
        dH += Math.pow(p_thermoData.dH,2);//svp
        dS += Math.pow(p_thermoData.dS,2);//svp
        dCp += Math.pow(p_thermoData.dCp,2);//svp


        if (p_thermoData.getName() != null) {
                if (name == null) name = p_thermoData.getName() + '\n';
                else name = name + p_thermoData.getName() + '\n';
        }




        //#]
    }

    //## operation toString()
    public String toString() {
        //#[ operation toString()
        return super.toString();
        //#]
    }

}
/*********************************************************************
        File Path	: RMG\RMG\jing\chem\ThermoData.java
*********************************************************************/

