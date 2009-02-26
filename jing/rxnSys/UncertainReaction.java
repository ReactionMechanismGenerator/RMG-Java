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



package jing.rxnSys;

import jing.rxn.*;
import jing.chem.*;
import java.util.*;
import jing.param.*;
import jing.chem.Species;
import jing.rxn.Reaction;
import jing.param.Pressure;
import jing.param.Temperature;
import jing.rxn.ReactionGenerator;

//## package jing::rxnSys

//----------------------------------------------------------------------------
// jing\rxnSys\UncertainReaction.java
//----------------------------------------------------------------------------



//## class UncertainReaction
class UncertainReaction {

  protected Reaction r;

  protected double sens;

  protected double unc;

  protected double product;

  //## operation UncertainReaction(Reaction, double, double)
  public UncertainReaction(Reaction p_reaction, double p_sens, double p_unc){
    //#[ operation UncertainReaction(Reaction, double, double)
    r = p_reaction;
    sens = p_sens;
    unc = p_unc;
    product = Math.abs(p_unc*p_sens);
    //#]
  }

  public Reaction getReaction(){
    return r;
  }

  public double getSensitivity(){
    return sens;
  }

  public double getUncertainty(){
    return unc;
  }

  public double getProduct(){
    return product;
  }
}
/*********************************************************************
        File Path	: RMG\RMG\jing\rxnSys\UncertainReaction.java
*********************************************************************/


