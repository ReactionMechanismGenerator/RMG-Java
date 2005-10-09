//!********************************************************************************
//!
//!    RMG: Reaction Mechanism Generator
//!
//!    Copyright: Jing Song, MIT, 2002, all rights reserved
//!
//!    Author's Contact: jingsong@mit.edu
//!
//!    Restrictions:
//!    (1) RMG is only for non-commercial distribution; commercial usage
//!        must require other written permission.
//!    (2) Redistributions of RMG must retain the above copyright
//!        notice, this list of conditions and the following disclaimer.
//!    (3) The end-user documentation included with the redistribution,
//!        if any, must include the following acknowledgment:
//!        "This product includes software RMG developed by Jing Song, MIT."
//!        Alternately, this acknowledgment may appear in the software itself,
//!        if and wherever such third-party acknowledgments normally appear.
//!
//!    RMG IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED
//!    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//!    OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//!    DISCLAIMED.  IN NO EVENT SHALL JING SONG BE LIABLE FOR
//!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
//!    OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
//!    OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//!    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
//!    THE USE OF RMG, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//!
//!******************************************************************************


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


