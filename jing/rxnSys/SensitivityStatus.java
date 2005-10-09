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


import java.util.*;
import jing.chem.Species;

//## package jing::rxnSys

//----------------------------------------------------------------------------
// jing\rxnSys\SensitivityStatus.java
//----------------------------------------------------------------------------

//## class SensitivityStatus
//svp
public class SensitivityStatus {

  protected double sensitivity;

  protected double sflux;

  protected int species_number;

  protected int reaction_number;

  //## operation SensitivityStatus(double,double,int,int)
  public SensitivityStatus(double p_sensitivity, double p_sflux, int p_snum, int p_rnum) {
    //#[ operation SensitivityStatus(double,double,int,int)
    sensitivity = p_sensitivity;
    sflux = p_sflux;
    species_number = p_snum;
    reaction_number = p_rnum;
    //#]
  }

  public double getSensitivity() {
    return sensitivity;
  }

  public double getSFlux(){
    return sflux;
  }
  public int getSpeciesNumber(){
    return species_number;
  }

  public int getReactionNumber(){
    return reaction_number;
  }

}
/*********************************************************************
        File Path	: RMG\RMG\jing\rxnSys\SensitivityStatus.java
*********************************************************************/

