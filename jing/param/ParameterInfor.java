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



package jing.param;


import java.util.*;

//## package jing::param

//----------------------------------------------------------------------------
// jing\param\ParameterInfor.java
//----------------------------------------------------------------------------

//## class ParameterInfor
public class ParameterInfor {

    protected int ID = -1;		//## attribute ID

    protected String name;		//## attribute name

    protected double range;		//## attribute range

    protected int type = -1;		//## attribute type


    // Constructors

    //## operation ParameterInfor(String,int,double)
    public  ParameterInfor(String p_name, int p_ID, double p_range) {
        //#[ operation ParameterInfor(String,int,double)
        name = p_name.toUpperCase();

        if (name.equals("C0")) type = 1;
        else if (name.equals("T")) type = 2;
         else if (name.equalsIgnoreCase("k")) type = 3;//svp
        else type = -1;

        ID = p_ID;
        range = p_range;


        //#]
    }
    public  ParameterInfor() {
    }

    public int getID() {
        return ID;
    }

    public String getName() {
        return name;
    }

    public double getRange() {
        return range;
    }

    public int getType() {
        return type;
    }

}
/*********************************************************************
	File Path	: RMG\RMG\jing\param\ParameterInfor.java
*********************************************************************/

