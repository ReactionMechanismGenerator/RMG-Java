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



package jing.chemUtil;


import jing.chem.*;
import java.util.*;
import jing.mathTool.*;
import jing.chem.FGAtom;
import jing.chemParser.*;
import jing.chem.FGElement;
import jing.chem.FreeElectron;

//## package jing::chemUtil

//----------------------------------------------------------------------------
//jing\chemUtil\Node.java
//----------------------------------------------------------------------------

/**
Node is the key graph component.  A node stores user-defined content.  Each node can have non-negative number of neighbors of arcs.  Each node has a special ID (integer) indicating the position number of a node in a graph.
*/
//## class Node
public class Node extends GraphComponent {

  /**
  The position identification of this node.
  */
  protected Integer ID;		//## attribute ID

  //protected Integer centralID = new Integer(-1);		//## attribute centralID

  protected FreeElectron feElement = null;		//## attribute feElement

  protected Object fgElement = null;		//## attribute fgElement


  // Constructors

  /**
  Requires:
  Effects: constructor:
  (1) store the p_element in this Node
  (2) set the ID of this node as p_position if p_position >=0. otherwise, throw invalidPositionIDException.
  Modifies:
  */
  //## operation Node(int,Object)
  public  Node(int p_position, Object p_element) throws InvalidNodeIDException {
      //#[ operation Node(int,Object)
      super(p_element);

      if (p_position < 0) throw new InvalidNodeIDException();
      else ID = new Integer(p_position);




      //#]
  }
  //## operation Node(Integer,Object)
  public  Node(Integer p_position, Object p_element) throws InvalidNodeIDException {
      //#[ operation Node(Integer,Object)
      super(p_element);

      if (p_position.intValue() < 0) throw new InvalidNodeIDException();
      else ID = p_position;



      //#]
  }
  public  Node() {
  }

  /**
  Requires:
  Effects: if pass-in GraphComponent is not an instance of Arc, throw InvalidNeighborException; otherwise, super.addNeighbor()
  Modifies: this.neighbor
  */
  //## operation addNeighbor(GraphComponent)
  public void addNeighbor(GraphComponent p_graphComponent) throws InvalidNeighborException {
      //#[ operation addNeighbor(GraphComponent)
      if (!(p_graphComponent instanceof Arc)) {
        throw new InvalidNeighborException("node");
      }
      else {
        super.addNeighbor(p_graphComponent);
        return;
      }
      //#]
  }

  /**
  Requires:
  Effects: super.addNeighbor()
  Modifies: this.neighbor
  */
  //## operation addNeighbor(Arc)
  public void addNeighbor(Arc p_arc) throws InvalidNeighborException {
      //#[ operation addNeighbor(Arc)
      super.addNeighbor(p_arc);



      //#]
  }

  //## operation changeChemNodeElement(int,Node)
  public Object changeChemNodeElement(int p_changedOrder, Node p_otherNode) {
      //#[ operation changeChemNodeElement(int,Node)
      Object element = getElement();
      if (element instanceof Collection) {
      	HashSet result = new HashSet();
      	Iterator iter = ((Collection)element).iterator();
      	while (iter.hasNext()) {
      		Object nodeElement = iter.next();
      		setElement(nodeElement);
      		Object newElement = changeChemNodeElement(p_changedOrder,p_otherNode);
      		if (newElement != null) {
      			if (newElement instanceof Collection) result.addAll((Collection)newElement);
      			else result.add(newElement);
      		}
      	}
      	setElement(result);
      	return result;
      }
      else if (element instanceof Atom) {
      	return element;
      }
      else if (element instanceof FGAtom) {
      	FGAtom fgAtom = (FGAtom)element;
         	Object newFgAtom = changeFgAtom(fgAtom,p_changedOrder,p_otherNode);
         	setElement(newFgAtom);
         	return newFgAtom;
      }
      else {
      	throw new InvalidChemNodeElementException();
      }
      //#]
  }

  //## operation changeFgAtom(FGAtom,int,Node)
  public Object changeFgAtom(FGAtom p_fgAtom, int p_changedOrder, Node p_otherNode) throws InvalidChangedOrderException {
      //#[ operation changeFgAtom(FGAtom,int,Node)
      FGElement Cs = FGElement.make("Cs");
      FGElement Cd = FGElement.make("Cd");
      FGElement Cdd = FGElement.make("Cdd");
      FGElement Ct = FGElement.make("Ct");
      FGElement CO = FGElement.make("CO");
      FGElement Cb = FGElement.make("Cb");
      FGElement Cbf = FGElement.make("Cbf");
      FGElement Os = FGElement.make("Os");
      FGElement Od = FGElement.make("Od");
      FGElement Oa = FGElement.make("Oa");
      FGElement H = FGElement.make("H");
      FGElement R1H = FGElement.make("R!H"); 
      FGElement R = FGElement.make("R");
      FGElement Sis = FGElement.make("Sis");
      FGElement Sid = FGElement.make("Sid");
      FGElement Sidd = FGElement.make("Sidd");
      FGElement Sit = FGElement.make("Sit");
      FGElement Ss = FGElement.make("Ss");
      FGElement Sd = FGElement.make("Sd");
      FGElement Sa = FGElement.make("Sa");

      FGElement fgElement = p_fgAtom.getFgElement();
      FreeElectron feElement = p_fgAtom.getFreeElectron();
      if (p_changedOrder == 1) {
      	// Cs and CO change into Cd
      	if (fgElement == Cs || fgElement == CO) {
      		return FGAtom.make(Cd,feElement);
      	}
      	// Os changes into Od
      	else if (fgElement == Os) {
      		return FGAtom.make(Od,feElement);
      	}
      	// no elemination can happen here
      	else if (fgElement == H || fgElement == Cdd || fgElement == Ct || fgElement == Cb ||fgElement == Cbf || fgElement == Oa || fgElement == Od ) {
      		return null;
      	}
      	// Cd, R, and R!H, return themselves
      	// addition to Cd, in our all reaction families, end up with Cd only
      	else if (fgElement == R || fgElement == R1H){
      		return p_fgAtom;
      	}
      	// Sis change into Sid
      	else if (fgElement == Sis) {
      		return FGAtom.make(Sid,feElement);
      	}
      	// Ss change into Sd
      	else if (fgElement == Ss) {
      		return FGAtom.make(Sd,feElement);
      	}
      	// Apparently Sidd, Sit, Sa, and Sd cannot do anything
      	else if (fgElement == Sidd || fgElement == Sit || fgElement == Sa || fgElement == Sd) {
      		return null;
      	}
      	// Nothing for Sid ... according to the above comment, all reaction familes
      	//	result in Cd staying as Cd, so I'll assume (for now) that Sid stays
      	//	as Sid
      	else{
      		// Why are we returning a Cdd?  If we encounter Cdd, we return null
      		return FGAtom.make(Cdd,feElement);
      	}
      }
      else if (p_changedOrder == -1) {
      	// Cd and CO change to Cs
      	if (fgElement == Cd || fgElement == CO) {
      		return FGAtom.make(Cs,feElement);
      	}
      	// Cdd and Ct change to Cd
       	else if (fgElement == Ct) {
      		return FGAtom.make(Cd,feElement);
      	}
      	// Cdd change to Cd or CO
      	else if (fgElement == Cdd) {
      		if (p_otherNode != null) {
      			Object otherElement = p_otherNode.getElement();
      			if (otherElement instanceof Collection) {
      				boolean carbon = false;
      				boolean oxygen = false;
      				Iterator iter = ((Collection)otherElement).iterator();
      				while (iter.hasNext()) {
      					ChemNodeElement cne = (ChemNodeElement)iter.next();
      					if (cne.isCarbon() || cne.isSilicon() || cne.isSulfur()) {
      						carbon = true;
      					}
      					else if (cne.isOxygen()) {
      						oxygen = true;
      					}
      					else if (cne.isAny() || cne.isNonH()) {
      						carbon = true;
      						oxygen = true;
      					}
      	     		}
      	     		if (carbon && oxygen) {
      					HashSet result = new HashSet();
      					result.add(FGAtom.make(CO,feElement));
      					result.add(FGAtom.make(Cd,feElement));
      					return result;
      				}
      				else if (carbon) {
      					return FGAtom.make(Cd,feElement);
      				}
      				else if (oxygen) {
      					return FGAtom.make(CO,feElement);
      				}
      				else {
      					throw new InvalidChangedOrderException();
      				}
      			}
      			else {
      				ChemNodeElement cne = (ChemNodeElement)otherElement;
      				if (cne.isCarbon() || cne.isSilicon() || cne.isSulfur()) {
      					return FGAtom.make(Cd,feElement);
      				}
      				else if (cne.isOxygen()) {
      					return FGAtom.make(CO,feElement);
      				}
      				else if (cne.isAny() || cne.isNonH()) {
      					HashSet result = new HashSet();
      					result.add(FGAtom.make(CO,feElement));
      					result.add(FGAtom.make(Cd,feElement));
      					return result;
      				}
      				else {
      					throw new InvalidChangedOrderException();
      				}
      			}
      		}
      		// Cdd + null other-end node, change to Cd and CO
      		else {
      			HashSet result = new HashSet();
      			result.add(FGAtom.make(CO,feElement));
      			result.add(FGAtom.make(Cd,feElement));
      			return result;
      		}
      	}
      	// Od change to Os
      	else if (fgElement == Od) {
      		return FGAtom.make(Os,feElement);
      	}
      	// no addition happens here
      	// now, i am not including Cb-Cb and Cbf-cb bond addition, should add later
      	else if (fgElement == H || fgElement == Cs || fgElement == Cb ||fgElement == Cbf || fgElement == Oa || fgElement == Os ) {
      		return null;
      	}
      	// Sid change to Sis
      	else if (fgElement == Sid) {
      		return FGAtom.make(Sis,feElement);
      	}
      	// Sit change to Sid
       	else if (fgElement == Sit) {
      		return FGAtom.make(Sid,feElement);
      	}
      	// Sidd change to Sid
      	else if (fgElement == Sidd) {
      		if (p_otherNode != null) {
      			Object otherElement = p_otherNode.getElement();
      			if (otherElement instanceof Collection) {
      				boolean carbon = false;
      				Iterator iter = ((Collection)otherElement).iterator();
      				while (iter.hasNext()) {
      					ChemNodeElement cne = (ChemNodeElement)iter.next();
      					if (cne.isCarbon() || cne.isSilicon() || cne.isSulfur() || cne.isOxygen()) {
      						carbon = true;
      					}
      					else if (cne.isAny() || cne.isNonH()) {
      						carbon = true;
      					}
      	     		}
      				if (carbon) {
      					return FGAtom.make(Sid,feElement);
      				}
      				else {
      					throw new InvalidChangedOrderException();
      				}
      			}
      			else {
      				ChemNodeElement cne = (ChemNodeElement)otherElement;
      				if (cne.isCarbon() || cne.isSilicon() || cne.isSulfur() || cne.isOxygen()) {
      					return FGAtom.make(Sid,feElement);
      				}
      				else if (cne.isAny() || cne.isNonH()) {
      					return FGAtom.make(Sid,feElement);
      				}
      				else {
      					throw new InvalidChangedOrderException();
      				}
      			}
      		}
      		// Sidd + null other-end node, change to Sid
      		else {
      			return FGAtom.make(Sid,feElement);
      		}
      	}
      	// Sd change to Ss
      	else if (fgElement == Sd) {
      		return FGAtom.make(Ss,feElement);
      	}
      	// no addition happens here
      	else if (fgElement == Sis || fgElement == Sa || fgElement == Ss ) {
      		return null;
      	}
      	// R and R!H, return themselves
      	else return p_fgAtom;
      }
      else {
      	throw new InvalidChangedOrderException();
      }
      //#]
  }

  //## operation checkFgElement(FGAtom)
  public boolean checkFgElement(FGAtom p_fgAtom) {
      //#[ operation checkFgElement(FGAtom)
      // add fgAtom checking here!
      return true;


      //#]
  }

  //## operation contentSub(GraphComponent)
  public boolean contentSub(GraphComponent p_graphComponent) throws InvalidChemNodeElementException {
      //#[ operation contentSub(GraphComponent)
      if (this == p_graphComponent) return false;
      if (!(p_graphComponent instanceof Node)) return false;

      Node node = (Node)p_graphComponent;

      // compare radical number
      FreeElectron fe1 = getFeElement();
      FreeElectron fe2 = node.getFeElement();
      if (fe1.getOrder() != fe2.getOrder()) return false;

      // compare spin of biradical
      if (fe1.getOrder() == 2) {
      	String spin1 = fe1.getSpin();
      	String spin2 = fe2.getSpin();
      	if (spin1 == null) {
      		if (spin2 != null) return false;
      	}
      	else {
      		if (spin2 != null) {
      			if (!spin1.equals(spin2)) return false;
      		}
      	}
      }

      // compare R, R!H
      Object fge1 = getFgElement();
      Object fge2 = node.getFgElement();

      FGElement any = FGElement.make("R");
      FGElement nonH = FGElement.make("R!H");
      FGElement H = FGElement.make("H");

      if (MathTool.isSub(any,fge2)) return true;
      if (MathTool.isSub(nonH,fge2) && !MathTool.isSub(H,fge1)) return true;

      return MathTool.isSub(fge1,fge2);



      //#]
  }

  //## operation generateFGElement(String,int,int,int,String,int,int)
  public HashSet generateFGElement(String p_atomType, int p_freeValency, int p_single, int p_double, String p_doubleEnd, int p_triple, int p_benzene) {
      //#[ operation generateFGElement(String,int,int,int,String,int,int)
      HashSet result = new HashSet();
      if (p_atomType.compareToIgnoreCase("H")==0) {
      	FGElement H = FGElement.make("H");
      	result.add(H);
      }
      else if (p_atomType.compareToIgnoreCase("C")==0) {
      	FGElement Cs = FGElement.make("Cs");
      	FGElement Cd = FGElement.make("Cd");
      	FGElement Cdd = FGElement.make("Cdd");
      	FGElement Ct = FGElement.make("Ct");
      	FGElement CO = FGElement.make("CO");
      	FGElement Cb = FGElement.make("Cb");
      	FGElement Cbf = FGElement.make("Cbf");

      	int benzene = 0;
      	if (p_benzene == 0) benzene = 0;
      	else if (p_benzene == 1 || p_benzene == 2) benzene = 3;
      	else if (p_benzene==3) benzene = 4;
      	else throw new InvalidConnectivityException();

      	int unocuppied = p_freeValency - p_single - 2*p_double - benzene - 3*p_triple;
      	if (unocuppied < -0.01) throw new InvalidConnectivityException();

      	if (p_triple != 0) {
      		if (p_triple == 1) result.add(Ct);
      		else throw new InvalidConnectivityException();
      	}
      	else if (p_double != 0) {
      		if (p_double == 1) {
      			if (unocuppied == 2) result.add(Cdd);
      			if (p_doubleEnd.compareToIgnoreCase("C")==0) {
      				result.add(Cd);
      			}
      			else if (p_doubleEnd.compareToIgnoreCase("O")==0) {
      				result.add(CO);
      			}
      			else if (p_doubleEnd.compareToIgnoreCase("R")==0 || p_doubleEnd.compareToIgnoreCase("R!H")==0) {
      					result.add(Cd);
      					result.add(CO);
      			}
      			else if (p_doubleEnd.compareToIgnoreCase("S")==0 || p_doubleEnd.compareToIgnoreCase("Si")==0) {
      				result.add(Cd);
      			}
      			else {
      				throw new InvalidConnectivityException();
      			}
      		}
      		else if (p_double == 2 && unocuppied == 0) {
      			result.add(Cdd);
      		}
      		else {
      			throw new InvalidConnectivityException();
      		}
      	}
      	else if (p_benzene != 0) {
      		if (p_benzene == 3) {
      			result.add(Cbf);
      		}
      		else if (p_benzene == 2) {
      			result.add(Cb);
      			if (unocuppied > 0) {
      				result.add(Cbf);
      			}
      		}
      		else if (p_benzene == 1) {
      			result.add(Cb);
      			if (unocuppied > 0) {
      				result.add(Cbf);
      			}
      		}
      		else {
      			throw new InvalidConnectivityException();
      		}
      	}
      	else {
      		result.add(Cs);
      		if (unocuppied==2) {
      			result.add(Cd);
      			result.add(CO);
      		}
      		else if (unocuppied==3) {
      			result.add(Cd);
      			result.add(CO);
      			result.add(Ct);
      			result.add(Cb);
      		}
      		else if (unocuppied==4) {
      			result.add(Cd);
      			result.add(CO);
      			result.add(Cdd);
      			result.add(Ct);
      			result.add(Cb);
      			result.add(Cbf);
      		}
      		else if (unocuppied > 4 || unocuppied < 0) {
      			throw new InvalidConnectivityException();
      		}
      	}
      }
      else if (p_atomType.compareToIgnoreCase("O")==0) {
      	FGElement Os = FGElement.make("Os");
      	FGElement Od = FGElement.make("Od");
      	FGElement Oa = FGElement.make("Oa");

      	int unocuppied = p_freeValency - p_single - 2*p_double;
      	if (unocuppied == 2) {
      		result.add(Os);
      		result.add(Od);
      	}
      	else if (unocuppied == 1) {
      		result.add(Os);
      	}
      	else if (unocuppied == 0) {
      		if (p_double == 1) {
      			result.add(Od);
      		}
      		else if (p_double == 0) {
      			if (p_single == 1 || p_single == 2) {
      				result.add(Os);
      			}
      			else if (p_single == 0) {
      				result.add(Oa);
      			}
      			else {
      				throw new InvalidConnectivityException();
      			}
      		}
      		else {
      			throw new InvalidConnectivityException();
      		}
      	}
      	else {
      		throw new InvalidConnectivityException();
      	}
      }
      else if (p_atomType.compareToIgnoreCase("Cl")==0) {//svp
              FGElement Cl = FGElement.make("Cl");
              result.add(Cl);
      }
      else if (p_atomType.compareToIgnoreCase("Si")==0) {
        	FGElement Sis = FGElement.make("Sis");
        	FGElement Sid = FGElement.make("Sid");
        	FGElement Sidd = FGElement.make("Sidd");
        	FGElement Sit = FGElement.make("Sit");

        	int unocuppied = p_freeValency - p_single - 2*p_double - 3*p_triple;
        	if (unocuppied < -0.01) throw new InvalidConnectivityException();

        	if (p_triple != 0) {
        		if (p_triple == 1) result.add(Sit);
        		else throw new InvalidConnectivityException();
        	}
        	else if (p_double != 0) {
        		if (p_double == 1) {
        			if (unocuppied == 2) result.add(Sidd);
        			if (p_doubleEnd.compareToIgnoreCase("C")==0 || p_doubleEnd.compareToIgnoreCase("O")==0 ||
        					p_doubleEnd.compareToIgnoreCase("S")==0 || p_doubleEnd.compareToIgnoreCase("Si")==0) {
        				result.add(Sid);
        			}
        			else if (p_doubleEnd.compareToIgnoreCase("R")==0 || p_doubleEnd.compareToIgnoreCase("R!H")==0) {
        					result.add(Sid);
        			}
        			else {
        				throw new InvalidConnectivityException();
        			}
        		}
        		else if (p_double == 2 && unocuppied == 0) {
        			result.add(Sidd);
        		}
        		else {
        			throw new InvalidConnectivityException();
        		}
        	}
        	else {
        		result.add(Sis);
        		if (unocuppied==2) {
        			result.add(Sid);
        		}
        		else if (unocuppied==3) {
        			result.add(Sid);
        			result.add(Sit);
        		}
        		else if (unocuppied==4) {
        			result.add(Sid);
        			result.add(Sidd);
        			result.add(Sit);
        		}
        		else if (unocuppied > 4 || unocuppied < 0) {
        			throw new InvalidConnectivityException();
        		}
        	}
      }
      else if (p_atomType.compareToIgnoreCase("S")==0) {
        	FGElement Ss = FGElement.make("Ss");
        	FGElement Sd = FGElement.make("Sd");
        	FGElement Sa = FGElement.make("Sa");

        	int unocuppied = p_freeValency - p_single - 2*p_double;
        	if (unocuppied == 2) {
        		result.add(Ss);
        		result.add(Sd);
        	}
        	else if (unocuppied == 1) {
        		result.add(Ss);
        	}
        	else if (unocuppied == 0) {
        		if (p_double == 1) {
        			result.add(Sd);
        		}
        		else if (p_double == 0) {
        			if (p_single == 1 || p_single == 2) {
        				result.add(Ss);
        			}
        			else if (p_single == 0) {
        				result.add(Sa);
        			}
        			else {
        				throw new InvalidConnectivityException();
        			}
        		}
        		else {
        			throw new InvalidConnectivityException();
        		}
        	}
        	else {
        		throw new InvalidConnectivityException();
        	}
        }
      


      return result;
      //#]
  }

  /**
  Requires:
  Effects: return a set of corresponding FGElement(s) of the atoms' element stored in this node.
  Return value: a collection of FGElement
  Modifies:
  */
  //## operation generateFgElement()
  public Collection generateFgElement() {
      //#[ operation generateFgElement()
      Object o = getElement();
      if (o == null) throw new NullGraphComponentException("node");
      HashSet result = new HashSet();
      if (o instanceof Collection) {
      	Iterator iter = ((Collection)o).iterator();
      	while (iter.hasNext()) {
      		Object thisObject = iter.next();
      		setElement(thisObject);
      		result.addAll(generateFgElement());
      	}
      	setElement(o);
      }
      // if o is an atom, identify its type as a FGAtom
      else if (o instanceof FGAtom) {
      	if (checkFgElement((FGAtom)o)) result.add(((FGAtom)o).getFgElement());
      }
      else if (o instanceof Atom) {
      	Atom a = (Atom)o;
      	int valency = (int)a.getValency();
      	int singleNum = 0;
      	int doubleNum = 0;
      	int tripleNum = 0;
      	int benzeneNum = 0;
      	Object [] otherEnd = new Object[2];
      	// if atom type is H
      	if (a.isHydrogen()) {
      		FGElement H = FGElement.make("H");
      		result.add(H);
      	}
              else if (a.isChlorine()){//svp
                FGElement Cl = FGElement.make("Cl");
                result.add(Cl);
              }
      	else if (a.isCarbon()) {
      		FGElement Cs = FGElement.make("Cs");
      		FGElement Cd = FGElement.make("Cd");
      		FGElement Cdd = FGElement.make("Cdd");
      		FGElement Ct = FGElement.make("Ct");
      		FGElement CO = FGElement.make("CO");
      		FGElement Cb = FGElement.make("Cb");
      		FGElement Cbf = FGElement.make("Cbf");

      		Iterator iter = getNeighbor();
      		while (iter.hasNext()) {
      			Arc arc = (Arc)iter.next();
      			Object b = arc.getElement();
      			if (!(b instanceof Bond)) {
      				System.out.println("encounter bond list during fgelement generation!!");
      				System.exit(0);
      			}
      			Bond bond = (Bond)b;
      			if (bond.isSingle()) singleNum++;
      			else if (bond.isDouble()) {
					try { otherEnd[doubleNum] = getOtherNode(arc).getElement(); }
					catch (NullPointerException e) {throw new InvalidBondException();}
      				doubleNum++;
      			}
      			else if (bond.isBenzene()) benzeneNum++;
      			else if (bond.isTriple()) tripleNum++;
      			else {
      				throw new InvalidBondException();
      			}
      		}

      		int benzene = 0;
      		if (benzeneNum == 0) benzene = 0;
      		else if (benzeneNum == 1 || benzeneNum == 2) benzene = 3;
      		else if (benzeneNum ==3) benzene = 4;
      		else throw new InvalidConnectivityException();

      		int unocuppied = valency - singleNum - 2*doubleNum - benzene - 3*tripleNum;
      		if (unocuppied < -0.01) throw new InvalidConnectivityException();

      		if (tripleNum != 0) {
      			result.add(Ct);
      		}
      		else if (doubleNum != 0) {
      			if (doubleNum == 1) {
      				Object cne = otherEnd[0];
      				if (unocuppied == 2) result.add(Cdd);
      				if (cne instanceof ChemNodeElement) {
      					ChemNodeElement thiscne = (ChemNodeElement)cne;
      					if (thiscne.isCarbon()) {
      						result.add(Cd);
      					}
      					else if (thiscne.isOxygen()) {
      						result.add(CO);
      					}
      					else if (thiscne.isAny() || thiscne.isNonH()) {
      						result.add(Cd);
      						result.add(CO);
      					}
      					else if (thiscne.isSilicon() || thiscne.isSulfur()) {
      						result.add(Cd);
      					}
      					else {
      						throw new InvalidConnectivityException();
      					}
      				}
      				else if (cne instanceof Collection) {
      					Iterator cne_iter = ((Collection)cne).iterator();
      					while (cne_iter.hasNext()) {
      						ChemNodeElement thiscne = (ChemNodeElement)cne_iter.next();
      						if (thiscne.isCarbon()) {
      							result.add(Cd);
      						}
      						else if (thiscne.isOxygen()) {
      							result.add(CO);
      						}
      						else if (thiscne.isAny() || thiscne.isNonH()) {
      							result.add(Cd);
      							result.add(CO);
      						}
      						else if (thiscne.isSilicon() || thiscne.isSulfur()) {
      							result.add(Cd);
      						}
      						else {
      							throw new InvalidConnectivityException();
      						}
      					}
      				}
      				else {
      					throw new InvalidConnectivityException();
      				}

      			}
      			else if (doubleNum == 2) {
      				result.add(Cdd);
      			}
      			else {
      				throw new InvalidConnectivityException();
      			}
      		}
      		else if (benzeneNum != 0) {
      			if (benzeneNum == 3) {
      				result.add(Cbf);
      			}
      			else if (benzeneNum == 2) {
      				result.add(Cb);
      				if (unocuppied > 0) {
      					result.add(Cbf);
      				}
      			}
      			else if (benzeneNum == 1) {
      				result.add(Cb);
      				if (unocuppied > 0) {
      					result.add(Cbf);
      				}
      			}
      			else {
      				throw new InvalidConnectivityException();
      			}
      		}
      		else {
      			result.add(Cs);
      			int freeValency = valency-singleNum;
      			if (freeValency==2) {
      				result.add(Cd);
      				result.add(CO);
      			}
      			else if (freeValency==3) {
      				result.add(Cd);
      				result.add(CO);
      				result.add(Ct);
      				result.add(Cb);
      			}
      			else if (freeValency==4) {
      				result.add(Cd);
      				result.add(CO);
      				result.add(Cdd);
      				result.add(Ct);
      				result.add(Cb);
      				result.add(Cbf);
      			}
      			else if (freeValency > 4 || freeValency < 0) {
      				throw new InvalidConnectivityException();
      			}
      		}
      	}
      	// if atom is O, identify if it is Os, or Od
      	else if (a.isOxygen()) {
      		FGElement Os = FGElement.make("Os");
      		FGElement Od = FGElement.make("Od");
      		FGElement Oa = FGElement.make("Oa");

      		Iterator iter = getNeighbor();
      		while (iter.hasNext()) {
      			Arc arc = (Arc)iter.next();
      			Bond bond = (Bond)arc.getElement();
      			if (bond.isSingle()) singleNum++;
      			else if (bond.isDouble()) doubleNum++;
      			else throw new InvalidBondException();
      		}
      		int freeValency = valency - singleNum - 2*doubleNum;
      		if (freeValency == 2) {
      			result.add(Os);
      			result.add(Od);
      		}
      		else if (freeValency == 1) {
      			result.add(Os);
      		}
      		else if (freeValency == 0) {
      			if (doubleNum == 1) {
      				result.add(Od);
      			}
      			else if (doubleNum == 0) {
      				if (singleNum == 1 || singleNum == 2) {
      					result.add(Os);
      				}
      				else if (singleNum == 0) {
      					result.add(Oa);
      				}
      				else {
      					throw new InvalidConnectivityException();
      				}
      			}
      			else {
      				throw new InvalidConnectivityException();
      			}
      		}
      		else {
      			throw new InvalidConnectivityException();
      		}
      	}
      	else if (a.isSilicon()) {
      		FGElement Sis = FGElement.make("Sis");
      		FGElement Sid = FGElement.make("Sid");
      		FGElement Sidd = FGElement.make("Sidd");
      		FGElement Sit = FGElement.make("Sit");

      		Iterator iter = getNeighbor();
      		while (iter.hasNext()) {
      			Arc arc = (Arc)iter.next();
      			Object b = arc.getElement();
      			if (!(b instanceof Bond)) {
      				System.out.println("encounter bond list during fgelement generation!!");
      				System.exit(0);
      			}
      			Bond bond = (Bond)b;
      			if (bond.isSingle()) singleNum++;
      			else if (bond.isDouble()) {
      				otherEnd[doubleNum] = getOtherNode(arc).getElement();
      				doubleNum++;
      			}
      			else if (bond.isTriple()) tripleNum++;
      			else {
      				throw new InvalidBondException();
      			}
      		}

      		int unocuppied = valency - singleNum - 2*doubleNum - 3*tripleNum;
      		if (unocuppied < -0.01) throw new InvalidConnectivityException();

      		if (tripleNum != 0) {
      			result.add(Sit);
      		}
      		else if (doubleNum != 0) {
      			if (doubleNum == 1) {
      				Object cne = otherEnd[0];
      				if (unocuppied == 2) result.add(Sidd);
      				if (cne instanceof ChemNodeElement) {
      					ChemNodeElement thiscne = (ChemNodeElement)cne;
      					if (thiscne.isCarbon() || thiscne.isOxygen() ||
      							thiscne.isSilicon() || thiscne.isSulfur()) {
      						result.add(Sid);
      					}
      					else if (thiscne.isAny() || thiscne.isNonH()) {
      						result.add(Sid);
      					}
      					else {
      						throw new InvalidConnectivityException();
      					}
      				}
      				else if (cne instanceof Collection) {
      					Iterator cne_iter = ((Collection)cne).iterator();
      					while (cne_iter.hasNext()) {
      						ChemNodeElement thiscne = (ChemNodeElement)cne_iter.next();
      						if (thiscne.isCarbon() || thiscne.isOxygen() ||
      								thiscne.isSilicon() || thiscne.isSulfur()) {
      							result.add(Sid);
      						}
      						else if (thiscne.isAny() || thiscne.isNonH()) {
      							result.add(Sid);
      						}
      						else {
      							throw new InvalidConnectivityException();
      						}
      					}
      				}
      				else {
      					throw new InvalidConnectivityException();
      				}

      			}
      			else if (doubleNum == 2) {
      				result.add(Sidd);
      			}
      			else {
      				throw new InvalidConnectivityException();
      			}
      		}
      		else {
      			result.add(Sis);
      			int freeValency = valency-singleNum;
      			if (freeValency==2) {
      				result.add(Sid);
      			}
      			else if (freeValency==3) {
      				result.add(Sid);
      				result.add(Sit);
      			}
      			else if (freeValency==4) {
      				result.add(Sid);
      				result.add(Sidd);
      				result.add(Sit);
      			}
      			else if (freeValency > 4 || freeValency < 0) {
      				throw new InvalidConnectivityException();
      			}
      		}
      	}
      	else if (a.isSulfur()) {
      		FGElement Ss = FGElement.make("Ss");
      		FGElement Sd = FGElement.make("Sd");
      		FGElement Sa = FGElement.make("Sa");

      		Iterator iter = getNeighbor();
      		while (iter.hasNext()) {
      			Arc arc = (Arc)iter.next();
      			Bond bond = (Bond)arc.getElement();
      			if (bond.isSingle()) singleNum++;
      			else if (bond.isDouble()) doubleNum++;
      			else throw new InvalidBondException();
      		}
      		int freeValency = valency - singleNum - 2*doubleNum;
      		if (freeValency == 2) {
      			result.add(Ss);
      			result.add(Sd);
      		}
      		else if (freeValency == 1) {
      			result.add(Ss);
      		}
      		else if (freeValency == 0) {
      			if (doubleNum == 1) {
      				result.add(Sd);
      			}
      			else if (doubleNum == 0) {
      				if (singleNum == 1 || singleNum == 2) {
      					result.add(Ss);
      				}
      				else if (singleNum == 0) {
      					result.add(Sa);
      				}
      				else {
      					throw new InvalidConnectivityException();
      				}
      			}
      			else {
      				throw new InvalidConnectivityException();
      			}
      		}
      		else {
      			throw new InvalidConnectivityException();
      		}
      	}
      	// if atom is not H, O, or C, throw InvalidAtomException
      	else {
      		throw new InvalidChemNodeElementException();
      	}
      }
      // if o is some other type, throw InvalidChemNodeElementException();
      else {
      	throw new InvalidChemNodeElementException();
      }

      return result;
      //#]
  }

  //## operation getFeElement()
  public FreeElectron getFeElement() {
      //#[ operation getFeElement()
      identifyFeElement();
      return feElement;
      //#]
  }

  //## operation getFgElement()
  public Object getFgElement() {
      //#[ operation getFgElement()
      identifyFgElement();
      return fgElement;
      //#]
  }

  //## operation getID()
  public Integer getID() {
      //#[ operation getID()
      return ID;
      //#]
  }

  //## operation getOtherArcs(Arc)
  public HashSet getOtherArcs(Arc p_arc) {
      //#[ operation getOtherArcs(Arc)
      HashSet otherArcs = new HashSet();

      Iterator iter = getNeighbor();
      while (iter.hasNext()) {
      	Arc arc = (Arc)iter.next();
      	if (!arc.equals(p_arc)) otherArcs.add(arc);
      }

      return otherArcs;
      //#]
  }
  
  public HashSet getNeighboringNodes(){
	  HashSet neighboringNodes = new HashSet();
	  Iterator iter = getNeighbor();
	  while (iter.hasNext()){
		  Arc arc = (Arc)iter.next();
		  neighboringNodes.add(getOtherNode(arc));
	  }
	  return neighboringNodes;
  }

  /**
  Requires:
  Effects: return the other node connected by p_arc.   if p_arc is not connecting this node with any other node, return null.
  Modifies:
  */
  //## operation getOtherNode(Arc)
  public Node getOtherNode(Arc p_arc) {
      //#[ operation getOtherNode(Arc)
      return p_arc.getOtherNode(this);



      //#]
  }

  //## operation identifyFeElement()
  public void identifyFeElement() {
      //#[ operation identifyFeElement()
      if (feElement == null) {
      	updateFeElement();
      }


      //#]
  }

  //## operation identifyFgElement()
  public void identifyFgElement() {
      //#[ operation identifyFgElement()
      if (fgElement == null) {
      	updateFgElement();
      }


      //#]
  }

  //## operation includeFgElement(FGElement)
  public void includeFgElement(FGElement p_fgElement) {
      //#[ operation includeFgElement(FGElement)
      //#]
  }
  //## operation setID(int)
  public void setID(int p_ID) {
      //#[ operation setID(int)
      ID = new Integer(p_ID);
      //#]
  }

  //## operation includeFgElementInChemNodeElement(FGElement)
  public boolean includeFgElementInChemNodeElement(FGElement p_fgElement) {
      //#[ operation includeFgElementInChemNodeElement(FGElement)
      Object element = getElement();
      if (element instanceof Collection) {
      	Iterator iter = ((Collection)element).iterator();
      	while (iter.hasNext()) {
      		ChemNodeElement cne = (ChemNodeElement)iter.next();
      		if (cne instanceof FGAtom) {
      			FGElement fge = ((FGAtom)cne).getFgElement();
      			if (fge.equals(p_fgElement)) return true;
      		}
      	}
      }
      else if (element instanceof FGAtom) {
      	FGElement fge = ((FGAtom)element).getFgElement();
      	if (fge.equals(p_fgElement)) return true;
      }

      return false;


      //#]
  }

 

  /**
  Requires:
  Effects: super.isConnected()
  Modifies:
  */
  //## operation isConnected(Arc)
  public boolean isConnected(Arc p_arc) {
      //#[ operation isConnected(Arc)
      return super.isConnected(p_arc);
      //#]
  }

  /**
  Requires:
  Effects: if pass-in GraphComponent is not an instance of Arc, return false; otherwise, return super.isConnected().
  Modifies:
  */
  //## operation isConnected(GraphComponent)
  public boolean isConnected(GraphComponent p_GraphComponent) {
      //#[ operation isConnected(GraphComponent)
      // node can only connect to arc
      if (!(p_GraphComponent instanceof Arc)) return false;
      else return super.isConnected(p_GraphComponent);



      //#]
  }

  //## operation isConnectedBy(Node)
  public Arc isConnectedBy(Node p_node) {
      //#[ operation isConnectedBy(Node)
      if (this == p_node) return null;

      Iterator iter = getNeighbor();
      while (iter.hasNext()) {
      	Arc arc = (Arc)iter.next();
      	if (p_node.neighbor.contains(arc)) return arc;
      }
      return null;
      //#]
  }

  /**
  Requires:
  Effects: if this node only has one or zero neighbor, return true;  otherwise, return false.
  Modifies:
  */
  //## operation isLeaf()
  public boolean isLeaf() {
      //#[ operation isLeaf()
      return (neighbor.size() <= 1);



      //#]
  }

  /**
  Requires:
  Effects: if all the neighbors of this node are instances of Arc, return true; otherwise, return false.
  Modifies:
  */
  //## operation neighborOk()
  public boolean neighborOk() {
      //#[ operation neighborOk()
      // check if all the neighbors are node
      Iterator iter = getNeighbor();
      while (iter.hasNext()) {
      	if (!(iter.next() instanceof Arc)) return false;
      }

      return true;
      //#]
  }

  //## operation printFgElement()
  public String printFgElement() {
      //#[ operation printFgElement()
      if (fgElement == null) return "";

      if (fgElement instanceof Collection) {
      	String s = "";
      	for (Iterator iter = ((Collection)fgElement).iterator(); iter.hasNext();) {
      		FGElement fge = (FGElement)iter.next();
      		s += fge.getName() + ",";
      	}
      	s = s.substring(0,s.length()-1);
      	return s;
      }
      else {
      	String s = ((FGElement)fgElement).getName();
      	return s;
      }


      //#]
  }

  //## operation printLable()
  public String printLable() {
      //#[ operation printLable()
      String s = "N" + String.valueOf(getID());

      return s;

      //#]
  }

  /**
  Requires:
  Effects:check rep in the following aspects
  (1) neighborOk()?
  Modifies:
  */
  //## operation repOk()
  public boolean repOk() {
      //#[ operation repOk()
      return (neighborOk());
      //#]
  }

  

 
  //## operation toString()
  public String toString() {
      //#[ operation toString()
      String s = "";
      s = s + getID() + " ";
      if (isCentralNode()) s = s + "*" + getCentralID() + " ";
      Object o = getElement();
      s = s + " " + ChemParser.writeChemNodeElement(o);
      FreeElectron fee = getFeElement();
      String fe = null;
      if (fee == null) fe = "0";
      else fe = fee.getName();

      s = s + " " + fe;
      Iterator iter = getNeighbor();
      while (iter.hasNext()) {
      	s = s + " {";
      	Arc arc = (Arc)iter.next();
      	Node node = getOtherNode(arc);
      	s = s + node.getID() + ",";
      	s = s + arc.toString() + "}";
      }
      return s;


      //#]
  }

  //## operation toStringWithoutCentralID()
  public String toStringWithoutCentralID() {
      //#[ operation toStringWithoutCentralID()
      String s = "";
      s = s + getID() + " ";

      Object o = getElement();
      s = s + " " + ChemParser.writeChemNodeElement(o);
      FreeElectron fee = getFeElement();
      String fe = null;
      if (fee == null) fe = "0";
      else fe = fee.getName();

      s = s + " " + fe;
      Iterator iter = getNeighbor();
      while (iter.hasNext()) {
      	s = s + " {";
      	Arc arc = (Arc)iter.next();
      	Node node = getOtherNode(arc);
      	s = s + node.getID() + ",";
      	s = s + arc.toString() + "}";
      }
      return s;


      //#]
  }

  //## operation toStringWithoutCentralIDAndH()
  public String toStringWithoutCentralIDAndH() {
      //#[ operation toStringWithoutCentralIDAndH()
      String s = "";
      s = s + getID() + " ";

      Atom atom = (Atom)getElement();
      s = s + " " + ChemParser.writeChemNodeElement(atom);
      FreeElectron fee = getFeElement();
      String fe = null;
      if (fee == null) fe = "0";
      else fe = fee.getName();

      s = s + " " + fe;
      Iterator iter = getNeighbor();
      while (iter.hasNext()) {
      	Arc arc = (Arc)iter.next();
      	Node node = getOtherNode(arc);
      	Atom neighborAtom = (Atom)node.getElement();
      	if (!neighborAtom.isHydrogen()) {
      		s = s + " {";
      		s = s + node.getID() + ",";
      		s = s + arc.toString() + "}";
      	}
      }
      return s;


      //#]
  }

  //## operation updateFeElement()
  public void updateFeElement() throws InvalidFreeElectronException {
      //#[ operation updateFeElement()
      Object o = getElement();
      FreeElectron fe = null;
      if (o instanceof ChemNodeElement) {
      	fe = ((ChemNodeElement)o).getFreeElectron();
      }
      else if (o instanceof Collection) {
      	Collection oc = (Collection)o;
      	Iterator iter = oc.iterator();
      	FreeElectron old = null;
      	while (iter.hasNext()) {
      		ChemNodeElement cne = (ChemNodeElement)iter.next();
      		fe = cne.getFreeElectron();
      		if (old == null) old = fe;
      		else {
      			if (old != fe) throw new InvalidFreeElectronException();
      		}
      	}
      }
      if (fe == null) throw new InvalidFreeElectronException();
      feElement = fe;


      //#]
  }

  //## operation updateFgElement()
  public void updateFgElement() {
      //#[ operation updateFgElement()
      Collection fgec = generateFgElement();
      if (fgec.size() == 1) {
      	FGElement fge = (FGElement)(fgec.iterator().next());
      	setFgElement(fge);
      }
      else {
      	setFgElement(fgec);
      }


      //#]
  }

  protected void setID(Integer p_ID) {
      ID = p_ID;
  }

  

  public void setFeElement(FreeElectron p_feElement) {
      feElement = p_feElement;
  }

  public void setFgElement(Object p_fgElement) {
      fgElement = p_fgElement;
  }

}
/*********************************************************************
	File Path	: RMG\RMG\jing\chemUtil\Node.java
*********************************************************************/

