package bondGroups;

import java.util.*;
import java.lang.*;
import jing.param.*;
import jing.chem.*;
import jing.chemUtil.*;
import jing.chem.FGElement;
import jing.chemUtil.Node;
import jing.chem.Atom;
import jing.chemUtil.Graph;

//## package bondGroups

/**
 *
 * <p>Title: </p>
 * <p>Description: Singleton thermo estimator for aromatic and non-aromatic (calls GATP) molecules</p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: </p>
 * @author Joanna Yu
 * @version 1.0
 *
 * BGATP is a Composed Class of jing.chem.GATP and also of bondGroups.BGThermoLibrary
 *       and implements the GeneralGAPP interface
 */

//## class BGATP
public class BGATP
    implements GeneralGAPP {

  // This is a singleton
  private static final BGATP INSTANCE = new BGATP();

  protected BGThermoLibrary bGLibrary;

  protected GATP gATP;

  // Constructor

  //## operation BGATP()
  private BGATP() {
    //#[ operation BGATP()
    initGATP();
    initBGLibrary();
//#]
  }

  /**
   * Requires:
   * Effects:
   * 1) Separates the ChemGraph into an aromatic and a non-aromatic part
   * 2) Calls GATP.getGAGroup to calculate the thermo of the non-aromatic part
   * 3) Calls BGATP.getBGAGroup to calculate the thermo of the aromatic part
   * Modifies:
   *
   */
  //## operation generateThermoData(ChemGraph)
  public ThermoData generateThermoData(ChemGraph p_chemGraph) {
//#[ operation generateThermoData(ChemGraph)
    ThermoData result = new ThermoData();

    int flag = determineArom(p_chemGraph);

    if (flag == 1) { // molecule has both aromatic and non-aromatic parts

      //=====ATTENTION!=====
      // I will have to make some corrections for the symmetry, because
      // GATP will return the entropy considering the symmetry of the
      // non-aromatic part only.
      // For example, for toluene, what would it do?
      // I have to think what I want to do with the symmetry of the
      // aromatic part!
      //====================

      // The first element of this list is a list of the separated chemgraphs
      // that is, each chemgraph is either purely aromatic or purely non-aromatic
      // The second element is a ThermoData with all the thermo from the connections
      LinkedList list = separateNonArom_Arom(p_chemGraph);

      LinkedList separatedCG = (LinkedList) list.getFirst();
      ThermoData thermoFromConnections = (ThermoData) list.getLast();

      Iterator separatedIter = separatedCG.listIterator();

      result.plus(thermoFromConnections);

      while (separatedIter.hasNext()) {
        Graph graph = (Graph) separatedIter.next();
        //          System.out.println(graph);
        try {
          ChemGraph chemGraph = ChemGraph.make(graph);
          // recursive function. But all the chemgraphs that come in here will next time
          // either have a flag=2 or flag=3
          result.plus(generateThermoData(chemGraph));
        }
        catch (ForbiddenStructureException e) {
          System.err.println(
              "Failed to generate make ChemGraph in BGATP.generateThermoData");
        }
      }
    }

    else if (flag == 2) { // molecule is purely aromatic
      result.plus(getBGAGroup(p_chemGraph));
    }

    else if (flag == 3) { // molecule is purely non-aromatic
      // call Jing's thermo calculator
      result.plus(gATP.generateThermoData(p_chemGraph));
    }

    else {
      System.out.println(
          "Something went wrong when determining whether is aromatic or not");
      return null;
    }

    //int numberK = p_chemGraph.getKekule();
    //System.out.println("Number Kekule Structures: " + numberK);
    //ThermoGAValue lnContribution = new ThermoGAValue(bGLibrary.lnValue *
        //Math.log(numberK), 0, 0, 0, 0, 0, 0, 0, 0, null);
    //result.plus(lnContribution);
    return result;
    //#]
  }

  /**
   * Requires:
   * Effects:
   * Modifies:
   * @param p_chemGraph ChemGraph
   * @return LinkedList
   */
  //## operation separateNonArom_Arom(ChemGraph)
  protected LinkedList separateNonArom_Arom(ChemGraph p_chemGraph) {
//#[ operation separateNonArom_Arom(ChemGraph)
    LinkedList list = new LinkedList();

    ThermoData thermoFromConnections = new ThermoData();

    Graph graph = p_chemGraph.getGraph();

    Iterator nodeIter = p_chemGraph.getNodeList();

    // Make a HashSet with all the arcs that will have to be broken
    HashSet separateArc = new HashSet();
    while (nodeIter.hasNext()) {
      Node node = (Node) nodeIter.next();

      FGElement fge = (FGElement) node.getFgElement();

      if ( (fge.getName()).equals("Cb")) {
        Iterator arcIter = node.getNeighbor();
        Arc arc = (Arc) arcIter.next();
        Bond bond = (Bond) arc.getElement();
        if (bond.isSingle()) {
          Node neighborNode = node.getOtherNode(arc);
          FGElement neighborFge = (FGElement) neighborNode.getFgElement();
          if (! (neighborFge.getName()).equals("H")) {
            separateArc.add(arc);
          }
        }
      }
    }

    // add the contributions of the bonds being broken to
    // thermoFromConnections



    // Break bonds between aromatic and non-aromatic parts,
    // make single bonds to H, and adjust thermoFromConnections accordingly
    Iterator sepArcIter = separateArc.iterator();
    while (sepArcIter.hasNext()) {

      // This piece of code is repeated. Probably can do this more
      // elegantly
      Arc removedArc = (Arc) sepArcIter.next();
      Iterator nodeIter2 = removedArc.getNeighbor();
      Node node1 = (Node) nodeIter2.next();
      Node node2 = removedArc.getOtherNode(node1);
      graph.clearCentralNode();
      graph.setCentralNode(1, node1);
      graph.setCentralNode(2, node2);
      ThermoGAValue thisGAValue2 = bGLibrary.findBondGroup(p_chemGraph);
      thermoFromConnections.plus(thisGAValue2);

      //   Arc removedArc = (Arc)sepArcIter.next();
      graph.removeArc(removedArc);

      Atom H = Atom.make(ChemElement.make("H"), FreeElectron.make("0"));
      Bond S = Bond.make("S");

      Iterator sepNodeIter = removedArc.getNeighbor();
      while (sepNodeIter.hasNext()) {
        Node sepNode = (Node) sepNodeIter.next();
        //       int nodeNumber = graph.getNodeNumber();
        Node newH = graph.addNode(H);
        graph.addArcBetween(sepNode, S, newH); // This also throws some Exceptions, but they are not specified in the declaration of the method. What does this mean?
        sepNode.updateFgElement(); // I am not sure if I need this
        p_chemGraph.resetThermoSite(sepNode);
        ThermoGAValue thisGAValue = gATP.thermoLibrary.findGAGroup(p_chemGraph); // also throws a bunch of exceptions!
        if (thisGAValue == null) System.err.println("Thermo group not found: " +
            sepNode.getID());
        else thermoFromConnections.minus(thisGAValue);
      }
    }

//    System.out.println("just before separation");
//    System.out.println(graph);
    LinkedList separatedCG = graph.partition();

    //   System.out.println(chemgraph);

    list.add(separatedCG);
    list.add(thermoFromConnections);

    return list;
    //#]
  }

  /**
   * Requires:
   * Effects: Calculate the thermo values for aromatic molecules
   * Modifies:
   * @param p_chemGraph ChemGraph
   * @return ThermoGAValue
   */
  //## operation getBGAGroup(ChemGraph)
  protected ThermoGAValue getBGAGroup(ChemGraph p_chemGraph) {
    //#[ operation getBGAGroup(ChemGraph)

    ThermoData result = new ThermoData();

    Graph graph = p_chemGraph.getGraph();

    // deal with radicals

    Iterator iterArc = p_chemGraph.getArcList();

    // somehow add here that you don't need to check this if the arc is for
    // a C-H bond.

    while (iterArc.hasNext()) {
      Arc arc = (Arc) iterArc.next();
      Bond bond = (Bond) arc.getElement();

      Iterator nodeIter = arc.getNeighbor();
      Node node1 = (Node) nodeIter.next();
      Node node2 = arc.getOtherNode(node1);
      Atom atom1 = (Atom) node1.getElement();
      Atom atom2 = (Atom) node2.getElement();

      //     System.out.println(graph.getCentralNode());

      if ( (!atom1.isHydrogen()) && (!atom2.isHydrogen())) { // if none of the atoms is H
        graph.clearCentralNode();
        graph.setCentralNode(1, node1);
        graph.setCentralNode(2, node2);
        // check symmetry! might have to change numbering of nodes
        //       p_chemGraph.refreshCentralNode();
//        System.out.println(graph.getCentralNode());
        //      System.out.println(p_chemGraph.toStringWithoutH());
        //      System.out.println(p_chemGraph.getCentralNode());
        ThermoGAValue thisGAValue = bGLibrary.findBondGroup(p_chemGraph);

        // I am not sure I need the following piece of code, test for an
        // assymetric bond!
        /*
                 if (thisGAValue==null){
          System.out.println("thisGAValue==null");
          graph.setCentralNode(1,node2);
          graph.setCentralNode(1,node1);
          System.out.println(p_chemGraph.getCentralNode());
          thisGAValue = bGLibrary.findBondGroup(p_chemGraph);
                 }
         */

        result.plus(thisGAValue);
      }

      /**     if (bond.isBenzene()){
             Iterator centralNodesIter = arc.getNeighbor();
             HashMap centralNodesMap = new HashMap();
             int count = 1;
             while(centralNodesIter.hasNext()){
               Node centralNode = (Node)centralNodesIter.next();
               centralNodesMap.put("count", centralNode);
               count = count + 1;
             }
             p_chemGraph.getGraph().setCentralNodes(centralNodesMap);
             ThermoGAValue thisGAValue = bGLibrary.findBondGroup(p_chemGraph);
             result.plus(thisGAValue);
           }
       */
    }

    return result;
    //#]
  }

  /**
   * Requires:
   * Effects: Returns 1 if the ChemGraph has an aromatic part and a non-aromatic part
   *                  2 if the ChemGraph is purely aromatic
   *                  3 if the ChemGraph is purely non-aromatic
   * Modifies:
   * @param p_chemGraph ChemGraph
   * @return int
   */
  //## operation determineArom(ChemGraph)
  protected int determineArom(ChemGraph p_chemGraph) {
//#[ operation determineArom(ChemGraph)
    int flag = 0;

    Iterator nodeIter = p_chemGraph.getNodeList();

    while (nodeIter.hasNext()) {
      Node node = (Node) nodeIter.next();
      Atom atom = (Atom) node.getElement();
      FGElement fge = (FGElement) node.getFgElement();
      if (!atom.isHydrogen()) {

        if (atom.isOxygen()) {
          // if O is bonded to two Cb or Cbf, treat it as aromatic
          // else as not aromatic
          int flagOxygen = 0;
          Iterator arcIter = node.getNeighbor();
          while (arcIter.hasNext()) {
            Arc arc = (Arc) arcIter.next();
            Node neighborNode = node.getOtherNode(arc);
            FGElement neighborFge = (FGElement) neighborNode.getFgElement();
            if (neighborFge.getName().equals("Cb") ||
                neighborFge.getName().equals("Cbf")) flagOxygen++;
          }
          if (flagOxygen == 2) { //treat this Oxygen as aromatic
            if (flag == 0) flag = 2;
            else if (flag == 3)return 1;
          }
          else {
            if (flag == 0) flag = 3;
            else if (flag == 2)return 1;
          }
        } // if (atom.isOxygen())
        else if (fge.getName().equals("Cb") || fge.getName().equals("Cbf")) {
          if (flag == 0) flag = 2;
          else if (flag == 3)return 1;
        }

        else {
          if (flag == 0) flag = 3;
          else if (flag == 2)return 1;
        }
      }
    }
    return flag;
    //#]
  }

  /**
   * OK, I will leave this part for later!
   * Requires:
   * Effects: Returns the ThermoGAValue related to the radical site
   * Modifies: Saturates the ChemGraph radical sites with H
   *           That's a problem! I probably can't change the ChemGraph, because
   *           it will be used for other things in RMG!!!! So I will have to
   *           recover the ChemGraph with its radical site.
   *           The generateThermoData function will have to know how to find
   *           groups even if there is a radical there.
   * @param p_chemGraph ChemGraph
   * @return ThermoGAValue
   */
  //## operation getHBI(ChemGraph)
  /*
     public ThermoGAValue getHBI(ChemGraph p_chemGraph){

    // 1) locate all the radical sites

    // 2) add H to saturate radical sites

    // 3) find the HBI for all the radical groups
    //    will probably have to call different libraries for aromatic and non-aromatic radical sites

    //

     } */

  /**
   * Requires:
   * Effects:
   *
   * Modifies:
   * @param p_chemGraph ChemGraph
   * @return ThermoGAValue
   */
  //## operation getBGAGroup(ChemGraph)
  /*  public ThermoGAValue getBGAGroup(ChemGraph p_chemGraph){
    }
   */

//## operation initGATP()
  protected void initGATP() {
    //#[ operation initGAGroupLibrary()
    gATP = GATP.getINSTANCE();
    //#]
  }

//## operation initBGLibrary()
  protected void initBGLibrary() {
    //#[ operation initBGLibrary()
    bGLibrary = BGThermoLibrary.getInstance();
    //#]
  }

//## operation getInstance()
  public static BGATP getInstance() {
    //#[ operation getInstance()
    return INSTANCE;
    //#]
  }

}
