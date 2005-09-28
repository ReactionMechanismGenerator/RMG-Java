package bondGroups;

import java.io.*;
import java.util.*;
import jing.chem.*;
import jing.chemUtil.*;
import jing.chemParser.*;
import jing.chemUtil.HierarchyTree;
import jing.chemUtil.HierarchyTreeNode;
import jing.chem.FunctionalGroup;

//## package bondGroups

/**
 *
 * <p>Title: </p>
 * <p>Description: Singleton library for bond-centered groups</p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: </p>
 * @author Joanna Yu
 * @version 1.0
 *
 * BGThermoLibrary is a Composed Class of jing.chem.ThermoGAGroupLibrary
 *
 */



//## class BGThermoLibrary
public class BGThermoLibrary {

  // This is a singleton
  private static final BGThermoLibrary INSTANCE = new BGThermoLibrary();      //## attribute INSTANCE

  protected HashMap bGroupDictionary;                                         //## attribute bGroupDictionary
  protected HierarchyTree bGroupTree;                                         //## attribute bGroupTree
  protected HashMap bGroupLibrary;                                            //## attribute bGroupLibrary

  protected ThermoGAGroupLibrary thermoLibrary;                               //## attribute BGroupLibrary
  protected double lnValue;

  // Constructor

  //## operation BcThermoLibrary()
  private BGThermoLibrary() {
    //#[ operation BcThermoLibrary()
    bGroupTree = new HierarchyTree();
    bGroupDictionary = new HashMap();
    bGroupLibrary = new HashMap();

    thermoLibrary = ThermoGAGroupLibrary.getINSTANCE();

    String directory = System.getProperty("bondGroup.BGThermoLibrary.pathName");
    if (directory == null) {
        System.out.println("undefined system property: bondGroup.BGThermoLibrary.pathName, exit!");
        System.exit(0);
      }

    String separator = System.getProperty("file.separator");
    if (!directory.endsWith(separator)) directory = directory + separator;

    String bGDictionary = directory + "BondGroup_Dictionary.txt";
    String bGTree = directory + "BondGroup_Tree.txt";
    String bGLibrary = directory + "BondGroup_Library.txt";
    String bGlnValue = directory + "BondGroup_lnValue.txt";

    readBGroup(bGDictionary,bGTree,bGLibrary, bGlnValue);
    //#]
  }

  //## operation readBGroup(String,String,String)
  public void readBGroup(String p_bGroupDictionary, String p_bGroupTree, String p_bGroupLibrary, String p_lnValue){
    //#[ operation readBGroup(String,String,String)
    readBGroupDictionary(p_bGroupDictionary);

    readBGroupTree(p_bGroupTree);

    readBGroupLibrary(p_bGroupLibrary);

    readBGroupLnValue(p_lnValue);
    //#]
  }

  //## operation readBGroupDictionary(String)
  public void readBGroupDictionary(String p_fileName){
    //#[ operation readBGroupDictionary(String)
    try {
        bGroupDictionary = thermoLibrary.readStandardDictionary(p_fileName);
        return;
      }
    catch (Exception e) {
        System.err.println("Error in readBGroupDictionary!");
        System.err.println("Error: " + e.getMessage());
        System.exit(0);
      }
      //#]
  }

  //## operation readBGroupTree(String)
  public void readBGroupTree(String p_fileName){
    //#[ operation readBGroupTree(String)
    try {
            bGroupTree = thermoLibrary.readStandardTree(p_fileName,bGroupDictionary,0);
    }
    catch (Exception e) {
            System.err.println("Error in readBGroupTree!");
            System.err.println("Error: " + e.getMessage());
            System.exit(0);
    }
    //#]
}

  //## operation readGroupLibrary(String)
  public void readBGroupLibrary(String p_fileName) {
    //#[ operation readGroupLibrary(String)
    try {
          bGroupLibrary = thermoLibrary.readStandardLibrary(p_fileName, bGroupDictionary);
          return;
    }
    catch (Exception e) {
          System.err.println("Error in readBGroupLibrary!");
          System.err.println("Error: " + e.getMessage());
          System.exit(0);
    }
    //#]
  }

  //## operation readGroupLnValue(String)
  public void readBGroupLnValue(String p_fileName){
    //##[ operation readGroupLnValue(String)
    try {
            FileReader in = new FileReader(p_fileName);
            BufferedReader data = new BufferedReader(in);

            String line = ChemParser.readMeaningfulLine(data);
            StringTokenizer token = new StringTokenizer(line);
            String tokenString = token.nextToken();
            lnValue = Double.valueOf(tokenString.trim()).doubleValue();

            in.close();
            return;

    }
    catch (FileNotFoundException e) {
      System.err.println("File BondGroup_lnValue.txt was not found");
      System.exit(0);
    }
    catch (IOException e) {
      System.err.println("Error in reading BondGroup_lnValue.txt!");
      System.exit(0);
    }
    //#]
  }

  /**
   * Requires:
   * Effects: find a matched bond-centered group for the chemGraph and returns its thermo value.
   * This is just like the ThermoGAGroupLibrary.findGAGroup, except that it looks into the
   * BondGroup files.
   * Modifies:
   * @param p_chemGraph ChemGraph
   * @return ThermoGAValue
   */
  //## operation findBondGroup(ChemGraph)
  public ThermoGAValue findBondGroup(ChemGraph p_chemGraph){
    //#[ operation findBondGroup(ChemGraph)
    if(p_chemGraph==null) return null;

    System.out.println(" ");
    System.out.println("stack1: (#1=node"+(p_chemGraph.getCentralNodeAt(1)).getID()+"), (#2=node"+(p_chemGraph.getCentralNodeAt(2)).getID()+")");

  //  Node centralNode1 = p_chemGraph.getCentralNodeAt(1);
   // System.out.println(centralNode1);
  //  System.out.println(centralNode1.getFgElement());
   // Matchable fgCentralNode1 = (Matchable)centralNode1.getElement();
   // System.out.println(fgCentralNode1.getName());

    Stack stack1 = bGroupTree.findMatchedPath(p_chemGraph);

//    System.out.println("stack 1 size:"+stack1.size());

    // test if find a better match changing the numbering of the central nodes
    Node node1 = p_chemGraph.getCentralNodeAt(1);
    Node node2 = p_chemGraph.getCentralNodeAt(2);

    Graph graph = p_chemGraph.getGraph();
    graph.clearCentralNode();
    graph.setCentralNode(1, node2);
    graph.setCentralNode(2, node1);

    Stack stack2 = bGroupTree.findMatchedPath(p_chemGraph);

    Stack stack = new Stack();

    // debugging
    Stack otherStack = new Stack();


 //   System.out.println("stack 2 size:"+stack2.size());
// Choosing the correct stack
 HierarchyTreeNode nodeOne = (HierarchyTreeNode)stack1.peek();
 FunctionalGroup fgOne = (FunctionalGroup) nodeOne.getElement();
 System.out.println("fgOne "+ fgOne.getName());

 HierarchyTreeNode nodeTwo = (HierarchyTreeNode)stack2.peek();
 FunctionalGroup fgTwo = (FunctionalGroup) nodeTwo.getElement();
 System.out.println("fgTwo "+fgTwo.getName());

// System.out.println("Depth of stack1: "+nodeOne.getDepth());
// System.out.println("Depth of stack2: " + nodeTwo.getDepth());
 int nodeOneDepth = nodeOne.getDepth();
 int nodeTwoDepth = nodeTwo.getDepth();
 int deeper = 0;
 if (nodeOneDepth > nodeTwoDepth) {
   deeper = indentifyStack(nodeOne, nodeTwo);
   if(deeper==1){
     stack=stack1;
     otherStack=stack2;
   }
   else if(deeper==-1){
     stack=stack2;
     otherStack=stack1;
   }
   else System.out.println("Something went wrong when choosing the correct stack");
 }
 else {
   deeper = indentifyStack(nodeTwo, nodeOne);
   if (deeper == 1) {
     stack = stack2;
     otherStack = stack1;
   }
   else if (deeper == -1) {
     stack = stack1;
     otherStack = stack2;
   }
   else System.out.println("Something went wrong when choosing the correct stack");
 }


    if (stack==null) return null;

    Iterator stackIter = stack.iterator();
    System.out.println(" ");
    System.out.println("Chosen stack:");
    while(stackIter.hasNext()){
      HierarchyTreeNode nodeI = (HierarchyTreeNode)stackIter.next();
//      System.out.println(nodeI.getDepth());
      Matchable fgI = (Matchable)nodeI.getElement();
      System.out.println(fgI.getName());

      ThermoGAValue gaI = (ThermoGAValue)bGroupLibrary.get(fgI);
      System.out.println(gaI);
    }

    while(!stack.empty()){
      HierarchyTreeNode node = (HierarchyTreeNode)stack.pop();
      Matchable fg = (Matchable)node.getElement();
      ThermoGAValue ga = (ThermoGAValue)bGroupLibrary.get(fg);

      HierarchyTreeNode otherNode = (HierarchyTreeNode)otherStack.pop();
      Matchable otherFG = (Matchable)otherNode.getElement();
      ThermoGAValue otherGA = (ThermoGAValue)bGroupLibrary.get(otherFG);

      if (ga != null) {
        System.out.println("REAL VALUE USED " + ga);
  //      System.out.println("otherStack " + otherGA);
        System.out.println(otherNode.getDepth());
        return ga;
      }
    }

    return null;
    //#[
  }


  /**
   * Requires:
   * Effects: Returns 1 if the deeper node should be chosen, -1 if the shallower node
   * should be chosen and 0 if there is an error
   * Modifies:
   * @param p_nodeDeep HierarchyTreeNode, p_nodeShallow HierarchyTreeNode
   * @return int
   */
  //## operation indentifyStack(HierarchyTreeNode, HierarchyTreeNode)
  public int indentifyStack(HierarchyTreeNode p_nodeDeep, HierarchyTreeNode p_nodeShallow) {
    //#[ operation indentifyStack(HierarchyTreeNode, HierarchyTreeNode)

    int deep = p_nodeDeep.getDepth();
    int shallow = p_nodeShallow.getDepth();

    while (deep > shallow) {
      p_nodeDeep = (HierarchyTreeNode) p_nodeDeep.getFather();
      deep--;
    }

    // Now both nodes are at the same level, and we have to choose according
    // to when the parents of these nodes are siblings, which sibling comes 1st
    while (p_nodeDeep.getFather() != p_nodeShallow.getFather()) {
      p_nodeDeep = (HierarchyTreeNode) p_nodeDeep.getFather();
      p_nodeShallow = (HierarchyTreeNode) p_nodeShallow.getFather();
    }
    TreeNode father = p_nodeDeep.getFather();
    Iterator siblings = father.getChildren();

    while (siblings.hasNext()) {
      HierarchyTreeNode thisSibling = (HierarchyTreeNode) siblings.next();
      if (thisSibling == p_nodeDeep) {
        return 1;
      }
      else if (thisSibling == p_nodeShallow) {
        return -1;
      }
    }
    return 0;
    //#]
  }


  //## operation getInstance()
  public static BGThermoLibrary getInstance(){
    //#[ operation getInstance()
    return INSTANCE;
    //#]
  }

  //## operation getBGroupDictionary()
  public HashMap getBGroupDictionary(){
    //#[ operation getBGroupDictionary()
    return bGroupDictionary;
    //#]
  }

  //## operation getBGroupTree()
  public HierarchyTree getBGroupTree(){
    //#[ operation getBGroupTree()
    return bGroupTree;
    //#]
  }

  //## operation getBGroupLibrary()
  public HashMap getBGroupLibrary(){
    //#[ operation getBGroupLibrary()
    return bGroupLibrary;
    //#]
  }
}
