


/**
 * @author User1
 *
 */
/* 
 * GUI.java creates a Graphical User Interface (GUI) for Prof.
 * Green's Group's Reaction Mechanism Generator (RMG).  More information
 * on RMG may be found at: http://www.sourceforge.net/projects/rmg
 * 
 * Michael Harper Jr.
 * MIT Chemical Engineering
 * Green Group
 * 
 * Last updated: February 10, 2009
 * 
 * Currently:
 * 	- creates & displays GUI
 *	- creates .txt file related to running RMG
 *		* condition.txt : file containing all inputs (species, tolerances, etc.) 
 *	- requires 2 folders + 1 program to compile / run properly
 *		* ACD ChemSketch: free program available at http://www.acdlabs.com/download/chemsk.html
 * 
 * Improvements:
 * 
 */

import jing.gui.*;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ListCellRenderer;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;

import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.text.DecimalFormat;
import java.text.SimpleDateFormat;

import java.io.*;
import java.util.*;

import jing.chem.Species;
import jing.chemParser.ChemParser;

public class GUI extends JPanel implements ActionListener {
    public static void main(String[] args) {
    	String workingDir = System.getenv("RMG");
    	System.setProperty("RMG.workingDirectory", workingDir);
    	theApp = new GUI();
    	theApp.createAndShowGUI();
    }
    	
    private void createAndShowGUI() {
    	GUIWindow frame = new GUIWindow("RMG", this);
    	frame.setContentPane(theApp.mainPanel);
    	Dimension wndSize = frame.getToolkit().getScreenSize();
    	frame.setBounds(wndSize.width*3/16, wndSize.height/16, wndSize.width*5/8, wndSize.height*7/8);
    	frame.setVisible(true);
    }

    public GUI() {
    	super(new GridLayout(1,1));
    	//	Create main panel + individual tabs
    	JTabbedPane tabbedPanel = new JTabbedPane();
    	JComponent tabInputs = createTabInputs();
    	JComponent tabInitialization = createTabInitialization();
    	JComponent tabTermination = createTabTermination();
    	JComponent tabSolver = createTabSolver();
    	JComponent tabSensitivity = createTabSensitivity();
    	//	Add the individual tabs to the main panel
    	tabbedPanel.addTab("Initial Conditions", null, tabInputs, "Specify Reactants, Inerts, & Temperature/Pressure Model");
    	tabbedPanel.addTab("Thermochemical Libraries", null, tabInitialization, "Specify Thermochemical Library");
    	tabbedPanel.addTab("Termination", null, tabTermination, "Specify Simulation Termination Conditions");
    	tabbedPanel.addTab("Dynamic Simulator", null, tabSolver, "Specify Solver Tolerances");
    	tabbedPanel.addTab("Sensitivity Analysis", null, tabSensitivity, "Specify Error/Sensitivity Analysis Criteria");
    	tabbedPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	//	Create the main panel to contain the main tabbed panel
    	mainPanel = new JPanel();
    	mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.PAGE_AXIS));
    	mainPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	//	Create & fill boxes to display panels
    	Box Main1 = Box.createHorizontalBox();
    	Box Main2 = Box.createHorizontalBox();
    	Box Main3 = Box.createVerticalBox();
    	Main1.add(tabbedPanel);
    	JPanel submitConditionFile = new JPanel();
    	Main2.add(submitConditionFile);
    	Main3.add(Main1);
    	Main3.add(Main2);
    	mainPanel.add(Main3);
    	// Change preferences on Submit button
    	JPanel comments = new JPanel();
    	comments.add(new JLabel("Comments to add to file's header:"));
    	headerComments = new JTextArea(3,20);
    	JScrollPane headerScroll = new JScrollPane(headerComments);
    	comments.add(headerScroll);
    	comments.add(Submit = new JButton("Save"));
    	
    	submitConditionFile.add(comments);
    	submitConditionFile.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder
    			("Create Initialization File: condition.txt"), BorderFactory.createEmptyBorder(5,5,5,5)));
    	Submit.addActionListener(this);
    	Submit.setActionCommand("Condition");
    	Submit.setToolTipText("Press to save file: condition.txt");
    }
    
    public void initializeJCB(JComboBox[] all) {
    	ListCellRenderer centerJComboBoxRenderer = new CenterJCBRenderer();
        //TableCellRenderer centerTableRenderer = new CenterTRenderer();
    	for (int i=0; i<all.length; i++) {
    		all[i].addActionListener(this);
            all[i].setRenderer(centerJComboBoxRenderer);
            all[i].setSelectedIndex(0);
    	}
    }
    
    /* TabInputs: Initial Conditions
 	This tab allows the user to input the initial conditions
 		(species, temperature, & pressure) of the system
	The user must input the following information:
		- Each species
			* Identity
			* Concentration
			* Units
		- Temperature & Pressure
			* Model: Constant
			* Value(s) and Units(s)
	Additionally, the user is asked to specify the reactivity of each species:
		- Reactive: Normal species that RMG will react based on library
		- Unreactive: Species that RMG will only react based on IncludeSpecies.txt (see GJB work for more details)
		- Inert Gas:
	The user is asked to supply an adjacency list for the species if it is reactive.
    */
    public JComponent createTabInputs() {
    	//	Create cell renderers for JComboBox/JTable - CENTER text
    	ListCellRenderer centerJComboBoxRenderer = new CenterJCBRenderer();
        TableCellRenderer centerRenderer = new CenterTRenderer();
    	
        /*	Create the "Species" panel
        		- This panel will allow the user to enter information about each species.
        		- Each species must be given a name, concentration, units, and reactivity
        		- If the reactivity is "Reactive," the user must supply an adjacency list for the molecule
        			This may be done by giving the location of the species .mol file or InChI
        		- If the reactivity is "Unreactive," the user must supply the location of the IncludeSpecies.txt file
        		- All of this information is stored in a table which the user can add/remove data	
         */
        JPanel Species = new JPanel();
        
        //	Create the "speciesData" subpanel
        //		This subpanel holds the data for each species
    	JPanel speciesData = new JPanel();
    	
    	//	Create labels
    	JLabel labelName = new JLabel("Name");
    	JLabel labelInChI = new JLabel("InChI");
    	JLabel labelConc = new JLabel("Concentration");
    	JLabel labelUnit = new JLabel("Units");
    	JLabel labelReact = new JLabel("Reactivity");
    	
    	//	Create boxes to store each of the pieces of data
    	Box boxName = Box.createVerticalBox();
    	Box boxInchi = Box.createVerticalBox();
    	Box boxConc = Box.createVerticalBox();
    	Box boxUnit = Box.createVerticalBox();
    	Box boxReact = Box.createVerticalBox();
    	
    	//	Add labels and data fields to boxes
    	boxName.add(labelName);
    	boxName.add(NameInput = new JTextField());
    	NameInput.setPreferredSize(new Dimension (75,25));
    	NameInput.setHorizontalAlignment(JTextField.CENTER);
        NameInput.addActionListener(this);
    	
    	boxInchi.add(labelInChI);
    	boxInchi.add(InChIInput = new JTextField());
    	InChIInput.setPreferredSize(new Dimension (150,25));
    	InChIInput.setHorizontalAlignment(JTextField.CENTER);
        InChIInput.addActionListener(this);
    	
    	boxConc.add(labelConc);
    	boxConc.add(ConcInput = new JTextField());
    	ConcInput.setPreferredSize(new Dimension(30,25));
    	ConcInput.setHorizontalAlignment(JTextField.CENTER);
    	ConcInput.addActionListener(this);
    	
    	boxUnit.add(labelUnit);
    	boxUnit.add(UnitsInput = new JComboBox(concUnits));
    	UnitsInput.setSelectedIndex(0);
    	UnitsInput.setRenderer(centerJComboBoxRenderer);
    	UnitsInput.addActionListener(this);
    	
    	boxReact.add(labelReact);
    	boxReact.add(ReactiveInput = new JComboBox(reactive));
    	ReactiveInput.setRenderer(centerJComboBoxRenderer);
    	ReactiveInput.setSelectedIndex(0);
    	ReactiveInput.addActionListener(this);
    	
    	//	Add the boxes to the "speciesData" subpanel
    	speciesData.add(boxName);
    	speciesData.add(boxInchi);
    	speciesData.add(boxConc);
    	speciesData.add(boxUnit);
    	speciesData.add(boxReact);
    	
    	//	Create the "speciesButton" panel
    	//		This subpanel holds the add/remove buttons
    	JPanel speciesButton = new JPanel();

    	//	Create the buttons: "Add" and "Remove"
    	JButton AddInput = new JButton("Add");
    	AddButtonListener addListenerInput = new AddButtonListener();
    	AddInput.addActionListener(addListenerInput);
    	AddInput.setActionCommand("AddInput");
    	AddInput.setToolTipText("Press to submit data for molecule");
    	
    	JButton DeleteInput = new JButton("Remove");
    	DeleteButtonListener deleteListenerInput = new DeleteButtonListener();
    	DeleteInput.addActionListener(deleteListenerInput);
    	DeleteInput.setActionCommand("DeleteInput");
    	DeleteInput.setToolTipText("Press to delete molecule's data");
    	
    	//	Add the box to the "speciesButton" subpanel
    	speciesButton.add(AddInput);
    	speciesButton.add(DeleteInput);

    	//	Create the "speciesTable" subpanel
    	JPanel speciesTable = new JPanel();
    	
    	//	Create table to hold the species data
    	tableInput = new JTable(tmodelInput = new MyTableModelInput());
    	tableInput.setPreferredScrollableViewportSize(new Dimension(500,50));
    	tableInput.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        tableInput.getColumnModel().getColumn(0).setPreferredWidth(75);
        tableInput.getColumnModel().getColumn(1).setPreferredWidth(100);
        tableInput.getColumnModel().getColumn(2).setPreferredWidth(75);
        tableInput.getColumnModel().getColumn(3).setPreferredWidth(100);
        tableInput.getColumnModel().getColumn(4).setPreferredWidth(150);
        
        //	Create a scroll panel and add it to "speciesTable" subpanel
    	JScrollPane scrollInput = new JScrollPane(tableInput);
        scrollInput.setBorder(BorderFactory.createLoweredBevelBorder());
        for (int i=0; i<tableInput.getColumnCount(); i++) {
        	TableColumn column = tableInput.getColumnModel().getColumn(i);
        	column.setCellRenderer(centerRenderer);
        }
        speciesTable.add(scrollInput);
        
        //	Create the "speciesAList" subpanel
        JPanel speciesAList = new JPanel();
        
        //	Create text field where the user can visualize the adjacency list
        speciesAdjList = new JTextArea(5,13);
        //	Add text field to a scroll panel
        JScrollPane scrollAdjList = new JScrollPane(speciesAdjList);
    	scrollAdjList.setBorder(BorderFactory.createLoweredBevelBorder());
    	
    	Box adj = Box.createVerticalBox();
    	adj.add(new JLabel("Adjacency List"));
    	adj.add(scrollAdjList);
    	
    	speciesAList.add(adj);
        
        //	Create the "speciesIS" subpanel
        //		If at least one of the species is "Unreactive," the user must specify the path of IncludeSpecies.txt
        JPanel speciesIS = new JPanel();
        
        //	Create the label for the speciesIS subpanel
        JLabel labelIS = new JLabel("Location of IncludeSpecies.txt: ");
        
        //	Populate the speicesIS subpanel
        speciesIS.add(labelIS);
        speciesIS.add(isPath = new JTextField(20));
        speciesIS.add(isButton = new JButton("Change"));
    	ChangeButtonListener addListenerIS = new ChangeButtonListener();
    	isButton.addActionListener(addListenerIS);
    	isButton.setActionCommand("isPath");
        //	Update the subpanel's properties
        isPath.addActionListener(this);
        isPath.setEditable(false);
        
        //	Create boxes to store all of the species subpanels
    	Box species1 = Box.createHorizontalBox();
    	Box species2 = Box.createHorizontalBox();
    	Box species3 = Box.createHorizontalBox();
    	Box species4 = Box.createHorizontalBox();
    	Box species5 = Box.createHorizontalBox();
    	Box speciesTotal = Box.createVerticalBox();
    	
    	//	Add the subpanels to the boxes
    	species1.add(speciesData);
    	species2.add(speciesButton);
    	species2.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
    	species3.add(speciesTable);
    	species4.add(speciesAList);
    	species4.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
    	species5.add(speciesIS);
        speciesTotal.add(species1);
        speciesTotal.add(species4);
        speciesTotal.add(species2);
        speciesTotal.add(species3);
        speciesTotal.add(species5);
        
        // Add the boxes to the main panel "Species"
        Species.add(speciesTotal);
        Species.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Species"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
        
        /*	Create the "Temperature" panel
        		The user must specify which temperature model to use.
        		Note: As of 10-Feb-2009, the only option is "Constant"
       			The user must also specify a list of temperatures to build the model
       	*/
        JPanel Temperature = new JPanel();
        
        //	Create the "tempModel" subpanel
        JPanel tempModel = new JPanel();
        tempModel.setBorder(BorderFactory.createEmptyBorder(0,0,10,0));
        
    	//	Create the labels
    	JLabel tempModelLabel = new JLabel("Select Model: ");
    	//	Add the labels to the panel
    	tempModel.add(tempModelLabel);
    	tempModel.add(tempModelCombo = new JComboBox(tpModel));
    	//	Adjust the labels properties
        tempModelCombo.setSelectedIndex(0);
    	tempModelCombo.addActionListener(this);
    	tempModelCombo.setRenderer(centerJComboBoxRenderer);
    	
    	//	Create the "tempValues" subpanel
    	JPanel tempValues = new JPanel();
    	
    	//	Create the labels
    	JLabel valueTLabel = new JLabel("List of Temperatures");
    	JLabel unitTLabel = new JLabel("Units");
    	
    	//	Create the data entries and edit their properties
    	tempConstant = new JTextArea(3,10);
    	
    	tempConstantUnit = new JComboBox(tempUnits);
    	tempConstantUnit.addActionListener(this);
    	tempConstantUnit.setRenderer(centerJComboBoxRenderer);
    	tempConstantUnit.setSelectedIndex(0);
    	
    	//	Create scroll panel for list of temperatures
    	JScrollPane tempScroll = new JScrollPane(tempConstant);
    	tempScroll.setBorder(BorderFactory.createLoweredBevelBorder());
    	
    	//	Create boxes to store labels and data
    	Box leftT = Box.createVerticalBox();
    	Box rightT = Box.createVerticalBox();
    	//	Add labels and data to boxes
    	leftT.add(valueTLabel);
    	leftT.add(tempScroll);
    	rightT.add(unitTLabel);
    	rightT.add(tempConstantUnit);
    	
    	//	Add boxes to "tempValues" subpanel
    	tempValues.add(leftT);
    	tempValues.add(rightT);
    	
    	//	Create boxes for Temperature panel
    	Box temp1 = Box.createHorizontalBox();
    	Box temp2 = Box.createHorizontalBox();
    	Box temp = Box.createVerticalBox();
    	//	Populate boxes
    	temp1.add(tempModel);
    	temp2.add(tempValues);
    	temp.add(temp1);
    	temp.add(temp2);
    	//	Add box to "Temperature" panel
    	Temperature.add(temp);
    	Temperature.setBorder(BorderFactory.createTitledBorder("Temperature"));

    	/*	Create the "Pressure" panel
		The user must specify which pressure model to use.
		Note: As of 10-Feb-2009, the only option is "Constant"
			The user must also specify a list of pressures to build the model
		*/
		JPanel Pressure = new JPanel();
		
		//	Create the "pressModel" subpanel
		JPanel pressModel = new JPanel();
		pressModel.setBorder(BorderFactory.createEmptyBorder(0,0,10,0));
		
		//	Create the labels
		JLabel pressModelLabel = new JLabel("Select Model: ");
		//	Add the labels to the panel
		pressModel.add(pressModelLabel);
		pressModel.add(pressModelCombo = new JComboBox(tpModel));
		//	Adjust the labels properties
		pressModelCombo.setSelectedIndex(0);
		pressModelCombo.addActionListener(this);
		pressModelCombo.setRenderer(centerJComboBoxRenderer);
		
		//	Create the "pressValues" subpanel
		JPanel pressValues = new JPanel();
		
		//	Create the labels
		JLabel valuePLabel = new JLabel("List of Pressures");
		JLabel unitPLabel = new JLabel("Units");
		
		//	Create the data entries and edit their properties
		pressConstant = new JTextArea(3,10);
		
		pressConstantUnit = new JComboBox(pressUnits);
		pressConstantUnit.addActionListener(this);
		pressConstantUnit.setRenderer(centerJComboBoxRenderer);
		pressConstantUnit.setSelectedIndex(0);
		
		//	Create scroll panel for list of pressures
		JScrollPane pressScroll = new JScrollPane(pressConstant);
		pressScroll.setBorder(BorderFactory.createLoweredBevelBorder());
		
		//	Create boxes to store labels and data
		Box leftP = Box.createVerticalBox();
		Box rightP = Box.createVerticalBox();
		//	Add labels and data to boxes
		leftP.add(valuePLabel);
		leftP.add(pressScroll);
		rightP.add(unitPLabel);
		rightP.add(pressConstantUnit);
		
		//	Add boxes to "pressValues" subpanel
		pressValues.add(leftP);
		pressValues.add(rightP);
		
		//	Create boxes for Pressure panel
		Box press1 = Box.createHorizontalBox();
		Box press2 = Box.createHorizontalBox();
		Box press = Box.createVerticalBox();
		//	Populate boxes
		press1.add(pressModel);
		press2.add(pressValues);
		press.add(press1);
		press.add(press2);
		//	Add box to "Pressure" panel
		Pressure.add(press);
		Pressure.setBorder(BorderFactory.createTitledBorder("Pressure"));
    	
    	//	Create boxes for the PressAndTemp combined panel
    	Box PressAndTemp = Box.createHorizontalBox();
    	PressAndTemp.add(Temperature);
    	PressAndTemp.add(Pressure);
    	
    	//	Create the PressAndTemp combined panel
    	JPanel tempandpress = new JPanel();
    	tempandpress.add(PressAndTemp);
    	tempandpress.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Temperature & Pressure"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	
    	//	Create boxes for the entire Inputs tab
    	Box TabTotal = Box.createVerticalBox();
    	TabTotal.add(Species);
    	TabTotal.add(tempandpress);
    	
    	//	Create the panel for the entire Input tab
    	JPanel input = new JPanel();
    	input.add(TabTotal);
        input.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Initial Conditions"),
        		BorderFactory.createEmptyBorder(5,5,5,5)));
        JScrollPane scrolltab2 = new JScrollPane(input);
        scrolltab2.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	return scrolltab2;
    }

    /* Tab0: Initialization
 	This tab allows the user to input the initialization specifications of the simulation.
	The user must input the following information:
		- Whether to restart the prior simulation
		- The database
		- The primary thermodynamic library
		- The reaction model enlarger
		- The dynamic simulator
		 	* The time steps for integration
		 	* The absolute tolerance
		 	* The relative tolerance
		- The primary reaction library
			* Name
			* Location
    */
    private JComponent createTabInitialization() {
    	//	Create cellrenderers for JComboBox/JTable - CENTER text
    	ListCellRenderer centerJComboBoxRenderer = new CenterJCBRenderer();
        TableCellRenderer centerTableRenderer = new CenterTRenderer();    	

    	//	Create the "Database" panel
    	JPanel Database = new JPanel();
    	Database.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Main database"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	
    	//	Populate the "Database" panel
    	JLabel databaseLabel = new JLabel("Choose database");
    	Database.add(databaseLabel);
    	databaseLabel.setToolTipText("Default = RMG/databases/RMG_database");
    	
    	Database.add(databasePath = new JTextField(25));
    	databasePath.setText(System.getProperty("RMG.workingDirectory") + 
			"/database/RMG_database");
    	databasePath.setEditable(false);
    	
    	Database.add(databaseButton = new JButton("Change"));
    	ChangeButtonListener addListenerDatabase = new ChangeButtonListener();
    	databaseButton.addActionListener(addListenerDatabase);
    	databaseButton.setActionCommand("databasePath");
    	
    	//	Craete the Primary Thermodynamic Library (PTL) panel
    	JPanel PTL = new JPanel();
    	PTL.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Primary Thermodynamic Library"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	
    	//	Populate the PTL panel
    	JLabel ptlLabel = new JLabel("Choose thermodynamic library");
    	PTL.add(ptlLabel);
    	ptlLabel.setToolTipText("Default = RMG/databases/RMG_database/thermo/primaryThermoLibrary");
    	
    	PTL.add(ptlPath = new JTextField(25));
    	ptlPath.setText(databasePath.getText() + "/thermo/PrimaryThermoLibrary");
    	ptlPath.setEditable(false);
    	
    	PTL.add(ptlButton = new JButton("Change"));
    	ChangeButtonListener addListenerPTL = new ChangeButtonListener();
    	ptlButton.addActionListener(addListenerPTL);
    	ptlButton.setActionCommand("ptlPath");
    	
    	//	Create the Primary Reaction Library (PRL) panel
    	JPanel PRL = new JPanel();
    	PRL.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Primary Reaction Library"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	
    	//	Create PRL On/Off switch label
    	JPanel prlOnOff = new JPanel();
    	prlOnOff.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	
    	//	Populate the On/Off label
    	JLabel libraryLabel = new JLabel("Enable reaction library");
    	prlOnOff.add(libraryLabel);
    	
    	prlOnOff.add(libraryCombo = new JComboBox(yesnoOptions));
        libraryCombo.setRenderer(centerJComboBoxRenderer);
        libraryCombo.setSelectedIndex(0);
        libraryCombo.setBorder(BorderFactory.createEmptyBorder(0,5,0,10));
        libraryCombo.setActionCommand("prlOnOff");
        libraryCombo.addActionListener(this);

    	//	Create PRL Name label
    	JPanel prlName = new JPanel();
    	prlName.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	
    	//	Populate the label
    	JLabel nameLabel = new JLabel("Name:");
    	prlName.add(nameLabel);
    	prlName.add(libName = new JTextField(20));
    	libName.setEnabled(false);

        //	Create PRL Location label
        JPanel prlLoc = new JPanel();
        prlLoc.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	
        //	Populate the label
        JLabel locationLabel = new JLabel("Location:");
        prlLoc.add(locationLabel);
    	locationLabel.setToolTipText("Default = " + 
			"RMG/databases/RMG_database/primaryReactionLibrary/combustion_core/version5");
    	
    	prlLoc.add(prlPath = new JTextField(20));
    	prlPath.setEnabled(false);
    	
    	prlLoc.add(prlButton = new JButton("Change"));
    	ChangeButtonListener addListenerLib = new ChangeButtonListener();
    	prlButton.addActionListener(addListenerLib);
    	prlButton.setActionCommand("prlPath");
    	prlButton.setEnabled(false);

        //	Create table and scroll panel to store PRL(s)
        tablePRL = new JTable(tmodelPRL = new MyTableModelPRL());
    	tablePRL.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    	tablePRL.setPreferredScrollableViewportSize(new Dimension(500,50));
        tablePRL.getColumnModel().getColumn(0).setPreferredWidth(150);
        tablePRL.getColumnModel().getColumn(1).setPreferredWidth(350);
        tablePRL.setEnabled(false);
        for (int i=0; i<tablePRL.getColumnCount(); i++) {
        	TableColumn column = tablePRL.getColumnModel().getColumn(i);
        	column.setCellRenderer(centerTableRenderer);
        }
        JScrollPane scrollPRL = new JScrollPane(tablePRL);
    	scrollPRL.setBorder(BorderFactory.createLoweredBevelBorder());

    	//	Create boxes to display the PRL table
    	Box PRLtable1 = Box.createHorizontalBox();
    	Box PRLtable2 = Box.createHorizontalBox();
    	Box PRLtable3 = Box.createHorizontalBox();
    	Box PRLtable4 = Box.createHorizontalBox();
    	Box PRLtable5 = Box.createVerticalBox();
    	
    	//	Fill the boxes with the appropriate components of the table
    	PRLtable1.add(AddPRL = new JButton("Add"));
    	AddButtonListener addListenerPRL = new AddButtonListener();
    	AddPRL.addActionListener(addListenerPRL);
    	AddPRL.setActionCommand("AddPRL");
    	AddPRL.setEnabled(false);
    	AddPRL.setToolTipText("Press to submit PRL Name & Location");
    	PRLtable1.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
    	
    	PRLtable2.add(DeletePRL = new JButton("Remove"));
    	DeleteButtonListener deleteListenerPRL = new DeleteButtonListener();
    	DeletePRL.addActionListener(deleteListenerPRL);
    	DeletePRL.setActionCommand("DeletePRL");
    	DeletePRL.setEnabled(false);
    	DeletePRL.setToolTipText("Press to remove PRL Name & Location");
    	
    	PRLtable3.add(PRLtable1);
    	PRLtable3.add(PRLtable2);
    	PRLtable4.add(scrollPRL);
    	PRLtable5.add(prlOnOff);
    	PRLtable5.add(prlName);
    	PRLtable5.add(prlLoc);
    	PRLtable5.add(PRLtable3);
    	PRLtable5.add(PRLtable4);
    	
    	PRL.add(PRLtable5);

    	//	Create & fill box for Initialization tab
    	Box TabTotal = Box.createVerticalBox();
    	TabTotal.add(Database);
    	TabTotal.add(PTL);
    	TabTotal.add(PRL);
        
    	//	Create the thermochemLibrary panel
    	JPanel thermochemLibrary = new JPanel();
    	thermochemLibrary.add(TabTotal);
    	thermochemLibrary.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder(
			"Thermochemistry Libraries"), BorderFactory.createEmptyBorder(10,10,10,10)));
        //	Embed the thermochemistry panel in a scroll panel
    	JScrollPane scrolltab = new JScrollPane(thermochemLibrary);
        scrolltab.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	return scrolltab;
	}
    
    protected JComponent createTabSolver() {
    	//	Create cellrenderers for JComboBox/JTable - CENTER text
    	ListCellRenderer centerJComboBoxRenderer = new CenterJCBRenderer();
    	
    	//	Create the Model Enlarger (ME) panel
    	JPanel ME = new JPanel();
    	ME.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Model Enlarger"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	//	Populate the ME panel
    	JLabel modelLabel = new JLabel("Choose reaction model enlarger");
    	ME.add(modelLabel);
    	modelLabel.setToolTipText("Default = RateBasedModelEnlarger");
    	ME.add(modelCombo = new JComboBox(meOptions));
    	modelCombo.setActionCommand("PDepME");
    	
    	//	Create the Kinetics Estimator (KE) panel
    	JPanel KE = new JPanel();
    	KE.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Kinetic Estimator"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	//	Populate the KE panel
    	JLabel estimatorLabel = new JLabel("Choose kinetics estimator");
    	KE.add(estimatorLabel);
    	KE.add(estimatorCombo = new JComboBox(keOptions));
    	estimatorCombo.setEnabled(false);
    	
    	//	Create the Dynamic Simulator (DS) panel
    	JPanel DS = new JPanel();
    	DS.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Dynamic Simulator"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	
    	//	Create the DS subpanel: Solver
    	JPanel dsSolver = new JPanel();
    	dsSolver.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	
    	//	Populate the Solver subpanel
    	JLabel solverLabel = new JLabel("Select Dynamic Simulator");
    	dsSolver.add(solverLabel);
    	solverLabel.setToolTipText("Default = DASSL");
    	
    	dsSolver.add(simulatorCombo = new JComboBox(simOptions));    	
        simulatorCombo.addActionListener(this);
        simulatorCombo.setRenderer(centerJComboBoxRenderer);
        simulatorCombo.setSelectedIndex(0);
        simulatorCombo.setActionCommand("saOnOff");
    	
    	//	Create the DS subpanel: aTolPanel
    	JPanel aTolPanel = new JPanel(); 
    	aTolPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	
    	//	Populate the aTolPanel subpanel
    	JLabel aToleranceLabel = new JLabel("Select absolute tolerance");
    	aTolPanel.add(aToleranceLabel);
    	aToleranceLabel.setToolTipText("Suggested value = 1E-12");
    	
    	aTolPanel.add(aTolerance = new JTextField());
    	aTolerance.setPreferredSize(new Dimension(100,25));
        aTolerance.addActionListener(this);
        aTolerance.setHorizontalAlignment(JTextField.CENTER);
        aTolerance.setText("1E-12");
    	
    	//	Create the DS subpanel: rTolPanel
        JPanel rTolPanel = new JPanel();
        rTolPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
        
        //	Populate the rTolPanel subpanel
        JLabel rToleranceLabel = new JLabel("Select relative tolerance");
        rTolPanel.add(rToleranceLabel);
        rToleranceLabel.setToolTipText("Suggested value = 1E-3");
    	
        rTolPanel.add(rTolerance = new JTextField());
    	rTolerance.setPreferredSize(new Dimension(100,25));
        rTolerance.addActionListener(this);
        rTolerance.setHorizontalAlignment(JTextField.CENTER);
        rTolerance.setText("1E-3");
    	
        //	Create the DS subpanel: interConv
        JPanel interConv = new JPanel();
        interConv.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
        //	Populate the subpanel
        JLabel multiConvLabel = new JLabel("Intermediate Conversions");
        interConv.add(multiConvLabel);
        multiConvLabel.setToolTipText("Default = AUTO");
        interConv.add(multiConv = new JTextArea(2,15));
        
    	//	Create the DS subpanel: interTS
    	JPanel interTS = new JPanel();
    	interTS.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	//	Populate the subpanel
        JLabel multiTSLabel = new JLabel("Intermediate Time Steps");
        interTS.add(multiTSLabel);
        multiTSLabel.setToolTipText("Default = AUTO");
        interTS.add(multiTS = new JTextArea(2,15));
    	multiTS.setEnabled(false);
    	
    	//	Create boxes for DS panel
        Box ds = Box.createVerticalBox();
        //	Populate the boxes
    	ds.add(dsSolver);
    	ds.add(aTolPanel);
    	ds.add(rTolPanel);
    	ds.add(interConv);
    	ds.add(interTS);
    	DS.add(ds);
    	
    	//	Create the Spectroscopic Data Estimator (SDE) subpanel
    	JPanel SDE = new JPanel();
    	SDE.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Spectroscopic Data Estimator"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	//	Populate the subpanel
    	JLabel sdeLabel = new JLabel("Spectroscopic Data Estimator");
    	SDE.add(sdeLabel);
    	SDE.add(sdeCombo = new JComboBox(sdeOptions));
    	
    	//	Create the Equation Of State (EOS) subpanel
    	JPanel EOS = new JPanel();
    	EOS.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Equation of State"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	//	Populate the subpanel
    	JLabel eosLabel = new JLabel("Equation of State");
    	EOS.add(eosLabel);
    	EOS.add(eosCombo = new JComboBox(eosOptions));
    	
    	//	Create the InChI generation subpanel
    	JPanel inchiPanel = new JPanel();
    	inchiPanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("InChI Generation"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	//	Populate the subpanel
    	JLabel inchiLabel = new JLabel("Generate InChIs");
    	inchiPanel.add(inchiLabel);
    	inchiPanel.add(inchiCombo = new JComboBox(yesnoOptions));
    	
        JComboBox[] allTab = {modelCombo, estimatorCombo, sdeCombo, eosCombo, inchiCombo};
        initializeJCB(allTab);
    	
    	Box TabTotal = Box.createVerticalBox();
    	TabTotal.add(ME);
    	TabTotal.add(KE);
    	TabTotal.add(DS);
    	TabTotal.add(SDE);
    	TabTotal.add(EOS);
    	TabTotal.add(inchiPanel);
    	
    	//	Create the dynamicSimulator panel
    	JPanel dynamicSimulator = new JPanel();
    	dynamicSimulator.add(TabTotal);
    	dynamicSimulator.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder(
			"Solver Options"), BorderFactory.createEmptyBorder(10,10,10,10)));
        //	Embed the dynamicSimulator panel in a scroll panel
    	JScrollPane scrolltab = new JScrollPane(dynamicSimulator);
        scrolltab.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	return scrolltab;
    }
    
    /* Tab1: Termination Sequence
 	This tab allows the user to input the terminating sequence of the simulation.
	The user must input the followig information:
		- The goal of the reaction
			* The limiting reactant & its associated conversion
			* The time of reaction
		- The error tolerance
    */
	protected JComponent createTabTermination() {
    	//	Create cellrenderers for JComboBox/JTable - CENTER text
    	ListCellRenderer centerJComboBoxRenderer = new CenterJCBRenderer();
    	
        //	Create "Goal" panel
        JPanel Goal = new JPanel();
    	Goal.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Termination Goal"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
        //	Populate the panel
    	JLabel controllerLabel = new JLabel("Goal of reaction");
    	Goal.add(controllerLabel);    	
    	Goal.add(controllerCombo = new JComboBox(goalOptions));
    	controllerCombo.setActionCommand("GoalofReaction");
        controllerCombo.addActionListener(this);
        controllerCombo.setRenderer(centerJComboBoxRenderer);

        // "If Goal = Conversion" subpanel
        JPanel Conversion = new JPanel();
    	Conversion.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Conversion"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	
    	//	Create conversionName subpanel
    	JPanel conversionName = new JPanel();
    	JLabel nameCLabel = new JLabel("Species");
    	conversionName.add(nameCLabel);
    	conversionName.add(SpeciesConvName = new JComboBox(species_blank));
    	
    	//	Create conversionValue subpanel
    	JPanel conversionValue = new JPanel();
    	JLabel convCLabel = new JLabel("Fractional Conversion");
    	conversionValue.add(convCLabel);
    	conversionValue.add(conversion = new JTextField());
    	conversion.setPreferredSize(new Dimension(50,25));
        conversion.addActionListener(this);
        conversion.setHorizontalAlignment(JTextField.CENTER);
    	
    	Box convBox = Box.createVerticalBox();
    	convBox.add(conversionName);
    	convBox.add(conversionValue);
    	Conversion.add(convBox);
    	
    	// "If Goal = ReactionTime" subpanel
    	JPanel ReactionTime = new JPanel();
    	ReactionTime.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Reaction Time"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));

    	//	Populate the subpanel
    	ReactionTime.add(new JLabel("Input reaction time"));
    	ReactionTime.add(time = new JTextField());
    	time.setPreferredSize(new Dimension(50,25));
    	time.setEnabled(false);
        time.addActionListener(this);
        time.setHorizontalAlignment(JTextField.CENTER);
    	ReactionTime.add(timeCombo = new JComboBox(timeUnits));
    	timeCombo.setEnabled(false);
        timeCombo.addActionListener(this);
        timeCombo.setRenderer(centerJComboBoxRenderer);
        timeCombo.setSelectedIndex(1);
        
        //	Finish controller panel
        JPanel Finish = new JPanel();
        Finish.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Select Finish Controller"),
        		BorderFactory.createEmptyBorder(5,5,5,5)));
        
        JLabel finishcLabel = new JLabel("Select finish controller");
        Finish.add(finishcLabel);
        
        Finish.add(finishcCombo = new JComboBox(fcOptions));
        finishcCombo.addActionListener(this);
        finishcCombo.setSelectedIndex(0);
        finishcCombo.setRenderer(centerJComboBoxRenderer);
    	
        //	Error tolerance subpanel
        JPanel eTolPanel = new JPanel();
    	eTolPanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Error tolerance"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	//	Populate subpanel
        JLabel eToleranceLabel = new JLabel("Select error tolerance");
    	eToleranceLabel.setToolTipText("Suggested value = 0.1");
    	eTolPanel.add(eToleranceLabel);
    	
    	eTolPanel.add(eTolerance = new JTextField());
    	eTolerance.setPreferredSize(new Dimension(100,25));
        eTolerance.addActionListener(this);
        eTolerance.setHorizontalAlignment(JTextField.CENTER);
        eTolerance.setText("0.1");
  
        //	Create boxes to store subpanels
    	Box term = Box.createVerticalBox();
    	Box goalTypes = Box.createHorizontalBox();
    	
    	goalTypes.add(Conversion);
    	goalTypes.add(ReactionTime);
    	term.add(Goal);
    	term.add(goalTypes);
    	term.add(Finish);
    	term.add(eTolPanel);
    	
    	JPanel terminationPanel = new JPanel();
    	terminationPanel.add(term);
    	terminationPanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Finish Controller"),
    			BorderFactory.createEmptyBorder(5,5,5,5))); 
        JScrollPane scrolltab = new JScrollPane(terminationPanel);
        scrolltab.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	return scrolltab;
    }
    
    /* Tab3: Error Analysis
 		This tab allows the user to choose whether error bars & sensitivity coefficients (and for which
 			speices) are available in the output .txt files
		The user must input the following information:
		- Whether to display error bars
		- Whether to display sensitivity coefficients
			* Which species to analyze
	*/ 
    protected JComponent createTabSensitivity() {
    	//	Create cellrenderers for JComboBox/JTable - CENTER text
    	ListCellRenderer centerJComboBoxRenderer = new CenterJCBRenderer();
        TableCellRenderer centerRenderer = new CenterTRenderer();
    	
        //	Create errorbars (EB) subpanel
        JPanel EB = new JPanel();
    	EB.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Error Bars"),
    			BorderFactory.createEmptyBorder(10,10,10,10)));
    	//	Populate the subpanel
        EB.add(new JLabel("Display error bars"));
        EB.add(eBarsCombo = new JComboBox(yesnoOptions));
        eBarsCombo.addActionListener(this);
        eBarsCombo.setRenderer(centerJComboBoxRenderer);
        eBarsCombo.setSelectedIndex(0);
        eBarsCombo.setEnabled(false);

    	//	Create Sensitivity Coefficients (SC) panel
        JPanel SC = new JPanel();
    	SC.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Sensitivity Coefficients"),
    			BorderFactory.createEmptyBorder(10,10,10,10)));
    	
    	//	Create SC on/off subpanel
    	JPanel scOnOff = new JPanel();
    	scOnOff.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	scOnOff.add(new JLabel("Display sensitivity coefficients"));
    	scOnOff.add(sensCoeffCombo = new JComboBox(yesnoOptions));
    	sensCoeffCombo.setSelectedIndex(0);
    	sensCoeffCombo.setActionCommand("scOnOff");
    	sensCoeffCombo.addActionListener(this);
    	sensCoeffCombo.setRenderer(centerJComboBoxRenderer);
    	sensCoeffCombo.setEnabled(false);
    	
    	//	Create SC list subpanel
    	JPanel scList = new JPanel();
    	scList.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	scList.add(new JLabel("Species Name"));
    	scList.add(scName = new JTextField());
    	scName.setEnabled(false);
    	scName.setHorizontalAlignment(JTextField.CENTER);
    	scName.setPreferredSize(new Dimension(60,25));
    	
    	//	Create sensitivity table, tablemodel, scroll panel
    	tableSens = new JTable(tmodelSens = new MyTableModelSensitivity());
    	tableSens.setPreferredScrollableViewportSize(new Dimension(30,50));
        tableSens.getColumnModel().getColumn(0).setPreferredWidth(30);
        tableSens.setEnabled(false);
        
    	JScrollPane scrollSens = new JScrollPane(tableSens);
    	scrollSens.setBorder(BorderFactory.createLoweredBevelBorder());
    	//			* Set the alignment for all cells in the table to CENTER
        for (int i=0; i<tableSens.getColumnCount(); i++) {
        	TableColumn column = tableSens.getColumnModel().getColumn(i);
        	column.setCellRenderer(centerRenderer);
        }

    	//	Create boxes to display panels
    	Box sens1 = Box.createHorizontalBox();
    	Box sens2 = Box.createHorizontalBox();
    	
    	//	Fill boxes
    	sens1.add(AddSens = new JButton("Add"));
    	sens1.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
    	AddButtonListener addListenerSens = new AddButtonListener();
    	AddSens.addActionListener(addListenerSens);
    	AddSens.setActionCommand("AddSens");
    	AddSens.setToolTipText("Press to submit molecule");
    	AddSens.setEnabled(false);
    	
    	sens2.add(DeleteSens = new JButton("Remove"));
    	DeleteButtonListener deleteListenerSens = new DeleteButtonListener();
    	DeleteSens.addActionListener(deleteListenerSens);
    	DeleteSens.setActionCommand("DeleteSens");
    	DeleteSens.setToolTipText("Press to delete molecule");
    	DeleteSens.setEnabled(false);
    	
    	JPanel scButtons = new JPanel();
    	scButtons.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	scButtons.add(sens1);
    	scButtons.add(sens2);
    	
    	Box scBox = Box.createVerticalBox();
    	scBox.add(scOnOff);
    	scBox.add(scList);
    	scBox.add(scButtons);
    	scBox.add(scrollSens);
    	SC.add(scBox);
    	
    	//	Total tab
    	Box Tab3Total = Box.createVerticalBox();
    	Tab3Total.add(EB);
    	Tab3Total.add(SC);
    	
    	JPanel SA = new JPanel();
    	SA.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Sensitivity Analysis"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	SA.add(Tab3Total);

    	JScrollPane scrolltab3 = new JScrollPane(SA);
    	scrolltab3.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
        return scrolltab3;
    }
      
    class DeleteButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			if ("DeletePRL".equals(e.getActionCommand())) {
				int highlightRow = tablePRL.getSelectedRow();
				if (highlightRow == -1) {
					System.out.println("Please select a row to delete");
				} else
					tmodelPRL.deleteRow(highlightRow);
			} else if ("DeleteInput".equals(e.getActionCommand())) {
				int highlightRow = tableInput.getSelectedRow();
				if (highlightRow == -1) {
					System.out.println("Please select a row to delete");
				} else {
					tmodelInput.deleteRow(highlightRow);
					limitingReactant();
				}
			} else if ("DeleteSens".equals(e.getActionCommand())) {
				int highlightRow = tableSens.getSelectedRow();
				if (highlightRow == -1) {
					System.out.println("Please select a row to delete");
				} else
					tmodelSens.deleteRow(highlightRow);
			}
		}
	}
    
    class AddButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			// IF the user adds a primary reaction library ...
			if ("AddPRL".equals(e.getActionCommand())) {
				// Extract the information
				int prlInput0 = tmodelPRL.nextEmptyRow+1;
				String prlInput1 = libName.getText();
				String prlInput2 = prlPath.getText();
				// Check that all information is present
				if (prlInput1.equals("") || prlInput2.equals("")) {
					System.out.println("Please input the Name and Location of a Primary Reaction Library");
				} else {
					PRLVector prlEntry = new PRLVector(prlInput0, prlInput1, prlInput2);
					tmodelPRL.updatePRL(prlEntry);
					libName.setText("");
    				prlPath.setText("");
				}
			// ELSE IF the user adds a species ...
			} else if ("AddInput".equals(e.getActionCommand())) {
				// Extract the information
				int icInput0 = tmodelInput.nextEmptyRow+1;
    			String icInput1 = NameInput.getText();
    			String icInput2 = ConcInput.getText();
    			String icInput3 = (String)UnitsInput.getSelectedItem();
    			String icInput4 = (String)ReactiveInput.getSelectedItem();
    			String icInput5 = "";
    			//	Ask user for adjacency list if species is marked "Reactive"
    			if (!icInput4.equals("Inert Gas")) {
	    			if (!speciesAdjList.getText().equals("")) {
	    				icInput5 = speciesAdjList.getText();
	    			} else if (!InChIInput.getText().equals("")) {
	    				String inchi = InChIInput.getText();
	    				icInput5 = Species.inchi2chemGraph(inchi);
	    			} else {
	    			    File molFile = null;
	    			    molFile = askUserForInput("Select corresponding .mol file", false);
	    			    if (molFile != null) icInput5 = Species.mol2chemGraph(molFile.getAbsolutePath());
	    			}
    			}
    			// If a species is marked "Unreactive," ask user for location of IncludeSpecies.txt file
    			if (icInput4.equals("Unreactive")) {
    				// If the path is not known
    				if (isPath.getText().equals("")) {
        				File isFile = null;
        				isFile = askUserForInput("Select IncludeSpecies.txt file", false);
        			    if (isFile != null) {
        			    	isPath.setText(isFile.getPath());
        			    	isPath.setEnabled(true);
        			    }
    				}
    			}
    			// Check that all information is present
				if (icInput1.equals("") || icInput2.equals("") || (icInput5.equals("") & !icInput4.equals("Inert Gas"))) {
        			System.out.println("Please input the molecule's Name, Concentration, and Adjacency List");
        		} else {
        			ICVector icEntry = new ICVector(icInput0, icInput1, icInput2, icInput3, icInput4, icInput5);
        			tmodelInput.updateIC(icEntry);
        			NameInput.setText("");
        			InChIInput.setText("");
        			ConcInput.setText("");
        			UnitsInput.setSelectedIndex(0);
        			ReactiveInput.setSelectedIndex(0);
        			speciesAdjList.setText("");
        			if (!icInput4.equals("Inert Gas")) limitingReactant();
        		}
			// ELSE IF the user adds a species to sensitivity analysis
			} else if ("AddSens".equals(e.getActionCommand())) {
				// Extract the information
				int sensInput0 = tmodelSens.nextEmptyRow+1;
    			String sensInput1 = scName.getText();
    			// Check if all information is present
    			if (sensInput1.equals("")) {
        			System.out.println("Please input the name of a Molecule");
        		} else {
        			SensVector sensEntry = new SensVector(sensInput0, sensInput1);
        			tmodelSens.updateSens(sensEntry);
        			scName.setText("");
        		}
			}
		}
	}
    
    class ChangeButtonListener implements ActionListener {
    	public void actionPerformed(ActionEvent e) {
    		// IF the user selects to switch the IncludeSpecies.txt path
			if ("isPath".equals(e.getActionCommand())) {
				File path = null;
				path = askUserForInput("Select IncludeSpecies.txt file", false);
				if (path != null) isPath.setText(path.getPath());
			// ELSE IF user selects to switch the database path
			}  else if ("databasePath".equals(e.getActionCommand())) {
				File path = null;
				path = askUserForInput("Select database folder", true);
				if (path != null) databasePath.setText(path.getPath());
			// If user selects to switch Primary Thermo Library path
			} else if ("ptlPath".equals(e.getActionCommand())) {
				File path = null;
				path = askUserForInput("Select PrimaryThermoLibrary folder", true);
				if (path != null) ptlPath.setText(path.getPath());
			// If user selects to switch Primary Reaction Library path
			} else if ("prlPath".equals(e.getActionCommand())) {
				File path = null;
				path = askUserForInput("Select PrimaryReactionLibrary folder", true);
				if (path != null) prlPath.setText(path.getPath());
			}
    	}
    }
    
    public File askUserForInput (String popupwindowTitle, boolean folder) {
    	String workingDirectory = System.getProperty("RMG.workingDirectory");
	    JFileChooser chooser = new JFileChooser(workingDirectory);
	    chooser.setApproveButtonText("Select");
	    chooser.setDialogTitle(popupwindowTitle);
	    chooser.rescanCurrentDirectory();
	    if (folder) chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
	    int result = chooser.showDialog(GUI.this, null);
	    File file = null;
	    if (result == chooser.APPROVE_OPTION) {
	    	file = chooser.getSelectedFile();
	    } else {
	    	chooser.cancelSelection();
	    }
	    return file;
    }
    
    /*
     * This function updates the Limiting Reactant combobox on the Termination
     * tab of the GUI - for a termination criteria of expected conversion, only
     * an input to the simulation is acceptable
     */
    public void limitingReactant() {
    	int permanentNumber = SpeciesConvName.getItemCount();
    	for (int i=1; i<permanentNumber; i++) {
    		SpeciesConvName.removeItemAt(1);
    	}
    	int nonInertCount = 0;
    	for (int j=0; j<tmodelInput.nextEmptyRow; j++) {
    		if (!tmodelInput.getValueAt(j,3).equals("Inert Gas")) {
    			++nonInertCount;
    			SpeciesConvName.insertItemAt(tmodelInput.getValueAt(j,0),nonInertCount);
    		}
    	}
    }
        
    public void createConditionFile() {
    	//	Let's create the condition.txt file!
    	String conditionFile = "";
    	
    	//	Add user supplied comments
    	conditionFile += "// User supplied comments:\r";
    	String[] header4Condition = headerComments.getText().split("[\n]");
    	for (int i=0; i<header4Condition.length; i++) {
    		conditionFile += "// " + header4Condition[i] + "\r";
    	}
    	//	Add date and time to the file
    	String DATE_FORMAT_NOW = "yyyy-MM-dd HH:mm:ss";
        SimpleDateFormat sdf = new SimpleDateFormat(DATE_FORMAT_NOW);
        Calendar cal = Calendar.getInstance();
        conditionFile += "//\r// File created at: " + sdf.format(cal.getTime()) + "\r\r";
        
        //	Add the name of the database and primary thermo library
        String[] tempString = databasePath.getText().split("[\\\\,/]");
        conditionFile += "Database: " + tempString[tempString.length-1] + "\r";
        
        tempString = ptlPath.getText().split("[\\\\,/]");
        conditionFile += "PrimaryThermoLibrary: " + tempString[tempString.length-1] + "\r";
        
        //	Add the temperature model, list of temperatures, and units
        conditionFile += "TemperatureModel: " + tempModelCombo.getSelectedItem() +
        	" (" + tempConstantUnit.getSelectedItem() + ")";
        String[] tempLines = tempConstant.getText().split("[\n]");
        StringTokenizer st = null;
        for (int i=0; i<tempLines.length; i++) {
        	st = new StringTokenizer(tempLines[i]);
        	while (st.hasMoreTokens()) {
        		String line = st.nextToken();
        		if (line.endsWith(",")) line = line.substring(0,line.length()-1);
        		conditionFile += " " + line;
        	}
        }
        conditionFile += "\r";
        
        //	Add the pressure model, list of pressures, and units
        conditionFile += "PressureModel: " + pressModelCombo.getSelectedItem() +
        	" (" + pressConstantUnit.getSelectedItem() + ")";
        String[] pressLines = pressConstant.getText().split("[\n]");
        for (int i=0; i<pressLines.length; i++) {
        	st = new StringTokenizer(pressLines[i]);
        	while (st.hasMoreTokens()) {
        		String line = st.nextToken();
        		if (line.endsWith(",")) line = line.substring(0,line.length()-1);
        		conditionFile += " " + line;
        	}
        }
        conditionFile += "\r\r";
        
        //	Add the equation of state (if equal to liquid)
        if (eosCombo.getSelectedItem().equals("Liquid")) {
        	conditionFile += "EquationOfState: Liquid\r";
        }
        
        //	Add the spectroscopic data estimator
        conditionFile += "SpectroscopicDataEstimator: ";
        if (sdeCombo.getSelectedItem().equals("Frequency Groups")) {
        	conditionFile += "frequencygroups\r";
        } else if (sdeCombo.getSelectedItem().equals("therfit")) {
        	conditionFile += "therfit\r";
        }
        
        //	Add whether to generate InChIs or not
        conditionFile += "InChIGeneration: ";
        if (inchiCombo.getSelectedItem().equals("Yes")) {
        	conditionFile += "on\r\r";
        } else if (inchiCombo.getSelectedItem().equals("No")) {
        	conditionFile += "off\r\r";
        }
        
        //	Determine if each species is reactive, unreactive, or inert
        int nonInertCount = 0;
        boolean unreactiveStatus = false;
        String reactiveSpecies = "";
        String inertSpecies = "";
        for (int i=0; i<tableInput.getRowCount(); i++) {
        	if (tableInput.getValueAt(i,3).equals("Inert Gas")) {
        		inertSpecies += tableInput.getValueAt(i,0) + " " +
        			tableInput.getValueAt(i,1) + " (" +
        			tableInput.getValueAt(i,2) + ")\r";
        	} else if (tableInput.getValueAt(i,3).equals("Reactive")) {
        		++nonInertCount;
        		reactiveSpecies += "(" + nonInertCount + ") " +
        		tableInput.getValueAt(i,0) + " " +
    			tableInput.getValueAt(i,1) + " (" +
    			tableInput.getValueAt(i,2) + ")\r" +
    			tableInput.getValueAt(i,4) + "\r";
        	} else if (tableInput.getValueAt(i,3).equals("Unreactive")) {
        		unreactiveStatus = true;
        		++nonInertCount;
        		reactiveSpecies += "(" + nonInertCount + ") " +
        		tableInput.getValueAt(i,0) + " " +
    			tableInput.getValueAt(i,1) + " (" +
    			tableInput.getValueAt(i,2) + ") unreactive\r" +
    			tableInput.getValueAt(i,4) + "\r";
        	}
        }
        //	Add the non-inert gas species
        conditionFile += "InitialStatus:\r\r" + reactiveSpecies + "\rEND\r\r";
        //	Add the inert gas species
        conditionFile += "InertGas:\r" + inertSpecies + "END\r\r";
        
        //	Add the model enlarger
        conditionFile += "ReactionModelEnlarger: " + modelCombo.getSelectedItem() + "\r";
        //	IF not PDepModelEnlarger, check if IncludeSpecies.txt file is given
        if (modelCombo.getSelectedItem().equals("RateBasedModelEnlarger")) {
	        if (unreactiveStatus) {
	        	conditionFile += "IncludeSpecies: " + isPath.getText() + "\r\r";
	        }
        } else if (modelCombo.getSelectedItem().equals("RateBasedPDepModelEnlarger")) {
        	conditionFile += "PDepKineticsEstimator: ";
        	if (estimatorCombo.getSelectedItem().equals("Modified Strong Collision")) {
        		conditionFile += "modifiedstrongcollision";
        	} else if (estimatorCombo.getSelectedItem().equals("Reservoir State")) {
        		conditionFile += "reservoirstate";
        	}
        	conditionFile += "\r\r";
        }
        
        //	Add the finish controller
        conditionFile += "FinishController: " + finishcCombo.getSelectedItem() + "\r" +
        	"(1) Goal ";
        //	If goal is conversion
        boolean conversionGoal = true;
        if (controllerCombo.getSelectedItem().equals("Conversion")) {
        	conditionFile += controllerCombo.getSelectedItem() + " " +
        		SpeciesConvName.getSelectedItem() + " " + conversion.getText() + "\r";
        } else if (controllerCombo.getSelectedItem().equals("ReactionTime")) {
        	conversionGoal = false;
        	conditionFile += controllerCombo.getSelectedItem() + " " +
        		time.getText() + " (" + timeCombo.getSelectedItem() + ")\r";
        }
        conditionFile += "(2) Error Tolerance: " + eTolerance.getText() + "\r\r";
        
        //	Add the dynamic simulator
        conditionFile += "DynamicSimulator: " + simulatorCombo.getSelectedItem() + "\r";
        if (conversionGoal) {
        	conditionFile += "Conversions:";
            String[] convLines = multiConv.getText().split("[\n]");
            for (int i=0; i<convLines.length; i++) {
            	st = new StringTokenizer(convLines[i]);
            	while (st.hasMoreTokens()) {
            		String line = st.nextToken();
            		if (line.endsWith(",")) line = line.substring(0,line.length()-1);
            		conditionFile += " " + line;
            	}
            }
        } else {
        	conditionFile += "TimeStep:";
        	String[] tsLines = multiTS.getText().split("[\n]");
            for (int i=0; i<tsLines.length; i++) {
            	st = new StringTokenizer(tsLines[i]);
            	while (st.hasMoreTokens()) {
            		String line = st.nextToken();
            		if (line.endsWith(",")) line = line.substring(0,line.length()-1);
            		conditionFile += " " + line;
            	}
            }
        }
    	conditionFile += "\rAtol: " + aTolerance.getText() + "\r" +
			"Rtol: " + rTolerance.getText() + "\r\r";
    	
    	//	Add error bars & sensitivity analysis information (if DASPK is selected)
    	if (simulatorCombo.getSelectedItem().equals("DASPK")) {
    		conditionFile += "Error bars: ";
    		if (eBarsCombo.getSelectedItem().equals("Yes")) {
    			conditionFile += "on\r";
    		} else if (eBarsCombo.getSelectedItem().equals("No")) {
    			conditionFile += "off\r";
    		}
    		conditionFile += "Display sensitivity coefficients: ";
    		if (sensCoeffCombo.getSelectedItem().equals("Yes")) {
    			conditionFile += "on\r";
    		} else if (sensCoeffCombo.getSelectedItem().equals("No")) {
    			conditionFile += "off\r";
    		}
    		conditionFile += "Display sensitivity information for:\r";
    		for (int j=0; j<tableSens.getRowCount(); j++) {
    			conditionFile += tableSens.getValueAt(j,0) + "\r";
    		}
    		conditionFile += "END\r\r";
    	}
    	
        //	Add the primary reaction library
    	String rmgEnvVar = System.getProperty("RMG.workingDirectory");
    	conditionFile += "PrimaryReactionLibrary: ";
    	if (libraryCombo.getSelectedItem().equals("Yes")) {
    		conditionFile += "on\r";
    		for (int k=0; k<tablePRL.getRowCount(); k++) {
    			conditionFile += "Name: " + tablePRL.getValueAt(k,0) + "\r" + "Location: ";
    			String prlDir = (String)tablePRL.getValueAt(k,1);
    	        if (prlDir.startsWith(rmgEnvVar)) {
    	        	int startIndex = rmgEnvVar.length() + 11;
    	        	conditionFile += prlDir.substring(startIndex) + "\r";
    	        } else {
    	        	conditionFile += prlDir + "\r";
    	        }
    		}
    		conditionFile += "END";
    	} else if (libraryCombo.getSelectedItem().equals("No")) {
    		conditionFile += "off\rEND";
    	}
    	
		File conditionPath = null;
		FileWriter fw = null;
		// Write the cTable to species.mol file
        try {
        	conditionPath = askUserForInput("Save file", false);
        	String[] rmgArgs = new String[1];
        	if (conditionPath != null) {
        		fw = new FileWriter(conditionPath);
        		fw.write(conditionFile);
            	fw.close();
            	rmgArgs[0] = conditionPath.getAbsolutePath();
            	RMG.main(rmgArgs);
        	}
        } catch (IOException e) {
        	String err = "Error writing condition file: ";
        	err += e.toString();
        	System.out.println(err);
        }
	}
    
    public void openConditionFile(File filePath) {
		//	Read in the .txt file
        FileReader in = null;
		try {
			in = new FileReader(filePath);
		} catch (FileNotFoundException e) {
			String err = "Error reading file: ";
			err += e.toString();
			System.out.println(err);
		}
		BufferedReader reader = new BufferedReader(in);
		String conditionFile = "";
        String line = ChemParser.readMeaningfulLine(reader);
        read: while (line != null) {
        	conditionFile += line + "\r";
        	line = ChemParser.readMeaningfulLine(reader);
        }
        String[] indivLines = conditionFile.split("[\r]");
        
        // Populate the GUI with the appropriate information
        int lineCounter = 0;
        StringTokenizer st = null;
        String tempString = null;
        
        //	Path of Database
        databasePath.setText(indivLines[lineCounter].substring(9));
        ++lineCounter;
        
        //	Path of PrimaryThermoLibrary
        ptlPath.setText(indivLines[lineCounter].substring(21));
        ++lineCounter;
        
        //	Temperature model
        st = new StringTokenizer(indivLines[lineCounter]);
        tempString = st.nextToken();	// Skip over "TemperatureModel:"
        tempModelCombo.setSelectedItem(st.nextToken());
        tempConstantUnit.setSelectedItem(st.nextToken().substring(1,2));
        String temps = st.nextToken();
        while (st.hasMoreTokens()) {
        	temps += " " + st.nextToken();
        }
        tempConstant.setText(temps);
        ++lineCounter;
        
        //	Pressure model
        st = new StringTokenizer(indivLines[lineCounter]);
        tempString = st.nextToken();	// Skip over "PressureModel:"
        pressModelCombo.setSelectedItem(st.nextToken());
        String pressUnit = st.nextToken();
        pressConstantUnit.setSelectedItem(pressUnit.substring(1,pressUnit.length()-1));
        String presses = st.nextToken();
        while (st.hasMoreTokens()) {
        	presses += " " + st.nextToken();
        }
        pressConstant.setText(presses);
        ++lineCounter;
        
        //	Equation of state (if present)
        st = new StringTokenizer(indivLines[lineCounter]);
        if (st.nextToken().startsWith("EquationOfState")) {
        	eosCombo.setSelectedItem(st.nextToken());
        	++lineCounter;
        }
        
        //	Spectroscopic data estimator
        st = new StringTokenizer(indivLines[lineCounter]);
        tempString = st.nextToken();	// Skip over "SpectroscopicDataEstimator:"
        String sde = st.nextToken();
        if (sde.equals("frequencygroups")) {
        	sdeCombo.setSelectedItem("Frequency Groups");
        } else if (sde.equals("therfit")) {
        	sdeCombo.setSelectedItem("therfit");
        }
        ++lineCounter;
        
        //	InChI generation
        st = new StringTokenizer(indivLines[lineCounter]);
        tempString = st.nextToken();	// Skip over "InChIGeneration:"
        String inchi = st.nextToken();
        if (inchi.equals("on")) {
        	inchiCombo.setSelectedItem("Yes");
        } else if (inchi.equals("off")) {
        	inchiCombo.setSelectedItem("No");
        }
        ++lineCounter;
        
        //	Non inert gas species
        ++lineCounter;	// Skip over "InitialStatus:"
        tmodelInput.clear();	// Clear the existing table
        String adjList = "";
        int numSpecies = 0;
        String name = "";
        String conc = "";
        String unit = "";
        String reactivity = "";
        while (!indivLines[lineCounter].startsWith("END")) {
        	if (indivLines[lineCounter].startsWith("(")) {
        		++numSpecies;
        		st = new StringTokenizer(indivLines[lineCounter]);
        		tempString = st.nextToken();	// Skip over (#)
        		name = st.nextToken();
        		conc = st.nextToken();
        		unit = st.nextToken();
        		reactivity = "Reactive";
        		if (st.hasMoreTokens()) {
        			reactivity = "Unreactive";
        		}
        	} else {
        		adjList += indivLines[lineCounter] + "\r";
        	}
        	
        	++lineCounter;
        	if (indivLines[lineCounter].startsWith("(") || indivLines[lineCounter].startsWith("END")) {
        		ICVector icEntry = new ICVector(numSpecies-1,name,conc,unit.substring(1,unit.length()-1),reactivity,adjList);
    			tmodelInput.updateIC(icEntry);
        		adjList = "";
        	}
        }
        limitingReactant();
        ++lineCounter;	// Skip over "END"
        ++lineCounter;	// Skip over "InertGas:"
        
        //	Inert gases
        while (!indivLines[lineCounter].startsWith("END")) {
    		++numSpecies;
    		st = new StringTokenizer(indivLines[lineCounter]);
    		name = st.nextToken();
    		conc = st.nextToken();
    		unit = st.nextToken();
    		reactivity = "Inert Gas";
    		ICVector icEntry = new ICVector(numSpecies-1,name,conc,unit.substring(1,unit.length()-1),reactivity,adjList);
			tmodelInput.updateIC(icEntry);
        	++lineCounter;
        }
        ++lineCounter;	// Skip over "END"
        
        //	Model enlarger
        st = new StringTokenizer(indivLines[lineCounter]);
        tempString = st.nextToken();	// Skip over "ModelEnlarger:"
        String me = st.nextToken();
        modelCombo.setSelectedItem(me);
        ++lineCounter;
        
        //	IncludeSpecies.txt file (if present)
        if (me.equals("RateBasedModelEnlarger")) {
        	if (indivLines[lineCounter].startsWith("IncludeSpecies")) {
        		isPath.setText(indivLines[lineCounter].substring(15));
        		++lineCounter;
        	}
        }
        //	P-dep kinetics estimator
        else if (me.equals("RateBasedPDepModelEnlarger")) {
        	st = new StringTokenizer(indivLines[lineCounter]);
        	tempString = st.nextToken();	// Skip over "PDepKineticsEstimator:"
        	String ke = st.nextToken();
        	if (ke.equals("reservoirstate")) {
        		estimatorCombo.setSelectedItem("Reservoir State");
        	} else if (ke.equals("modifiedstrongcollision")) {
        		estimatorCombo.setSelectedItem("Modified Strong Collision");
        	}
        	++lineCounter;
        }
        
        //	Finish controller
        st = new StringTokenizer(indivLines[lineCounter]);
        tempString = st.nextToken();	// Skip over "FinishController:"
        finishcCombo.setSelectedItem(st.nextToken());
        ++lineCounter;
        
        //	Goal
        st = new StringTokenizer(indivLines[lineCounter]);
        tempString = st.nextToken();	// Skip over "(1)"
        tempString = st.nextToken();	// Skip over "Goal"
        String goal = st.nextToken();
        if (goal.startsWith("Conversion")) {
        	controllerCombo.setSelectedItem("Conversion");
        	SpeciesConvName.setSelectedItem(st.nextToken());
        	conversion.setText(st.nextToken());
        } else if (goal.startsWith("Reaction")) {
        	controllerCombo.setSelectedItem("ReactionTime");
        	time.setText(st.nextToken());
        	String units = st.nextToken();
        	timeCombo.setSelectedItem(units.substring(1,units.length()-1));
        }
        ++lineCounter;
        
        //	Error tolerance
        st = new StringTokenizer(indivLines[lineCounter]);
        tempString = st.nextToken();	// Skip over "(2)"
        tempString = st.nextToken();	// Skip over "Error"
        tempString = st.nextToken();	// Skip over "Tolerance:"
        eTolerance.setText(st.nextToken());
        ++lineCounter;
        
        //	Dynamic simulator
        st = new StringTokenizer(indivLines[lineCounter]);
        tempString = st.nextToken();	// Skip over "DynamicSimulator:"
        simulatorCombo.setSelectedItem(st.nextToken());
        ++lineCounter;
        
        //	Intermediate integration step
        st = new StringTokenizer(indivLines[lineCounter]);
        String stepName = st.nextToken();
        String steps = "";
        while (st.hasMoreTokens()) {
        	steps += st.nextToken() + " ";
        }
        if (stepName.startsWith("Conversions")) {
        	multiConv.setText(steps);
        } else if (stepName.startsWith("TimeStep")) {
        	multiTS.setText(steps);
        }
        ++lineCounter;
        
        //	Absolute and Relative tolerance
        st = new StringTokenizer(indivLines[lineCounter]);
        tempString = st.nextToken();	// Skip over "Atol:"
        aTolerance.setText(st.nextToken());
        ++lineCounter;
        
        st = new StringTokenizer(indivLines[lineCounter]);
        tempString = st.nextToken();	// Skip over "Rtol:"
        rTolerance.setText(st.nextToken());
        ++lineCounter;

        //	Error bars and Sensitivity analysis (if present)
        if (indivLines[lineCounter].startsWith("Error")) {
        	st = new StringTokenizer(indivLines[lineCounter]);
        	tempString = st.nextToken();	// Skip over "Error"
        	tempString = st.nextToken();	// Skip over "bars:"
        	String onoff = st.nextToken();
        	if (onoff.equals("on")) {
        		eBarsCombo.setSelectedItem("Yes");
        	} else if (onoff.equals("off")) {
        		eBarsCombo.setSelectedItem("No");
        	}
        	++lineCounter;
        	
            st = new StringTokenizer(indivLines[lineCounter]);
            tempString = st.nextToken();	// Skip over "Display"
            tempString = st.nextToken();	// Skip over "sensitivity"
            tempString = st.nextToken();	// Skip over "coefficients:"
            onoff = st.nextToken();
            if (onoff.equals("on")) {
        		sensCoeffCombo.setSelectedItem("Yes");
        	} else if (onoff.equals("off")) {
        		sensCoeffCombo.setSelectedItem("No");
        	}
        	++lineCounter;        	
        }
        
        //	Sensitivity coefficients (if present)
        tmodelSens.clear();	// Clear the existing table
        if (indivLines[lineCounter].startsWith("Display sensitivity information for")) {
        	++lineCounter;	// Skip over this entire line
        	
        	numSpecies = 0;
        	while (!indivLines[lineCounter].startsWith("END")) {
        		++numSpecies;
        		SensVector sensEntry = new SensVector(numSpecies-1,indivLines[lineCounter]);
    			tmodelSens.updateSens(sensEntry);
    			++lineCounter;
        	}
        	++lineCounter;	// Skip over "END"
        }
        
        //	Primary reaction libraries
        tmodelPRL.clear();	// Clear the existing table
        st = new StringTokenizer(indivLines[lineCounter]);
        tempString = st.nextToken();	// Skip "PrimaryReactionLibrary:"
        String onoff = st.nextToken();
        if (onoff.equals("on")) {
    		libraryCombo.setSelectedItem("Yes");
    		++lineCounter;
    		
    		String libName = "";
    		String libLoc = "";
    		int numLib = 0;
    		while (!indivLines[lineCounter].startsWith("END")) {
    			++numLib;
    			libName = indivLines[lineCounter].substring(5);
    			++lineCounter;
    			libLoc = indivLines[lineCounter].substring(9);
    			++lineCounter;
    			
    			PRLVector prlEntry = new PRLVector(numLib-1,libName,libLoc);
				tmodelPRL.updatePRL(prlEntry);
    		}
    	} else if (onoff.equals("off")) {
    		libraryCombo.setSelectedItem("No");
    	}
    }
    
    public void enableComponents (JComponent[] allComponents) {
    	for (int i=0; i<allComponents.length; i++)
    		allComponents[i].setEnabled(true);
    }
    
    public void disableComponents (JComponent[] allComponents) {
    	for (int i=0; i<allComponents.length; i++)
    		allComponents[i].setEnabled(false);
    }
    
    public void actionPerformed(ActionEvent event) {    	
    	//	- If Termination sequence is specified, highlight appropriate boxes
    	if ("GoalofReaction".equals(event.getActionCommand())) {
    		JComponent[] timeComps = {timeCombo, time, multiTS};
    		JComponent[] convComps = {SpeciesConvName, conversion, multiConv};
    		if (controllerCombo.getSelectedItem().equals("ReactionTime")) {
    			enableComponents(timeComps);
    			disableComponents(convComps);
    		} else if (controllerCombo.getSelectedItem().equals("Conversion")) {
    			enableComponents(convComps);
    			disableComponents(timeComps);
    		}
    	}
    	//	- If Sensitivity analysis is turned on/off, highlight the appropriate boxes
    	else if ("scOnOff".equals(event.getActionCommand())) {
    		JComponent[] scComps = {AddSens, DeleteSens, tableSens, scName};
    		if (sensCoeffCombo.getSelectedItem().equals("Yes")) {
    			enableComponents(scComps);
    		} else if (sensCoeffCombo.getSelectedItem().equals("No")){
    			disableComponents(scComps);
    		}
    	}
    	//	If the user enables/disables a primary reaction library
    	else if ("prlOnOff".equals(event.getActionCommand())) {
    		JComponent[] prlComps = {libName, prlPath, prlButton, AddPRL, DeletePRL, tablePRL};
    		if (libraryCombo.getSelectedItem().equals("No")) {
    			disableComponents(prlComps);
    		} else if (libraryCombo.getSelectedItem().equals("Yes")) {
    			enableComponents(prlComps);
    		}
    	}
    	//	If the user enables/disables the p-dep model enlarger
    	else if ("PDepME".equals(event.getActionCommand())) {
    		JComponent[] meComps = {estimatorCombo};
    		if (modelCombo.getSelectedItem().equals("RateBasedModelEnlarger")) {
    			disableComponents(meComps);
    		} else if (modelCombo.getSelectedItem().equals("RateBasedPDepModelEnlarger")) {
    			enableComponents(meComps);
    		}
    	}
    	//	If the user enables/disables DASPK
    	else if ("saOnOff".equals(event.getActionCommand())) {
    		JComponent[] saComps =  {eBarsCombo, sensCoeffCombo};
    		JComponent[] otherComps = {AddSens, DeleteSens, tableSens, scName};
    		if (simulatorCombo.getSelectedItem().equals("DASPK")) {
    			enableComponents(saComps);
    			if (sensCoeffCombo.getSelectedItem().equals("Yes")) enableComponents(otherComps);
    			else if (sensCoeffCombo.getSelectedItem().equals("No")) disableComponents(otherComps);
    		} else if (simulatorCombo.getSelectedItem().equals("DASSL")) {
    			disableComponents(saComps);
    			disableComponents(otherComps);
    		}
    			
    	}
    	//	- Create the condition.txt file
    	else if ("Condition".equals(event.getActionCommand())) {
    		createConditionFile();
    	}
    	
	}
    
	//	Main panel
    private static GUI theApp;
	
	//	Tab0: Initialization
	String[] timeUnits = {"sec", "min", "hr", "day"};
	//	Tab1: Termination Sequence
	String[] species_blank = {" "};
	String[] goalOptions = {"Conversion", "ReactionTime"};
	String[] fcOptions = {"RateBasedFinishController", "RateBasedPDepFinishController"};
	//	Tab2: Initial Condition
	String[] tpModel = {"Constant"};
	String[] concUnits = {"mol/cm3", "mol/m3", "mol/liter"};
	String[] tempUnits = {"K", "C", "F"};
	String[] pressUnits = {"atm", "bar", "Pa", "torr"};
	String[] reactive = {"Reactive", "Unreactive", "Inert Gas"};
	//	Tab3: Error Analysis
	String[] yesnoOptions = {"No", "Yes"};
	//	Tab4: Dynamic Simulator
	String[] meOptions = {"RateBasedModelEnlarger", "RateBasedPDepModelEnlarger"};
	String[] keOptions = {"Reservoir State", "Modified Strong Collision"};
	String[] simOptions = {"DASSL", "DASPK"}; //"CHEMKIN"
	String[] sdeOptions = {"Frequency Groups", "therfit"};
	String[] eosOptions = {"Ideal Gas", "Liquid"};
	
	JPanel
	//	Main Panel
	mainPanel;
	
    JComboBox
    //	Tab0: Initialization
    simulatorCombo, timeStepCombo, libraryCombo,
    //	Tab1: Termination Sequence
    SpeciesConvName, controllerCombo, timeCombo, finishcCombo,
    //	Tab2: Initial Condition
    UnitsInput, ReactiveInput, tempModelCombo, tempConstantUnit,
    	pressModelCombo, pressConstantUnit,
    //	Tab3: Error Analysis
    SpeciesSensName, eBarsCombo, sensCoeffCombo,
    //	Tab4: Dynamic Simulator
    modelCombo, estimatorCombo, sdeCombo, eosCombo, inchiCombo;
    
    JTextField
    //	Tab0: Initialization
    nameFiller, locationFiller, databasePath, prlPath, ptlPath, libName,
    	timeStep, aTolerance, rTolerance,
    //	Tab1: Termination Sequence
    conversion, time, eTolerance,
    //	Tab2: Initial Condition
    NameInput, InChIInput, isPath, ConcInput,
    //	Tab3: Error Analysis
    scName;
    
    JTextArea
    //	Main Panel
    headerComments,
    //	Tab2: Initial Condition
    speciesAdjList, tempConstant, pressConstant,
    //	Tab4: Dynamic Simulator
    multiTS, multiConv;
    
    JButton 
    //	Main panel
    Submit, 
    //	Tab0: Initialization
    AddPRL, DeletePRL, databaseButton, prlButton, ptlButton,
    //  Tab1: Termination
    //	Tab2: Initial Condition
    AddInput, DeleteInput, isButton,
    //	Tab3: Error Analysis
    AddSens, DeleteSens;
    
    JTable
    //	Tab0: Initialization
    tablePRL,
    //  Tab1: Termination
    //	Tab2: Initial Condition
    tableInput,
    //	Tab3: Error Analysis
    tableSens;
    
    //	Tab0: Initialization
    MyTableModelPRL tmodelPRL;
    //  Tab1: Termination
    //	Tab2: Initial Condition
    MyTableModelInput tmodelInput;
    //	Tab3: Error Analysis
    MyTableModelSensitivity tmodelSens;
}
