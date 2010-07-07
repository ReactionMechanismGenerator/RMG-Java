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


/* 
 * GUI.java creates a Graphical User Interface (GUI) for Prof.
 * Green's Group's Reaction Mechanism Generator (RMG).  More information
 * on RMG may be found at: http://www.sourceforge.net/projects/rmg
 * 
 * Michael Harper Jr.
 * MIT Chemical Engineering
 * Green Group
 * 
 * Last updated: March 3, 2009
 *  - Modified reactants & inert gas input table (Initial Conditions tab)
 *  	-- If user adds "Inert Gas", the entry must be "N2", "Ar", or "Ne"; otherwise, an error message is displayed to the user (but the GUI's entries remain untouched).
 *  - Modified how RMG-GUI writes the "condition.txt" file
 *  	-- If user does not add all of the necessary information for the condition.txt file, an error message is displayed (the GUI's entries remain untouched) and no condition.txt file will be written.
 *  	-- If user does not add any "Comments" to be placed in the header of the condition.txt file, a warning if displayed
 *  - Modified how RMG-GUI reads the "condition.txt" file
 *  	-- If condition.txt file is missing a field, GUI instructs user which field is missing (and where it should be located in file)
 *  	-- If condition.txt file contains unknown entries (e.g. units of Temperature not equal to K, F, or C), GUI informs user that entry is incorrect and gives user available options
 * Last updated: February 17, 2009
 *  - Adapted GUI to handle new condition.txt file format
 *      -- Different location of SpectroscopicDataEstimator
 *      -- PressureDependence replacing ModelEnlarger + PDepKineticsEstimator
 *      -- No FinishController options
 * 
 * Currently:
 * 	- creates & displays GUI
 *	- creates .txt file related to running RMG
 *		* condition.txt : file containing all inputs (species, tolerances, etc.) 
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
import javax.swing.JOptionPane;
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

import jing.chem.PrimaryThermoLibrary;
import jing.chem.Species;
import jing.chemParser.ChemParser;

public class GUI extends JPanel implements ActionListener {
    public static void main(String[] args) {
    	RMG.globalInitializeSystemProperties();
//    	String workingDir = System.getenv("RMG");
//    	System.setProperty("RMG.workingDirectory", workingDir);
    	theApp = new GUI();
    	theApp.createAndShowGUI();
    }
    	
    private void createAndShowGUI() {
    	GUIWindow frame = new GUIWindow("RMG", this);
    	frame.setContentPane(theApp.mainPanel);
    	Dimension wndSize = frame.getToolkit().getScreenSize();
    	frame.setBounds(wndSize.width*1/8, 0, wndSize.width*3/4, wndSize.height*15/16);
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
    	JComponent tabOptions = createTabOptions();
    	JComponent tabSensitivity = createTabSensitivity();
    	//	Add the individual tabs to the main panel
    	tabbedPanel.addTab("Initial Conditions", null, tabInputs, "Specify Reactants, Inerts, & Temperature/Pressure Model");
    	tabbedPanel.addTab("Species/Reaction Libraries", null, tabInitialization, "Specify thermochemistry/transport properties for individual species & rates for individual reactions");
    	tabbedPanel.addTab("Termination Criteria", null, tabTermination, "Specify Simulation Termination Conditions");
    	tabbedPanel.addTab("Dynamic Simulator", null, tabSolver, "Specify Solver Tolerances");
    	tabbedPanel.addTab("Additional Options", null, tabOptions, "Specify Other Options");
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
    	comments.add(save = new JButton("Save"));
        save.addActionListener(this);
    	save.setActionCommand("saveCondition");
    	save.setToolTipText("Press to save file");
        comments.add(saveAndRun = new JButton("Save and Run"));
        saveAndRun.addActionListener(this);
    	saveAndRun.setActionCommand("runCondition");
    	saveAndRun.setToolTipText("Press to save and run file");
    	
    	submitConditionFile.add(comments);
    	submitConditionFile.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder
    			("Create Initialization File: condition.txt"), BorderFactory.createEmptyBorder(5,5,5,5)));
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
    	JLabel labelName = new JLabel("          Name");
    	JLabel labelInChI = new JLabel("                 InChI");
    	JLabel labelConc = new JLabel("   Concentration");
    	JLabel labelUnit = new JLabel("Units             ");
    	JLabel labelReact = new JLabel("Reactivity                   ");
    	
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
    	tableInput.setPreferredScrollableViewportSize(new Dimension(500,150));
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
    	adj.add(new JLabel("             Adjacency List"));
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
    	JLabel valueTLabel = new JLabel("            List of Temperatures");
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
		JLabel valuePLabel = new JLabel("           List of Pressures");
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
		- The primary kinetic library
			* Name
			* Location
    */
    public JComponent createTabInitialization() {
    	//	Create cellrenderers for JComboBox/JTable - CENTER text
        TableCellRenderer centerTableRenderer = new CenterTRenderer();    	

    	//	Create the "Database" panel
    	JPanel Database = new JPanel();
    	Database.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Main database"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	
    	//	Populate the "Database" panel
    	JLabel databaseLabel = new JLabel("Choose database");
    	Database.add(databaseLabel);
    	databaseLabel.setToolTipText("Default = $RMG/databases/RMG_database");
    	
    	Database.add(databasePath = new JTextField(25));
    	databasePath.setText(System.getProperty("RMG.workingDirectory") + 
			"/databases/RMG_database");
    	databasePath.setEditable(false);
    	
    	Database.add(databaseButton = new JButton("Select"));
    	ChangeButtonListener addListenerDatabase = new ChangeButtonListener();
    	databaseButton.addActionListener(addListenerDatabase);
    	databaseButton.setActionCommand("databasePath");
    	
    	//	Create the Primary Thermodynamic Library (PTL) panel
    	JPanel PTL = new JPanel();
    	PTL.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Primary Thermo Library"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	PTL.setToolTipText("Override RMG's group-additivity algorithm to compute species' thermochemistry");

    	//	Create PTL Name label
    	JPanel ptlName = new JPanel();
    	ptlName.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	
    	//	Populate the label
    	JLabel ptlNameLabel = new JLabel("Name:");
    	ptlName.add(ptlNameLabel);
    	ptlName.add(ptlLibName = new JTextField(20));

        //	Create PTL Location label
        JPanel ptlLoc = new JPanel();
        ptlLoc.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	
        //	Populate the label
        JLabel ptlLocationLabel = new JLabel("Location:");
        ptlLoc.add(ptlLocationLabel);
    	ptlLocationLabel.setToolTipText("Default = " + 
			"RMG/databases/RMG_database/thermo_libraries/primaryThermoLibrary");
    	
    	ptlLoc.add(ptlPath = new JTextField(20));
    	
    	ptlLoc.add(ptlButton = new JButton("Select"));
    	ChangeButtonListener ptlAddListenerLib = new ChangeButtonListener();
    	ptlButton.addActionListener(ptlAddListenerLib);
    	ptlButton.setActionCommand("ptlPath");

        //	Create table and scroll panel to store PTL(s)
        tablePTL = new JTable(tmodelPTL = new MyTableModelPRL());
    	tablePTL.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    	tablePTL.setPreferredScrollableViewportSize(new Dimension(700,50));
        tablePTL.getColumnModel().getColumn(0).setPreferredWidth(100);
        tablePTL.getColumnModel().getColumn(1).setPreferredWidth(600);
        for (int i=0; i<tablePTL.getColumnCount(); i++) {
        	TableColumn column = tablePTL.getColumnModel().getColumn(i);
        	column.setCellRenderer(centerTableRenderer);
        }
        JScrollPane scrollPTL = new JScrollPane(tablePTL);
    	scrollPTL.setBorder(BorderFactory.createLoweredBevelBorder());

    	//	Create boxes to display the PTL table
    	Box PTLtable1 = Box.createHorizontalBox();
    	Box PTLtable2 = Box.createHorizontalBox();
    	Box PTLtable3 = Box.createHorizontalBox();
    	Box PTLtable4 = Box.createHorizontalBox();
    	Box PTLtable5 = Box.createVerticalBox();
    	Box PTLtable6 = Box.createHorizontalBox();
    	
    	//	Fill the boxes with the appropriate components of the table
    	PTLtable1.add(AddPTL = new JButton("Add"));
    	AddButtonListener addListenerPTL = new AddButtonListener();
    	AddPTL.addActionListener(addListenerPTL);
    	AddPTL.setActionCommand("AddPTL");
    	AddPTL.setToolTipText("Press to submit PTL Name & Location");
    	PTLtable1.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
    	
    	PTLtable2.add(DeletePTL = new JButton("Remove"));
    	DeleteButtonListener deleteListenerPTL = new DeleteButtonListener();
    	DeletePTL.addActionListener(deleteListenerPTL);
    	DeletePTL.setActionCommand("DeletePTL");
    	DeletePTL.setToolTipText("Press to remove PTL Name & Location");
    	
    	PTLtable3.add(PTLtable1);
    	PTLtable3.add(PTLtable2);
    	PTLtable4.add(scrollPTL);
    	PTLtable6.add(ptlName);
    	PTLtable6.add(ptlLoc);
    	PTLtable5.add(PTLtable6);
    	PTLtable5.add(PTLtable3);
    	PTLtable5.add(PTLtable4);
    	
    	PTL.add(PTLtable5);
    	
    	//	Initialize the PTL with the RMG default library
    	//		This library contains H and H2, which cannot be estimated using Benson's group additivity scheme
		PRLVector initialPTL = new PRLVector(0, "Default_H_H2", databasePath.getText()+"/thermo_libraries/primaryThermoLibrary");
		tmodelPTL.updatePRL(initialPTL);
    	
		/*
		 * Primary Transport Library field
		 */
    	//	Create the Primary Transport Library (PTransL) panel
    	JPanel PTransL = new JPanel();
    	PTransL.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Primary Transport Library"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	PTransL.setToolTipText("Override RMG's group-additivity algorithm to compute species' Lennard-Jones parameters");

    	//	Create PTransL Name label
    	JPanel ptranslName = new JPanel();
    	ptranslName.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	
    	//	Populate the label
    	JLabel ptranslNameLabel = new JLabel("Name:");
    	ptranslName.add(ptranslNameLabel);
    	ptranslName.add(ptranslLibName = new JTextField(20));

        //	Create PTransL Location label
        JPanel ptranslLoc = new JPanel();
        ptranslLoc.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	
        //	Populate the label
        JLabel ptranslLocationLabel = new JLabel("Location:");
        ptranslLoc.add(ptranslLocationLabel);
    	ptranslLocationLabel.setToolTipText("Default = " + 
			"RMG/databases/RMG_database/transport_libraries/GRI-Mech3.0");
    	
    	ptranslLoc.add(ptranslPath = new JTextField(20));
    	
    	ptranslLoc.add(ptranslButton = new JButton("Select"));
    	ChangeButtonListener ptranslAddListenerLib = new ChangeButtonListener();
    	ptranslButton.addActionListener(ptranslAddListenerLib);
    	ptranslButton.setActionCommand("ptranslPath");

        //	Create table and scroll panel to store PTransL(s)
        tablePTransL = new JTable(tmodelPTransL = new MyTableModelPRL());
    	tablePTransL.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    	tablePTransL.setPreferredScrollableViewportSize(new Dimension(700,50));
        tablePTransL.getColumnModel().getColumn(0).setPreferredWidth(100);
        tablePTransL.getColumnModel().getColumn(1).setPreferredWidth(600);
        for (int i=0; i<tablePTransL.getColumnCount(); i++) {
        	TableColumn column = tablePTransL.getColumnModel().getColumn(i);
        	column.setCellRenderer(centerTableRenderer);
        }
        JScrollPane scrollPTransL = new JScrollPane(tablePTransL);
    	scrollPTransL.setBorder(BorderFactory.createLoweredBevelBorder());

    	//	Create boxes to display the PTransL table
    	Box PTransLtable1 = Box.createHorizontalBox();
    	Box PTransLtable2 = Box.createHorizontalBox();
    	Box PTransLtable3 = Box.createHorizontalBox();
    	Box PTransLtable4 = Box.createHorizontalBox();
    	Box PTransLtable5 = Box.createVerticalBox();
    	Box PTransLtable6 = Box.createHorizontalBox();
    	
    	//	Fill the boxes with the appropriate components of the table
    	PTransLtable1.add(AddPTransL = new JButton("Add"));
    	AddButtonListener addListenerPTransL = new AddButtonListener();
    	AddPTransL.addActionListener(addListenerPTransL);
    	AddPTransL.setActionCommand("AddPTransL");
    	AddPTransL.setToolTipText("Press to submit PTransportL Name & Location");
    	PTransLtable1.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
    	
    	PTransLtable2.add(DeletePTransL = new JButton("Remove"));
    	DeleteButtonListener deleteListenerPTransL = new DeleteButtonListener();
    	DeletePTransL.addActionListener(deleteListenerPTransL);
    	DeletePTransL.setActionCommand("DeletePTransL");
    	DeletePTransL.setToolTipText("Press to remove PTransportL Name & Location");
    	
    	PTransLtable3.add(PTransLtable1);
    	PTransLtable3.add(PTransLtable2);
    	PTransLtable4.add(scrollPTransL);
    	PTransLtable6.add(ptranslName);
    	PTransLtable6.add(ptranslLoc);
    	PTransLtable5.add(PTransLtable6);
    	PTransLtable5.add(PTransLtable3);
    	PTransLtable5.add(PTransLtable4);
    	
    	PTransL.add(PTransLtable5);
    	
    	//	Initialize the PTransL with the GRI-Mech 3.0 library
    	//		This library contains small-molecule chemistry, which cannot be estimated by the Joback group additivity scheme
		PRLVector initialPTransL = new PRLVector(0, "GRI-Mech 3.0", databasePath.getText()+"/transport_libraries/GRI-Mech3.0");
		tmodelPTransL.updatePRL(initialPTransL);
		
		
    	//	Create the Primary Kinetic Library (PKL) panel
    	JPanel PRL = new JPanel();
    	PRL.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Primary Kinetic Library"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	PRL.setToolTipText("Override RMG's algorithm to compute a reaction's rate coefficient");

    	//	Create PKL Name label
    	JPanel prlName = new JPanel();
    	prlName.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	
    	//	Populate the label
    	JLabel prlNameLabel = new JLabel("Name:");
    	prlName.add(prlNameLabel);
    	prlName.add(prlLibName = new JTextField(20));

        //	Create PRL Location label
        JPanel prlLoc = new JPanel();
        prlLoc.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	
        //	Populate the label
        JLabel prlLocationLabel = new JLabel("Location:");
        prlLoc.add(prlLocationLabel);
    	prlLocationLabel.setToolTipText("Default = " + 
			"RMG/databases/RMG_database/kinetics_libraries/primaryKineticLibrary"); // doesn't exist!
    	
    	prlLoc.add(prlPath = new JTextField(20));
    	
    	prlLoc.add(prlButton = new JButton("Select"));
    	ChangeButtonListener prlAddListenerLib = new ChangeButtonListener();
    	prlButton.addActionListener(prlAddListenerLib);
    	prlButton.setActionCommand("prlPath");

        //	Create table and scroll panel to store PRL(s)
        tablePRL = new JTable(tmodelPRL = new MyTableModelPRL());
    	tablePRL.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    	tablePRL.setPreferredScrollableViewportSize(new Dimension(700,50));
        tablePRL.getColumnModel().getColumn(0).setPreferredWidth(100);
        tablePRL.getColumnModel().getColumn(1).setPreferredWidth(600);
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
    	Box PRLtable6 = Box.createHorizontalBox();
    	
    	//	Fill the boxes with the appropriate components of the table
    	PRLtable1.add(AddPRL = new JButton("Add"));
    	AddButtonListener addListenerPRL = new AddButtonListener();
    	AddPRL.addActionListener(addListenerPRL);
    	AddPRL.setActionCommand("AddPRL");
    	AddPRL.setToolTipText("Press to submit PRL Name & Location");
    	PRLtable1.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
    	
    	PRLtable2.add(DeletePRL = new JButton("Remove"));
    	DeleteButtonListener deleteListenerPRL = new DeleteButtonListener();
    	DeletePRL.addActionListener(deleteListenerPRL);
    	DeletePRL.setActionCommand("DeletePRL");
    	DeletePRL.setToolTipText("Press to remove PRL Name & Location");
    	
    	PRLtable3.add(PRLtable1);
    	PRLtable3.add(PRLtable2);
    	PRLtable4.add(scrollPRL);
    	PRLtable6.add(prlName);
    	PRLtable6.add(prlLoc);
    	PRLtable5.add(PRLtable6);
    	PRLtable5.add(PRLtable3);
    	PRLtable5.add(PRLtable4);
    	
    	PRL.add(PRLtable5);
    	
    	//	Create the Seed Mechanism (SM) panel
    	JPanel SM = new JPanel();
    	SM.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Seed Mechanism"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	SM.setToolTipText("Populate the mechanism with a list of species / reaction, without using RMG's rate-based algorithm");

    	//	Create SM Name label
    	JPanel smName = new JPanel();
    	smName.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	
    	//	Populate the label
    	JLabel smNameLabel = new JLabel("Name:");
    	smName.add(smNameLabel);
    	smName.add(smLibName = new JTextField(20));

        //	Create SM Location label
        JPanel smLoc = new JPanel();
        smLoc.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	
        //	Populate the label
        JLabel smLocationLabel = new JLabel("Location:");
        smLoc.add(smLocationLabel);
    	smLocationLabel.setToolTipText("Default = " + 
			"RMG/databases/RMG_database/kinetics_libraries/combustion_core/version5");
    	
    	smLoc.add(smPath = new JTextField(20));
    	
    	smLoc.add(smButton = new JButton("Select"));
    	ChangeButtonListener smAddListenerLib = new ChangeButtonListener();
    	smButton.addActionListener(smAddListenerLib);
    	smButton.setActionCommand("smPath");
    	
    	JPanel smReact = new JPanel();
    	smReact.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	smReact.add(new JLabel("Generate reactions between Seed Mechanism species?"));
    	smReact.add(smCombo = new JComboBox(yesnoOptions));

        //	Create table and scroll panel to store SM(s)
        tableSM = new JTable(tmodelSM = new MyTableModelSM());
    	tableSM.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    	tableSM.setPreferredScrollableViewportSize(new Dimension(700,50));
        tableSM.getColumnModel().getColumn(0).setPreferredWidth(100);
        tableSM.getColumnModel().getColumn(1).setPreferredWidth(150);
        tableSM.getColumnModel().getColumn(2).setPreferredWidth(450);
        for (int i=0; i<tableSM.getColumnCount(); i++) {
        	TableColumn column = tableSM.getColumnModel().getColumn(i);
        	column.setCellRenderer(centerTableRenderer);
        }
        JScrollPane scrollSM = new JScrollPane(tableSM);
    	scrollSM.setBorder(BorderFactory.createLoweredBevelBorder());

    	//	Create boxes to display the PRL table
    	Box SMtable1 = Box.createHorizontalBox();
    	Box SMtable2 = Box.createHorizontalBox();
    	Box SMtable3 = Box.createHorizontalBox();
    	Box SMtable4 = Box.createHorizontalBox();
    	Box SMtable5 = Box.createVerticalBox();
    	Box SMtable6 = Box.createHorizontalBox();
    	
    	//	Fill the boxes with the appropriate components of the table
    	SMtable1.add(AddSM = new JButton("Add"));
    	AddButtonListener addListenerSM = new AddButtonListener();
    	AddSM.addActionListener(addListenerSM);
    	AddSM.setActionCommand("AddSM");
    	AddSM.setToolTipText("Press to submit Seed Mechanism Name & Location");
    	SMtable1.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
    	
    	SMtable2.add(DeleteSM = new JButton("Remove"));
    	DeleteButtonListener deleteListenerSM = new DeleteButtonListener();
    	DeleteSM.addActionListener(deleteListenerSM);
    	DeleteSM.setActionCommand("DeleteSM");
    	DeleteSM.setToolTipText("Press to remove Seed Mechanism Name & Location");
    	
    	SMtable3.add(SMtable1);
    	SMtable3.add(SMtable2);
    	SMtable4.add(scrollSM);
    	SMtable6.add(smName);
    	SMtable6.add(smLoc);
    	SMtable5.add(SMtable6);
    	SMtable5.add(smReact);
    	SMtable5.add(SMtable3);
    	SMtable5.add(SMtable4);
    	
    	SM.add(SMtable5);

    	//	Create & fill box for Initialization tab
    	Box TabTotal = Box.createVerticalBox();
    	TabTotal.add(Database);
    	TabTotal.add(PTL);
    	TabTotal.add(PTransL);
    	TabTotal.add(PRL);
    	TabTotal.add(SM);
    	
        JComboBox[] allTab = {smCombo};
        initializeJCB(allTab);
        
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
    
    public JComponent createTabSolver() {
    	//	Create cellrenderers for JComboBox/JTable - CENTER text
    	ListCellRenderer centerJComboBoxRenderer = new CenterJCBRenderer();
    	
    	//	Create the Dynamic Simulator (DS) panel
    	JPanel DS = new JPanel();
    	DS.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Dynamic Simulator"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	
    	JPanel Prune = new JPanel();
    	Prune.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Pruning Options"),
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
    	aToleranceLabel.setToolTipText("Suggested value = 1E-15");
    	
    	aTolPanel.add(aTolerance = new JTextField());
    	aTolerance.setPreferredSize(new Dimension(100,25));
        aTolerance.addActionListener(this);
        aTolerance.setHorizontalAlignment(JTextField.CENTER);
        aTolerance.setText("1E-15");
    	
    	//	Create the DS subpanel: rTolPanel
        JPanel rTolPanel = new JPanel();
        rTolPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
        
        //	Populate the rTolPanel subpanel
        JLabel rToleranceLabel = new JLabel("Select relative tolerance");
        rTolPanel.add(rToleranceLabel);
        rToleranceLabel.setToolTipText("Suggested value = 1E-6");
    	
        rTolPanel.add(rTolerance = new JTextField());
    	rTolerance.setPreferredSize(new Dimension(100,25));
        rTolerance.addActionListener(this);
        rTolerance.setHorizontalAlignment(JTextField.CENTER);
        rTolerance.setText("1E-6");
    	
        //	Create the DS subpanel: interConv
        JPanel interConv = new JPanel();
        interConv.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
        //	Populate the subpanel
        JLabel multiConvLabel = new JLabel("Intermediate Conversions/TimeSteps");
        interConv.add(multiConvLabel);
        interConv.add(interCombo = new JComboBox(interOptions));
        interCombo.addActionListener(this);
        interCombo.setRenderer(centerJComboBoxRenderer);
        interCombo.setActionCommand("Pruning");
        
    	/*
    	 * Pruning options
    	 */
        JPanel indivSteps = new JPanel();
        indivSteps.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
        indivSteps.add(new JLabel("Specific time/conversion steps"));
        indivSteps.add(textIndivSteps = new JTextField());
        textIndivSteps.setPreferredSize(new Dimension(100,25));
        textIndivSteps.setHorizontalAlignment(JTextField.CENTER);
        textIndivSteps.setEnabled(false);
        
        JPanel terminTol = new JPanel();
        terminTol.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
        terminTol.add(new JLabel("Termination Tolerance"));
        terminTol.add(textTerminTol = new JTextField());
        textTerminTol.setPreferredSize(new Dimension(100,25));
        textTerminTol.setHorizontalAlignment(JTextField.CENTER);
        textTerminTol.setText("1.0E4");
        textTerminTol.setEnabled(false);
        
        JPanel pruneTol = new JPanel();
        pruneTol.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
        pruneTol.add(new JLabel("Pruning tolerance"));
        pruneTol.add(textPruneTol = new JTextField());
        textPruneTol.setPreferredSize(new Dimension(100,25));
        textPruneTol.setHorizontalAlignment(JTextField.CENTER);
        textPruneTol.setText("1.0E-15");
        textPruneTol.setEnabled(false);
        
        JPanel minSpc4Prune = new JPanel();
        minSpc4Prune.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
        minSpc4Prune.add(new JLabel("Minimum # of species before initiating pruning"));
        minSpc4Prune.add(textMinSpc = new JTextField());
        textMinSpc.setPreferredSize(new Dimension(100,25));
        textMinSpc.setHorizontalAlignment(JTextField.CENTER);
        textMinSpc.setText("1000");
        textMinSpc.setEnabled(false);
        
        JPanel maxSpcNEdge = new JPanel();
        maxSpcNEdge.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
        maxSpcNEdge.add(new JLabel("Maximum # of species in edge after pruning"));
        maxSpcNEdge.add(textMaxSpc = new JTextField());
        textMaxSpc.setPreferredSize(new Dimension(100,25));
        textMaxSpc.setHorizontalAlignment(JTextField.CENTER);
        textMaxSpc.setText("1000");
        textMaxSpc.setEnabled(false);
        
    	Box pruningOptionsBox = Box.createVerticalBox();
    	pruningOptionsBox.add(terminTol);
    	pruningOptionsBox.add(pruneTol);
    	pruningOptionsBox.add(minSpc4Prune);
    	pruningOptionsBox.add(maxSpcNEdge);
    	
    	//	Create boxes for DS panel
        Box ds = Box.createVerticalBox();
        //	Populate the boxes
    	ds.add(dsSolver);
    	ds.add(aTolPanel);
    	ds.add(rTolPanel);
    	//ds.add(interConv);
    	ds.add(interConv);
    	ds.add(indivSteps);
    	DS.add(ds);
    	
    	Prune.add(pruningOptionsBox);
    	
    	Box dsNprune = Box.createHorizontalBox();
    	dsNprune.add(DS);
    	dsNprune.add(Prune);
    	
    	// Create the total pressure dependence panel
    	JPanel totalPDep = new JPanel();
    	totalPDep.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Pressure dependence"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	
    	//	Create the Pressure dependence (Pdep) model panel
    	JPanel Pdep = new JPanel();
    	Pdep.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	//	Populate the ME panel
    	JLabel modelLabel = new JLabel("Choose pressure dependence model");
    	Pdep.add(modelLabel);
    	modelLabel.setToolTipText("Default = Off");
    	Pdep.add(pdepCombo = new JComboBox(pdepOptions));
    	pdepCombo.setActionCommand("pDep");
    	
    	//	Create the Spectroscopic Data Estimator (SDE) subpanel
    	JPanel SDE = new JPanel();
    	SDE.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	//	Populate the subpanel
    	JLabel sdeLabel = new JLabel("Spectroscopic Data Estimator");
    	SDE.add(sdeLabel);
    	SDE.add(sdeCombo = new JComboBox(sdeOptions));
    	sdeCombo.setEnabled(false);
    	
    	//	Create the PDep Kinetics Model (PDKM) subpanel
    	JPanel PDKM = new JPanel();
    	PDKM.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	//	Populate the subpanel
    	JLabel pdkmLabel = new JLabel("PDep Kinetics Model");
    	PDKM.add(pdkmLabel);
    	PDKM.add(pdkmCombo = new JComboBox(pdkmOptions));
    	pdkmCombo.setEnabled(false);
    	pdkmCombo.setActionCommand("pdkm");
    	
    	//	Create the "decrease grain size" subpanel
    	JPanel GrainSize = new JPanel();
    	GrainSize.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	JLabel grainsizeLabel = new JLabel("Decrease Grain Size");
    	GrainSize.add(grainsizeLabel);
    	GrainSize.add(grainsizeCombo = new JComboBox(yesnoOptions));
    	grainsizeCombo.setEnabled(false);
    	grainsizeLabel.setToolTipText("If the pressure-dependent rate > 2*(High-P-Limit rate):" +
    			" Increase the number of energy steps used in solving the Master Equation" +
    			" and re-run the fame.exe calculation.");
    	
    	JPanel Cheby = new JPanel();
    	Cheby.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	
    	//	Create boxes to store each of the pieces of data
    	Box boxTP = Box.createVerticalBox();
    	Box boxMin = Box.createVerticalBox();
    	Box boxMax = Box.createVerticalBox();
    	Box boxUnits = Box.createVerticalBox();
    	Box boxNumGenerate = Box.createVerticalBox();
    	Box boxNumReport = Box.createVerticalBox();
    	Box boxAll = Box.createVerticalBox();
    	Box boxListorMaxMin = Box.createVerticalBox();
    	
    	//	Add labels and data fields to boxes
    	boxTP.add(new JLabel(" "));
    	boxTP.add(new JLabel("Temperature"));
    	boxTP.add(new JLabel(" "));
    	boxTP.add(new JLabel(" "));
    	boxTP.add(new JLabel("Pressure"));
    	
    	boxMin.add(new JLabel("  Min. Value"));
    	boxMin.add(chebyTMin = new JTextField());
    	boxMin.add(new JLabel(" "));
    	boxMin.add(chebyPMin = new JTextField());
    	chebyTMin.setPreferredSize(new Dimension (25,25));
    	chebyTMin.setHorizontalAlignment(JTextField.CENTER);
        chebyTMin.addActionListener(this);
    	chebyPMin.setPreferredSize(new Dimension (25,25));
    	chebyPMin.setHorizontalAlignment(JTextField.CENTER);
        chebyPMin.addActionListener(this);
    	
    	boxMax.add(new JLabel("  Max. Value"));
    	boxMax.add(chebyTMax = new JTextField());
    	boxMax.add(new JLabel(" "));
    	boxMax.add(chebyPMax = new JTextField());
    	chebyTMax.setPreferredSize(new Dimension (25,25));
    	chebyTMax.setHorizontalAlignment(JTextField.CENTER);
        chebyTMax.addActionListener(this);
    	chebyPMax.setPreferredSize(new Dimension (25,25));
    	chebyPMax.setHorizontalAlignment(JTextField.CENTER);
        chebyPMax.addActionListener(this);

    	boxUnits.add(new JLabel("Units    "));
    	boxUnits.add(chebyTUnits = new JComboBox(tempUnits));
    	boxUnits.add(new JLabel(" "));
    	boxUnits.add(chebyPUnits = new JComboBox(pressUnits));
    	
    	boxNumGenerate.add(new JLabel("  # to Generate"));
    	boxNumGenerate.add(chebyTGen = new JTextField());
    	boxNumGenerate.add(new JLabel(" "));
    	boxNumGenerate.add(chebyPGen = new JTextField());
    	chebyTGen.setPreferredSize(new Dimension (25,25));
    	chebyTGen.setHorizontalAlignment(JTextField.CENTER);
        chebyTGen.addActionListener(this);
    	chebyPGen.setPreferredSize(new Dimension (25,25));
    	chebyPGen.setHorizontalAlignment(JTextField.CENTER);
        chebyPGen.addActionListener(this);
    	
    	boxNumReport.add(new JLabel("  # to Report"));
    	boxNumReport.add(chebyTRep = new JTextField());
    	boxNumReport.add(new JLabel(" "));
    	boxNumReport.add(chebyPRep = new JTextField());
    	chebyTRep.setPreferredSize(new Dimension(25,25));
    	chebyTRep.setHorizontalAlignment(JTextField.CENTER);
        chebyTRep.addActionListener(this);
    	chebyPRep.setPreferredSize(new Dimension(25,25));
    	chebyPRep.setHorizontalAlignment(JTextField.CENTER);
        chebyPRep.addActionListener(this);
        
        boxAll.add(new JLabel("                    Values"));
        boxAll.add(chebyT = new JTextField());
        boxAll.add(new JLabel(" "));
        boxAll.add(chebyP = new JTextField());
        chebyT.setPreferredSize(new Dimension(150,25));
        chebyT.setHorizontalAlignment(JTextField.CENTER);
        chebyT.addActionListener(this);
        chebyP.setPreferredSize(new Dimension(150,25));
        chebyP.setHorizontalAlignment(JTextField.CENTER);
        chebyP.addActionListener(this);
        
        boxListorMaxMin.add(new JLabel("List/Limits              "));
        boxListorMaxMin.add(chebyTCombo = new JComboBox(listormaxmin));
    	chebyTCombo.setActionCommand("tcheby");
        boxListorMaxMin.add(new JLabel(" "));
        boxListorMaxMin.add(chebyPCombo = new JComboBox(listormaxmin));
    	chebyPCombo.setActionCommand("pcheby");
    	
    	//	Add the boxes to the "Cheby" subpanel
    	Cheby.add(boxTP);
    	Cheby.add(boxNumReport);
    	Cheby.add(boxNumGenerate);
    	Cheby.add(boxListorMaxMin);
    	Cheby.add(boxMin);
    	Cheby.add(boxMax);
    	Cheby.add(boxAll);
    	Cheby.add(boxUnits);
    	
    	Box pDepBox = Box.createVerticalBox();
    	pDepBox.add(Pdep);
    	pDepBox.add(SDE);
    	pDepBox.add(PDKM);
    	pDepBox.add(GrainSize);
    	pDepBox.add(Cheby);
    	totalPDep.add(pDepBox);
    	
        JComboBox[] allTab = {pdepCombo, sdeCombo, pdkmCombo, chebyTUnits, chebyPUnits,
        		chebyTCombo, chebyPCombo, grainsizeCombo, interCombo};
        initializeJCB(allTab);
        
		JComponent[] pdkmComps = {chebyTMin, chebyTMax, chebyTUnits, chebyTGen, chebyTRep,
				chebyPMin, chebyPMax, chebyPUnits, chebyPGen, chebyPRep,
				chebyT, chebyTCombo, chebyP, chebyPCombo};
		disableComponents(pdkmComps);
    	
    	Box TabTotal = Box.createVerticalBox();
    	TabTotal.add(dsNprune);
    	TabTotal.add(totalPDep);    	
    	
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
	public JComponent createTabTermination() {
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
    public JComponent createTabSensitivity() {
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
			} else if ("DeletePTL".equals(e.getActionCommand())) {
				int highlightRow = tablePTL.getSelectedRow();
				if (highlightRow == -1) {
					System.out.println("Please select a row to delete");
				} else
					tmodelPTL.deleteRow(highlightRow);
			}
			else if ("DeleteSM".equals(e.getActionCommand())) {
				int highlightRow = tableSM.getSelectedRow();
				if (highlightRow == -1) {
					System.out.println("Please select a row to delete");
				} else
					tmodelSM.deleteRow(highlightRow);
			}
			else if ("DeleteFS".equals(e.getActionCommand())) {
				int highlightRow = tableFS.getSelectedRow();
				if (highlightRow == -1)
					System.out.println("Please select a row to delete");
				else
					tmodelFS.deleteRow(highlightRow);
			}
			else if ("DeletePTransL".equals(e.getActionCommand())) {
				int highlightRow = tablePTransL.getSelectedRow();
				if (highlightRow == -1) System.out.println("Please select a row to delete");
				else tmodelPTransL.deleteRow(highlightRow);
			}
		}
	}
    
    class AddButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			// IF the user adds a primary kinetic library ...
			if ("AddPRL".equals(e.getActionCommand())) {
				// Extract the information
				int prlInput0 = tmodelPRL.nextEmptyRow+1;
				String prlInput1 = prlLibName.getText();
				String prlInput2 = prlPath.getText();
				// Check that all information is present
				if (prlInput1.equals("") || prlInput2.equals("")) {
					System.out.println("Please input the Name and Location of a Primary Kinetic Library");
				} else {
					PRLVector prlEntry = new PRLVector(prlInput0, prlInput1, prlInput2);
					tmodelPRL.updatePRL(prlEntry);
					prlLibName.setText("");
    				prlPath.setText("");
				}
			} else if ("AddPTL".equals(e.getActionCommand())) {
				// Extract the information
				int ptlInput0 = tmodelPTL.nextEmptyRow+1;
				String ptlInput1 = ptlLibName.getText();
				String ptlInput2 = ptlPath.getText();
				// Check that all information is present
				if (ptlInput1.equals("") || ptlInput2.equals("")) {
					System.out.println("Please input the Name and Location of a Primary Thermo Library");
				} else {
					PRLVector ptlEntry = new PRLVector(ptlInput0, ptlInput1, ptlInput2);
					tmodelPTL.updatePRL(ptlEntry);
					ptlLibName.setText("");
    				ptlPath.setText("");
				}
			} else if ("AddSM".equals(e.getActionCommand())) {
				// Extract the information
				int smInput0 = tmodelSM.nextEmptyRow+1;
				String smInput1 = smLibName.getText();
				String smInput2 = (String)smCombo.getSelectedItem();
				String smInput3 = smPath.getText();
				// Check that all information is present
				if (smInput1.equals("") || smInput2.equals("") || smInput3.equals("")) {
					System.out.println("Please input the Name and Location of a Seed Mechanism");
				} else {
					SMVector smEntry = new SMVector(smInput0, smInput1, smInput2, smInput3);
					tmodelSM.updateSM(smEntry);
					smLibName.setText("");
					smCombo.setSelectedIndex(0);
    				smPath.setText("");
				}
			// ELSE IF the user adds a species ...
			} else if ("AddInput".equals(e.getActionCommand())) {
				// Extract the information
				int icInput0 = tmodelInput.nextEmptyRow+1;
    			String icInput1 = NameInput.getText();
    			
    			boolean renameName = false;
    			try {
    				int doesNameBeginWithNumber = Integer.parseInt(icInput1.substring(0,1));
    				renameName = true;
    			} catch (NumberFormatException e1) {
    				// We're good
    			}    			
    			
    			String icInput2 = ConcInput.getText();
    			String icInput3 = (String)UnitsInput.getSelectedItem();
    			String icInput4 = (String)ReactiveInput.getSelectedItem();
    			String icInput5 = "";
    			//	Ask user for adjacency list if species is marked "Reactive"
    			if (!icInput4.equals("Inert")) {
	    			if (!speciesAdjList.getText().equals("")) {
	    				icInput5 = speciesAdjList.getText();
	    			} else if (!InChIInput.getText().equals("")) {
	    				String inchi = InChIInput.getText();
	    				icInput5 = Species.inchi2AdjList(inchi);
	    			} else {
	    			    File molFile = null;
	    			    molFile = askUserForInput("Select corresponding .mol file", false, null);
	    			    if (molFile != null) icInput5 = Species.mol2AdjList(molFile.getAbsolutePath(), molFile.getAbsolutePath());
	    			}
    			} else {
    				String[] inerts = {"N2", "Ar", "Ne", "He"};
    				boolean correctInertFormat = false;
    				for (int i=0; i<inerts.length; i++) {
    					if (icInput1.equals(inerts[i])) correctInertFormat = true;
    				}
    				if (!correctInertFormat) {
    					System.out.println("ERROR: Adding Inert Gas species to model - RMG recognizes the following species as Inert Gas (case sensitive): N2, Ar, Ne, He");
    					return;
    				}
    			}
    			// If a species is marked "Reactive-User", ask user for location of IncludeSpecies.txt file
    			if (icInput4.equals("Reactive-User")) {
    				// If the path is not known
    				if (isPath.getText().equals("")) {
        				File isFile = null;
        				isFile = askUserForInput("Select IncludeSpecies.txt file", false, null);
        			    if (isFile != null) {
        			    	isPath.setText(isFile.getPath());
        			    	isPath.setEnabled(true);
        			    }
    				}
    			}
    			// Check that all information is present
				if (icInput1.equals("") || icInput2.equals("") || (icInput5.equals("") & !icInput4.equals("Inert"))) {
        			System.out.println("Please input the molecule's Name, Concentration, and Adjacency List");
        		} else if (renameName) {
    				System.out.println("\nA species name should not begin with a number." +
    						" Please rename species: " + icInput1 + "\n");
        		}
				else {
        			ICVector icEntry = new ICVector(icInput0, icInput1, icInput2, icInput3, icInput4, icInput5);
        			tmodelInput.updateIC(icEntry);
        			NameInput.setText("");
        			InChIInput.setText("");
        			ConcInput.setText("");
        			UnitsInput.setSelectedIndex(0);
        			ReactiveInput.setSelectedIndex(0);
        			speciesAdjList.setText("");
        			if (!icInput4.equals("Inert")) limitingReactant();
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
			// ELSE IF the user adds a forbidden structure
			else if ("AddFS".equals(e.getActionCommand())) {
				int fsInput0 = tmodelFS.nextEmptyRow+1;
				String fsInput1 = FSName.getText();
				String fsInput2 = FSAdjList.getText();
				if (fsInput1.equals("") || fsInput2.equals("")) {
					System.out.println("Please input the name and adjacency list of the forbidden structure");
				} else {
					FSVector fsEntry = new FSVector(fsInput0, fsInput1, fsInput2);
					tmodelFS.updateFS(fsEntry);
					FSName.setText("");
					FSAdjList.setText("");
				}
			}
			else if ("AddPTransL".equals(e.getActionCommand())) {
				// Extract the information
				int ptranslInput0 = tmodelPTransL.nextEmptyRow+1;
				String ptranslInput1 = ptranslLibName.getText();
				String ptranslInput2 = ptranslPath.getText();
				// Check that all information is present
				if (ptranslInput1.equals("") || ptranslInput2.equals("")) {
					System.out.println("Please input the Name and Location of a Primary Transport Library");
				} else {
					PRLVector ptranslEntry = new PRLVector(ptranslInput0, ptranslInput1, ptranslInput2);
					tmodelPTransL.updatePRL(ptranslEntry);
					ptranslLibName.setText("");
    				ptranslPath.setText("");
				}
			}
		}
	}
    
    class ChangeButtonListener implements ActionListener {
    	public void actionPerformed(ActionEvent e) {
    		// IF the user selects to switch the IncludeSpecies.txt path
			if ("isPath".equals(e.getActionCommand())) {
				File path = null;
				path = askUserForInput("Select IncludeSpecies.txt file", false, null);
				if (path != null) isPath.setText(path.getPath());
			// ELSE IF user selects to switch the database path
			}  else if ("databasePath".equals(e.getActionCommand())) {
				File path = null;
				String workingDirectory = System.getProperty("RMG.workingDirectory");
				path = askUserForInput("Select database folder", true, workingDirectory + "/databases");
				if (path != null) databasePath.setText(path.getPath());
			// If user selects to switch Primary Thermo Library path
			} else if ("ptlPath".equals(e.getActionCommand())) {
				File path = null;
				path = askUserForInput("Select PrimaryThermoLibrary folder", true, databasePath.getText()+"/thermo_libraries");
				if (path != null) ptlPath.setText(path.getPath());
			// If user selects to switch Primary Kinetic Library path
			} else if ("prlPath".equals(e.getActionCommand())) {
				File path = null;
				path = askUserForInput("Select PrimaryKineticLibrary folder", true, databasePath.getText()+"/kinetics_libraries");
				if (path != null) prlPath.setText(path.getPath());
			} else if ("smPath".equals(e.getActionCommand())) {
				File path = null;
				path = askUserForInput("Select SeedMechanism folder", true, databasePath.getText()+"/kinetics_libraries");
				if (path != null) smPath.setText(path.getPath());
			} else if ("ptranslPath".equals(e.getActionCommand())) {
				File path = null;
				path = askUserForInput("Select PrimaryTransportLibrary folder", true, databasePath.getText()+"/transport_libraries");
				if (path != null) ptranslPath.setText(path.getPath());
			}
    	}
    }
    
    public File askUserForInput (String popupwindowTitle, boolean folder, String additionalPath) {
    	//String workingDirectory = System.getProperty("RMG.workingDirectory");
    	String workingDirectory = System.getProperty("user.dir");
    	String completePath ="";
    	if (additionalPath != null)
    		completePath = additionalPath;
    	else
    		completePath = workingDirectory;
    	JFileChooser chooser = new JFileChooser(completePath);
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
    		if (!tmodelInput.getValueAt(j,3).equals("Inert")) {
    			++nonInertCount;
    			SpeciesConvName.insertItemAt(tmodelInput.getValueAt(j,0),nonInertCount);
    		}
    	}
    }
        
    public void createConditionFile(boolean saveAndRun) {
    	//String workingDir = System.
    	
    	//	Let's create the condition.txt file!
    	String conditionFile = "";
    	String[] tempString;
    	
    	//	Add user supplied comments
    	conditionFile += "// User supplied comments:\r";
    	if (headerComments.getText().equals("")) System.out.println("Warning: Writing condition.txt file: No comments being read in.");
    	String[] header4Condition = headerComments.getText().split("[\n]");
    	for (int i=0; i<header4Condition.length; i++) {
    		conditionFile += "// " + header4Condition[i] + "\r";
    	}
    	//	Add date and time to the file
    	String DATE_FORMAT_NOW = "yyyy-MM-dd HH:mm:ss";
        SimpleDateFormat sdf = new SimpleDateFormat(DATE_FORMAT_NOW);
        Calendar cal = Calendar.getInstance();
        conditionFile += "//\r// File created at: " + sdf.format(cal.getTime()) + "\r\r";
        
        // String storing the folder of the Main Database (w.r.t. $RMG/databases)
        String stringMainDatabase = "";
        
        //	Add the name of the database and primary thermo library
        if (databasePath.getText().equals("")) {
        	System.out.println("ERROR: Writing condition.txt file: Could not read Main Database (Thermochemical Libraries tab)");
        	return;
        }
        else {
        	tempString = databasePath.getText().split("[\\\\,/]");
        	stringMainDatabase = tempString[tempString.length-1];
        	conditionFile += "Database: " + stringMainDatabase + "\r";
        }
        
        //	Add the user-specified maximum number of carbons, oxygens, and radicals / species
        if (!maxCarbonNum.getText().equals(""))
        	conditionFile += "\rMaxCarbonNumberPerSpecies: " + maxCarbonNum.getText();
        if (!maxOxygenNum.getText().equals(""))
        	conditionFile += "\rMaxOxygenNumberPerSpecies: " + maxOxygenNum.getText();
        if (!maxRadicalNum.getText().equals(""))
        	conditionFile += "\rMaxRadicalNumberPerSpecies: " + maxRadicalNum.getText();
        if (!maxSulfurNum.getText().equals(""))
        	conditionFile += "\rMaxSulfurNumberPerSpecies: " + maxSulfurNum.getText();
        if (!maxSiliconNum.getText().equals(""))
        	conditionFile += "\rMaxSiliconNumberPerSpecies: " + maxSiliconNum.getText();
        if (!maxHeavyAtom.getText().equals(""))
        	conditionFile += "\rMaxHeavyAtomPerSpecies: " + maxHeavyAtom.getText();
                
    	//	Add the name(s)/location(s) of the primary thermo library
        String rmgEnvVar = System.getProperty("RMG.workingDirectory");
        String ptlReferenceDirectory_linux = rmgEnvVar + "/databases/" + stringMainDatabase + "/thermo_libraries/";
        String ptlReferenceDirectory_windows = rmgEnvVar + "\\databases\\" + stringMainDatabase + "\\thermo_libraries\\";
    	conditionFile += "\r\rPrimaryThermoLibrary:\r";

		if (tablePTL.getRowCount()==0) {
        	System.out.println("Warning: Writing condition.txt file: Could not read Primary Thermo Library (Thermochemical Libraries tab)");
		} else {
    		for (int k=0; k<tablePTL.getRowCount(); k++) {
    			conditionFile += "Name: " + tablePTL.getValueAt(k,0) + "\r" + "Location: ";
    			String ptlDir = (String)tablePTL.getValueAt(k,1);
    	        if (ptlDir.toLowerCase().startsWith(ptlReferenceDirectory_linux.toLowerCase()) ||
    	        		ptlDir.toLowerCase().startsWith(ptlReferenceDirectory_windows.toLowerCase())) {
    	        	int startIndex = ptlReferenceDirectory_linux.length();
    	        	conditionFile += ptlDir.substring(startIndex) + "\r";
    	        } else {
    	        	conditionFile += ptlDir + "\r";
    	        }
    		}
		}
		conditionFile += "END\r\r";
		
    	//	Add the name(s)/location(s) of the primary transport library
        String ptranslReferenceDirectory_linux = rmgEnvVar + "/databases/" + stringMainDatabase + "/transport_libraries/";
        String ptranslReferenceDirectory_windows = rmgEnvVar + "\\databases\\" + stringMainDatabase + "\\transport_libraries\\";
    	conditionFile += "PrimaryTransportLibrary:\r";

		if (tablePTransL.getRowCount()==0) {
        	System.out.println("Warning: Writing condition.txt file: Could not read Primary Transport Library (Thermochemical Libraries tab)");
		} else {
    		for (int k=0; k<tablePTransL.getRowCount(); k++) {
    			conditionFile += "Name: " + tablePTransL.getValueAt(k,0) + "\r" + "Location: ";
    			String ptranslDir = (String)tablePTransL.getValueAt(k,1);
    	        if (ptranslDir.toLowerCase().startsWith(ptranslReferenceDirectory_linux.toLowerCase()) ||
    	        		ptranslDir.toLowerCase().startsWith(ptranslReferenceDirectory_windows.toLowerCase())) {
    	        	int startIndex = ptranslReferenceDirectory_linux.length();
    	        	conditionFile += ptranslDir.substring(startIndex) + "\r";
    	        } else {
    	        	conditionFile += ptranslDir + "\r";
    	        }
    		}
		}
		conditionFile += "END\r\r";
		
		/*
		 * 4APR2010
		 */
		// Forbidden structures
		if (tableFS.getRowCount()==0) {
			System.out.println("Writing condition.txt file: No forbidden structures specified by user.");
		}
        else {
        	conditionFile += "ForbiddenStructures:\r";
	        for (int i=0; i<tableFS.getRowCount(); i++) {
	        	conditionFile += 
	        		tableFS.getValueAt(i,0) + "\r" +
	        		tableFS.getValueAt(i,1) + "\r\r";
	        }
	        conditionFile += "END\r\r";
        }
		
		conditionFile += "ReadRestart: " + readRestartOnOff.getSelectedItem() + "\r";
		conditionFile += "WriteRestart: " + writeRestartOnOff.getSelectedItem() + "\r\r";
        
        //	Add the temperature model, list of temperatures, and units
        conditionFile += "TemperatureModel: " + tempModelCombo.getSelectedItem() +
        	" (" + tempConstantUnit.getSelectedItem() + ")";
        if (tempConstant.getText().equals("")) {
        	System.out.println("ERROR: Writing condition.txt file: Could not read list of temperatures (Initial Conditions tab)");
        	return;
        }
        else {
        	String[] tempLines = tempConstant.getText().split("[\n]");
        	StringTokenizer st = null;
        	for (int i=0; i<tempLines.length; i++) {
	        	st = new StringTokenizer(tempLines[i]);
	        	while (st.hasMoreTokens()) {
	        		String line = st.nextToken();
	        		if (line.endsWith(",")) line = line.substring(0,line.length()-1).trim();
	        		conditionFile += " " + line;
	        	}
	        }
        }
        conditionFile += "\r";
        
        //	Add the pressure model, list of pressures, and units
        conditionFile += "PressureModel: " + pressModelCombo.getSelectedItem() +
        	" (" + pressConstantUnit.getSelectedItem() + ")";
        if (pressConstant.getText().equals("")) {
        	System.out.println("ERROR: Writing condition.txt file: Could not read list of pressures (Initial Conditions tab)");
        	return;
        }
        else {
	        String[] pressLines = pressConstant.getText().split("[\n]");
	        StringTokenizer st = null;
	        for (int i=0; i<pressLines.length; i++) {
	        	st = new StringTokenizer(pressLines[i]);
	        	while (st.hasMoreTokens()) {
	        		String line = st.nextToken();
	        		if (line.endsWith(",")) line = line.substring(0,line.length()-1).trim();
	        		conditionFile += " " + line;
	        	}
	        }
        }
        conditionFile += "\r\r";
        
        //	Add the equation of state (if equal to liquid)
        if (eosCombo.getSelectedItem().equals("Liquid")) {
        	conditionFile += "EquationOfState: Liquid\r";
        }
        
        //	Add whether to generate InChIs or not
        if (inchiCombo.getSelectedItem().equals("Yes"))
        	conditionFile += "InChIGeneration: on\r";
        
        //	Add whether to turn solvation on or not
        
        //	Add whether to run on-the-fly quantum jobs or not
        
        //	Determine if each species is reactive, unreactive, or inert
        int nonInertCount = 0;
        boolean unreactiveStatus = false;
        String reactiveSpecies = "";
        String inertSpecies = "";
        if (tableInput.getRowCount()==0) {
        	System.out.println("ERROR: Writing condition.txt file: Could not read table of reactants (Initial Conditions tab)");
        	return;
        }
        else {
	        for (int i=0; i<tableInput.getRowCount(); i++) {
	        	if (tableInput.getValueAt(i,3).equals("Inert")) {
	        		inertSpecies += tableInput.getValueAt(i,0) + " " +
	        			tableInput.getValueAt(i,1) + " (" +
	        			tableInput.getValueAt(i,2) + ")\r";
	        	} else if (tableInput.getValueAt(i,3).equals("Reactive-RMG")) {
	        		++nonInertCount;
	        		reactiveSpecies += "(" + nonInertCount + ") " +
	        		tableInput.getValueAt(i,0) + " " +
	    			tableInput.getValueAt(i,1) + " (" +
	    			tableInput.getValueAt(i,2) + ")\r" +
	    			tableInput.getValueAt(i,4) + "\r\r";
	        	} else if (tableInput.getValueAt(i,3).equals("Reactive-User")) {
	        		unreactiveStatus = true;
	        		++nonInertCount;
	        		reactiveSpecies += "(" + nonInertCount + ") " +
	        		tableInput.getValueAt(i,0) + " " +
	    			tableInput.getValueAt(i,1) + " (" +
	    			tableInput.getValueAt(i,2) + ") unreactive\r" +
	    			tableInput.getValueAt(i,4) + "\r\r";
	        	}
	        }
	        //	Add the non-inert gas species
	        conditionFile += "\rInitialStatus:\r\r" + reactiveSpecies + "END\r\r";
	        //	Add the inert gas species
	        conditionFile += "InertGas:\r" + inertSpecies + "END\r\r";
        }
        
        //	Add the spectroscopic data estimator
        conditionFile += "SpectroscopicDataEstimator: ";
        if (sdeCombo.getSelectedItem().equals("Frequency Groups")) {
        	conditionFile += "FrequencyGroups\r";
//        } else if (sdeCombo.getSelectedItem().equals("Three Frequency Model")) {
//        	conditionFile += "ThreeFrequencyModel\r";
        } else if (sdeCombo.getSelectedItem().equals("off")) {
        	conditionFile += "off\r";
        }
        
        //	Add the pressure dependence model
        conditionFile += "PressureDependence: ";
    	if (pdepCombo.getSelectedItem().equals("Modified Strong Collision")) {
    		conditionFile += "modifiedstrongcollision\r";
    	} else if (pdepCombo.getSelectedItem().equals("Reservoir State")) {
    		conditionFile += "reservoirstate\r";
    	} else if (pdepCombo.getSelectedItem().equals("off")) {
    		conditionFile += "off\r";
    	}
    	
    	// Add the pdep kinetics model
    	if (!pdepCombo.getSelectedItem().equals("off")) {
    		conditionFile += "PDepKineticsModel: ";
    		if (pdkmCombo.getSelectedItem().equals("CHEB")) {
    			conditionFile += "Chebyshev";
    			if (chebyTRep.getText().equals("") || chebyPRep.getText().equals("")) {
    				System.out.println("Number of chebyshev basis functions not supplied: Using default 4x4.");
    				conditionFile += "\r";
    			}
    			else {
    				conditionFile += " " + chebyTRep.getText() + " " + chebyPRep.getText() + "\r";
    			}
    		} else if (pdkmCombo.getSelectedItem().equals("PLOG"))
    			conditionFile += "PDepArrhenius\r";
    		else if (pdkmCombo.getSelectedItem().equals("Rate"))
    			conditionFile += "Rate\r";
    	}
    	
    	// Add the temperatures/pressure OR trange/prange
    	if (!chebyTMin.getText().equals("") && !chebyTMax.getText().equals("") &&
    			!chebyTGen.getText().equals("")) {
    		conditionFile += "TRange: (" + chebyTUnits.getSelectedItem() + 
				") " + chebyTMin.getText() + " " + chebyTMax.getText() + " " + 
				chebyTGen.getText();
    	}
    	else if (!chebyT.getText().equals("") && !chebyTGen.getText().equals("")) {
    		conditionFile += "Temperatures: " + chebyTGen.getText() + " (" +
    			chebyTUnits.getSelectedItem() + ")";
        	String[] tempLines = chebyT.getText().split("[\n]");
        	StringTokenizer st = null;
        	for (int i=0; i<tempLines.length; i++) {
	        	st = new StringTokenizer(tempLines[i]);
	        	while (st.hasMoreTokens()) {
	        		String line = st.nextToken();
	        		if (line.endsWith(",")) line = line.substring(0,line.length()-1).trim();
	        		conditionFile += " " + line;
	        	}
	        }
    	}
    	else
    		System.out.println("Temperature information for PDepKineticsModel not supplied:  Using default values.");
    	conditionFile += "\r";
    	
    	if (!chebyPMin.getText().equals("") && !chebyPMax.getText().equals("") &&
    			!chebyPGen.getText().equals("")) {
    		conditionFile += "PRange: (" + chebyPUnits.getSelectedItem() + 
				") " + chebyPMin.getText() + " " + chebyPMax.getText() + " " + 
				chebyPGen.getText();
    	}
    	else if (!chebyP.getText().equals("") && !chebyPGen.getText().equals("")) {
    		conditionFile += "Pressures: " + chebyPGen.getText() + " (" +
    			chebyPUnits.getSelectedItem() + ")";
        	String[] tempLines = chebyP.getText().split("[\n]");
        	StringTokenizer st = null;
        	for (int i=0; i<tempLines.length; i++) {
	        	st = new StringTokenizer(tempLines[i]);
	        	while (st.hasMoreTokens()) {
	        		String line = st.nextToken();
	        		if (line.endsWith(",")) line = line.substring(0,line.length()-1).trim();
	        		conditionFile += " " + line;
	        	}
	        }
    	}
    	else
    		System.out.println("Pressure information for PDepKineticsModel not supplied:  Using default values.");
    	conditionFile += "\r";
    	
    	// Add the "DecreaseGrainSize" field
    	if (grainsizeCombo.getSelectedItem().equals("Yes")) {
    		conditionFile += "DecreaseGrainSize: yes\r";
    	}
    	
    	conditionFile += "\r";
        //	Add the path of the IncludeSpecies.txt file (if present)
        if (unreactiveStatus) {
        	conditionFile += "IncludeSpecies: ";
        	if (isPath.getText().equals("")) {
            	System.out.println("ERROR: Writing condition.txt file: Could not read location of IncludeSpecies.txt (Initial Conditions tab)");
            	return;
        	}
        	else {
    	        if (isPath.getText().startsWith(rmgEnvVar)) {
    	        	int startIndex = rmgEnvVar.length()+1;
    	        	conditionFile += isPath.getText().substring(startIndex) + "\r\r";
    	        } else {
    	        	conditionFile += isPath.getText() + "\r\r";
    	        }
        	}
        }
        
        //	Add the finish controller
        conditionFile += "FinishController:\r" + "(1) Goal ";
        //	If goal is conversion
        boolean conversionGoal = true;
        if (controllerCombo.getSelectedItem().equals("Conversion")) {
        	if (SpeciesConvName.getSelectedItem().equals(" ") || conversion.getText().equals("")) {
            	System.out.println("ERROR: Writing condition.txt file: Could not read limiting reactant and/or fractional conversion (Termination tab)");
            	return;
        	}
        	else {
	        	conditionFile += controllerCombo.getSelectedItem() + " " +
	        		SpeciesConvName.getSelectedItem() + " " + conversion.getText() + "\r";
        	}
        } else if (controllerCombo.getSelectedItem().equals("ReactionTime")) {
        	if (time.getText().equals("")) {
            	System.out.println("ERROR: Writing condition.txt file: Could not read time of reaction (Termination tab)");
            	return;
        	}
        	else {
	        	conversionGoal = false;
	        	conditionFile += controllerCombo.getSelectedItem() + " " +
	        		time.getText() + " (" + timeCombo.getSelectedItem() + ")\r";
        	}
        }
        
        if (eTolerance.getText().equals("")) {
        	System.out.println("ERROR: Writing condition.txt file: Could not read model's error tolerance (Termination tab)");
        	return;
        } else {
        	conditionFile += "(2) Error Tolerance: " + eTolerance.getText() + "\r\r";
        }
        
        //	Add the dynamic simulator
        conditionFile += "DynamicSimulator: " + simulatorCombo.getSelectedItem() + "\r";
        if (conversionGoal) {
        	conditionFile += "Conversions:";
        	if (interCombo.getSelectedItem().equals("AUTO")) {
        		conditionFile += " AUTO";
        	}
        	else if (interCombo.getSelectedItem().equals("AUTOPRUNE")) {
        		conditionFile += " AUTOPRUNE\r";
        		conditionFile += "TerminationTolerance: " + textTerminTol.getText() + "\r";
        		conditionFile += "PruningTolerance: " + textPruneTol.getText() + "\r";
        		conditionFile += "MinSpeciesForPruning: " + textMinSpc.getText() + "\r";
        		conditionFile += "MaxEdgeSpeciesAfterPruning: " + textMaxSpc.getText();
        	}
        	else {
	            String[] convLines = textIndivSteps.getText().split("[\n]");
	            for (int i=0; i<convLines.length; i++) {
	            	StringTokenizer st = null;
	            	st = new StringTokenizer(convLines[i]);
	            	while (st.hasMoreTokens()) {
	            		String line = st.nextToken();
	            		if (line.endsWith(",")) line = line.substring(0,line.length()-1);
	            		conditionFile += " " + line;
	            	}
	            }
        	}   
        } else {
        	conditionFile += "TimeStep:";
        	if (interCombo.getSelectedItem().equals("AUTO")) {
        		conditionFile += " AUTO";
        	}
        	else if (interCombo.getSelectedItem().equals("AUTOPRUNE")) {
        		conditionFile += " AUTOPRUNE\r";
        		conditionFile += "TerminationTolerance: " + textTerminTol.getText() + "\r";
        		conditionFile += "PruningTolerance: " + textPruneTol.getText() + "\r";
        		conditionFile += "MinSpeciesForPruning: " + textMinSpc.getText() + "\r";
        		conditionFile += "MaxEdgeSpeciesAfterPruning: " + textMaxSpc.getText();
        	}
        	else {
	        	String[] tsLines = textIndivSteps.getText().split("[\n]");
	            for (int i=0; i<tsLines.length; i++) {
	            	StringTokenizer st = null;
	            	st = new StringTokenizer(tsLines[i]);
	            	while (st.hasMoreTokens()) {
	            		String line = st.nextToken();
	            		if (line.endsWith(",")) line = line.substring(0,line.length()-1);
	            		conditionFile += " " + line;
	            	}
	            }
        	}
        }
        
        if (aTolerance.getText().equals("") || rTolerance.getText().equals("")) {
        	System.out.println("ERROR: Writing condition.txt file: Could not read absolute and/or relative tolerance (Dynamic Simulator tab)");
        	return;
        }
        else {
        	conditionFile += "\rAtol: " + aTolerance.getText() + "\r" + "Rtol: " + rTolerance.getText() + "\r\r";
        }
    	
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
    			conditionFile += "on\rDisplay sensitivity information for:\r";
    			if (tableSens.getRowCount()==0) {
    	        	System.out.println("ERROR: Writing condition.txt file: Could not read table of species for which to perform sensitivity analysis (Sensitivity Analysis tab)");
    	        	return;
    			}
        		for (int j=0; j<tableSens.getRowCount(); j++) {
        			conditionFile += tableSens.getValueAt(j,0) + "\r";
        		}
        		conditionFile += "END\r\r";
    		} else if (sensCoeffCombo.getSelectedItem().equals("No")) {
    			conditionFile += "off\r";
    		}
    	}
    	
        //	Add the name(s)/location(s) of the primary kinetic library
        String prlReferenceDirectory_linux = rmgEnvVar + "/databases/" + stringMainDatabase + "/kinetics_libraries/";
        String prlReferenceDirectory_windows = rmgEnvVar + "\\databases\\" + stringMainDatabase + "\\kinetics_libraries\\";
    	conditionFile += "PrimaryKineticLibrary:\r";

		if (tablePRL.getRowCount()==0) {
        	System.out.println("Warning: Writing condition.txt file: Could not read Primary Kinetic Library (Thermochemical Libraries tab)");
		} else {
    		for (int k=0; k<tablePRL.getRowCount(); k++) {
    			conditionFile += "Name: " + tablePRL.getValueAt(k,0) + "\r" + "Location: ";
    			String prlDir = (String)tablePRL.getValueAt(k,1);
    	        if (prlDir.toLowerCase().startsWith(prlReferenceDirectory_linux.toLowerCase()) ||
    	        		prlDir.toLowerCase().startsWith(prlReferenceDirectory_windows.toLowerCase())) {
    	        	int startIndex = prlReferenceDirectory_linux.length();
    	        	conditionFile += prlDir.substring(startIndex) + "\r";
    	        } else {
    	        	conditionFile += prlDir + "\r";
    	        }
    		}
		}
		conditionFile += "END\r\r";
		
        //	Add the name(s)/location(s) of the primary kinetic library
    	conditionFile += "SeedMechanism:\r";

		if (tableSM.getRowCount()==0) {
        	System.out.println("Warning: Writing condition.txt file: Could not read Seed Mechanism (Thermochemical Libraries tab)");
		} else {
    		for (int k=0; k<tableSM.getRowCount(); k++) {
    			conditionFile += "Name: " + tableSM.getValueAt(k,0) + "\r" + "Location: ";
    			String smDir = (String)tableSM.getValueAt(k,2);
    	        if (smDir.toLowerCase().startsWith(prlReferenceDirectory_linux.toLowerCase()) ||
    	        		smDir.toLowerCase().startsWith(prlReferenceDirectory_windows.toLowerCase())) {
    	        	int startIndex = prlReferenceDirectory_linux.length();
    	        	conditionFile += smDir.substring(startIndex) + "\r";
    	        } else {
    	        	conditionFile += smDir + "\r";
    	        }
    	        conditionFile += "GenerateReactions: " + tableSM.getValueAt(k,1) + "\r";
    		}
		}
		conditionFile += "END\r\r";
		
		//	Add the Chemkin chem.inp file options
		conditionFile += "ChemkinUnits:\r";
		if (chemkinVerbosity.getSelectedItem().equals("Yes"))
			conditionFile += "Verbose: on\r";
		if (chemkinSMILES.getSelectedItem().equals("Yes"))
			conditionFile += "SMILES: on\r";
		conditionFile += "A: " + chemkinAUnits.getSelectedItem() + "\r";
		conditionFile += "Ea: " + chemkinEaUnits.getSelectedItem() + "\r";		

        File conditionPath = null;
		FileWriter fw = null;
		// Write the conditionFile string to user-specified file
        try {
        	
        	conditionPath = askUserForInput("Save file", false, null);
        	if (conditionPath != null) {
        		fw = new FileWriter(conditionPath);
        		fw.write(conditionFile);
            	fw.close();
        	}
        } catch (IOException e) {
        	String err = "Error writing condition file: ";
        	err += e.toString();
        	System.out.println(err);
        }

        // Run RMG if "Save and Run" selected
        if (saveAndRun && conditionPath != null) {
            runConditionFile(conditionPath.getAbsolutePath());
        }
	}
    
    public void openConditionFile(File filePath) {
    	int numSpecies = 0;
    	tmodelPTL.clear();
    	tmodelPRL.clear();
    	tmodelSM.clear();
    	tmodelInput.clear();
    	tmodelSens.clear();
    	tmodelFS.clear();
    	tmodelPTransL.clear();
        System.out.println("GUI cannot read in file's header (comments)");
    	
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
        String line = ChemParser.readMeaningfulLine(reader);
        read: while (line != null) {
	        // Populate the GUI with as much information as possible
	        
	        StringTokenizer st = null;
	        String tempString = "";
	        String[] tempStringVector = null;
	        
	        if (line.startsWith("Database")) {
		        //	Path of Database
		        st = new StringTokenizer(line);
		        tempString = st.nextToken();	// Skip over "Database:"
				tempString = st.nextToken();   // tempString now contains the database directory name
				String joinedPath = new File(System.getProperty("RMG.workingDirectory"), "databases").toString(); // find the folder containing databases
				joinedPath = new File(joinedPath, tempString).toString(); // find the path to the database specified
		        databasePath.setText(joinedPath);  // set the full path, not the shorthand (so that other fields can find it)
	        }
	        
	        else if (line.startsWith("PrimaryThermoLibrary")) {
		        //	Name(s)/Path(s) of PrimaryThermoLibrary
             	line = ChemParser.readMeaningfulLine(reader);
             	int ptlCounter = 0;
             	while (!line.equals("END")) {
             		++ptlCounter;
             		tempStringVector = line.split("Name: ");
             		String name = tempStringVector[tempStringVector.length-1].trim();
             		line = ChemParser.readMeaningfulLine(reader);
             		tempStringVector = line.split("Location: ");
             		String path = tempStringVector[tempStringVector.length-1].trim();
             		PRLVector ptlEntry = new PRLVector(ptlCounter-1,name,path);
					tmodelPTL.updatePRL(ptlEntry);
					line = ChemParser.readMeaningfulLine(reader);
             	}
	        }
	        
	        else if (line.startsWith("PrimaryTransportLibrary")) {
		        //	Name(s)/Path(s) of PrimaryTransportLibrary
             	line = ChemParser.readMeaningfulLine(reader);
             	int ptranslCounter = 0;
             	while (!line.equals("END")) {
             		++ptranslCounter;
             		tempStringVector = line.split("Name: ");
             		String name = tempStringVector[tempStringVector.length-1].trim();
             		line = ChemParser.readMeaningfulLine(reader);
             		tempStringVector = line.split("Location: ");
             		String path = tempStringVector[tempStringVector.length-1].trim();
             		PRLVector ptranslEntry = new PRLVector(ptranslCounter-1,name,path);
					tmodelPTransL.updatePRL(ptranslEntry);
					line = ChemParser.readMeaningfulLine(reader);
             	}
	        }
	        
	        else if (line.startsWith("TemperatureModel")) {
		        //	Temperature model
		        st = new StringTokenizer(line);
		        tempString = st.nextToken();	// Skip over "TemperatureModel:"
		        
		        String modelTemp = st.nextToken();
		        if (!modelTemp.equals("Constant")) {
		        	System.out.println("ERROR: Reading in condition.txt file - Invalid TemperatureModel.  RMG recognizes 'Constant' TemperatureModel, but not " + modelTemp + ".  GUI does not contain all information present in condition.txt file.");
		        }
		        tempModelCombo.setSelectedItem(modelTemp);
		        
		        String unitTemp = st.nextToken().substring(1,2);
		        if (unitTemp.equals("K") || unitTemp.equals("F") || unitTemp.equals("C")) {
		        	tempConstantUnit.setSelectedItem(unitTemp);
		        }
		        else {
		        	System.out.println("ERROR: Reading in condition.txt file - Invalid temperature unit.  RMG recognizes temperature units of 'K', 'C', or 'F', but not " + unitTemp + ".  GUI does not contain all information present in condition.txt file.");
		        }
		        
		        String temps = st.nextToken();
		        while (st.hasMoreTokens()) {
		        	temps += " " + st.nextToken();
		        }
		        tempConstant.setText(temps);
	        }
	        
	        else if (line.startsWith("PressureModel")) {
		        //	Pressure model
		        st = new StringTokenizer(line);
		        tempString = st.nextToken();	// Skip over "PressureModel:"
		        
		        String modelPress = st.nextToken();
		        if (!modelPress.equals("Constant")) {
		        	System.out.println("ERROR: Reading in condition.txt file - Invalid PressureModel.  RMG recognizes 'Constant' PressureModel, but not " + modelPress + ".  GUI does not contain all information present in condition.txt file.");
		        }
		        else {
		        	pressModelCombo.setSelectedItem(modelPress);
		        }
		        
		        String temp = st.nextToken();
		        String unitPress = temp.substring(1,temp.length()-1);
		        if (unitPress.equals("atm") || unitPress.equals("Pa") || unitPress.equals("bar") || unitPress.equals("torr")) {
		        	pressConstantUnit.setSelectedItem(unitPress);
		        }
		        else {
		        	System.out.println("ERROR: Reading in condition.txt file - Invalid pressure unit.  RMG recognizes pressure units of 'atm', 'Pa', 'bar' or 'torr', but not " + unitPress + ".  GUI does not contain all information present in condition.txt file.");
		        }
		        
		        String presses = st.nextToken();
		        while (st.hasMoreTokens()) {
		        	presses += " " + st.nextToken();
		        }
		        pressConstant.setText(presses);
	        }
	        
	        // GUI only recognized Ideal Gas EOS (30Jul2009)
//	        else if (line.startsWith("EquationOfState")) {
//		        //	Equation of state (if present)
//		        st = new StringTokenizer(line);
//		        tempString = st.nextToken();	// Skip over "EquationOfState"
//		        tempString = st.nextToken();
//		        if (tempString.equals("Liquid")) {
//		        	eosCombo.setSelectedItem(tempString);
//		        } else {
//		        	System.out.println("ERROR: Reading in condition.txt file - Invalid equation of state.  RMG recognizes 'Liquid', but not " + tempString + ".  GUI does not contain all information present in condition.txt file.");
//		        }
//	        }
	        
	        else if (line.startsWith("MaxCarbonNumber")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over "MaxCarbonNumberPerSpecies"\
	        	maxCarbonNum.setText(st.nextToken());
	        }
	        
	        else if (line.startsWith("MaxOxygenNumber")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over "MaxOxygenNumberPerSpecies"\
	        	maxOxygenNum.setText(st.nextToken());
	        }
	        
	        else if (line.startsWith("MaxRadicalNumber")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over "MaxRadicalNumberPerSpecies"\
	        	maxRadicalNum.setText(st.nextToken());
	        }
	        
	        else if (line.startsWith("MaxSulfurNumber")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over "MaxSulfurNumberPerSpecies"\
	        	maxSulfurNum.setText(st.nextToken());
	        }
	        
	        else if (line.startsWith("MaxSiliconNumber")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over "MaxSiliconNumberPerSpecies"\
	        	maxSiliconNum.setText(st.nextToken());
	        }
	        
	        else if (line.startsWith("MaxHeavyAtom")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over "MaxHeavyAtomPerSpecies"\
	        	maxHeavyAtom.setText(st.nextToken());
	        }
	        
	        else if (line.startsWith("InChIGeneration")) {
		        //	InChI generation
		        st = new StringTokenizer(line);
		        tempString = st.nextToken();	// Skip over "InChIGeneration:"
		        String inchi = st.nextToken();
		        if (inchi.toLowerCase().equals("on")) {
		        	inchiCombo.setSelectedItem("Yes");
		        } else if (inchi.toLowerCase().equals("off")) {
		        	inchiCombo.setSelectedItem("No");
		        } else {
		        	System.out.println("ERROR: Reading in condition.txt file - Invalid argument for InChIGeneration.  RMG recognizes an argument of 'On' or 'Off', but not " + inchi + ".  GUI does not contain all information present in condition.txt file.");
		        }
	        }
	        
	        else if (line.toLowerCase().startsWith("readrestart")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over "ReadRestart:"
	        	String rRestart = st.nextToken();
	        	if (rRestart.toLowerCase().equals("yes")) {
		        	readRestartOnOff.setSelectedItem("Yes");
		        } else if (rRestart.toLowerCase().equals("no")) {
		        	readRestartOnOff.setSelectedItem("No");
		        } else {
		        	System.out.println("ERROR: Reading in condition.txt file - Invalid argument for ReadRestart.  RMG recognizes an argument of 'Yes' or 'No', but not " + rRestart + ".  GUI does not contain all information present in condition.txt file.");
		        }
	        }
	        
	        else if (line.toLowerCase().startsWith("writerestart")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over "WriteRestart:"
	        	String wRestart = st.nextToken();
	        	if (wRestart.toLowerCase().equals("yes")) {
		        	writeRestartOnOff.setSelectedItem("Yes");
		        } else if (wRestart.toLowerCase().equals("no")) {
		        	writeRestartOnOff.setSelectedItem("No");
		        } else {
		        	System.out.println("ERROR: Reading in condition.txt file - Invalid argument for WriteRestart.  RMG recognizes an argument of 'Yes' or 'No', but not " + wRestart + ".  GUI does not contain all information present in condition.txt file.");
		        }
	        }
	        
	        // GUI does not handle Solvation (30Jul2009)
//	        else if (line.startsWith("Solvation")) {
//		        //	Solvation
//		        st = new StringTokenizer(line);
//		        tempString = st.nextToken();	// Skip over "Solvation:"
//		        String solvation = st.nextToken();
//		        if (solvation.toLowerCase().equals("on")) {
//		        	solvationCombo.setSelectedItem("Yes");
//		        } else if (solvation.toLowerCase().equals("off")) {
//		        	solvationCombo.setSelectedItem("No");
//		        } else {
//		        	System.out.println("ERROR: Reading in condition.txt file - Invalid argument for Solvation.  RMG recognizes an argument of 'On' or 'Off', but not " + solvation + ".  GUI does not contain all information present in condition.txt file.");
//		        }
//	        }
	        
	        else if (line.startsWith("InitialStatus")) {
		        String adjList = "";
		        String name = "";
		        String conc = "";
		        String unit = "";
		        String reactivity = "";
		        line = ChemParser.readMeaningfulLine(reader); 
		        while (!line.startsWith("END")) {
		        	if (line.startsWith("(")) {
		        		++numSpecies;
		        		st = new StringTokenizer(line);
		        		tempString = st.nextToken();	// Skip over (#)
		        		name = st.nextToken();
		        		conc = st.nextToken();
		        		unit = st.nextToken();
		        		reactivity = "Reactive-RMG";
		        		if (st.hasMoreTokens()) {
		        			reactivity = "Reactive-User";
		        		}
		        	} else {
		        		adjList += line + "\r";
		        	}
		        	
		        	line = ChemParser.readMeaningfulLine(reader);

		        	if (line.startsWith("(") || line.startsWith("END")) {
		        		ICVector icEntry = new ICVector(numSpecies-1,name,conc,unit.substring(1,unit.length()-1),reactivity,adjList);
		    			tmodelInput.updateIC(icEntry);
		        		adjList = "";
		        	}
		        }
		        limitingReactant();
	        }
	        
	        else if (line.startsWith("InertGas")) {		        
		        //	Inert gases
		        String name = "";
		        String conc = "";
		        String unit = "";
		        String reactivity = "";
		        String adjList = "";
	        	line = ChemParser.readMeaningfulLine(reader);
		        while (!line.startsWith("END")) {
		    		++numSpecies;
		    		st = new StringTokenizer(line);
		    		name = st.nextToken();
		    		conc = st.nextToken();
		    		unit = st.nextToken();
		    		reactivity = "Inert";
		    		if (name.equals("Ar") || name.equals("Ne") || name.equals("N2") || name.equals("He")) {
		    			ICVector icEntry = new ICVector(numSpecies-1,name,conc,unit.substring(1,unit.length()-1),reactivity,adjList);
		    			tmodelInput.updateIC(icEntry);
		    		} else
		    			System.out.println("Warning: Could not read in one of the inert gas species: " + name + ".  The Inert Gas options are 'N2', 'Ar', 'He', and 'Ne' (case sensitive)");
		    		line = ChemParser.readMeaningfulLine(reader);
		        }
	        }
	        
	        else if (line.startsWith("Spectroscopic")) {
		        //	Spectroscopic data estimator
		        st = new StringTokenizer(line);
		        tempString = st.nextToken();	// Skip over SpectroscopicDataEstimator
		        String sde = st.nextToken();
		        if (sde.toLowerCase().equals("frequencygroups")) {
		        	sdeCombo.setSelectedItem("Frequency Groups");
//		        } else if (sde.toLowerCase().equals("threefrequencymodel")) {
//		        	sdeCombo.setSelectedItem("Three Frequency Model");
		        } else if (sde.toLowerCase().equals("off")) {
		        	sdeCombo.setSelectedItem("off");
		        } else {
		        	System.out.println("ERROR: Reading in condition.txt file - Invalid argument for SpectroscopicDataEstimator.  RMG recognizes an argument of 'frequencygroups' or 'off' but not " + sde + ".  GUI does not contain all information present in condition.txt file.");
		        }
	        }
	        
	        else if (line.startsWith("PressureDep")) {	        
		        //	Pressure dependence model
		        st = new StringTokenizer(line);
		        tempString = st.nextToken();	// Skip over PressureDependence
		        String pdep = st.nextToken();
		        if (pdep.toLowerCase().equals("off")) {
		        	pdepCombo.setSelectedItem("off");
		        } else if (pdep.toLowerCase().equals("reservoirstate")) {
		        	pdepCombo.setSelectedItem("Reservoir State");
		        } else if (pdep.toLowerCase().equals("modifiedstrongcollision")) {
		        	pdepCombo.setSelectedItem("Modified Strong Collision");
		        } else {
		        	System.out.println("ERROR: Reading in condition.txt file - Invalid argument for PressureDependence.  RMG recognizes an argument of 'reservoirstate', 'modifiedstrongcollision', or 'off' but not " + pdep + ".  GUI does not contain all information present in condition.txt file.");
		        }
	        }
	        
	        else if (line.startsWith("PDepKinetics")) {
	        	//	Pressure dependent kinetics model
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over PDepKineticsModel
	        	String pdkm = st.nextToken();
		        if (pdkm.toLowerCase().equals("chebyshev")) {
		        	pdkmCombo.setSelectedItem("CHEB");
		        	if (st.hasMoreTokens()) {
		        		chebyTRep.setText(st.nextToken());
		        		chebyPRep.setText(st.nextToken());
		        	}
		        } else if (pdkm.toLowerCase().equals("pdeparrhenius")) {
		        	pdkmCombo.setSelectedItem("PLOG");
		        } else if (pdkm.toLowerCase().equals("rate")) {
		        	pdkmCombo.setSelectedItem("Rate");
		        } else {
		        	System.out.println("ERROR: Reading in condition.txt file - Invalid argument for PDepKineticsModel.  RMG recognizes an argument of 'chebyshev' or 'pdeparrhenius' but not " + pdkm + ".  GUI does not contain all information present in condition.txt file.");
		        }
	        }
	        
	        else if (line.startsWith("TRange")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over TRange
	        	String units = st.nextToken();
	        	String unit = units.substring(1,units.length()-1);
	        	if (unit.equals("K") || unit.equals("F") || unit.equals("C"))
	        		chebyTUnits.setSelectedItem(unit);
	        	else
	        		System.out.println("ERROR: Reading in condition.txt file - Invalid chebyshev temperature unit.  RMG recognizes temperature units of 'K', 'C', or 'F', but not " + unit + ".  GUI does not contain all information present in condition.txt file.");
	        	chebyTMin.setText(st.nextToken());
	        	chebyTMax.setText(st.nextToken());
	        	chebyTGen.setText(st.nextToken());
	        	chebyTCombo.setSelectedItem("Max/Min");
	        }
	        
	        else if (line.startsWith("PRange")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over PRange
	        	String units = st.nextToken();
	        	String unit = units.substring(1,units.length()-1);
	        	if (unit.equals("atm") || unit.equals("Pa") || unit.equals("bar") || unit.equals("torr"))
	        		chebyPUnits.setSelectedItem(unit);
	        	else
	        		System.out.println("ERROR: Reading in condition.txt file - Invalid chebyshev pressure unit.  RMG recognizes pressure units of 'atm', 'bar', 'torr', or 'Pa', but not " + unit + ".  GUI does not contain all information present in condition.txt file.");
	        	chebyPMin.setText(st.nextToken());
	        	chebyPMax.setText(st.nextToken());
	        	chebyPGen.setText(st.nextToken());
	        	chebyPCombo.setSelectedItem("Max/Min");
	        }
	        
	        else if (line.startsWith("Temperatures")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over Temperatures:
	        	chebyTGen.setText(st.nextToken());
	        	String units = st.nextToken();
	        	String unit = units.substring(1,units.length()-1);
	        	if (unit.equals("K") || unit.equals("F") || unit.equals("C"))
	        		chebyTUnits.setSelectedItem(unit);
	        	else
	        		System.out.println("ERROR: Reading in condition.txt file - Invalid chebyshev temperature unit.  RMG recognizes temperature units of 'K', 'C', or 'F', but not " + unit + ".  GUI does not contain all information present in condition.txt file.");
	        	String Tvalues = "";
	        	while (st.hasMoreTokens()) {
	        		Tvalues += st.nextToken() + " ";
	        	}
	        	chebyT.setText(Tvalues);
	        	chebyTCombo.setSelectedItem("List");
	        }
	        
	        else if (line.startsWith("Pressures")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over Pressures:
	        	chebyPGen.setText(st.nextToken());
	        	String units = st.nextToken();
	        	String unit = units.substring(1,units.length()-1);
	        	if (unit.equals("atm") || unit.equals("Pa") || unit.equals("bar") || unit.equals("torr"))
	        		chebyPUnits.setSelectedItem(unit);
	        	else
	        		System.out.println("ERROR: Reading in condition.txt file - Invalid chebyshev pressure unit.  RMG recognizes pressure units of 'atm', 'bar', 'torr', or 'Pa', but not " + unit + ".  GUI does not contain all information present in condition.txt file.");
	        	String Pvalues = "";
	        	while (st.hasMoreTokens()) {
	        		Pvalues += st.nextToken() + " ";
	        	}
	        	chebyP.setText(Pvalues);
	        	chebyPCombo.setSelectedItem("List");
	        }
	        
	        else if (line.toLowerCase().startsWith("decreasegrainsize")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over "DecreaseGrainSize:"
	        	String yesno = st.nextToken();
	        	if (yesno.toLowerCase().equals("yes")) {
		        	grainsizeCombo.setSelectedItem("Yes");
		        } else if (yesno.toLowerCase().equals("no")) {
		        	grainsizeCombo.setSelectedItem("No");
		        } else {
		        	System.out.println("ERROR: Reading in condition.txt file - Invalid argument for DecreaseGrainSize.  RMG recognizes an argument of 'Yes' or 'No', but not " + yesno + ".  GUI does not contain all information present in condition.txt file.");
		        }
	        }
	        
	        else if (line.startsWith("IncludeSpecies")) {
		        //	IncludeSpecies.txt file (if present)
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over IncludeSpecies
	        	isPath.setText(st.nextToken());
	        }
	        
	        else if (line.startsWith("(1) Goal")) {
		        //	Goal
		        st = new StringTokenizer(line);
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
		        } else {
		        	System.out.println("ERROR: Reading in condition.txt file - Invalid argument for Goal.  RMG recognizes an argument of 'Conversion' or 'ReactionTime' but not " + goal + ".  GUI does not contain all information present in condition.txt file.");
		        }
	        }
	        
	        else if (line.startsWith("(2) Error")) {
		        //	Error tolerance
		        st = new StringTokenizer(line);
		        tempString = st.nextToken();	// Skip over "(2)"
		        tempString = st.nextToken();	// Skip over "Error"
		        tempString = st.nextToken();	// Skip over "Tolerance:"
		        eTolerance.setText(st.nextToken());
	        }
	        
	        else if (line.startsWith("DynamicSim")) {
	        //	Dynamic simulator
		        st = new StringTokenizer(line);
		        tempString = st.nextToken();
		        String sim = st.nextToken();
		        if (sim.equals("DASPK") || sim.equals("DASSL")) {
		        	simulatorCombo.setSelectedItem(sim);
		        } else {
		        	System.out.println("ERROR: Reading in condition.txt file - Invalid argument for DynamicSimulator.  RMG recognizes an argument of 'DASPK' or 'DASSL' but not " + sim + ".  GUI does not contain all information present in condition.txt file.");
		        }
	        }
	        
	        else if (line.startsWith("Conversions")) {
		        //	Intermediate integration step
		        st = new StringTokenizer(line);
		        String stepName = st.nextToken();	// Skip over Conversions
		        String steps = st.nextToken();
		        if (steps.toLowerCase().equals("auto")) {
		        	interCombo.setSelectedItem("AUTO");
		        }
		        else if (steps.toLowerCase().equals("autoprune")) {
		        	interCombo.setSelectedItem("AUTOPRUNE");
		        }
		        else {
			        while (st.hasMoreTokens()) {
			        	steps += " " + st.nextToken();
			        }
			        textIndivSteps.setText(steps);
			        interCombo.setSelectedIndex(0);
		        }
	        }
	        
	        else if (line.startsWith("TimeStep")) {
		        //	Intermediate integration step
		        st = new StringTokenizer(line);
		        String stepName = st.nextToken();	// Skip over TimeStep
		        String steps = st.nextToken();
		        if (steps.toLowerCase().equals("auto")) {
		        	interCombo.setSelectedItem("AUTO");
		        }
		        else if (steps.toLowerCase().equals("autoprune")) {
		        	interCombo.setSelectedItem("AUTOPRUNE");
		        }
		        else {
			        while (st.hasMoreTokens()) {
			        	steps += " " + st.nextToken();
			        }
			        textIndivSteps.setText(steps);
			        interCombo.setSelectedIndex(0);
		        }
	        }

	        /*
	         * Pruning options
	         */
	        else if (line.toLowerCase().startsWith("terminationtolerance")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over TerminationTolerance:
	        	textTerminTol.setText(st.nextToken());
	        }
	        
	        else if (line.toLowerCase().startsWith("pruningtolerance")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over PruningTolerance:
	        	textPruneTol.setText(st.nextToken());
	        }
	        
	        else if (line.toLowerCase().startsWith("minspeciesforpruning")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over MinSpeciesForPruning:
	        	textMinSpc.setText(st.nextToken());
	        }
	        
	        else if (line.toLowerCase().startsWith("maxedgespeciesafterpruning")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over MaxEdgeSpeciesAfterPruning:
	        	textMaxSpc.setText(st.nextToken());
	        }
	        	        
	        else if (line.startsWith("Atol")) {
		        //	Absolute tolerance
		        st = new StringTokenizer(line);
		        tempString = st.nextToken();	// Skip over Atol
		        aTolerance.setText(st.nextToken());
	        }
	        
	        else if (line.startsWith("Rtol")) {
		        //	Relative tolerance
		        st = new StringTokenizer(line);
		        tempString = st.nextToken();	// Skip over Rtol
		        rTolerance.setText(st.nextToken());
	        }
	        
	        else if (line.startsWith("Error bars")) {
		        //	Error bars and Sensitivity analysis (if present)
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over "Error"
	        	tempString = st.nextToken();	// Skip over "bars:"
	        	String onoff = st.nextToken();
	        	if (onoff.equals("on")) {
	        		eBarsCombo.setSelectedItem("Yes");
	        	} else if (onoff.equals("off")) {
	        		eBarsCombo.setSelectedItem("No");
	        	} else {
	        		System.out.println("ERROR: Reading in condition.txt file - Invalid argument for Error bars.  RMG recognizes an argument of 'Off' or 'On' but not " + onoff + ".  GUI does not contain all information present in condition.txt file.");
	        	}
	        }
	        
	        else if (line.startsWith("Display sensitivity coefficients")) {
	            st = new StringTokenizer(line);
	            tempString = st.nextToken();	// Skip over "Display"
	            tempString = st.nextToken();	// Skip over "sensitivity"
	            tempString = st.nextToken();	// Skip over "coefficients:"
	            String onoff = st.nextToken();
	            if (onoff.equals("on")) {
	        		sensCoeffCombo.setSelectedItem("Yes");
	        	} else if (onoff.equals("off")) {
	        		sensCoeffCombo.setSelectedItem("No");
	        	} else {
	        		System.out.println("ERROR: Reading in condition.txt file - Invalid argument for Display sensitivity coefficients.  RMG recognizes an argument of 'Off' or 'On' but not " + onoff + ".  GUI does not contain all information present in condition.txt file.");
	        	}     	
	        }
	        
	        else if (line.startsWith("Display sensitivity information")) {
		        //	Sensitivity coefficients (if present)
	        	line = ChemParser.readMeaningfulLine(reader);
	        	numSpecies = 0;
	        	while (!line.startsWith("END")) {
	        		++numSpecies;
	        		SensVector sensEntry = new SensVector(numSpecies-1,line);
	    			tmodelSens.updateSens(sensEntry);
	    			line = ChemParser.readMeaningfulLine(reader);
	        	}
	        }
	        
	        else if (line.startsWith("PrimaryKineticLibrary")) {
		        
             	line = ChemParser.readMeaningfulLine(reader);
             	int prlCounter = 0;
             	while (!line.equals("END")) {
             		++prlCounter;
             		tempStringVector = line.split("Name: ");
             		String name = tempStringVector[tempStringVector.length-1].trim();
             		line = ChemParser.readMeaningfulLine(reader);
             		tempStringVector = line.split("Location: ");
             		String path = tempStringVector[tempStringVector.length-1].trim();
             		PRLVector prlEntry = new PRLVector(prlCounter-1,name,path);
					tmodelPRL.updatePRL(prlEntry);
					line = ChemParser.readMeaningfulLine(reader);
             	}
	        }
	        
	        else if (line.startsWith("SeedMechanism")) {
		        //	Name(s)/Path(s) of SeedMechanism
             	line = ChemParser.readMeaningfulLine(reader);
             	int smCounter = 0;
             	while (!line.equals("END")) {
             		++smCounter;
             		tempStringVector = line.split("Name: ");
             		String name = tempStringVector[tempStringVector.length-1].trim();
             		line = ChemParser.readMeaningfulLine(reader);
             		tempStringVector = line.split("Location: ");
             		String path = tempStringVector[tempStringVector.length-1].trim();
             		line = ChemParser.readMeaningfulLine(reader);
             		tempStringVector = line.split("GenerateReactions: ");
             		String react = tempStringVector[tempStringVector.length-1].trim();
             		SMVector smEntry = new SMVector(smCounter-1,name,react,path);
					// need to read in "GenerateReactions:" line
					tmodelSM.updateSM(smEntry);
					line = ChemParser.readMeaningfulLine(reader);
             	}
	        }
	        
	        else if (line.startsWith("Verbose")) {
        		st = new StringTokenizer(line);
        		tempString = st.nextToken();	// Skip over Verbose
	            String onoff = st.nextToken();
	            if (onoff.equals("on")) {
	        		chemkinVerbosity.setSelectedItem("Yes");
	        	} else if (onoff.equals("off")) {
	        		chemkinVerbosity.setSelectedItem("No");
	        	} else {
	        		System.out.println("ERROR: Reading in condition.txt file - Invalid argument for Chemkin Verbosity.  RMG recognizes an argument of 'Off' or 'On' but not " + onoff + ".  GUI does not contain all information present in condition.txt file.");
	        	}
	        }
	        
	        else if (line.toLowerCase().startsWith("smiles")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over SMILES
	        	String onoff = st.nextToken();
	            if (onoff.equals("on")) {
	        		chemkinSMILES.setSelectedItem("Yes");
	        	} else if (onoff.equals("off")) {
	        		chemkinSMILES.setSelectedItem("No");
	        	} else {
	        		System.out.println("ERROR: Reading in condition.txt file - Invalid argument for Chemkin SMILES.  RMG recognizes an argument of 'Off' or 'On' but not " + onoff + ".  GUI does not contain all information present in condition.txt file.");
	        	}
	        }
	        
	        else if (line.startsWith("A:")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over A
	        	String units = st.nextToken();
	        	if (units.equals("moles") || units.equals("molecules"))
	        		chemkinAUnits.setSelectedItem(units);
	        	else
	        		System.out.println("ERROR: Reading in condition.txt file - Invalid argument for Chemkin A units.  RMG recognizes an argument of 'moles' or 'molecules' but not " + units + ".  GUI does not contain all information present in condition.txt file.");		
	        }
	        
	        else if (line.startsWith("Ea:")) {
	        	st = new StringTokenizer(line);
	        	tempString = st.nextToken();	// Skip over Ea
	        	String units = st.nextToken();
	        	if (units.equals("cal/mol") || units.equals("kcal/mol") || units.equals("kJ/mol") ||
	        			units.equals("J/mol") || units.equals("Kelvins"))
	        		chemkinEaUnits.setSelectedItem(units);
	        	else
	        		System.out.println("ERROR: Reading in condition.txt file - Invalid argument for Chemkin Ea units.  RMG recognizes an argument of 'kcal/mol', 'cal/mol', 'kJ/mol', 'J/mol', or 'Kelvins' but not " + units + ".  GUI does not contain all information present in condition.txt file.");		
	        }
	        
	        else if (line.startsWith("ForbiddenStruct")) {
		        String adjList = "";
		        String name = "";
		        int numForbidStructs = 0;
		        line = ChemParser.readMeaningfulLine(reader);
		        while (!line.startsWith("END")) {
        			try {
        				int doesLineBeginWithNumber = Integer.parseInt(line.substring(0,1));
        				adjList += line + "\r";
        			} catch (NumberFormatException e) {
        				name = line;
        				++numForbidStructs;
        			}
		        	
		        	line = ChemParser.readMeaningfulLine(reader);

		        	if (line.startsWith("END")) {
		        		FSVector fsEntry = new FSVector(numForbidStructs,name,adjList);
		    			tmodelFS.updateFS(fsEntry);
		        		adjList = "";
		        	}
		        	else {
			        	try {
				        	Integer.parseInt(line.substring(0,1));
			        	} catch (NumberFormatException e) {
			        		FSVector fsEntry = new FSVector(numForbidStructs,name,adjList);
			    			tmodelFS.updateFS(fsEntry);
			        		adjList = "";
			        	}
		        	}
		        }
	        }
	        
	        else
	        	if (!line.startsWith("FinishController") & !line.startsWith("ChemkinUnits"))
	        		System.out.println("Warning: Reading in condition.txt file - Interpreter does not recognize the following line: " + line);

	    	line = ChemParser.readMeaningfulLine(reader);
        }
    }
    
    public JComponent createTabOptions() {
    	//	Create cellrenderers for JComboBox/JTable - CENTER text
    	ListCellRenderer centerJComboBoxRenderer = new CenterJCBRenderer();
        TableCellRenderer centerRenderer = new CenterTRenderer();
    	
        //	Create "Options" panel
        JPanel Options = new JPanel();
    	Options.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Additional Options"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	
    	//	Create the InChI generation subpanel
    	JPanel inchiPanel = new JPanel();
    	inchiPanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("InChI Generation"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	//	Populate the subpanel
    	JLabel inchiLabel = new JLabel("Generate InChIs");
    	inchiPanel.add(inchiLabel);
    	inchiPanel.add(inchiCombo = new JComboBox(yesnoOptions));
    	
    	//	Create the Max carbon/oxygen/radical subpanel
    	JPanel maxPanel = new JPanel();
    	maxPanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Species properties"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	//	Populate the subpanel
    	JPanel carbonPanel = new JPanel();
    	carbonPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	JLabel maxCarbonLabel = new JLabel("Maximum Carbon # per species");
    	carbonPanel.add(maxCarbonLabel);
    	carbonPanel.add(maxCarbonNum = new JTextField());
    	maxCarbonNum.setPreferredSize(new Dimension(50,25));
    	maxCarbonNum.addActionListener(this);
    	maxCarbonNum.setHorizontalAlignment(JTextField.CENTER);
    	    	
    	JPanel oxygenPanel = new JPanel();
    	oxygenPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	JLabel maxOxygenLabel = new JLabel("Maximum Oxygen # per species");
    	oxygenPanel.add(maxOxygenLabel);
    	oxygenPanel.add(maxOxygenNum = new JTextField());
    	maxOxygenNum.setPreferredSize(new Dimension(50,25));
    	maxOxygenNum.addActionListener(this);
    	maxOxygenNum.setHorizontalAlignment(JTextField.CENTER);
    	
    	JPanel radicalPanel = new JPanel();
    	radicalPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	JLabel maxRadicalLabel = new JLabel("Maximum Radical # per species");
    	radicalPanel.add(maxRadicalLabel);
    	radicalPanel.add(maxRadicalNum = new JTextField());
    	maxRadicalNum.setPreferredSize(new Dimension(50,25));
    	maxRadicalNum.addActionListener(this);
    	maxRadicalNum.setHorizontalAlignment(JTextField.CENTER);
    	
    	JPanel sulfurPanel = new JPanel();
    	sulfurPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	JLabel maxSulfurLabel = new JLabel("Maximum Sulfur # per species");
    	sulfurPanel.add(maxSulfurLabel);
    	sulfurPanel.add(maxSulfurNum = new JTextField());
    	maxSulfurNum.setPreferredSize(new Dimension(50,25));
    	maxSulfurNum.addActionListener(this);
    	maxSulfurNum.setHorizontalAlignment(JTextField.CENTER);
    	
    	JPanel siliconPanel = new JPanel();
    	siliconPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	JLabel maxSiliconLabel = new JLabel("Maximum Silicon # per species");
    	siliconPanel.add(maxSiliconLabel);
    	siliconPanel.add(maxSiliconNum = new JTextField());
    	maxSiliconNum.setPreferredSize(new Dimension(50,25));
    	maxSiliconNum.addActionListener(this);
    	maxSiliconNum.setHorizontalAlignment(JTextField.CENTER);
    	
    	JPanel heavyatomPanel = new JPanel();
    	heavyatomPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	JLabel maxheavyatomLabel = new JLabel("Maximum Heavy Atoms per species");
    	heavyatomPanel.add(maxheavyatomLabel);
    	heavyatomPanel.add(maxHeavyAtom = new JTextField());
    	maxHeavyAtom.setPreferredSize(new Dimension(50,25));
    	maxHeavyAtom.addActionListener(this);
    	maxHeavyAtom.setHorizontalAlignment(JTextField.CENTER);
    	
    	Box maxSpeciesProps1 = Box.createVerticalBox();
    	maxSpeciesProps1.add(carbonPanel);
    	maxSpeciesProps1.add(oxygenPanel);
    	maxSpeciesProps1.add(radicalPanel);
    	
    	Box maxSpeciesProps2 = Box.createVerticalBox();
    	maxSpeciesProps2.add(sulfurPanel);
    	maxSpeciesProps2.add(siliconPanel);
    	maxSpeciesProps2.add(heavyatomPanel);
    	
    	Box maxSpeciesProps = Box.createHorizontalBox();
    	maxSpeciesProps.add(maxSpeciesProps1);
    	maxSpeciesProps.add(maxSpeciesProps2);
    	
    	maxPanel.add(maxSpeciesProps);
    	
    	//	Create the restart options subpanel
    	JPanel restartPanel = new JPanel();
    	restartPanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Restart"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	//	Populate the subpanel
    	JPanel readRestartPanel = new JPanel();
    	readRestartPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	JLabel readRestartLabel = new JLabel("Read restart files?");
    	readRestartPanel.add(readRestartLabel);
    	readRestartPanel.add(readRestartOnOff = new JComboBox(yesnoOptions));
    	
    	JPanel writeRestartPanel = new JPanel();
    	writeRestartPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	JLabel writeRestartLabel = new JLabel("Write restart files?");
    	writeRestartPanel.add(writeRestartLabel);
    	writeRestartPanel.add(writeRestartOnOff = new JComboBox(yesnoOptions));
    	
    	Box restartProps = Box.createHorizontalBox();
    	restartProps.add(readRestartPanel);
    	restartProps.add(writeRestartPanel);
    	
    	restartPanel.add(restartProps);
    	
    	//	Create the chemkin options subpanel
    	JPanel chemkinPanel = new JPanel();
    	chemkinPanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("CHEMKIN chem.inp properties"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	//	Populate the subpanel
    	JPanel AUnitsPanel = new JPanel();
    	AUnitsPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	JLabel ALabel = new JLabel("Units of A in chem.inp file");
    	AUnitsPanel.add(ALabel);
    	AUnitsPanel.add(chemkinAUnits = new JComboBox(AUnitsOptions));
    	
    	JPanel EaUnitsPanel = new JPanel();
    	EaUnitsPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	JLabel EaLabel = new JLabel("Units of Ea in chem.inp file");
    	EaUnitsPanel.add(EaLabel);
    	EaUnitsPanel.add(chemkinEaUnits = new JComboBox(EaUnitsOptions));
    	
    	JPanel verbosePanel = new JPanel();
    	verbosePanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	JLabel verboseLabel = new JLabel("Detailed kinetics comments");
    	verbosePanel.add(verboseLabel);
    	verbosePanel.add(chemkinVerbosity = new JComboBox(yesnoOptions));
    	
    	JPanel smilesPanel = new JPanel();
    	smilesPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	JLabel smilesLabel = new JLabel("Write InChI in chem.inp");
    	smilesPanel.add(smilesLabel);
    	smilesPanel.add(chemkinSMILES = new JComboBox(yesnoOptions));
    	
    	Box chemkinProps = Box.createVerticalBox();
    	chemkinProps.add(AUnitsPanel);
    	chemkinProps.add(EaUnitsPanel);
    	chemkinProps.add(verbosePanel);
    	chemkinProps.add(smilesPanel);
    	
    	chemkinPanel.add(chemkinProps);
    	
    	//	Create the Equation Of State (EOS) subpanel
    	JPanel EOS = new JPanel();
    	EOS.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Equation of State"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));
    	//	Populate the subpanel
    	JLabel eosLabel = new JLabel("Equation of State");
    	EOS.add(eosLabel);
    	EOS.add(eosCombo = new JComboBox(eosOptions));
    	
    	//	Create the "forbiddenStructures" panel
    	//		This subpanel holds the add/remove buttons
    	JPanel forbiddenStructures = new JPanel();
    	forbiddenStructures.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Forbidden Structures"),
    			BorderFactory.createEmptyBorder(5,5,5,5)));

        //	Create the "speciesAList" subpanel
        JPanel fsAdjList = new JPanel();
        
        //	Create text field where the user can visualize the adjacency list
        FSAdjList = new JTextArea(5,13);
        //	Add text field to a scroll panel
        JScrollPane scrollfsAdjList = new JScrollPane(FSAdjList);
    	scrollfsAdjList.setBorder(BorderFactory.createLoweredBevelBorder());
    	
    	Box fsadj = Box.createVerticalBox();
    	fsadj.add(new JLabel("           Name"));
    	fsadj.add(FSName = new JTextField());
    	fsadj.add(new JLabel("        Structure"));
    	fsadj.add(scrollfsAdjList);
    	
    	fsAdjList.add(fsadj);
    	
    	//	Create the buttons: "Add" and "Remove"
    	JButton AddFS = new JButton("Add");
    	AddButtonListener addListenerFS = new AddButtonListener();
    	AddFS.addActionListener(addListenerFS);
    	AddFS.setActionCommand("AddFS");
    	AddFS.setToolTipText("Press to submit forbidden structure");
    	
    	JButton DeleteFS = new JButton("Remove");
    	DeleteButtonListener deleteListenerFS = new DeleteButtonListener();
    	DeleteFS.addActionListener(deleteListenerFS);
    	DeleteFS.setActionCommand("DeleteFS");
    	DeleteFS.setToolTipText("Press to delete forbidden structure");
    	
    	//	Add the box to the "fsButton" subpanel
    	JPanel fsButton = new JPanel();
    	fsButton.add(AddFS);
    	fsButton.add(DeleteFS);

    	//	Create the "fsTable" subpanel
    	JPanel fsTable = new JPanel();
    	
    	//	Create table to hold the species data
    	tableFS = new JTable(tmodelFS = new MyTableModelFS());
    	tableFS.setPreferredScrollableViewportSize(new Dimension(300,50));
    	tableFS.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        tableFS.getColumnModel().getColumn(0).setPreferredWidth(100);
        tableFS.getColumnModel().getColumn(1).setPreferredWidth(200);
        
        //	Create a scroll panel and add it to "speciesTable" subpanel
    	JScrollPane scrollFS = new JScrollPane(tableFS);
        scrollFS.setBorder(BorderFactory.createLoweredBevelBorder());
        for (int i=0; i<tableFS.getColumnCount(); i++) {
        	TableColumn column = tableFS.getColumnModel().getColumn(i);
        	column.setCellRenderer(centerRenderer);
        }
        fsTable.add(scrollFS);
        
        Box fsBox = Box.createVerticalBox();
        fsBox.add(fsAdjList);
        fsBox.add(fsTable);
        fsBox.add(fsButton);
        forbiddenStructures.add(fsBox);
    	
    	Box totalOptionsBox = Box.createVerticalBox();
    	Box horizontalBox1 = Box.createHorizontalBox();
    	Box verticalBox1 = Box.createVerticalBox();
    	
    	totalOptionsBox.add(restartPanel);
    	totalOptionsBox.add(maxPanel);
    	horizontalBox1.add(chemkinPanel);
    	verticalBox1.add(inchiPanel);
    	verticalBox1.add(EOS);
    	horizontalBox1.add(verticalBox1);
    	totalOptionsBox.add(horizontalBox1);
    	totalOptionsBox.add(forbiddenStructures);
    	
    	Options.add(totalOptionsBox);
    	
        JComboBox[] allTab = {eosCombo, inchiCombo, chemkinAUnits, 
        		chemkinEaUnits, chemkinVerbosity, chemkinSMILES};
        initializeJCB(allTab);
    	
    	JScrollPane scrolltab = new JScrollPane(Options);
        scrolltab.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    	return scrolltab;
    }

    public void runConditionFile(String c_path) {
        String[] rmgArgs = {c_path};
        RMG.main(rmgArgs);
        System.exit(0);
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
    		JComponent[] timeComps = {timeCombo, time};
    		JComponent[] convComps = {SpeciesConvName, conversion};
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
    	else if ("pDep".equals(event.getActionCommand())) {
    		JComponent[] pDepComps = {pdkmCombo, grainsizeCombo};
    		JComponent[] pdkmComps = {chebyTMin, chebyTMax, chebyTUnits, chebyTGen, chebyTRep,
    				chebyPMin, chebyPMax, chebyPUnits, chebyPGen, chebyPRep,
    				chebyPCombo, chebyTCombo, chebyP, chebyT};
    		if (pdepCombo.getSelectedItem().equals("off")) {
    			sdeCombo.setSelectedItem("off");
    			disableComponents(pDepComps);
    			disableComponents(pdkmComps);
    		} else {
    			sdeCombo.setSelectedItem("Frequency Groups");
    			enableComponents(pDepComps);
    			enableComponents(pdkmComps);
    			JComponent[] tlistComps = {chebyT};
    			JComponent[] plistComps = {chebyP};
    			JComponent[] tmaxminComps = {chebyTMin, chebyTMax};
    			JComponent[] pmaxminComps = {chebyPMin, chebyPMax};
    			if (chebyTCombo.getSelectedItem().equals("List")) {
    				enableComponents(tlistComps);
    				disableComponents(tmaxminComps);
    			} else if (chebyTCombo.getSelectedItem().equals("Max/Min")) {
    				disableComponents(tlistComps);
    				enableComponents(tmaxminComps);
    			}
    			if (chebyPCombo.getSelectedItem().equals("List")) {
    				enableComponents(plistComps);
    				disableComponents(pmaxminComps);
    			} else if (chebyPCombo.getSelectedItem().equals("Max/Min")) {
    				disableComponents(plistComps);
    				enableComponents(pmaxminComps);
    			}
    		}
    	}
    	else if ("pdkm".equals(event.getActionCommand())) {
    		JComponent[] maxminComps = {chebyTMin, chebyTMax, chebyPMin, chebyPMax};
    		JComponent[] listComps = {chebyT, chebyP};
    		JComponent[] nonRateComps = {chebyTGen, chebyPGen, chebyTCombo, chebyPCombo};
    		JComponent[] chebyComps = {chebyTRep, chebyPRep};

    		if (pdkmCombo.getSelectedItem().equals("Rate")) {
    			chebyTGen.setText("1");
    			chebyPGen.setText("1");
    			chebyTCombo.setSelectedIndex(0);
    			chebyPCombo.setSelectedIndex(0);
    			disableComponents(chebyComps);
    			disableComponents(nonRateComps);
    			disableComponents(maxminComps);
    			enableComponents(listComps);
    		} else {
    			enableComponents(nonRateComps);
    			enableComponents(maxminComps);
    			enableComponents(listComps);
        		if (pdkmCombo.getSelectedItem().equals("CHEB")) enableComponents(chebyComps);
        		else disableComponents(chebyComps);
				JComponent[] comps1 = {chebyT};
				JComponent[] comps2 = {chebyTMin, chebyTMax};
				JComponent[] comps3 = {chebyP};
				JComponent[] comps4 = {chebyPMin, chebyPMax};
    			if (chebyTCombo.getSelectedItem().equals("List")) {
    				enableComponents(comps1);
    				disableComponents(comps2);
    			}
    			else {
    				enableComponents(comps2);
    				disableComponents(comps1);
    			}
    			if (chebyPCombo.getSelectedItem().equals("List")) {
    				enableComponents(comps3);
    				disableComponents(comps4);
    			}
    			else {
    				enableComponents(comps4);
    				disableComponents(comps3);
    			}
    		}
    	}
    	else if ("pcheby".equals(event.getActionCommand())) {
			JComponent[] plistComps = {chebyP};
			JComponent[] pmaxminComps = {chebyPMin, chebyPMax};
			if (chebyPCombo.getSelectedItem().equals("List")) {
				enableComponents(plistComps);
				disableComponents(pmaxminComps);
			} else if (chebyPCombo.getSelectedItem().equals("Max/Min")) {
				disableComponents(plistComps);
				enableComponents(pmaxminComps);
			}
    	}
    	else if ("tcheby".equals(event.getActionCommand())) {
			JComponent[] tlistComps = {chebyT};
			JComponent[] tmaxminComps = {chebyTMin, chebyTMax};
			if (chebyTCombo.getSelectedItem().equals("List")) {
				enableComponents(tlistComps);
				disableComponents(tmaxminComps);
			} else if (chebyTCombo.getSelectedItem().equals("Max/Min")) {
				disableComponents(tlistComps);
				enableComponents(tmaxminComps);
			}
    	}
    	// Pruning options
    	else if ("Pruning".equals(event.getActionCommand())) {
    		JComponent[] indivComps = {textIndivSteps};
    		JComponent[] pruneComps = {textTerminTol, textPruneTol, textMinSpc, textMaxSpc};
    		if (interCombo.getSelectedItem().equals("AUTO")) {
    			disableComponents(indivComps);
    			disableComponents(pruneComps);
    		} else if (interCombo.getSelectedItem().equals("AUTOPRUNE")) {
    			disableComponents(indivComps);
    			enableComponents(pruneComps);
    		} else {
    			enableComponents(indivComps);
    			disableComponents(pruneComps);
    		}
    	}
    	//	Save the condition.txt file
    	else if ("saveCondition".equals(event.getActionCommand())) {
    		createConditionFile(false);
    	}
        //  Save and run the condition.txt file
        else if ("runCondition".equals(event.getActionCommand())) {
            createConditionFile(true);
        }
    	
	}
    
	//	Main panel
    private static GUI theApp;
	
	//	Tab0: Initialization
	String[] timeUnits = {"sec", "min", "hr", "day"};
	//	Tab1: Termination Sequence
	String[] species_blank = {" "};
	String[] goalOptions = {"Conversion", "ReactionTime"};
	//	Tab2: Initial Condition
	String[] tpModel = {"Constant"};
	String[] concUnits = {"mol/cm3", "mol/m3", "mol/liter"};
	String[] tempUnits = {"K", "C", "F"};
	String[] pressUnits = {"atm", "bar", "Pa", "torr"};
	String[] reactive = {"Reactive-RMG", "Reactive-User", "Inert"};
	//	Tab3: Error Analysis
	String[] yesnoOptions = {"No", "Yes"};
	//	Tab4: Dynamic Simulator
	String[] pdepOptions = {"off", "Reservoir State", "Modified Strong Collision"};
	String[] simOptions = {"DASSL", "DASPK"}; //"CHEMKIN"
	String[] sdeOptions = {"off","Frequency Groups"}; //"Three Frequency Model"
	String[] pdkmOptions = {"CHEB", "PLOG", "Rate"};
	String[] listormaxmin = {"List", "Max/Min"};
	String[] interOptions = {"AUTO","AUTOPRUNE","Indiv. time steps"};
	//	Tab : Additional Options
	String[] AUnitsOptions = {"moles", "molecules"};
	String[] EaUnitsOptions = {"kcal/mol", "cal/mol", "kJ/mol", "J/mol", "Kelvins"};
	String[] eosOptions = {"Ideal Gas"}; // "Liquid"
		
	JPanel
	//	Main Panel
	mainPanel;
	
    JComboBox
    //	Tab0: Initialization
    simulatorCombo, timeStepCombo, libraryCombo, smCombo,
    //	Tab1: Termination Sequence
    SpeciesConvName, controllerCombo, timeCombo,
    //	Tab2: Initial Condition
    UnitsInput, ReactiveInput, tempModelCombo, tempConstantUnit,
    	pressModelCombo, pressConstantUnit,
    //	Tab3: Error Analysis
    SpeciesSensName, eBarsCombo, sensCoeffCombo,
    //	Tab4: Dynamic Simulator
    pdepCombo, sdeCombo, pdkmCombo, chebyTCombo, chebyPCombo, grainsizeCombo,
    	interCombo,
    //	Tab : Additional Options
    eosCombo, inchiCombo, chemkinAUnits, chemkinEaUnits, chemkinVerbosity,
    	readRestartOnOff, writeRestartOnOff, chemkinSMILES,
    //
    chebyTUnits, chebyPUnits;
    
    JTextField
    //	Tab0: Initialization
    nameFiller, locationFiller, databasePath, prlPath, ptlPath,
    	timeStep, aTolerance, rTolerance, ptlLibName, prlLibName,
    	smLibName, smPath, ptranslLibName, ptranslPath,
    //	Tab1: Termination Sequence
    conversion, time, eTolerance,
    //	Tab2: Initial Condition
    NameInput, InChIInput, isPath, ConcInput,
    //	Tab3: Error Analysis
    scName,
    //	Tab : Additional Options
    maxCarbonNum, maxOxygenNum, maxRadicalNum, maxSulfurNum, maxSiliconNum,
    maxHeavyAtom, FSName,
    //
    chebyTMin, chebyTMax, chebyPMin, chebyPMax, chebyTGen, chebyTRep,
    	chebyPGen, chebyPRep, chebyT, chebyP,
    //
    textIndivSteps, textTerminTol, textPruneTol, textMinSpc, textMaxSpc;
    
    JTextArea
    //	Main Panel
    headerComments,
    //	Tab2: Initial Condition
    speciesAdjList, tempConstant, pressConstant,
    //	Tab4: Dynamic Simulator
    
    //	Tab : Additional Options
    FSAdjList;
    
    JButton 
    //	Main panel
    save, saveAndRun,
    //	Tab0: Initialization
    AddPRL, DeletePRL, databaseButton, prlButton, ptlButton,
    	AddPTL, DeletePTL, smButton, AddSM, DeleteSM,
    	ptranslButton, AddPTransL, DeletePTransL,
    //  Tab1: Termination
    //	Tab2: Initial Condition
    AddInput, DeleteInput, isButton,
    //	Tab3: Error Analysis
    AddSens, DeleteSens,
    //	Tab : Additional Options
    AddFS, DeleteFS;
    
    JTable
    //	Tab0: Initialization
    tablePRL, tablePTL, tableSM, tablePTransL,
    //  Tab1: Termination
    //	Tab2: Initial Condition
    tableInput,
    //	Tab3: Error Analysis
    tableSens,
    //	Tab : Additional Options
    tableFS;
    
    //	Tab0: Initialization
    MyTableModelPRL tmodelPRL, tmodelPTL, tmodelPTransL;
    MyTableModelSM tmodelSM;
    //  Tab1: Termination
    //	Tab2: Initial Condition
    MyTableModelInput tmodelInput;
    //	Tab3: Error Analysis
    MyTableModelSensitivity tmodelSens;
    //	Tab : Additional Options
    MyTableModelFS tmodelFS;
}
