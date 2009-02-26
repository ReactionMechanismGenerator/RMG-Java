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



/**
 * @author User1
 *
 */

import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.KeyStroke;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import static java.awt.event.InputEvent.*;

public class GUIWindow extends JFrame {
	// Constructor
	public GUIWindow(String title, GUI theApp) {
		setTitle(title);								// Set the window title

		this.theApp = theApp;
		setDefaultCloseOperation(EXIT_ON_CLOSE);						// Default is exit the application
		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);											// Add the menu bar to the window
		
		JMenu fileMenu = new JMenu("File");			// Create File menu
        JMenu runMenu = new JMenu("Run");           // Create Run menu
		JMenu helpMenu = new JMenu("Help");			// Create Help menu
		fileMenu.setMnemonic('F');					// Create shortcut
        runMenu.setMnemonic('R');                   //  Create shortcut
		helpMenu.setMnemonic('H');					// Create shortcut
		
		// Construct the File drop-down menu
		JMenuItem openItem = fileMenu.add("Open");								// Add Open item
		openItem.addActionListener(new MenuListener("Open"));
		openItem.setAccelerator(KeyStroke.getKeyStroke('O',CTRL_DOWN_MASK));
		JMenuItem saveItem = fileMenu.add("Save");
		saveItem.addActionListener(new MenuListener("Save"));
		saveItem.setAccelerator(KeyStroke.getKeyStroke('S',CTRL_DOWN_MASK));
		JMenuItem closeItem = fileMenu.add("Close");							// Add Close item
		closeItem.addActionListener(new MenuListener("Close"));

        //  Construct the Run drop-down menu
        JMenuItem runItem = runMenu.add("Run RMG");                             // Add Run item
        runItem.addActionListener(new MenuListener("Run"));
        runItem.setAccelerator(KeyStroke.getKeyStroke('R',CTRL_DOWN_MASK));

		// Construct the Help drop-down menu
		JMenuItem aboutItem = helpMenu.add("About RMG");
		aboutItem.addActionListener(new MenuListener("Help"));
		aboutItem.setAccelerator(KeyStroke.getKeyStroke('H',CTRL_DOWN_MASK));
		
		menuBar.add(fileMenu);			// Add the File menu
        menuBar.add(runMenu);           // Add the Run menu
		menuBar.add(helpMenu);			// Add the Help menu
	}
	
    class MenuListener implements ActionListener {
    	MenuListener(String selection) {
    		this.selection = selection;
    	}

		public void actionPerformed(ActionEvent event) {
			if (event.getActionCommand().equals("Open")) {
				File openFile = null;
				openFile = theApp.askUserForInput("Open file", false);
				if (openFile != null) theApp.openConditionFile(openFile);
			}
			else if (event.getActionCommand().equals("Save")) {
				theApp.createConditionFile(false);
			}
			else if (event.getActionCommand().equals("Close")) {
				System.exit(0);
			}
            else if (event.getActionCommand().equals("Run RMG")) {
                File runFile = null;
                runFile = theApp.askUserForInput("Run file", false);
                if (runFile != null) theApp.runConditionFile(runFile.getAbsolutePath());
            }
			else if (event.getActionCommand().equals("About RMG")) {
				JOptionPane.showMessageDialog(theApp, 
						"RMG v3.0\n\n(c) Copyright Prof. William H. Green and RMG developers 2009.\n" +
						"http://sourceforge.net/projects/rmg/\n" +
						"http://rmg.mit.edu/\n" + 
						"rmg_dev@mit.edu",
						"About RMG", JOptionPane.INFORMATION_MESSAGE);
			}
		}
		
		String selection;
    }
	
	GUI theApp;
}