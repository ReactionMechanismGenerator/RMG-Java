package jing.gui;

import java.util.Vector;

import javax.swing.table.AbstractTableModel;

public class MyTableModelInput extends AbstractTableModel {
	protected int NUM_COLUMNS = 5;
    protected int START_NUM_ROWS = 0;
    public int nextEmptyRow = 0;
    protected int numRows = 0;

    public String getColumnName(int column) {
    	switch (column) {
    	case 0:
    		return "Molecule";
    	case 1:
    		return "Concentration";
    	case 2:
    		return "Units";
    	case 3:
    		return "Status";
    	case 4:
    		return "Adjacency List";
    	}
    	return "";
    }

    public int getColumnCount() {
        return NUM_COLUMNS;
    }

    public int getRowCount() {
        return numRows;
    }

    public Object getValueAt(int row, int column) {
    	try {
            ICVector icvect = (ICVector)data2.elementAt(row);
            switch (column) {
            case 0:
            	return icvect.Name;
            case 1:
            	return icvect.Concentration;
            case 2:
            	return icvect.Unit;
            case 3:
            	return icvect.React;
            case 4:
            	return icvect.AdjList;
            }
    	} catch (Exception e) {
    		//System.out.println(e);
    	}
    	return "";
    }

    public void updateIC(ICVector icvector) {
        int rowIndex = icvector.Index-1;
        ICVector icvect = null;
        int index = -1; 
        boolean found = false;
        boolean addedRow = false;
        boolean tableChanged = false;
        
        int i = 0;
        while (!found && (i < nextEmptyRow)) {
        	icvect = (ICVector)data2.elementAt(i);
            if (icvect.Index == rowIndex) {
                found = true;
                index = i;
            } else {
                i++;
            }
        }

        if (found) {
        	data2.addElement(icvector);
        	numRows++;
        	addedRow = true;
        } else {
        	if (numRows <= nextEmptyRow) {
        		numRows++;
        		addedRow = true;
            }
            index = nextEmptyRow;
            if (data2 == null)	data2 = new Vector<Object>();
            data2.addElement(icvector);
        }
    
        nextEmptyRow++;
        
		if (numRows > START_NUM_ROWS) {
			++START_NUM_ROWS;
			tableChanged = true;
		}
        
        //Notify listeners that the data changed.
        if (addedRow) {
        	fireTableRowsInserted(index,index);
        } else if (tableChanged) {
        	fireTableStructureChanged();
        } else {
        	fireTableRowsUpdated(index,index);
        }
    }

    public void deleteRow(int row) {
		String dNameInput = "";
		String dConcInput = "";
		String dUnitInput = "";
		String dReactInput = "";
		String dAdjInput = "";
    	for (int i=0; i<nextEmptyRow; i++) {
    		if (i!=row) {
    			dNameInput += getValueAt(i,0) + "\r";
    			dConcInput += getValueAt(i,1) + "\r";
    			dUnitInput += getValueAt(i,2) + "\r";
    			dReactInput += getValueAt(i,3) + "\r";
    			if (getValueAt(i,4).equals("")) {
    				dAdjInput += " " + "\t";
    			} else
    				dAdjInput += getValueAt(i,4) + "\t";
    		}
		}
    	clear();
    	String[] dNameInputIndiv = dNameInput.split("[\r]",0);
    	String[] dConcInputIndiv = dConcInput.split("[\r]",0);
    	String[] dUnitInputIndiv = dUnitInput.split("[\r]",0);
    	String[] dReactInputIndiv = dReactInput.split("[\r]",0);
    	String[] dAdjInputIndiv = dAdjInput.split("[\t]",0);
    	if (dNameInputIndiv.length >= START_NUM_ROWS) {
    		START_NUM_ROWS = dNameInputIndiv.length;
    	} else
    		START_NUM_ROWS = 1;
    	if (!dNameInput.equals("")) {
    		for (int j=0; j<dNameInputIndiv.length; j++) {
    			ICVector Molecule = new ICVector(j+1, dNameInputIndiv[j], dConcInputIndiv[j], dUnitInputIndiv[j], dReactInputIndiv[j], dAdjInputIndiv[j]);
    			updateIC(Molecule);
    		}
    	}
    }

    public void clear() {
    	int oldNumRows = numRows;
    	if (data2 != null) data2.removeAllElements();
        nextEmptyRow = 0;
        numRows = 0;

        if (oldNumRows > START_NUM_ROWS) {
        	fireTableRowsDeleted(START_NUM_ROWS, oldNumRows - 1);
        }
        fireTableRowsUpdated(0, START_NUM_ROWS - 1);
    }
    
    Vector<Object> data2;
}