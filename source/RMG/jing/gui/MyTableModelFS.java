package jing.gui;

import java.util.Vector;

import javax.swing.table.AbstractTableModel;

public class MyTableModelFS extends AbstractTableModel {
	protected int NUM_COLUMNS = 2;
    protected int START_NUM_ROWS = 0;
    public int nextEmptyRow = 0;
    protected int numRows = 0;

    public String getColumnName(int column) {
    	switch (column) {
    	case 0:
    		return "Name";
    	case 1:
    		return "Adjacency List";
    	}
    	return "";
    }

    public int getColumnCount() {
        return NUM_COLUMNS;
    }

    public int getRowCount() {
        if (numRows < START_NUM_ROWS) {
            return START_NUM_ROWS;
        } else
            return numRows;
    }

    public Object getValueAt(int row, int column) {
    	try {
            FSVector fsvect = (FSVector)data.elementAt(row);
            switch (column) {
            case 0:
            	return fsvect.Name;
            case 1:
            	return fsvect.AdjList;
            }
    	} catch (Exception e) {
    		//System.out.println(e);
    	}
    	return "";
    }

    public void updateFS(FSVector fsvector) {
        int rowIndex = fsvector.Index-1;
        FSVector fsvect = null;
        int index = -1; 
        boolean found = false;
        boolean addedRow = false;
        boolean tableChanged = false;
        
        int i = 0;
        while (!found && (i < nextEmptyRow)) {
        	fsvect = (FSVector)data.elementAt(i);
            if (fsvect.Index == rowIndex) {
                found = true;
                index = i;
            } else {
                i++;
            }
        }

        if (found) {
    		data.addElement(fsvector);
        	numRows++;
        	addedRow = true;            	
        	
        } else {
        	if (numRows <= nextEmptyRow) {
        		numRows++;
        		addedRow = true;            		
            }
        	
            index = nextEmptyRow;
            if (data == null) data = new Vector<Object>();
            data.addElement(fsvector);
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
		String dNameFS = "";
		String dAdjListFS = "";
    	for (int i=0; i<nextEmptyRow; i++) {
    		if (i!=row) {
    			dNameFS += getValueAt(i,0) + "\r";
    			if (getValueAt(i,1).equals("")) {
    				dAdjListFS += " " + "\t";
    			} else
    				dAdjListFS += getValueAt(i,1) + "\t";
    		}
		}
    	clear();
    	String[] dNameFSIndiv = dNameFS.split("[\r]",0);
    	String[] dAdjListFSIndiv = dAdjListFS.split("[\t]",0);
    	if (dNameFSIndiv.length >= START_NUM_ROWS) {
    		START_NUM_ROWS = dNameFSIndiv.length;
    	} else
    		START_NUM_ROWS = 1;
    	if (!dNameFS.equals("")) {
    		for (int j=0; j<dNameFSIndiv.length; j++) {
    			FSVector forbidStruct = new FSVector(j+1, dNameFSIndiv[j], dAdjListFSIndiv[j]);
    			updateFS(forbidStruct);
    		}
    	}
    }

    public void clear() {
    	int oldNumRows = numRows;
        if (data != null) data.removeAllElements();
        nextEmptyRow = 0;
        numRows = 0;

        if (oldNumRows > START_NUM_ROWS) {
        	fireTableRowsDeleted(START_NUM_ROWS, oldNumRows - 1);
        }
        fireTableRowsUpdated(0, START_NUM_ROWS - 1);
    }
    
    Vector<Object> data;
}