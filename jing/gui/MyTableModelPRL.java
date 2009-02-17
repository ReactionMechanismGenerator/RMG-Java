package jing.gui;

import java.util.Vector;

import javax.swing.table.AbstractTableModel;

public class MyTableModelPRL extends AbstractTableModel {
	protected int NUM_COLUMNS = 2;
    protected int START_NUM_ROWS = 0;
    public int nextEmptyRow = 0;
    protected int numRows = 0;

    public String getColumnName(int column) {
    	switch (column) {
    	case 0:
    		return "Name";
    	case 1:
    		return "Location";
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
            PRLVector prlvect = (PRLVector)data1.elementAt(row);
            switch (column) {
            case 0:
            	return prlvect.Name;
            case 1:
            	return prlvect.Location;
            }
    	} catch (Exception e) {
    		//System.out.println(e);
    	}
    	return "";
    }

    public void updatePRL(PRLVector prlvector) {
        int rowIndex = prlvector.Index-1;
        PRLVector prlvect = null;
        int index = -1; 
        boolean found = false;
        boolean addedRow = false;
        boolean tableChanged = false;
        
        int i = 0;
        while (!found && (i < nextEmptyRow)) {
        	prlvect = (PRLVector)data1.elementAt(i);
            if (prlvect.Index == rowIndex) {
                found = true;
                index = i;
            } else {
                i++;
            }
        }

        if (found) {        		
        	data1.addElement(prlvector);
        	numRows++;
        	addedRow = true;
        } else {
        	if (numRows <= nextEmptyRow) {
        		numRows++;
        		addedRow = true;
            }
            index = nextEmptyRow;
            if (data1 == null) data1 = new Vector<Object>();
            data1.addElement(prlvector);
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
		String dNamePRL = "";
		String dLocPRL = "";
    	for (int i=0; i<nextEmptyRow; i++) {
    		if (i!=row) {
    			dNamePRL += getValueAt(i,0) + "\r";
    			dLocPRL += getValueAt(i,1) + "\r";
    		}
		}
    	clear();
    	String[] dNamePRLIndiv = dNamePRL.split("[\r]",0);
    	String[] dLocPRLIndiv = dLocPRL.split("[\r]",0);
    	if (dNamePRLIndiv.length >= START_NUM_ROWS) {
    		START_NUM_ROWS = dNamePRLIndiv.length;
    	} else
    		START_NUM_ROWS = 1;
    	if (!dNamePRL.equals("")) { 
    		for (int j=0; j<dNamePRLIndiv.length; j++) {
    			PRLVector PRL = new PRLVector(j+1, dNamePRLIndiv[j], dLocPRLIndiv[j]);
    			updatePRL(PRL);
    		}
    	}
    }

    public void clear() {
    	int oldNumRows = numRows;
        if (data1 != null) data1.removeAllElements();
        nextEmptyRow = 0;
        numRows = 0;

        if (oldNumRows > START_NUM_ROWS) {
        	fireTableRowsDeleted(START_NUM_ROWS, oldNumRows - 1);
        }
        fireTableRowsUpdated(0, START_NUM_ROWS - 1);
    }
    
    Vector<Object> data1;
}