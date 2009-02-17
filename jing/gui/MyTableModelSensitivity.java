package jing.gui;

import java.util.Vector;

import javax.swing.table.AbstractTableModel;

public class MyTableModelSensitivity extends AbstractTableModel {
	protected int NUM_COLUMNS = 1;
    protected int START_NUM_ROWS = 0;
    public int nextEmptyRow = 0;
    protected int numRows = 0;

    public String getColumnName(int column) {
    	switch (column) {
    	case 0:
    		return "Molecule";
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
            SensVector sensvect = (SensVector)data3.elementAt(row);
            switch (column) {
            case 0:
            	return sensvect.Name;
            }
    	} catch (Exception e) {
    		//System.out.println(e);
    	}
    	return "";
    }

    public void updateSens(SensVector sensvector) {
        int rowIndex = sensvector.Index-1;
        SensVector sensvect = null;
        int index = -1; 
        boolean found = false;
        boolean addedRow = false;
        boolean tableChanged = false;
        
        int i = 0;
        while (!found && (i < nextEmptyRow)) {
        	sensvect = (SensVector)data3.elementAt(i);
            if (sensvect.Index == rowIndex) {
                found = true;
                index = i;
            } else {
                i++;
            }
        }

        if (found) {
    		data3.addElement(sensvector);
        	numRows++;
        	addedRow = true;            	
        	
        } else {
        	if (numRows <= nextEmptyRow) {
        		numRows++;
        		addedRow = true;            		
            }
        	
            index = nextEmptyRow;
            if (data3 == null) data3 = new Vector<Object>();
            data3.addElement(sensvector);
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
		String dNameSens = "";
    	for (int i=0; i<nextEmptyRow; i++) {
    		if (i!=row) {
    			dNameSens += getValueAt(i,0) + "\r";
    		}
		}
    	clear();
    	String[] dNameSensIndiv = dNameSens.split("[\r]",0);
    	if (dNameSensIndiv.length >= START_NUM_ROWS) {
    		START_NUM_ROWS = dNameSensIndiv.length;
    	} else
    		START_NUM_ROWS = 1;
    	if (!dNameSens.equals("")) {
    		for (int j=0; j<dNameSensIndiv.length; j++) {
    			SensVector Sensitivity = new SensVector(j+1, dNameSensIndiv[j]);
    			updateSens(Sensitivity);
    		}
    	}
    }

    public void clear() {
    	int oldNumRows = numRows;
        if (data3 != null) data3.removeAllElements();
        nextEmptyRow = 0;
        numRows = 0;

        if (oldNumRows > START_NUM_ROWS) {
        	fireTableRowsDeleted(START_NUM_ROWS, oldNumRows - 1);
        }
        fireTableRowsUpdated(0, START_NUM_ROWS - 1);
    }
    
    Vector<Object> data3;
}