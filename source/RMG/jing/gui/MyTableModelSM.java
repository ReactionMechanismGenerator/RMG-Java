package jing.gui;

import java.util.Vector;

import javax.swing.table.AbstractTableModel;

public class MyTableModelSM extends AbstractTableModel {
	protected int NUM_COLUMNS = 3;
    protected int START_NUM_ROWS = 0;
    public int nextEmptyRow = 0;
    protected int numRows = 0;

    public String getColumnName(int column) {
    	switch (column) {
    	case 0:
    		return "Name";
    	case 1:
    		return "Generate Reactions";
    	case 2:
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
            SMVector smvect = (SMVector)data1.elementAt(row);
            switch (column) {
            case 0:
            	return smvect.Name;
            case 1:
            	return smvect.React;
            case 2:
            	return smvect.Location;
            }
    	} catch (Exception e) {
    		//System.out.println(e);
    	}
    	return "";
    }

    public void updateSM(SMVector smvector) {
        int rowIndex = smvector.Index-1;
        SMVector smvect = null;
        int index = -1; 
        boolean found = false;
        boolean addedRow = false;
        boolean tableChanged = false;
        
        int i = 0;
        while (!found && (i < nextEmptyRow)) {
        	smvect = (SMVector)data1.elementAt(i);
            if (smvect.Index == rowIndex) {
                found = true;
                index = i;
            } else {
                i++;
            }
        }

        if (found) {        		
        	data1.addElement(smvector);
        	numRows++;
        	addedRow = true;
        } else {
        	if (numRows <= nextEmptyRow) {
        		numRows++;
        		addedRow = true;
            }
            index = nextEmptyRow;
            if (data1 == null) data1 = new Vector<Object>();
            data1.addElement(smvector);
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
		String dNameSM = "";
		String dLocSM = "";
		String dReactSM = "";
    	for (int i=0; i<nextEmptyRow; i++) {
    		if (i!=row) {
    			dNameSM += getValueAt(i,0) + "\r";
    			dReactSM += getValueAt(i,1) + "\r";
    			dLocSM += getValueAt(i,2) + "\r";
    		}
		}
    	clear();
    	String[] dNameSMIndiv = dNameSM.split("[\r]",0);
    	String[] dReactSMIndiv = dReactSM.split("[\r]",0);
    	String[] dLocSMIndiv = dLocSM.split("[\r]",0);
    	if (dNameSMIndiv.length >= START_NUM_ROWS) {
    		START_NUM_ROWS = dNameSMIndiv.length;
    	} else
    		START_NUM_ROWS = 1;
    	if (!dNameSM.equals("")) { 
    		for (int j=0; j<dNameSMIndiv.length; j++) {
    			SMVector SM = new SMVector(j+1, dNameSMIndiv[j], dReactSMIndiv[j], dLocSMIndiv[j]);
    			updateSM(SM);
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