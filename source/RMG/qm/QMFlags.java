package qm;

import jing.chem.QMTP;

/**
 * Container type with some user input flags for QM applications like
 * {@link ThermoDataEstimator}
 *
 */
public class QMFlags {
	public Boolean qmActive;
	
	/**
	 * TODO method attribute should become Enum subtype
	 * However, for now, this will break {@link QMTP} and others 
	 */
	public String method;
	
	public Boolean qmOnCyclicsOnly;
	
	public Integer maxRadNumForQM;

	public Integer connectivityCheck;

	public Boolean qmOnNonAromaticFusedCyclicsOnly;
	
	public QMFlags(){
		qmActive = null;
		method = null;
		qmOnCyclicsOnly = null;
		qmOnNonAromaticFusedCyclicsOnly = null;
		maxRadNumForQM = null;
		connectivityCheck = null;
	}
}
