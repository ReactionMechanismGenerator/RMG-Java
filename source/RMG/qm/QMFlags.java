package qm;

import jing.chem.QMTP;

/**
 * Container type with some user input flags for QM applications like {@link ThermoDataEstimator}
 */
public class QMFlags {
    public String TDSTRATEGY;
    /**
     * TODO method attribute should become Enum subtype However, for now, this will break {@link QMTP} and others
     */
    public String method;
    public Integer maxRadNumForQM;
    public Integer connectivityCheck;

    public QMFlags() {
        TDSTRATEGY = null;
        method = null;
        maxRadNumForQM = null;
        connectivityCheck = null;
    }
}
