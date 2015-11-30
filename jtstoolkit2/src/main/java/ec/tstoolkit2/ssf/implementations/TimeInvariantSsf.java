/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.implementations;

import ec.tstoolkit2.ssf.multivariate.IMultivariateSsf;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import ec.tstoolkit2.ssf.multivariate.ISsfMeasurements;
import ec.tstoolkit2.ssf.multivariate.MultivariateSsf;
import ec.tstoolkit2.ssf.univariate.Ssf;

/**
 *
 * @author Jean Palate
 */
public class TimeInvariantSsf extends Ssf{
    public static ISsf of(ISsf ssf){
        return of(ssf, StateInfo.Forecast);
    }

    public static ISsf of(ISsf ssf, StateInfo info){
        TimeInvariantDynamics td=TimeInvariantDynamics.of(ssf.getDynamics(), info);
        if (td == null)
            return null;
        TimeInvariantMeasurement tm=TimeInvariantMeasurement.of(ssf.getStateDim(), ssf.getMeasurement());
        return new TimeInvariantSsf(td, tm);
    }
    
    private TimeInvariantSsf(final ISsfDynamics dynamics, ISsfMeasurement measurement) {
        super(dynamics, measurement);
    }
}
