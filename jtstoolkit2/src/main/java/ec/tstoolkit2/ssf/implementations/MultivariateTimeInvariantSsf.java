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
import ec.tstoolkit2.ssf.multivariate.ISsfMeasurements;
import ec.tstoolkit2.ssf.multivariate.MultivariateSsf;

/**
 *
 * @author Jean Palate
 */
public class MultivariateTimeInvariantSsf extends MultivariateSsf{
    public static IMultivariateSsf of(IMultivariateSsf ssf){
        return of(ssf, StateInfo.Forecast);
    }

    public static IMultivariateSsf of(IMultivariateSsf ssf, StateInfo info){
        TimeInvariantDynamics td=TimeInvariantDynamics.of(ssf.getDynamics(), info);
        if (td == null)
            return null;
        TimeInvariantMeasurements tm=TimeInvariantMeasurements.of(ssf.getStateDim(), ssf.getMeasurements());
        return new MultivariateTimeInvariantSsf(td, tm);
    }
    
    public static IMultivariateSsf of(ISsf ssf){
        return of(ssf, StateInfo.Forecast);
    }
    
    public static IMultivariateSsf of(ISsf ssf, StateInfo info){
        TimeInvariantDynamics td=TimeInvariantDynamics.of(ssf.getDynamics(), info);
        if (td == null)
            return null;
        TimeInvariantMeasurements tm=TimeInvariantMeasurements.of(ssf.getStateDim(), ssf.getMeasurement());
        return new MultivariateTimeInvariantSsf(td, tm);
    }

    private MultivariateTimeInvariantSsf(final ISsfDynamics dynamics, ISsfMeasurements measurement) {
        super(dynamics, measurement);
    }
}
