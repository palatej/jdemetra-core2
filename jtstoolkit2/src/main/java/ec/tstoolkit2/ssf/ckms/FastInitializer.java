/*
 * Copyright 2013 National Bank of Belgium
 *
 * Licensed under the EUPL, Version 1.1 or – as soon they will be approved 
 * by the European Commission - subsequent versions of the EUPL (the "Licence");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 *
 * http://ec.europa.eu/idabc/eupl
 *
 * Unless required by applicable law or agreed to in writing, software 
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and 
 * limitations under the Licence.
 */
package ec.tstoolkit2.ssf.ckms;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.design.Development;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.SsfException;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;

/**
 *
 * @param <S>
 * @author Jean Palate
 */
@Development(status = Development.Status.Alpha)
public class FastInitializer<S extends ISsf> implements IFastInitializer<S> {

    /**
     * K = TPZ', L=K, F=ZPZ'+H
     * @param ssf
     * @param fstate
     * @return
     */
    @Override
    public boolean initialize(final S ssf, final FastState fstate) {
        if (!ssf.isTimeInvariant()) {
            return false;
        }
        ISsfDynamics dynamics = ssf.getDynamics();
        ISsfMeasurement measurement = ssf.getMeasurement();
        if (dynamics.isDiffuse()) {
            return false;
        }
        SubMatrix P0 = Matrix.square(dynamics.getStateDim()).subMatrix();
        dynamics.Pf0(P0, StateInfo.Forecast);
        measurement.ZM(0, P0, fstate.k);
        dynamics.TX(0, fstate.k);
        fstate.l.copy(fstate.k);
        fstate.f = measurement.ZX(0, fstate.k);
        if (measurement.hasErrors()) {
            fstate.f += measurement.errorVariance(0);
        }
        return true;
    }
}
