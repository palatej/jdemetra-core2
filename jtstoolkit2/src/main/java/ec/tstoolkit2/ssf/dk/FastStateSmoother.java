/*
 * Copyright 2013-2014 National Bank of Belgium
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
package ec.tstoolkit2.ssf.dk;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.DataBlockStorage;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.univariate.DefaultDisturbanceSmoothingResults;
import ec.tstoolkit2.ssf.univariate.DisturbanceSmoother;
import ec.tstoolkit2.ssf.univariate.IDisturbanceSmoothingResults;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfData;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;

/**
 *
 * @author Jean Palate
 */
public class FastStateSmoother {

    private ISsfDynamics dynamics;
    private ISsfMeasurement measurement;
    private Matrix S;

    public DataBlockStorage process(ISsf ssf, ISsfData data) {
        initSsf(ssf);
        int dim = dynamics.getStateDim();
        int n = data.getCount();
        DataBlockStorage storage = new DataBlockStorage(dim, n);
        DefaultDisturbanceSmoothingResults srslts = DefaultDisturbanceSmoothingResults.light(measurement.hasErrors());
        srslts.prepare(ssf, 0, n);
        DataBlock a = initialState(ssf, data, srslts);
        storage.save(0, a);
        int pos = 1;
        do {
            // next: a(t+1) = T a(t) + S*r(t)
            dynamics.TX(pos, a);
            loadInfo(pos);
            a.addProduct(srslts.u(pos), S.rows());
            // T
            storage.save(pos++, a);
        } while (pos < n);

        return storage;
    }

    private void initSsf(ISsf ssf) {
        dynamics = ssf.getDynamics();
        measurement = ssf.getMeasurement();
        int dim = dynamics.getStateDim(), resdim = dynamics.getInnovationsDim();
        S = new Matrix(dim, resdim);
        if (dynamics.isTimeInvariant()) {
            dynamics.S(0, S.subMatrix());
        }
    }

    private DataBlock initialState(ISsf ssf, ISsfData data, IDisturbanceSmoothingResults srslts) {
        if (dynamics.isDiffuse()) {
            DiffuseDisturbanceSmoother sm = new DiffuseDisturbanceSmoother();
            sm.setCalcVariances(false);
            sm.process(ssf, data, srslts);
            return sm.firstSmoothedState();
        } else {
            DisturbanceSmoother sm = new DisturbanceSmoother();
            sm.setCalcVariances(false);
            sm.process(ssf, data, srslts, 0);
            return sm.firstSmoothedState();
        }
    }

    private void loadInfo(int pos) {
        if (!dynamics.isTimeInvariant()) {
            dynamics.S(pos, S.subMatrix());
        }
    }

}
