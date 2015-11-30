/*
 * Copyright 2015 National Bank of Belgium
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
package ec.tstoolkit2.ssf.array;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit2.ssf.multivariate.IMultivariateSsf;
import ec.tstoolkit2.ssf.multivariate.ISsfMeasurements;
import ec.tstoolkit2.ssf.multivariate.IMultivariateSsfData;
import ec.tstoolkit2.ssf.multivariate.IMultivariateFilteringResults;
import ec.tstoolkit.maths.matrices.ElementaryTransformations;
import ec.tstoolkit.maths.matrices.LowerTriangularMatrix;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.multivariate.TransformedPredictionErrors;

/**
 * Array form of the Kalman filter
 *
 * @author Jean Palate
 */
public class MultivariateArrayFilter {

    private LState state;
    private TransformedPredictionErrors perrors;
    private ISsfMeasurements measurements;
    private ISsfDynamics dynamics;
    private IMultivariateSsfData data;
    private int pos, end, nm, dim, nres;
    private Matrix A;

    /**
     *
     */
    public MultivariateArrayFilter() {
    }

    /**
     */
    protected void error() {
        DataBlock U=perrors.getTransformedPredictionErrors();
        Matrix L=perrors.getCholeskyFactor();
        U.set(0);
        for (int i = 0; i < nm; ++i) {
            double y = data.get(pos, i);
            U.set(i, y - measurements.ZX(pos, i, state.a));
        }
        LowerTriangularMatrix.rsolve(L, U, State.ZERO);
    }

    private boolean initFilter() {
        pos = 0;
        end = data.getCount();
        nm = measurements.getMaxCount();
        nres = dynamics.getInnovationsDim();
        dim = dynamics.getStateDim();
        A = new Matrix(dim + nm, dim + nm + nres);
        return true;
    }

    private int initState() {
        state = new LState(L());
        perrors = new TransformedPredictionErrors(dim, nm);
        if (!dynamics.a0(state.a, StateInfo.Forecast)) {
            return -1;
        }
        Matrix P0 = new Matrix(dim, dim);
        if (!dynamics.Pf0(P0.subMatrix(), StateInfo.Forecast)) {
            return -1;
        }
        SymmetricMatrix.lcholesky(P0, State.ZERO);
        state.L.copy(P0.subMatrix());

        return 0;
    }

    /**
     *
     * @param ssf
     * @param data
     * @param rslts
     * @return
     */
    public boolean process(final IMultivariateSsf ssf, final IMultivariateSsfData data, final IMultivariateFilteringResults rslts) {
        measurements = ssf.getMeasurements();
        if (!measurements.isHomogeneous()) {
            return false;
        }
        dynamics = ssf.getDynamics();
        this.data = data;
        if (!initFilter()) {
            return false;
        }
        pos = initState();
        if (pos < 0) {
            return false;
        }
        rslts.open(ssf, this.data);
        do {
            preArray();
            ElementaryTransformations.fastGivensTriangularize(A.subMatrix());
            postArray();
            error();
            rslts.save(pos, perrors);
            nextState();
        } while (++pos < end);
        rslts.close();
        return true;
    }

    private void preArray() {
        measurements.ZM(pos, L(), ZL());
        dynamics.TM(pos, L());
        R().set(0);
        measurements.R(pos, R());
        U().set(0);
        dynamics.U(pos, U());
        K().set(0);
    }

    private void postArray() {
        perrors.getCholeskyFactor().subMatrix().copy(R());
        perrors.getK().subMatrix().copy(K());
    }

    private void nextState() {
        dynamics.TX(pos, state.a);
        for (int i = 0; i < nm; ++i) {
            state.a.addAY(perrors.getTransformedPredictionErrors().get(i), perrors.getK().column(i));
        }

    }

    private SubMatrix R() {
        return A.subMatrix(0, nm, 0, nm);
    }

    private SubMatrix K() {
        return A.subMatrix(nm, nm + dim, 0, nm);
    }

    private SubMatrix ZL() {
        return A.subMatrix(0, nm, nm, nm + dim);
    }

    private SubMatrix L() {
        return A.subMatrix(nm, nm + dim, nm, nm + dim);
    }

    private SubMatrix U() {
        return A.subMatrix(nm, nm + dim, nm + dim, -1);
    }
}
