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
package ec.tstoolkit2.ssf.multivariate;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.DataBlockIterator;
import ec.tstoolkit.maths.matrices.LowerTriangularMatrix;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;

/**
 *
 * @author Jean Palate
 */
public class MultivariateOrdinaryFilter {

    public static interface Initializer {

        int initialize(State state, IMultivariateSsf ssf, IMultivariateSsfData data);
    }

    private final Initializer initializer;
    private State state;
    private TransformedPredictionErrors perrors;
    private ISsfMeasurements measurements;
    private ISsfDynamics dynamics;
    private IMultivariateSsfData data;
    private int pos, end;

    /**
     *
     */
    public MultivariateOrdinaryFilter() {
        initializer = null;
    }

    /**
     *
     * @param initializer
     */
    public MultivariateOrdinaryFilter(final Initializer initializer) {
        this.initializer = initializer;
    }

    private int countMeasurements() {
        int n = 0;
        for (int i = 0; i < measurements.getCount(pos); ++i) {
            if (!data.isMissing(pos, i)) {
                ++n;
            }
        }
        return n;
    }

    /**
     * Computes zm = Z * M
     *
     * @param M
     * @param zm
     */
    private void ZM(SubMatrix M, SubMatrix zm) {
        DataBlockIterator zrows = zm.rows();
        DataBlock zr = zrows.getData();
        int imax = measurements.getCount(pos);
        for (int i = 0; i < imax; ++i) {
            if (!data.isMissing(pos, i)) {
                measurements.ZM(pos, i, M, zr);
                if (!zrows.next()) {
                    return;
                }
            }
        }
    }

    private void addH(Matrix P) {
        int nm = measurements.getCount(pos);
        Matrix H = Matrix.square(nm);
        measurements.H(pos, H.subMatrix());
        for (int i = 0, r = 0; i < nm; ++i) {
            if (!data.isMissing(pos, i)) {
                for (int j = 0, c = 0; j < i; ++j) {
                    if (!data.isMissing(pos, j)) {
                        double h = H.get(i, j);
                        P.add(r, c, h);
                        P.add(c, r, h);
                        ++c;
                    }
                }
                P.add(r, r, H.get(i, i));
                ++r;
            }
        }
    }

    /**
     * Computes a(t+1|t), P(t+1|t) from a(t|t), P(t|t) a(t+1|t) = T(t)a(t|t)
     * P(t+1|t) = T(t)P(t|t)T'(t)
     */
    protected void pred() {
        if (state.getInfo() != StateInfo.Forecast) {
            state.setInfo(StateInfo.Forecast);
            SubMatrix P = state.P().subMatrix();
            DataBlock a = state.a();
            dynamics.TX(pos, a);
            dynamics.TVT(pos, P);
            dynamics.addV(pos, P);
        }
    }

    /**
     * Computes: e(t)=y(t) - Z(t)a(t|t-1)) F(t)=Z(t)P(t|t-1)Z'(t)+H(t) F(t) =
     * L(t)L'(t) E(t) = e(t)L'(t)^-1 K(t)= P(t|t-1)Z'(t)L'(t)^-1
     *
     * Not computed for missing values
     */
    protected void error() {
        int nobs = countMeasurements();
        if (nobs == 0) {
            perrors = null;
        } else {
            perrors = new TransformedPredictionErrors(dynamics.getStateDim(), nobs);
            Matrix L = perrors.getCholeskyFactor();
            // K = PZ'(ZPZ'+H)^-1/2
            // computes (ZP)' in K'. Missing values are set to 0 
            // Z~v x r, P~r x r, K~r x v
            SubMatrix F = L.subMatrix(), K = perrors.getK().subMatrix();
            ZM(state.P().subMatrix(), K.transpose());
            // computes ZPZ'; results in pe_.L
            ZM(K, F);
            addH(L);
            SymmetricMatrix.reinforceSymmetry(L);

            // pe_L contains the Cholesky factor !!!
            SymmetricMatrix.lcholesky(L, State.ZERO);

            // We put in K  PZ'*(ZPZ'+H)^-1/2 = PZ'* F^-1 = PZ'*(LL')^-1/2 = PZ'(L')^-1
            // K L' = PZ' or L K' = ZP
            LowerTriangularMatrix.rsolve(L, K.transpose(), State.ZERO);
            DataBlock U = perrors.getTransformedPredictionErrors();
            for (int i = 0, j = 0; i < measurements.getCount(pos); ++i) {
                if (!data.isMissing(pos, i)) {
                    double y = data.get(pos, i);
                    U.set(j, y - measurements.ZX(pos, i, state.a()));
                    ++j;
                }
            }
            // E = e*L'^-1 or E L' = e or L*E' = e'
            LowerTriangularMatrix.rsolve(L, U, State.ZERO);
        }
    }

    /**
     * Updates the state vector and its covariance a(t|t) = a(t|t-1) + e(t)
     */
    protected void update() {
        state.setInfo(StateInfo.Concurrent);
        if (perrors == null) {
            return;
        }
        int n = perrors.getK().getColumnsCount();
        // P = P - (M)* F^-1 *(M)' --> Symmetric
        // PZ'(LL')^-1 ZP' =PZ'L'^-1*L^-1*ZP'
        // A = a + (M)* F^-1 * v
//        for (int i = 0; i < n; ++i) {
//            state.P().addXaXt(-1, perrors.getK().column(i));//, state_.K.column(i));
//            state.a().addAY(perrors.getTransformedPredictionErrors().get(i), perrors.getK().column(i));
//        }
        Matrix P = state.P();
        Matrix K = perrors.getK();
        DataBlock U = perrors.getTransformedPredictionErrors();
        for (int i = 0; i < n; ++i) {
            P.addXaXt(-1, K.column(i));//, state_.K.column(i));
            state.a().addAY(U.get(i), K.column(i));
        }
    }

    /**
     *
     * @return
     */
    public State getState() {
        return state;
    }

    private boolean initialize(IMultivariateSsf ssf, IMultivariateSsfData data) {
        this.data = data;
        measurements = ssf.getMeasurements();
        dynamics = ssf.getDynamics();
        pos = 0;
        end = data.getCount();
        if (initializer == null) {
            state = State.of(dynamics, StateInfo.Forecast);
        } else {
            State initial = new State(dynamics.getStateDim());
            pos = initializer.initialize(initial, ssf, data);
            if (pos >= 0) {
                state = initial;
            }
        }
        return state != null;
    }

    /**
     *
     * @param ssf
     * @param data
     * @param rslts
     * @return
     */
    public boolean process(final IMultivariateSsf ssf, final IMultivariateSsfData data, final IMultivariateFilteringResults rslts) {
        if (!initialize(ssf, data)) {
            return false;
        }
        if (rslts != null) {
            rslts.open(ssf, this.data);
        }
            pred();
        while (pos < end) {
            if (rslts != null) {
                rslts.save(pos, state);
            }
            error();
            if (rslts != null) {
                rslts.save(pos, perrors);
            }
            update();
            if (rslts != null) {
                rslts.save(pos, state);
            }
            pred();
             ++pos;
        }
        if (rslts != null) {
            rslts.close();
        }
        return true;
    }

}
