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
/*
 */
package ec.tstoolkit2.ssf.akf;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.DataBlockIterator;
import ec.tstoolkit.maths.matrices.LowerTriangularMatrix;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.multivariate.IMultivariateSsf;
import ec.tstoolkit2.ssf.multivariate.IMultivariateSsfData;
import ec.tstoolkit2.ssf.multivariate.ISsfMeasurements;

/**
 *
 * @author Jean Palate
 */
public class MultivariateAugmentedFilter {

    private AugmentedState state;
    private AugmentedPredictionErrors perrors;
    private ISsfMeasurements measurements;
    private ISsfDynamics dynamics;
    private IMultivariateSsfData data;
    private int pos, end;
    private final boolean collapsing;
    private int collapsingPos = -1;

    /**
     *
     */
    public MultivariateAugmentedFilter() {
        collapsing = false;
    }

    public MultivariateAugmentedFilter(final boolean collapsing) {
        this.collapsing = collapsing;
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
     * P(t+1|t) = T(t)P(t|t)T'(t) The same transformation is applied on state.B
     * (diffuse constraints)
     */
    protected void pred() {
        if (state.getInfo() != StateInfo.Forecast) {
            state.setInfo(StateInfo.Forecast);
            SubMatrix P = state.P().subMatrix();
            DataBlock a = state.a();
            dynamics.TX(pos, a);
            dynamics.TM(pos, state.B());
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
    private boolean error() {
        int nobs = countMeasurements();
        if (nobs == 0) {
            perrors = null;
            return false;
        } else {
            perrors = new AugmentedPredictionErrors(dynamics.getStateDim(), nobs, dynamics.getNonStationaryDim());
            Matrix L = perrors.getCholeskyFactor();
            // K = PZ'(ZPZ'+H)^-1/2
            // computes (ZP)' in K'. 
            // Z~v x r, P~r x r, K~r x v
            // F = ZPZ'+H ~ v x v
            // L = F^-1/2 ~ v x v
            SubMatrix F = L.subMatrix(), K = perrors.getK().subMatrix();
            // K' = ZP or K = PZ'
            ZM(state.P().subMatrix(), K.transpose());
            // computes ZPZ'; results in L
            ZM(K, F);
            addH(L);
            // to avoid numerical problems
            SymmetricMatrix.reinforceSymmetry(L);

            // L contains now the Cholesky factor !!!
            SymmetricMatrix.lcholesky(L, State.ZERO);

            // We put in K  PZ'*(ZPZ'+H)^-1/2 = PZ'* L'^-1
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
            // U = e*L'^-1 or U L' = e or L*U' = e'
            LowerTriangularMatrix.rsolve(L, U, State.ZERO);
            Matrix E = perrors.E();
            // E is ndiffuse x nobs. Each column contains the diffuse effects
            // on the corresponding variable
            ZM(state.B(), E.subMatrix().transpose());
            E.chs();
            DataBlockIterator erows = E.rows();
            DataBlock erow = erows.getData();
//            Matrix B = state.B();
//            do {
//                for (int i = 0, j = 0; i < measurements.getCount(pos); ++i) {
//                    if (!data.isMissing(pos, i)) {
//                        erow.set(j, -measurements.ZX(pos, i, B.column(erows.getPosition())));
//                        ++j;
//                    }
//                }
//            } while (erows.next());
//            erows.begin();
            do {
                LowerTriangularMatrix.rsolve(L, erow, State.ZERO);
            } while (erows.next());
            return true;
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
        Matrix P = state.P();
        DataBlock U = perrors.getTransformedPredictionErrors();
        Matrix K = perrors.getK();
        for (int i = 0; i < n; ++i) {
//            for (int j = 0; j < n; ++j) {
//                P.addXaYt(-1, K.column(i), K.column(j));//, state_.K.column(i));
//            }
            state.a().addAY(U.get(i), K.column(i));
            P.addXaXt(-1, K.column(i));//, state_.K.column(i));
        }
        DataBlockIterator acols = state.B().columns();
        DataBlock acol = acols.getData();
        do {
            DataBlock row = perrors.E().row(acols.getPosition());
            for (int i = 0; i < n; ++i) {
                acol.addAY(row.get(i), K.column(i));
            }
        } while (acols.next());
    }

    /**
     *
     * @return
     */
    public AugmentedState getState() {
        return state;
    }

    public int getCollapsingPosition() {
        return collapsingPos;
    }

    private boolean initFilter() {
        pos = 0;
        end = data.getCount();
        return true;
    }

    private boolean initState() {
        pos = 0;
        state = AugmentedState.of(dynamics, StateInfo.Forecast);
        return (state != null);
    }

    /**
     *
     * @param ssf
     * @param data
     * @param rslts
     * @return
     */
    public boolean process(final IMultivariateSsf ssf, final IMultivariateSsfData data, final IMultivariateAugmentedFilteringResults rslts) {
        measurements = ssf.getMeasurements();
        dynamics = ssf.getDynamics();
        this.data = data;
        if (!initFilter()) {
            return false;
        }
        if (!initState()) {
            return false;
        }
        rslts.open(ssf, data);
        while (pos < end) {
            pred();
            if (collapse(rslts)) {
                break;
            }
            if (error()) {
                rslts.save(pos, perrors);
                update();
            } else {
                state.setInfo(StateInfo.Concurrent);
            }
            ++pos;
        }
        rslts.close();
        return true;
    }

    protected boolean collapse(IMultivariateAugmentedFilteringResults decomp) {
        if (!collapsing) {
            return false;
        }
        // update the state vector
        if (!decomp.collapse(state)) {
            return false;
        }
        collapsingPos = pos;
        return true;
    }

}
