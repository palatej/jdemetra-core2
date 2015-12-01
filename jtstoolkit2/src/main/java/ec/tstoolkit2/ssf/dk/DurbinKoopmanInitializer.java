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
package ec.tstoolkit2.ssf.dk;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.design.Development;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.SsfException;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfData;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import ec.tstoolkit2.ssf.univariate.OrdinaryFilter;

/**
 *
 * @author Jean Palate
 */
@Development(status = Development.Status.Alpha)
public class DurbinKoopmanInitializer implements OrdinaryFilter.Initializer {

    private final IDiffuseFilteringResults results;
    private DiffuseState state;
    private DiffusePredictionError pe;
    private ISsfMeasurement measurement;
    private ISsfDynamics dynamics;
    private ISsfData data;
    private int pos;
    private double norm = 0;

    public DurbinKoopmanInitializer() {
        this.results = null;
    }

    /**
     *
     * @param results
     */
    public DurbinKoopmanInitializer(IDiffuseFilteringResults results) {
        this.results = results;
    }

    /**
     * Computes a(t+1|t), P(t+1|t) from a(t|t), P(t|t) a(t+1|t) = T(t)a(t|t)
     * P(t+1|t) = T(t)P(t|t)T'(t) + V(t)
     */
    protected void pred() {
        if (state.getInfo() != StateInfo.Forecast) {
            state.setInfo(StateInfo.Forecast);
            SubMatrix P = state.P().subMatrix();
            DataBlock a = state.a();
            dynamics.TX(pos, a);
            dynamics.TVT(pos, P);
            dynamics.addV(pos, P);
            dynamics.TVT(pos, state.Pi().subMatrix());
        }
    }

    /**
     * Computes: e(t)=y(t) - Z(t)a(t|t-1)) F(t)=Z(t)P(t|t-1)Z'(t)+H(t) F(t) =
     * L(t)L'(t) E(t) = e(t)L'(t)^-1 K(t)= P(t|t-1)Z'(t)L'(t)^-1
     *
     * Not computed for missing values
     *
     * @return false if it has not been computed (missing value), true otherwise
     */
    protected boolean error() {
        if (state.getInfo() != StateInfo.Forecast) {
            throw new SsfException(SsfException.STATUS);
        }
        // calc f and fi
        // fi = Z Pi Z' , f = Z P Z' + H
        double fi = measurement.ZVZ(pos, state.Pi().subMatrix());
        if (Math.abs(fi) < State.ZERO) {
            fi = 0;
        }
        pe.setDiffuseNorm2(fi);
        double f = measurement.ZVZ(pos, state.P().subMatrix());
        if (measurement.hasErrors()) {
            f += measurement.errorVariance(pos);
        }
        if (Math.abs(f) / norm < State.ZERO) {
            f = 0;
        }
        pe.setVariance(f);
        if (data.hasData()) {
            double y = data.get(pos);
            if (Double.isNaN(y)) {
                pe.setMissing();
                return false;
            } else {
                pe.set(y - measurement.ZX(pos, state.a()));
            }
        }
        measurement.ZM(pos, state.P().subMatrix(), pe.M());
        if (pe.isDiffuse()) {
            measurement.ZM(pos, state.Pi().subMatrix(), pe.Mi());
        }
        return true;
    }

    /**
     *
     * @param fstate
     * @param ssf
     * @param data
     * @return
     */
    @Override
    public int initialize(final State fstate, final ISsf ssf, final ISsfData data) {
        measurement = ssf.getMeasurement();
        dynamics = ssf.getDynamics();
        this.data = data;
        int r = ssf.getStateDim();
        pos = 0;
        int end = data.getCount();
        if (!initState()) {
            return -1;
        }
        pred();
        while (pos < end) {
            if (isZero(this.state.Pi())) {
                break;
            }
            if (results != null) {
                results.save(pos, state);
            }
            if (error()) {
                if (results != null) {
                    results.save(pos, pe);
                }
                update();
            } else {
                this.state.setInfo(StateInfo.Concurrent);
            }
            if (results != null) {
                results.save(pos, state);
            }
            pred();
            ++pos;
        }
        if (results != null) {
            results.close(pos);
        }
        fstate.setInfo(StateInfo.Forecast);
        fstate.P().copy(state.P());
        fstate.a().copy(state.a());
        return pos;
    }

    private boolean initState() {
        state = DiffuseState.of(dynamics, StateInfo.Forecast);
        if (state == null) {
            return false;
        }
        norm = state.Pi().nrm2();
        pe = new DiffusePredictionError(dynamics.getStateDim());
        return true;
    }

    private boolean isZero(final Matrix P) {
        return P.isZero(1e-9 * norm);
    }

    private void update() {
        if (pe.isDiffuse()) {
            update1();
        } else {
            update0();
        }
        state.setInfo(StateInfo.Concurrent);
    }

    private void update0() {
        // variance

        double f = pe.getVariance(), e = pe.get();
        DataBlock C = pe.M();
        SymmetricMatrix.addXaXt(state.P(), -1 / f, C);

        // state
        // a0 = Ta0 + f1*TMi*v0. Reuse Mf as temporary buffer
        if (data.hasData()) {
            // prod(n, m_T, m_a0, m_tmp);
            double c = e / f;
//            for (int i = 0; i < m_r; ++i)
//                state.A.set(i, state.A.get(i) + state.C.get(i) * c);
            state.a().addAY(c, C);
        }
    }

    private void update1() {
        // calc f0, f1, f2
//        double f1 = 1 / pe.fi;
//        double f2 = -pe.f * f1 * f1;
        double f = pe.getVariance(), e = pe.get(), fi = pe.getDiffuseNorm2();
        DataBlock C = pe.M(), Ci = pe.Mi();

        // Pi = Pi - f1* (Ci)(Ci)'
        SymmetricMatrix.addXaXt(state.Pi(), -1 / fi, Ci);

        // P = P - f2*(Ci)(Ci)'-f1(Ci*Cf' + Cf*Ci')
        // = P + f/(fi*fi)(Ci)(Ci)' - 1/fi(Ci*Cf' + Cf*Ci')
        // = P - 1/f (Cf)(Cf') + f/(fi*fi)(Ci)(Ci)'- 1/fi(Ci*Cf' + Cf*Ci')+ 1/f (Cf)(Cf')
        // = P  - 1/f (Cf)(Cf') + (1/f)(Cf - (f/fi)Ci)(Cf - (f/fi)Ci)'
        SymmetricMatrix.addXaXt(state.P(), -1 / f, C);
        DataBlock tmp = C.deepClone();
        tmp.addAY(-f / fi, Ci);
        SymmetricMatrix.addXaXt(state.P(), 1 / f, tmp);

        if (data.hasData()) {
            // a0 = Ta0 + f1*TMi*v0. 
            state.a().addAY(e / fi, Ci);
        }
    }

}
