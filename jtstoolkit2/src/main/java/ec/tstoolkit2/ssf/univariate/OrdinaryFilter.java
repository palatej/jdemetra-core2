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
package ec.tstoolkit2.ssf.univariate;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;

/**
 * Ordinary Kalman filter for univariate time series
 *
 * @author Jean Palate
 */
public class OrdinaryFilter {

    public static interface Initializer {

        int initialize(State state, ISsf ssf, ISsfData data);
    }

    private final Initializer initializer;
    private State state;
    private PredictionError pe;
    private ISsfMeasurement measurement;
    private ISsfDynamics dynamics;
    private int pos, end;
    private boolean missing;

    /**
     *
     * @param initializer
     */
    public OrdinaryFilter(Initializer initializer) {
        this.initializer = initializer;
    }

    public OrdinaryFilter() {
        this.initializer = null;
    }

    /**
     * Computes a(t+1|t), P(t+1|t) from a(t|t), P(t|t)
     */
    protected void pred() {
        if (state.getInfo() != StateInfo.Forecast) {
            SubMatrix P = state.P().subMatrix();
            DataBlock a = state.a();
            dynamics.TX(pos, a);
            dynamics.TVT(pos, P);
            dynamics.addV(pos, P);
            state.setInfo(StateInfo.Forecast);
        }
    }


    protected boolean error(ISsfData data) {
        missing = data.isMissing(pos);
        if (missing) {
            // pe_ = null;
            pe.setMissing();
            return false;
        } else {
            // pe_ = new PredictionError(ssf_.getStateDim(), 1);
            // K = PZ'/f
            // computes (ZP)' in K'. Missing values are set to 0 
            // Z~v x r, P~r x r, K~r x v
            DataBlock C = pe.M();
            // computes ZPZ'; results in pe_.L
            //measurement.ZVZ(pos_, state_.P.subMatrix(), F);
            measurement.ZM(pos, state.P().subMatrix(), C);
            double v = measurement.ZX(pos, C);
            if (measurement.hasErrors()) {
                v += measurement.errorVariance(pos);
            }
            pe.setVariance(v);
            // We put in K  PZ'*(ZPZ'+H)^-1 = PZ'* F^-1 = PZ'*(LL')^-1/2 = PZ'(L')^-1
            // K L' = PZ' or L K' = ZP

            double y = data.get(pos);
            pe.set(y - measurement.ZX(pos, state.a()));
            return true;
        }
    }

    protected void update() {
        state.setInfo(StateInfo.Concurrent);
        if (pe == null) {
            return;
        }
        double e = pe.get();
        double v = pe.getVariance();
        DataBlock C = pe.M();

        // P = P - (M)* F^-1 *(M)' --> Symmetric
        // PZ'(LL')^-1 ZP' =PZ'L'^-1*L^-1*ZP'
        // A = a + (M)* F^-1 * v
        state.a().addAY(e / v, C);
        update(state.P(), v, C);//, state_.K.column(i));
    }

    /**
     * Retrieves the final state (which is a(N|N-1))
     *
     * @return
     */
    public State getState() {
        return state;
    }

    private boolean initialize(ISsf ssf, ISsfData data) {
        measurement = ssf.getMeasurement();
        dynamics = ssf.getDynamics();
        pos = 0;
        end = data.getCount();
        pe = new PredictionError(dynamics.getStateDim());
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
    public boolean process(final ISsf ssf, final ISsfData data, final IFilteringResults rslts) {
        // intialize the state with a(0|-1)
        if (!initialize(ssf, data)) {
            return false;
        }
        pred();
        while (pos < end) {
            rslts.save(pos, state);
            if (error(data)) {
                rslts.save(pos, pe);
                update();
            } else {
                rslts.save(pos, pe);
                state.setInfo(StateInfo.Concurrent);
            }
            rslts.save(pos, state);
            pred();
            ++pos;
        }
        return true;
    }

    // P -= c*r
    private void update(Matrix P, double v, DataBlock C) {//, DataBlock r) {
        SymmetricMatrix.addXaXt(P, -1 / v, C);
    }

}
