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
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfData;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;

/**
 *
 * @author Jean Palate
 */
public class AugmentedFilter {

    private AugmentedState state;
    private AugmentedPredictionError pe;
    private ISsfMeasurement measurement;
    private ISsfDynamics dynamics;
    private ISsfData data;
    private int pos, end;
    private boolean missing;
    private final boolean collapsing;
    private int collapsingPos = -1;

    /**
     *
     */
    public AugmentedFilter() {
        collapsing = false;
    }

    public AugmentedFilter(final boolean collapsing) {
        this.collapsing = collapsing;
    }

    /**
     * Computes a(t+1|t), P(t+1|t) from a(t|t), P(t|t)
     */
    private void pred() {
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

    protected boolean error() {
        missing = data.isMissing(pos);
        if (missing) {
            pe.E().set(0);
            pe.M().set(0);
            // pe_ = null;
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
            if (v < State.ZERO) {
                v = 0;
            }
            pe.setVariance(v);
            // We put in K  PZ'*(ZPZ'+H)^-1 = PZ'* F^-1 = PZ'*(LL')^-1/2 = PZ'(L')^-1
            // K L' = PZ' or L K' = ZP

            double y = data.get(pos);
            pe.set(y - measurement.ZX(pos, state.a()));
            measurement.ZM(pos, state.B(), pe.E());
            pe.E().chs();
            return true;
        }
    }

    protected void update() {
        state.setInfo(StateInfo.Concurrent);
        double v = pe.getVariance(), e = pe.get();
        // P = P - (M)* F^-1 *(M)' --> Symmetric
        // PZ'(LL')^-1 ZP' =PZ'L'^-1*L^-1*ZP'
        // a = a + (M)* F^-1 * v
        state.a().addAY(e / v, pe.M());
        DataBlockIterator acols = state.B().columns();
        DataBlock acol = acols.getData();
        do {
            acol.addAY(pe.E().get(acols.getPosition()) / v, pe.M());
        } while (acols.next());
        update(state.P(), v, pe.M());
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
        if (state == null) {
            return false;
        }
        pe = new AugmentedPredictionError(dynamics.getStateDim(), dynamics.getNonStationaryDim());
        return true;
    }

    /**
     *
     * @param ssf
     * @param data
     * @param rslts
     * @return
     */
    public boolean process(final ISsf ssf, final ISsfData data, final IAugmentedFilteringResults rslts) {
        measurement = ssf.getMeasurement();
        dynamics = ssf.getDynamics();
        this.data = data;
        if (!initFilter()) {
            return false;
        }
        if (!initState()) {
            return false;
        }
        pred();
        while (pos < end) {
            if (rslts != null) {
                rslts.save(pos, state);
            }
            if (collapse(rslts)) {
                break;
            }
            if (error()) {
                if (rslts != null) {
                    rslts.save(pos, pe);
                }
                update();
            } else {
                this.state.setInfo(StateInfo.Concurrent);
            }
            if (rslts != null) {
                rslts.save(pos, state);
            }
            pred();
            ++pos;
        }
        return true;
    }

    // P -= c*r
    private void update(Matrix P, double v, DataBlock C) {//, DataBlock r) {
        SymmetricMatrix.addXaXt(P, -1 / v, C);
    }

    protected boolean collapse(IAugmentedFilteringResults decomp) {
        if (!collapsing) {
            return false;
        }
        // update the state vector
        if (!decomp.collapse(pos, state)) {
            return false;
        }
        collapsingPos = pos;
        return true;
    }

}
