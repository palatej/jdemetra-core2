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
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.univariate.IFilteringResults;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfData;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import ec.tstoolkit2.ssf.univariate.PredictionError;

/**
 * Chandrasekhar recursions
 *
 * @param <F>
 * @author Jean Palate
 */
@Development(status = Development.Status.Alpha)
public class FastFilter<F extends ISsf> {

    private int pos_, end_;
    private double eps_ = 1e-12;
    private double neps_;

    private PredictionError pe_;

    private final IFastInitializer<F> initializer_;
    private ISsfMeasurement m_;
    private ISsfDynamics dyn_;

    private ISsfData data_;

    private DataBlock a, l, k;
    private double f;
    private double[] L, K;
    private int steadypos_;

    /**
     *
     */
    public FastFilter() {
        initializer_ = new FastInitializer();
    }

    public FastFilter(IFastInitializer<F> initializer) {
        initializer_ = initializer;
    }
    
    public void setEpsilon(double eps){
        eps_=eps;
    }
    
    public double getEpsilon(){
        return eps_;
    }
    
    public int getSteadyStatePosition(){
        return steadypos_;
    }

    private boolean initialize(F ssf) {
        steadypos_ = -1;
        dyn_ = ssf.getDynamics();
        m_ = ssf.getMeasurement();
        FastState state = new FastState(dyn_.getStateDim());
        pos_ = 0;
        end_ = data_.getLength();
        if (!initializer_.initialize(ssf, state)) {
            return false;
        }
        int dim = dyn_.getStateDim();
        a = new DataBlock(dim);
        k = state.k.deepClone();
        l = state.l.deepClone();
        f = state.f;
        K = k.getData();
        L = l.getData();
        dyn_.a0(a, StateInfo.Forecast);
        neps_ = eps_ * f;
        pe_ = new PredictionError(dyn_.getStateDim());
        return true;
    }

    /**
     *
     * @param ssf
     * @param data
     * @param rslts
     * @return
     */
    public boolean process(final F ssf, final ISsfData data,
            final IFilteringResults rslts) {
        data_ = data;
        if (!initialize(ssf)) {
            return false;
        }
        while (pos_ < end_) {
            if (pos_ > 0) {
                next();
            }
            error();
            if (rslts != null) {
                rslts.save(pos_, pe_);
            }
            update();
            //
            ++pos_;
        }
        return true;
    }

    private void next() {
        if (steadypos_ >= 0) {
            return;
        }
	// K(i+1) = K(i) - T L(i) * (Z*L(i))/V(i)
        // L(i+1) = T L(i) - K(i) * (Z*L(i))/V(i)
        // F(i+1) = F(i) - (Z*L(i))^2/V(i)

        // ZLi, V(i+1)
        double zl = m_.ZX(pos_, l);

        // TL(i)
        dyn_.TX(pos_, l);

        if (Math.abs(zl) > neps_) {
            // C, L
            double zlv = zl / f;
            f -= zl * zlv;
            for (int i = 0; i < L.length; ++i) {
                double tl = L[i];
                L[i] -= K[i] * zlv;
                K[i] -= tl * zlv;
            }
        } else if (l.nrm2() < eps_) {
            steadypos_ = pos_;
        }

    }

    private void update() {
        dyn_.TX(pos_, a);
        a.addAY(pe_.get() / pe_.getVariance(), k);
    }

    private void error() {
        double y = data_.get(pos_);
        pe_.set(y - m_.ZX(pos_, a));
        pe_.setVariance(f);

    }
}
