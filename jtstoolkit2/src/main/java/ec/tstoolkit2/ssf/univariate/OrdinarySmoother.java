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
package ec.tstoolkit2.ssf.univariate;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.DataBlockIterator;
import ec.tstoolkit.data.DescriptiveStatistics;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.ResultsRange;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;

/**
 *
 * @author Jean Palate
 */
public class OrdinarySmoother {

    private State state;
    private ISsfDynamics dynamics;
    private ISsfMeasurement measurement;
    private ISmoothingResults srslts;
    private IFilteringResults frslts;

    private double err, errVariance;
    private DataBlock M, R;
    private Matrix N;
    private boolean missing, calcvar = true;
    private int pos, stop;

    public boolean process(ISsf ssf, ISsfData data) {
        if (ssf.getDynamics().isDiffuse()) {
            return false;
        }
        OrdinaryFilter filter = new OrdinaryFilter();
        DefaultFilteringResults fresults = DefaultFilteringResults.full();
        if (!filter.process(ssf, data, fresults)) {
            return false;
        }
        return process(ssf, 0, data.getCount(), fresults);
    }

    public boolean process(ISsf ssf, DefaultFilteringResults results) {
        if (ssf.getDynamics().isDiffuse()) {
            return false;
        }
        ResultsRange range = results.getRange();
        return process(ssf, range.getStart(), range.getEnd(), results);
    }

    public boolean process(ISsf ssf, int start, int end, IFilteringResults results) {
        ISmoothingResults sresults;
        if (calcvar) {
            sresults = DefaultSmoothingResults.full();
        } else {
            sresults = DefaultSmoothingResults.light();
        }

        return process(ssf, start, end, results, sresults);
    }

    public boolean process(ISsf ssf, final int start, final int end, IFilteringResults results, ISmoothingResults sresults) {
        frslts = results;
        srslts = sresults;
        stop = start;
        pos = end;
        initFilter(ssf);
        initSmoother(ssf);
        while (--pos >= stop) {
            loadInfo();
            if (iterate()) {
                srslts.save(pos, state);
            }
        }

        return true;
    }

    public boolean resume(final int start) {
        stop = start;
        while (pos >= stop) {
            loadInfo();
            if (iterate()) {
                srslts.save(pos, state);
            }
            pos--;
        }
        return true;
    }

    public ISmoothingResults getResults() {
        return srslts;
    }

    public DataBlock getFinalR() {
        return R;
    }

    public Matrix getFinalN() {
        return N;
    }

    private void initSmoother(ISsf ssf) {
        int dim = ssf.getStateDim();
        state = new State(dim);
        state.setInfo(StateInfo.Smoothed);

        R = new DataBlock(dim);
        M = new DataBlock(dim);

        if (calcvar) {
            N = Matrix.square(dim);
        }
    }

    private void loadInfo() {
        err = frslts.error(pos);
        errVariance = frslts.errorVariance(pos);
        M.copy(frslts.M(pos));
        missing = !DescriptiveStatistics.isFinite(err);
    }

    private boolean iterate() {
        iterateR();
        if (calcvar) {
            iterateN();
        }
        DataBlock fa = frslts.a(pos);
        SubMatrix fP = frslts.P(pos);
        if (fP == null) {
            return false;
        }
        // a = a + r*P
        DataBlock a = state.a();
        a.copy(fa);
        a.addProduct(R, fP.columns());
        if (calcvar) {
            // P = P-PNP
            Matrix P = state.P();
            P.subMatrix().copy(fP);
            Matrix V = SymmetricMatrix.quadraticForm(N, P);
            P.sub(V);
        }
        return true;
    }
    // 

    private void xL(DataBlock x) {
        // xL = x(T-KZ) = x(T-Tc/f*Z) = xT - ((xT)*c)/f * Z
        // compute xT
        dynamics.XT(pos, x);
        // compute q=xT*c
        double q = x.dot(M);
        // remove q/f*Z
        measurement.XpZd(pos, x, -q / errVariance);
    }

    private void XL(DataBlockIterator X) {
        DataBlock x = X.getData();
        do {
            xL(x);
        } while (X.next());
    }

    /**
     *
     */
    private void iterateN() {
        if (!missing && errVariance != 0) {
            // N(t-1) = Z'(t)*Z(t)/f(t) + L'(t)*N(t)*L(t)
            XL(N.rows());
            XL(N.columns());

            // Compute V = C'U
            DataBlock v = new DataBlock(M.getLength());
            v.product(N.columns(), M);

            DataBlockIterator columns = N.columns();
            DataBlock col = columns.getData();
            DataBlockIterator rows = N.rows();
            DataBlock row = rows.getData();
            int i = 0;
            do {
                double k = v.get(i++);
                if (k != 0) {
                    measurement.XpZd(pos, row, -k);
                    measurement.XpZd(pos, col, -k);
                }
            } while (rows.next() && columns.next());

            measurement.VpZdZ(pos, N.subMatrix(), 1 / errVariance);
            SymmetricMatrix.reinforceSymmetry(N);
        } else {
            //T'*N(t)*T
            DataBlockIterator columns = N.columns();
            DataBlock col = columns.getData();
            do {
                dynamics.XT(pos, col);
            } while (columns.next());
            DataBlockIterator rows = N.rows();
            DataBlock row = rows.getData();
            do {
                dynamics.XT(pos, row);
            } while (rows.next());
            SymmetricMatrix.reinforceSymmetry(N);
        }
    }

    /**
     *
     */
    private void iterateR() {
        // R(t-1)=v(t)/f(t)Z(t)+R(t)L(t)
        //   = v/f*Z + R*(T-TC/f*Z)
        //  = (v - RT*C)/f*Z + RT
        dynamics.XT(pos, R);
        if (!missing && errVariance != 0) {
            // RT
            double c = (err - R.dot(M)) / errVariance;
            measurement.XpZd(pos, R, c);
        }
    }

    private void initFilter(ISsf ssf) {
        dynamics = ssf.getDynamics();
        measurement = ssf.getMeasurement();
    }

    public void setCalcVariances(boolean b) {
        calcvar = b;
    }

    public boolean isCalcVariances() {
        return calcvar;
    }

}
