/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
public class DisturbanceSmoother {

    private ISsfDynamics dynamics;
    private ISsfMeasurement measurement;
    private IDisturbanceSmoothingResults srslts;
    private IFilteringResults frslts;

    private double err, errVariance, esm, esmVariance, h;
    private DataBlock K, R, U;
    private Matrix N, UVar, S, Q;
    private boolean missing, res, calcvar = true;
    private int pos, stop;
    // temporary
    private DataBlock tmp;
    private double c, v;

    public boolean process(ISsf ssf, ISsfData data) {
        if (ssf.getDynamics().isDiffuse()) {
            return false;
        }
        OrdinaryFilter filter = new OrdinaryFilter();
        DefaultFilteringResults fresults = DefaultFilteringResults.light();
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
        IDisturbanceSmoothingResults sresults;
        boolean hasErrors = ssf.getMeasurement().hasErrors();
        if (calcvar) {
            sresults = DefaultDisturbanceSmoothingResults.full(hasErrors);
        } else {
            sresults = DefaultDisturbanceSmoothingResults.light(hasErrors);
        }

        return process(ssf, start, end, results, sresults);
    }

    public boolean process(ISsf ssf, final int start, final int end, IFilteringResults results, IDisturbanceSmoothingResults sresults) {
        frslts = results;
        srslts = sresults;
        stop = start;
        pos = end;
        initFilter(ssf);
        initSmoother(ssf);
        while (--pos >= stop) {
            loadInfo();
            if (iterate()) {
                srslts.saveSmoothedTransitionDisturbances(pos, U, UVar.subMatrix());
                if (res) {
                    srslts.saveSmoothedMeasurementDisturbance(pos, esm, esmVariance);
                }
            }
        }
        return true;
    }

    public boolean resume(final int start) {
        stop = start;
        while (pos >= stop) {
            loadInfo();
            if (iterate()) {
                srslts.saveSmoothedTransitionDisturbances(pos, U, UVar.subMatrix());
                if (res) {
                    srslts.saveSmoothedMeasurementDisturbance(pos, esm, esmVariance);
                }
            }
            pos--;
        }
        return true;
    }

    public IDisturbanceSmoothingResults getResults() {
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
        int resdim = dynamics.getInnovationsDim();

        R = new DataBlock(dim);
        K = new DataBlock(dim);
        U = new DataBlock(resdim);
        S = new Matrix(dim, resdim);
        if (calcvar) {
            N = Matrix.square(dim);
            tmp = new DataBlock(dim);
            UVar = Matrix.square(resdim);
            Q = Matrix.square(resdim);
            if (measurement.isTimeInvariant()) {
                h = measurement.errorVariance(0);
            }
        }
        if (dynamics.isTimeInvariant()) {
            dynamics.S(0, S.subMatrix());
            if (Q != null) {
                dynamics.Q(0, Q.subMatrix());
            }
        }
    }

    private void loadInfo() {
        err = frslts.error(pos);
        errVariance = frslts.errorVariance(pos);
        K.setAY(1 / errVariance, frslts.M(pos));
        dynamics.TX(pos, K);
        missing = !DescriptiveStatistics.isFinite(err);
        if (!dynamics.isTimeInvariant()) {
            dynamics.S(pos, S.subMatrix());
            if (Q != null) {
                dynamics.Q(pos, Q.subMatrix());
            }
        }
        if (!measurement.isTimeInvariant()) {
            h = measurement.errorVariance(pos);
        }
    }

    private boolean iterate() {
        iterateR();
        if (calcvar) {
            iterateN();
        }
        if (!missing) {
            // updates the smoothed disturbances
            if (res) {
                esm = c * h;
            }
            U.product(R, S.columns());
            if (calcvar) {
                if (res) {
                    esmVariance = h - h * h * v;
                }
                // v(U) = Q-S'NS
                UVar.copy(Q);
                Matrix V = SymmetricMatrix.quadraticForm(N, S);
                UVar.sub(V);
            }
            return true;
        } else {
            return false;
        }
    }
    // 

    /**
     *
     */
    private void iterateN() {
        if (!missing && errVariance != 0) {
            // N(t-1) = Z'(t)*Z(t)/f(t) + L'(t)*N(t)*L(t)
            // = Z'(t)*Z(t)/f(t) + (T' - Z'K')N(T - KZ)
            // =  Z'(t)*Z(t)(1/f(t) + K'NK) + T'NT - <T'NKZ>
            // 1. NK 
            tmp.product(N.rows(), K);
            // 2. v
            v = 1 / errVariance + tmp.dot(K);
            // 3. T'NK
            dynamics.XT(pos, tmp);
            // TNT
            tvt(N);
            measurement.VpZdZ(pos, N.subMatrix(), v);
            subZ(N.rows(), tmp);
            subZ(N.columns(), tmp);
        } else {
            tvt(N);
        }
        SymmetricMatrix.reinforceSymmetry(N);
    }

    /**
     *
     */
    private void iterateR() {
        // R(t-1)=(v/f + R(t)*K)Z + R(t)*T
        // R(t-1)=esm*Z +  R(t)*T
        if (!missing && errVariance != 0) {
            // RT
            c = (err / errVariance - R.dot(K));
            dynamics.XT(pos, R);
            measurement.XpZd(pos, R, c);
        } else {
            dynamics.XT(pos, R);
            esm = Double.NaN;
        }
    }

    private void initFilter(ISsf ssf) {
        dynamics = ssf.getDynamics();
        measurement = ssf.getMeasurement();
        res = measurement.hasErrors();
    }

    public void setCalcVariances(boolean b) {
        calcvar = b;
    }

    public boolean isCalcVariances() {
        return calcvar;
    }

    private void tvt(Matrix N) {
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

    }

    private void subZ(DataBlockIterator rows, DataBlock b) {
        DataBlock row = rows.getData();
        do {
            double cur = b.get(rows.getPosition());
            if (cur != 0) {
                measurement.XpZd(pos, row, -cur);
            }
        } while (rows.next());
    }

}
