/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.dk;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.DataBlockIterator;
import ec.tstoolkit.data.DescriptiveStatistics;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.univariate.ISmoothingResults;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfData;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import ec.tstoolkit2.ssf.univariate.OrdinarySmoother;

/**
 *
 * @author Jean Palate
 */
public class DiffuseSmoother {

    private DiffuseState state;
    private ISsfDynamics dynamics;
    private ISsfMeasurement measurement;
    private ISmoothingResults srslts;
    private IDiffuseFilteringResults frslts;

    private double e, f, fi;
    private DataBlock C, Ci, Rf, Ri;
    private Matrix N0, N1, N2;
    private boolean missing, hasinfo, calcvar = true;
    private int pos;

    public boolean process(final ISsf ssf, final ISsfData data, ISmoothingResults sresults) {
        IDiffuseFilteringResults fresults = DkToolkit.filter(ssf, data, true);
        return process(ssf, data.getCount(), fresults, sresults);
    }

    public boolean process(ISsf ssf, final int endpos, IDiffuseFilteringResults results, ISmoothingResults sresults) {
        frslts = results;
        srslts = sresults;
        initFilter(ssf);
        initSmoother(ssf);
        ordinarySmoothing(ssf, endpos);
        pos = frslts.getEndDiffusePosition();
        while (--pos >= 0) {
            loadInfo();
            iterate();
            if (hasinfo) {
                srslts.save(pos, state);
            }
        }

        return true;
    }

    public ISmoothingResults getResults() {
        return srslts;
    }

    private void initSmoother(ISsf ssf) {
        int dim = ssf.getStateDim();
        state = new DiffuseState(dim);
        state.setInfo(StateInfo.Smoothed);

        Rf = new DataBlock(dim);
        C = new DataBlock(dim);
        Ri = new DataBlock(dim);
        Ci = new DataBlock(dim);

        if (calcvar) {
            N0 = Matrix.square(dim);
            N1 = Matrix.square(dim);
            N2 = Matrix.square(dim);
        }
    }

    private void loadInfo() {
        e = frslts.error(pos);
        f = frslts.errorVariance(pos);
        fi = frslts.diffuseNorm2(pos);
        C.copy(frslts.M(pos));
        Ci.copy(frslts.Mi(pos));
        if (fi != 0) {
            C.addAY(-f / fi, Ci);
        }
        missing = !DescriptiveStatistics.isFinite(e);
        DataBlock fa = frslts.a(pos);
        hasinfo = fa != null;
        if (!hasinfo) {
            return;
        }
        state.a().copy(fa);
        if (calcvar) {
            state.P().subMatrix().copy(frslts.P(pos));
            state.Pi().subMatrix().copy(frslts.Pi(pos));
        }
    }

    private void iterate() {
        iterateR();
        updateA();
        if (calcvar) {
            // P = P-PNP
            iterateN();
            updateP();
        }
    }
    // 

    private void updateA() {
        DataBlock a = state.a();
        if (calcvar) {
            a.addProduct(Rf, state.P().columns());
            a.addProduct(Ri, state.Pi().columns());
        } else { // to avoid unnecessary copies
            a.addProduct(Rf, frslts.P(pos).columns());
            a.addProduct(Ri, frslts.Pi(pos).columns());
        }
    }

    private void updateP() {

    }

    private void xL(DataBlock x) {
        // xL = x(T-KZ) = x(T-Tc/f*Z) = xT - ((xT)*c)/f * Z
        // compute xT
        dynamics.XT(pos, x);
        // compute q=xT*c
        double q = x.dot(C);
        // remove q/f*Z
        measurement.XpZd(pos, x, -q / f);
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
        if (fi == 0) {
            iterateRegularN();
        }
    }

    private void iterateRegularN() {
        if (!missing && f != 0) {
            // N(t-1) = Z'(t)*Z(t)/f(t) + L'(t)*N(t)*L(t)
            XL(N0.rows());
            XL(N0.columns());

            double cuc = SymmetricMatrix.quadraticForm(N0, C);

            // Compute V = C'U
            DataBlock v = new DataBlock(C.getLength());
            v.product(N0.columns(), C);

            DataBlockIterator columns = N0.columns();
            DataBlock col = columns.getData();
            DataBlockIterator rows = N0.rows();
            DataBlock row = rows.getData();
            int i = 0;
            do {
                double k = v.get(i++);
                if (k != 0) {
                    measurement.XpZd(pos, row, -k);
                    measurement.XpZd(pos, col, -k);
                }
            } while (rows.next() && columns.next());

            measurement.VpZdZ(pos, N0.subMatrix(), 1 / f);
            SymmetricMatrix.reinforceSymmetry(N0);
        } else {
            //T'*N(t)*T
            tvt(N0);
            tvt(N1);
            tvt(N2);
        }
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
        SymmetricMatrix.reinforceSymmetry(N);

    }

    /**
     *
     */
    private void iterateR() {
        if (fi == 0) {
            iterateRegularR();
        } else {
            iterateDiffuseR();
        }
    }

    private void iterateRegularR() {
        // R(t-1)=v(t)/f(t)Z(t)+R(t)L(t)
        //   = v/f*Z + R*(T-TC/f*Z)
        //  = (v - RT*C)/f*Z + RT
        dynamics.XT(pos, Rf);
        dynamics.XT(pos, Ri);
        if (!missing && f != 0) {
            // RT
            double c = (e - Rf.dot(C)) / f;
            measurement.XpZd(pos, Rf, c);
        }
    }

    private void iterateDiffuseR() {
        dynamics.XT(pos, Rf);
        dynamics.XT(pos, Ri);
        if (!missing && f != 0) {
            // Ri(t-1)=c*Z(t) +Ri(t)*T(t)
            // c = e/fi-(Ri(t)*T(t)*Ci(t))/fi-(Rf(t)*T(t)*Cf(t))/f
            double ci = (e - Ri.dot(Ci) - Rf.dot(C)) / fi;
            measurement.XpZd(pos, Ri, ci);
            // Rf(t-1)=c*Z(t)+Rf(t)*T(t)
            // c =  - Rf(t)T(t)*Ci/fi
            double cf = -Rf.dot(Ci) / fi;
            measurement.XpZd(pos, Rf, cf);
        }
    }

    private void initFilter(ISsf ssf) {
        dynamics = ssf.getDynamics();
        measurement = ssf.getMeasurement();
    }

    public void setCalcVariances(boolean b) {
        calcvar = false;
    }

    public boolean isCalcVariances() {
        return calcvar;
    }

    private void ordinarySmoothing(ISsf ssf, final int endpos) {
        OrdinarySmoother smoother = new OrdinarySmoother();
        smoother.setCalcVariances(calcvar);
        smoother.process(ssf, frslts.getEndDiffusePosition(), endpos, frslts, srslts);
        // updates R, N
        Rf.copy(smoother.getFinalR());
        if (calcvar) {
            N0.copy(smoother.getFinalN());
        }
    }

}
