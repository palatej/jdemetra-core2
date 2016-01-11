/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.dk;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.DataBlockIterator;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.univariate.ISmoothingResults;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;

/**
 *
 * @author Jean Palate
 */
public abstract class BaseDiffuseSmoother {

    protected ISsfDynamics dynamics;
    protected ISsfMeasurement measurement;
    protected ISmoothingResults srslts;

    protected double e, f, fi;
    protected DataBlock C, Ci, Rf, Ri, tmp0, tmp1;
    protected Matrix N0, N1, N2;
    protected boolean missing, hasinfo, calcvar = true;
    protected int pos;

    public ISmoothingResults getResults() {
        return srslts;
    }

    protected void iterate() {
        iterateR();
        updateA();
        if (calcvar) {
            // P = P-PNP
            iterateN();
            updateP();
        }
    }
    // 

    protected abstract void updateA();

    protected abstract void updateP();
    /**
     * Computes in place x = x-c/f*z
     *
     * @param x
     * @param k
     */
    private void xQ(DataBlock x) {
        measurement.XpZd(pos, x, -x.dot(C));
    }

    private void XQ(DataBlockIterator X) {
        DataBlock x = X.getData();
        do {
            xQ(x);
        } while (X.next());
    }

    private void xQi(DataBlock x) {
        measurement.XpZd(pos, x, -x.dot(Ci));
    }

    private void XQi(DataBlockIterator X) {
        DataBlock x = X.getData();
        do {
            xQi(x);
        } while (X.next());
    }

    /**
     *
     */
    private void iterateN() {
        if (missing || (f == 0 && fi == 0)) {
            iterateMissingN();
        } else if (fi == 0) {
            iterateRegularN();
        } else {
            iterateDiffuseN();
        }
    }

    private void iterateMissingN() {
        tvt(N0);
        tvt(N1);
        tvt(N2);
        // reinforceSymmetry();
    }

    private void iterateRegularN() {
        // N(t-1) = Z'(t)*Z(t)/f(t) + L'(t)*N(t)*L(t)
        tvt(N0);
        XQ(N0.rows());
        XQ(N0.columns());
        measurement.VpZdZ(pos, N0.subMatrix(), 1 / f);
        tvt(N1);
        XQ(N1.columns());
        tvt(N2);
    }

    private void iterateDiffuseN() {
        // Nf = Li'*Nf*Li
        // N1 = Z'Z/Fi + Li'*N1*Li - < Z'Kf'*Nf'*Li >
        // N2 = Z'Z * c + Li'*N2*Li - < Z'Kf'*N1'*Li >, c= Kf'*Nf*Kf-Ff/(Fi*Fi)
        // compute first N2 then N1 and finally Nf
        tvt(N0);
        tvt(N1);
        tvt(N2);

        tmp0.product(C, N0.columns());
        tmp1.product(C, N1.columns());

        double kn0k = tmp0.dot(C);

        XQi(N0.rows());
        XQi(N0.columns());
        XQi(N1.rows());
        XQi(N2.columns());
        XQi(N2.rows());
        XQi(N1.columns());
        xQi(tmp0);
        xQi(tmp1);

        measurement.VpZdZ(pos, N1.subMatrix(), 1 / fi); //
        measurement.VpZdZ(pos, N2.subMatrix(), kn0k - f / (fi * fi));

        subZ(N1.rows(), tmp0);
        subZ(N1.columns(), tmp0);
        subZ(N2.rows(), tmp1);
        subZ(N2.columns(), tmp1);
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
            double c = e / f - Rf.dot(C);
            measurement.XpZd(pos, Rf, c);
        }
    }

    private void iterateDiffuseR() {
        dynamics.XT(pos, Rf);
        dynamics.XT(pos, Ri);
        if (!missing && f != 0) {
            // Ri(t-1)=c*Z(t) +Ri(t)*T(t)
            // c = e/fi-(Ri(t)*T(t)*Ci(t))/fi-(Rf(t)*T(t)*Cf(t))/f
            double ci = e / fi - Ri.dot(Ci) - Rf.dot(C);
            measurement.XpZd(pos, Ri, ci);
            // Rf(t-1)=c*Z(t)+Rf(t)*T(t)
            // c =  - Rf(t)T(t)*Ci/fi
            double cf = -Rf.dot(Ci);
            measurement.XpZd(pos, Rf, cf);
        }
    }

    protected void initFilter(ISsf ssf) {
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
