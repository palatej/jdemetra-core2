/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.dk.sqrt;

import ec.tstoolkit2.ssf.dk.*;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.DataBlockIterator;
import ec.tstoolkit.data.DescriptiveStatistics;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.akf.AugmentedState;
import ec.tstoolkit2.ssf.univariate.ISmoothingResults;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfData;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import ec.tstoolkit2.ssf.univariate.OrdinarySmoother;

/**
 *
 * @author Jean Palate
 */
public class DiffuseSquareRootSmoother {

    private AugmentedState state;
    private ISsfDynamics dynamics;
    private ISsfMeasurement measurement;
    private ISmoothingResults srslts;
    private IDiffuseSquareRootFilteringResults frslts;

    private double e, f, fi;
    private DataBlock C, Ci, Rf, Ri, tmp0, tmp1;
    private Matrix N0, N1, N2;
    private boolean missing, calcvar = true;
    private int pos;

    public boolean process(final ISsf ssf, final ISsfData data, ISmoothingResults sresults) {
        IDiffuseSquareRootFilteringResults fresults = DkToolkit.sqrtFilter(ssf, data, true);
        return process(ssf, data.getCount(), fresults, sresults);
    }

    public boolean process(ISsf ssf, final int endpos, IDiffuseSquareRootFilteringResults results, ISmoothingResults sresults) {
        frslts = results;
        srslts = sresults;
        initFilter(ssf);
        initSmoother(ssf);
        ordinarySmoothing(ssf, endpos);
        pos = frslts.getEndDiffusePosition();
        while (--pos >= 0) {
            loadInfo();
            iterate();
            srslts.save(pos, state);
        }

        return true;
    }

    public ISmoothingResults getResults() {
        return srslts;
    }

    private void initSmoother(ISsf ssf) {
        int dim = ssf.getStateDim();
        state = new AugmentedState(dim, ssf.getDynamics().getNonStationaryDim());
        state.setInfo(StateInfo.Smoothed);

        Rf = new DataBlock(dim);
        C = new DataBlock(dim);
        Ri = new DataBlock(dim);
        Ci = new DataBlock(dim);

        if (calcvar) {
            tmp0 = new DataBlock(dim);
            tmp1 = new DataBlock(dim);
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
        if (fi != 0) {
            Ci.copy(frslts.Mi(pos));
            Ci.mul(1 / fi);
            C.addAY(-f, Ci);
            C.mul(1 / fi);
        } else {
            C.mul(1 / f);
            Ci.set(0);
        }
        missing = !DescriptiveStatistics.isFinite(e);
        state.a().copy(frslts.a(pos));
        if (calcvar) {
            state.P().subMatrix().copy(frslts.P(pos));
            SubMatrix B = frslts.B(pos);
            state.restoreB(B);
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
            // Pi=B*B'
            SubMatrix B = state.B();
            DataBlock tmp = new DataBlock(state.getDiffuseDim());
            tmp.product(Ri, B.columns());
            a.addProduct(tmp, B.rows());
        } else { // to avoid unnecessary copies
            a.addProduct(Rf, frslts.P(pos).columns());
            SubMatrix B = frslts.B(pos);
            DataBlock tmp = new DataBlock(B.getColumnsCount());
            tmp.product(Ri, B.columns());
            a.addProduct(tmp, B.rows());
        }
    }

    private void updateP() {
        Matrix P = state.P();
        Matrix PN0P = SymmetricMatrix.quadraticForm(N0, P);
        SubMatrix B = state.B();
        Matrix BN2B = SymmetricMatrix.quadraticForm(N2.subMatrix(), B);
        Matrix PN2P=SymmetricMatrix.quadraticFormT(BN2B.subMatrix(), B);
        Matrix N1B= new Matrix(N1.getRowsCount(), B.getColumnsCount());
        N1B.subMatrix().product(N1.subMatrix(), B);
        Matrix PN1B=P.times(N1B);
        Matrix PN1P=Matrix.square(P.getRowsCount());
        PN1P.subMatrix().product(PN1B.subMatrix(), B.transpose());
        P.sub(PN0P);
        P.sub(PN2P);
        P.sub(PN1P);
        P.subMatrix().sub(PN1P.subMatrix().transpose());
        SymmetricMatrix.reinforceSymmetry(P);
     }

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
        XQi(N1.columns());
        XQi(N2.rows());
        XQi(N2.columns());
        xQi(tmp0);
        xQi(tmp1);

        measurement.VpZdZ(pos, N1.subMatrix(), 1 / fi);
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
            double c = e / f - Rf.dot(C);;
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
