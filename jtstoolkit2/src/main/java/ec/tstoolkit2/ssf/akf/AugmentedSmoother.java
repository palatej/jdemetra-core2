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
import ec.tstoolkit.data.DescriptiveStatistics;
import ec.tstoolkit.maths.matrices.LowerTriangularMatrix;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
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
public class AugmentedSmoother {

    private AugmentedState state;
    private ISsfDynamics dynamics;
    private ISsfMeasurement measurement;
    private ISmoothingResults srslts;
    private IAugmentedFilteringResults frslts;

    private double e, f;
    private DataBlock C, E, R;
    private Matrix N, Rd, U, V, RNA, S;
    private Matrix Psi;
    private DataBlock delta;
    private boolean missing, hasinfo, calcvar = true;
    private int pos;

    public boolean process(final ISsf ssf, final ISsfData data, ISmoothingResults sresults) {
        IAugmentedFilteringResults fresults = AkfToolkit.filter(ssf, data, true);
        return process(ssf, data.getLength(), fresults, sresults);
    }

    public boolean process(ISsf ssf, final int endpos, IAugmentedFilteringResults results, ISmoothingResults sresults) {
        frslts = results;
        srslts = sresults;
        initFilter(ssf);
        initSmoother(ssf);
        ordinarySmoothing(ssf, endpos);
        pos = frslts.getCollapsingPosition();
        calcSmoothedDiffuseEffects();
        while (--pos >= 0) {
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
        int nd = ssf.getDynamics().getNonStationaryDim();
        state = new AugmentedState(dim, nd);
        state.setInfo(StateInfo.Smoothed);

        R = new DataBlock(dim);
        C = new DataBlock(dim);
        E = new DataBlock(nd);
        Rd = new Matrix(dim, nd);
        U = new Matrix(dim, nd);

        if (calcvar) {
            N = Matrix.square(dim);
            V = new Matrix(dim, nd);
            RNA = new Matrix(dim, nd);
        }
    }

    private void loadInfo() {
        e = frslts.error(pos);
        f = frslts.errorVariance(pos);
        E.copy(frslts.E(pos));
        C.copy(frslts.M(pos));
        missing = !DescriptiveStatistics.isFinite(e);
        DataBlock fa = frslts.a(pos);
        hasinfo = fa != null;
        if (!hasinfo) {
            return;
        }
        state.a().copy(fa);
        if (calcvar) {
            state.restoreB(frslts.B(pos));
            state.P().subMatrix().copy(frslts.P(pos));
        }
    }

    private void iterate() {
        loadInfo();
        iterateR();
        calcU();
        updateA();
        if (calcvar) {
            // P = P-PNP
            iterateN();
            calcV();
            updateP();
        }
    }

    // 
    private void calcU() {
        // U = A + PR
        DataBlockIterator columns = U.columns();
        DataBlock uc = columns.getData();
        DataBlockIterator rcolumns = Rd.columns();
        DataBlock rc = rcolumns.getData();
        DataBlockIterator acolumns = calcvar ? state.B().columns() : frslts.B(pos).columns();
        DataBlock ac = acolumns.getData();
        DataBlockIterator prows = calcvar ? state.P().rows() : frslts.P(pos).rows();
        do {
            prows.begin();
            uc.product(prows, rc);
            uc.add(ac);
        } while (columns.next() && rcolumns.next() && acolumns.next());
    }

    private void calcV() {
        // V =  PR + PNA = P(R+NA)

        // RNA = R + NA
        DataBlockIterator rnacolumns = RNA.columns();
        DataBlock rnac = rnacolumns.getData();
        DataBlockIterator rcolumns = Rd.columns();
        DataBlock rc = rcolumns.getData();
        DataBlockIterator acolumns = state.B().columns();
        DataBlock ac = acolumns.getData();
        do {
            rnac.product(N.rows(), ac);
            rnac.add(rc);
        } while (rnacolumns.next() && rcolumns.next() && acolumns.next());

        DataBlockIterator columns = V.columns();
        rnacolumns.begin();
        DataBlock vc = columns.getData();
        do {
            vc.product(state.P().rows(), rnac);
        } while (columns.next() && rnacolumns.next());

    }

    private void updateA() {
        DataBlock a = state.a();
        // normal iteration
        a.addProduct(R, calcvar ? state.P().columns() : frslts.P(pos).columns());
        // diffuse correction
        a.addProduct(U.rows(), delta);
    }

    private void updateP() {
        Matrix P = state.P();
        // normal iteration
        Matrix PNP = SymmetricMatrix.quadraticForm(N, P);
        P.sub(PNP);
        // diffuse correction
        Matrix UPsiU = SymmetricMatrix.quadraticFormT(Psi, U);
        P.add(UPsiU);
        SubMatrix u = U.subMatrix(), vt=V.subMatrix().transpose();
        LowerTriangularMatrix.rsolve(S, u.transpose());
        LowerTriangularMatrix.rsolve(S, vt);
        // compute U*V'
        Matrix UV=Matrix.square(U.getRowsCount());
        SubMatrix uv=UV.subMatrix();
        uv.product(u, vt);
        P.subMatrix().sub(uv);
        P.subMatrix().sub(uv.transpose());
        SymmetricMatrix.reinforceSymmetry(P);
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
        if (!missing && f != 0) {
            // N(t-1) = Z'(t)*Z(t)/f(t) + L'(t)*N(t)*L(t)
            XL(N.rows());
            XL(N.columns());
            measurement.VpZdZ(pos, N.subMatrix(), 1 / f);
            SymmetricMatrix.reinforceSymmetry(N);
        } else {
            //T'*N(t)*T
            tvt(N);
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
        // R(t-1)=v(t)/f(t)Z(t)+R(t)L(t)
        //   = v/f*Z + R*(T-TC/f*Z)
        //  = (v - RT*C)/f*Z + RT
        dynamics.XT(pos, R);
        dynamics.MT(pos, Rd.subMatrix().transpose());
        if (!missing && f != 0) {
            // RT
            double c = (e - R.dot(C)) / f;
            measurement.XpZd(pos, R, c);
            // apply the same to the colums of Rd
            DataBlockIterator rcols = Rd.columns();
            DataBlock rcol = rcols.getData();
            do {
                c = (E.get(rcols.getPosition()) - rcol.dot(C)) / f;
                measurement.XpZd(pos, rcol, c);
            } while (rcols.next());
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
        smoother.process(ssf, frslts.getCollapsingPosition(), endpos, frslts, srslts);
        // updates R, N
        R.copy(smoother.getFinalR());
        if (calcvar) {
            N.copy(smoother.getFinalN());
        }
    }

    private void calcSmoothedDiffuseEffects() {
        // computes the smoothed diffuse effects and their covariance...

        QAugmentation q = frslts.getAugmentation();
        // delta = S(s+B*R), psi = Psi= S - S*B*N*B*S 
        // delta = a'^-1*a^-1(-a*b' + B*R)
        // delta = - (b * a^-1)' + a'^-1*a^-1*B*r = a'^-1 * (a^-1*B*r - b)
        // Psi = = a'^-1*(I - a^-1*B'*N*B*a'^-1)* a^-1
        SubMatrix B = frslts.B(pos);
        S = new Matrix(q.a());
        // computes B*r
        delta = new DataBlock(B.getColumnsCount());
        delta.product(B.columns(), R);
        // t1 = - b*a^-1 <-> -t1*a=b
        LowerTriangularMatrix.rsolve(S, delta);
        delta.sub(q.b());
        LowerTriangularMatrix.lsolve(S, delta);
        // B'NB 
        if (N != null) {
            // we have to make a copy of B
            Matrix A = new Matrix(B);
            // a^-1*B' =C <-> B'=aC
            LowerTriangularMatrix.rsolve(S, A.subMatrix().transpose());
            Psi = SymmetricMatrix.quadraticForm(N, A);
            Psi.chs();
            Psi.diagonal().add(1);
            // B*a^-1* =C <->B =Ca
            LowerTriangularMatrix.lsolve(S, Psi.subMatrix());
            // a'^-1*B = C <-> B' = C'a
            LowerTriangularMatrix.lsolve(S, Psi.subMatrix().transpose());
        }
    }

    public IAugmentedFilteringResults getFilteringResults() {
        return frslts; 
    }
}
