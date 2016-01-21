/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.dk.sqrt;

import ec.tstoolkit2.ssf.dk.*;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.DescriptiveStatistics;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.akf.AugmentedState;
import ec.tstoolkit2.ssf.univariate.ISmoothingResults;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfData;
import ec.tstoolkit2.ssf.univariate.OrdinarySmoother;

/**
 *
 * @author Jean Palate
 */
public class DiffuseSquareRootSmoother extends BaseDiffuseSmoother{

    private AugmentedState state;
    private IDiffuseSquareRootFilteringResults frslts;

    public boolean process(final ISsf ssf, final ISsfData data, ISmoothingResults sresults) {
        IDiffuseSquareRootFilteringResults fresults = DkToolkit.sqrtFilter(ssf, data, true);
        return process(ssf, data.getLength(), fresults, sresults);
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

    // 

    @Override
    protected void updateA() {
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

    @Override
    protected void updateP() {
        Matrix P = state.P();
        Matrix PN0P = SymmetricMatrix.quadraticForm(N0, P);
        SubMatrix B = state.B();
//        Matrix Pi=Matrix.square(B.getRowsCount());
//        SymmetricMatrix.XXt(B, Pi.subMatrix());
        Matrix BN2B = SymmetricMatrix.quadraticForm(N2.subMatrix(), B);
        Matrix PN2P=SymmetricMatrix.quadraticFormT(BN2B.subMatrix(), B);
        Matrix N1B= new Matrix(N1.getRowsCount(), B.getColumnsCount());
        N1B.subMatrix().product(N1.subMatrix(), B);
        Matrix PN1B=P.times(N1B);
        Matrix PN1Pi=Matrix.square(P.getRowsCount());
        PN1Pi.subMatrix().product(PN1B.subMatrix(), B.transpose());
//        Matrix PN2P = SymmetricMatrix.quadraticForm(N2, Pi);
//        Matrix PN1 = P.times(N1);
//        Matrix PN1Pi = PN1.times(Pi);
        P.sub(PN0P);
        P.sub(PN2P);
        P.sub(PN1Pi);
        P.subMatrix().sub(PN1Pi.subMatrix().transpose());
        SymmetricMatrix.reinforceSymmetry(P);
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

    public IDiffuseSquareRootFilteringResults getFilteringResults() {
        return frslts;
      }

}
