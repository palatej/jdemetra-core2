/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.akf;

import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.LogSign;
import ec.tstoolkit.eco.Determinant;
import ec.tstoolkit.maths.matrices.LowerTriangularMatrix;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;

/**
 *
 * @author Admin
 */
public class QAugmentation {

    // Q is the cholesky factor of the usual "Q matrix" of De Jong.
    // Q(dj) = |S   -s|
    //         |-s'  q|
    // Q = |a 0|
    //     |b c|
    // so that we have:
    // q = b * b' + c * c
    // S = a * a' 
    // -s = a * b'
    // s' * S^-1 * s = b * a' * S^-1 * a * b' = b * b'
    // q - s' * S^-1 * s = c * c
    // s' * S^-1 = b * a' * S^-1 = b * a^-1 
    private Matrix Q;
    private int n, nd;
    private Determinant det = new Determinant();

    public void prepare(final int nd, final int nvars) {
        clear();
        this.nd = nd;
        Q = new Matrix(nd + 1, nd + 1 + nvars);
    }

    public void clear() {
        n = 0;
        Q = null;
        det.clear();
    }
    
    public int getDegreesofFreedom(){
        return n-nd;
    }

//    public void update(Matrix E, DataBlock U) {
//        Q.subMatrix(0, nd, nd + 1, nd + 1 + nvars).copy(E.subMatrix());
//        Q.row(nd).range(nd + 1, nd + 1 + nvars).copy(U);
//        ec.tstoolkit.maths.matrices.ElementaryTransformations.fastGivensTriangularize(Q.subMatrix());
//    }
//
    public void update(AugmentedPredictionError pe) {
        ++n;
        double v = pe.getVariance();
        double e = pe.get();
        det.add(v);
        DataBlock col = Q.column(nd + 1);
        double se = Math.sqrt(v);
        col.range(0, nd).setAY(1 / se, pe.E());
        col.set(nd, e / se);
        ec.tstoolkit.maths.matrices.ElementaryTransformations.fastGivensTriangularize(Q.subMatrix());
    }

    public SubMatrix a() {
        return Q.subMatrix(0, nd, 0, nd);
    }

    public DataBlock b() {
        return Q.row(nd).range(0, nd);
    }

    public double c() {
        return Q.get(nd, nd);
    }

    public AkfDiffuseLikelihood likelihood() {
        AkfDiffuseLikelihood ll = new AkfDiffuseLikelihood();
        double cc = c();
        cc *= cc;
        LogSign dsl = a().diagonal().sumLog();
        double dcorr = 2 * dsl.value;
        ll.set(cc, det.getLogDeterminant(), dcorr, n, nd);
        return ll;
    }

    public boolean canCollapse() {
        return isPositive(Q.diagonal().drop(0, 1));
    }

    public boolean collapse(AugmentedState state) {
        if (state.getInfo() != StateInfo.Forecast) {
            return false;
        }
        if (!isPositive(Q.diagonal().drop(0, 1))) {
            return false;
        }

        // update the state vector
        Matrix A = new Matrix(state.B());
        int d = A.getColumnsCount();
        Matrix S = new Matrix(a());
        LowerTriangularMatrix.rsolve(S, A.subMatrix().transpose());
        for (int i = 0; i < d; ++i) {
            DataBlock col = A.column(i);
            state.a().addAY(-Q.get(d, i), col);
            state.P().addXaXt(1, col);
        }
        state.dropAllConstraints();
        return true;
    }
    // TODO Update with Java 8

    public static boolean isPositive(DataBlock q) {
        for (int i = 0; i < q.getLength(); ++i) {
            if (q.get(i) < State.ZERO) {
                return false;
            }
        }
        return true;
    }
}
