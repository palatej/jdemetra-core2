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
package ec.tstoolkit2.ssf.akf;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.LogSign;
import ec.tstoolkit.design.Development;
import ec.tstoolkit.eco.Determinant;
import ec.tstoolkit.maths.matrices.LowerTriangularMatrix;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit2.ssf.IPredictionErrorDecomposition;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.multivariate.IMultivariateSsf;
import ec.tstoolkit2.ssf.multivariate.IMultivariateSsfData;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfData;

/**
 *
 * @author Jean Palate
 */
@Development(status = Development.Status.Alpha)
public class AugmentedPredictionErrorDecomposition implements IPredictionErrorDecomposition, IAugmentedFilteringResults {

    private final Determinant det = new Determinant();
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
    private Matrix Q;
    private int n, nd;

    /**
     *
     */
    public AugmentedPredictionErrorDecomposition() {
    }

    /**
     *
     */
    @Override
    public void clear() {
        det.clear();
        Q=null;
        n=0;
        nd=0;
    }

    /**
     *
     */
    @Override
    public void close() {
    }
    // TODO Update with Java 8
    private static boolean isPositive(DataBlock q) {
        for (int i = 0; i < q.getLength(); ++i) {
            if (q.get(i) < State.ZERO) {
                return false;
            }
        }
        return true;
    }

    @Override
    public boolean canCollapse() {
        return !isPositive(Q.diagonal().drop(0, 1));
    }

    @Override
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
        DataBlock D = b().deepClone();
        LowerTriangularMatrix.lsolve(S, D);
        for (int i = 0; i < d; ++i) {
            DataBlock col = A.column(i);
            state.a().addAY(-Q.get(d, i), col);
            state.P().addXaXt(1, col);
        }
        state.dropAllConstraints();
        return true;
    }

    private SubMatrix a() {
        return Q.subMatrix(0, nd, 0, nd);
    }

    private DataBlock b() {
        return Q.row(nd).range(0, nd);
    }

    private double c() {
        return Q.get(nd, nd);
    }

    /**
     *
     * @param nd
     */
    private void prepare(final int nd) {
        det.clear();
        this.n = 0;
        this.nd = nd;
        Q = new Matrix(nd + 1, nd + 2);
    }

    @Override
    public void save(final int t, final AugmentedPredictionError pe) {
        if (pe.isMissing()) {
            return;
        }
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

    public Matrix getFinalQ() {
        return Q;
    }

    @Override
    public AkfDiffuseLikelihood likelihood() {
        AkfDiffuseLikelihood ll = new AkfDiffuseLikelihood();
        double cc = c();
        cc *= cc;
        LogSign dsl = a().diagonal().sumLog();
        double dcorr = 2 * dsl.value;
        ll.set(cc, det.getLogDeterminant(), dcorr, n, nd);
        return ll;
    }

    @Override
    public void open(ISsf ssf, ISsfData data) {
        prepare(ssf.getDynamics().getNonStationaryDim());
    }

//    @Override
    public void open(IMultivariateSsf ssf, IMultivariateSsfData data) {
        prepare(ssf.getDynamics().getNonStationaryDim());
    }

    @Override
    public void save(int t, AugmentedState state) {
        // nothing to do. We are just interested by the prediction error...
    }

}