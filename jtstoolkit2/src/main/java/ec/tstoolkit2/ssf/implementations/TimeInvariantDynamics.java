/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.implementations;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.DataBlockIterator;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.StateInfo;

/**
 *
 * @author Jean Palate
 */
public class TimeInvariantDynamics implements ISsfDynamics {

    public static class Innovations {

        private static Innovations of(ISsfDynamics sd) {
            int n = sd.getStateDim();
            int ne = sd.getInnovationsDim();
            Matrix V = Matrix.square(n);
            sd.V(0, V.subMatrix());
            Matrix S = new Matrix(n, ne);
            sd.S(0, S.subMatrix());
            return new Innovations(V, S);
        }

        public Innovations(final Matrix V) {
            this.V = V;
            S = null;
        }

        public Innovations(final Matrix V, final Matrix S) {
            this.S = S;
            if (V == null && S != null) {
                this.V = SymmetricMatrix.XXt(S);
            } else {
                this.V = V;
            }
        }

        public final Matrix S, V;
    }

    public static class Initialization {

        public static Initialization of(ISsfDynamics sd, StateInfo info) {
            int n = sd.getStateDim();
            Matrix P0 = Matrix.square(n);
            DataBlock a0 = new DataBlock(n);
            SubMatrix p0 = P0.subMatrix();
            sd.Pf0(p0, info);
            sd.a0(a0, info);
            if (!sd.isDiffuse()) {
                return new Initialization(P0, a0, info);
            }
            int nd = sd.getNonStationaryDim();
            Matrix B0 = new Matrix(n, nd);
            Matrix Pi0 = Matrix.square(n);
            sd.diffuseConstraints(B0.subMatrix());
            sd.Pi0(Pi0.subMatrix());
            return new Initialization(P0, B0, Pi0, a0, info);
        }

        public Initialization(final Matrix P0, StateInfo info) {
            this.P0 = P0;
            Pi0 = null;
            B0 = null;
            a0 = null;
            initial = info;
        }

        public Initialization(final Matrix P0, final Matrix B0, final StateInfo info) {
            initial = info;
            this.P0 = P0;
            this.B0 = B0;
            Pi0 = SymmetricMatrix.XXt(B0);
            a0 = null;
        }

        public Initialization(final Matrix P0, final Matrix Pi0, final Matrix B0, final StateInfo info) {
            initial = info;
            this.P0 = P0;
            this.Pi0 = Pi0;
            this.B0 = B0;
            a0 = null;
        }

        public Initialization(final Matrix P0, final DataBlock a0, final StateInfo info) {
            initial = info;
            this.P0 = P0;
            Pi0 = null;
            B0 = null;
            this.a0 = a0;
        }

        public Initialization(final Matrix P0, final Matrix B0, final DataBlock a0, final StateInfo info) {
            initial = info;
            this.P0 = P0;
            this.B0 = B0;
            Pi0 = SymmetricMatrix.XXt(B0);
            this.a0 = a0;
        }

        public Initialization(final Matrix P0, final Matrix B0, final Matrix Pi0, final DataBlock a0, final StateInfo info) {
            initial = info;
            this.P0 = P0;
            this.B0 = B0;
            this.Pi0 = Pi0;
            this.a0 = a0;
        }

        public final Matrix P0, B0, Pi0;
        public final DataBlock a0;
        public final StateInfo initial;
    }

    private final Matrix T;
    private final Matrix V;
    private transient Matrix S;

    private final Matrix Pf0, B0;
    private final DataBlock a0;
    private final StateInfo initial;

    public TimeInvariantDynamics(Matrix T, Innovations E, Initialization I) {
        this.T = T;
        this.B0 = I.B0;
        this.Pf0 = I.P0;
        this.a0 = I.a0;
        this.initial = I.initial;
        this.S = E.S;
        this.V = E.V;

    }

    public static TimeInvariantDynamics of(ISsfDynamics sd, StateInfo info) {
        if (!sd.isTimeInvariant()) {
            return null;
        }
        int n = sd.getStateDim();
        Matrix t = Matrix.square(n);
        sd.T(0, t.subMatrix());
        Innovations e = Innovations.of(sd);
        if (e == null) {
            return null;
        }
        Initialization i = Initialization.of(sd, info);
        if (i == null) {
            return null;
        }
        return new TimeInvariantDynamics(t, e, i);
    }

    private synchronized void checkS() {
        if (S == null) {
            S = V.clone();
            SymmetricMatrix.lcholesky(S);
        }
    }

    @Override
    public int getStateDim() {
        return T.getColumnsCount();
    }

    @Override
    public boolean isTimeInvariant() {
        return true;
    }

    @Override
    public boolean isValid() {
        return true;
    }

    @Override
    public int getInnovationsDim() {
        return S == null ? getStateDim() : S.getColumnsCount();
    }

    @Override
    public void V(int pos, SubMatrix qm) {
        qm.copy(V.subMatrix());
    }

    @Override
    public boolean hasInnovations(int pos) {
        return V != null;
    }

    @Override
    public void S(int pos, SubMatrix sm) {
        checkS();
        sm.copy(S.subMatrix());
    }

    @Override
    public void addSU(int pos, DataBlock x, DataBlock u) {
        checkS();
        x.addProduct(S.rows(), u);
    }

    @Override
    public void XS(int pos, DataBlock x, DataBlock xs) {
        checkS();
        xs.product(x, S.columns());
    }

    @Override
    public void T(int pos, SubMatrix tr) {
        tr.copy(T.subMatrix());
    }

    @Override
    public boolean isDiffuse() {
        return B0 != null;
    }

    @Override
    public int getNonStationaryDim() {
        return B0 == null ? 0 : B0.getColumnsCount();
    }

    @Override
    public void diffuseConstraints(SubMatrix b) {
        if (B0 != null) {
            b.copy(B0.subMatrix());
        }
    }

    @Override
    public boolean a0(DataBlock a0, StateInfo info) {
        if (this.a0 == null) {
            return true;
        }
        if (info == StateInfo.Concurrent && initial == StateInfo.Forecast) {
            return false;
        }
        a0.copy(this.a0);
        if (info != initial) {
            TX(-1, a0);
        }
        return true;
    }

    @Override
    public boolean Pf0(SubMatrix pf0, StateInfo info) {
        if (info != initial) {
            return false;
        }
        pf0.copy(this.Pf0.subMatrix());
        return true;
    }

    @Override
    public void TM(int pos, SubMatrix tm) {
        DataBlock tx = new DataBlock(T.getColumnsCount());
        DataBlockIterator cols = tm.columns();
        DataBlock col = cols.getData();
        do {
            tx.product(T.rows(), col);
            col.copy(tx);
        } while (cols.next());
    }

    @Override
    public void TVT(int pos, SubMatrix tvt) {
        Matrix V = new Matrix(tvt);
        SymmetricMatrix.quadraticFormT(V.subMatrix(), T.subMatrix(), tvt);
    }

    @Override
    public void TX(int pos, DataBlock x) {
        DataBlock tx = new DataBlock(x.getLength());
        tx.product(T.rows(), x);
        x.copy(tx);
    }

    @Override
    public void XT(int pos, DataBlock x) {
        DataBlock tx = new DataBlock(x.getLength());
        tx.product(x, T.columns());
        x.copy(tx);
    }

    @Override
    public void addV(int pos, SubMatrix p) {
        p.add(V.subMatrix());
    }

}
