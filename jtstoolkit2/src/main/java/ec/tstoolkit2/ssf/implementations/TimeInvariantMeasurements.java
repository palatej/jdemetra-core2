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
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import ec.tstoolkit2.ssf.multivariate.ISsfMeasurements;

/**
 *
 * @author Jean Palate
 */
public class TimeInvariantMeasurements implements ISsfMeasurements {

    private final Matrix Z, H, R;

    public static TimeInvariantMeasurements of(int dim, ISsfMeasurement measurement) {
        return of(dim, Measurements.proxy(measurement));
    }

    public static TimeInvariantMeasurements of(int dim, ISsfMeasurements measurements) {
        if (!measurements.isTimeInvariant()) {
            return null;
        }
        int m = measurements.getMaxCount();
        Matrix Z = new Matrix(m, dim);
        measurements.Z(0, Z.subMatrix());
        if (!measurements.hasErrors()) {
            return new TimeInvariantMeasurements(Z, null);
        }
        Matrix H = Matrix.square(m), R = Matrix.square(m);
        measurements.H(0, H.subMatrix());
        measurements.R(0, R.subMatrix());
        return new TimeInvariantMeasurements(Z, H, R);

    }

    public TimeInvariantMeasurements(Matrix Z, Matrix H) {
        this.Z = Z;
        this.H = H;
        if (H != null) {
            R = H.clone();
            SymmetricMatrix.lcholesky(H, State.ZERO);
        } else {
            R = null;
        }
    }

    public TimeInvariantMeasurements(Matrix Z, Matrix H, Matrix R) {
        this.Z = Z;
        this.H = H;
        this.R = R;
    }

    @Override
    public boolean isTimeInvariant() {
        return true;
    }

    @Override
    public int getCount(int pos) {
        return Z.getRowsCount();
    }

    @Override
    public int getMaxCount() {
        return Z.getRowsCount();
    }

    @Override
    public boolean isHomogeneous() {
        return true;
    }

    @Override
    public void Z(int pos, int var, DataBlock z) {
        z.copy(Z.row(var));
    }

    @Override
    public void Z(int pos, SubMatrix z) {
        z.copy(Z.subMatrix());
    }

    @Override
    public boolean hasErrors() {
        return H != null;
    }

    @Override
    public boolean hasError(int pos) {
        return H != null;
    }

    @Override
    public boolean hasIndependentErrors() {
        return H == null || H.isDiagonal();
    }

    @Override
    public void H(int pos, SubMatrix h) {
        if (H != null) {
            h.copy(H.subMatrix());
        }
    }

    @Override
    public void R(int pos, SubMatrix r) {
        if (R != null) {
            r.copy(R.subMatrix());
        }
    }

    @Override
    public double ZX(int pos, int var, DataBlock m) {
        return Z.row(var).dot(m);
    }

    @Override
    public double ZVZ(int pos, int ivar, int jvar, SubMatrix V) {
        DataBlock zv = new DataBlock(Z.getColumnsCount());
        zv.product(Z.row(ivar), V.columns());
        return zv.dot(Z.row(jvar));
    }

    @Override
    public void ZVZ(int pos, SubMatrix V, SubMatrix zvz) {
        SymmetricMatrix.quadraticFormT(V, Z.subMatrix(), zvz);
    }

    @Override
    public void addH(int pos, SubMatrix V) {
        if (H != null) {
            V.add(H.subMatrix());
        }
    }

    @Override
    public void VpZdZ(int pos, int ivar, int jvar, SubMatrix V, double d) {
        DataBlockIterator cols = V.columns();
        DataBlock col = cols.getData();
        DataBlock zi = Z.row(ivar), zj = Z.row(jvar);
        int i = 0;
        do {
            col.addAY(d * zj.get(i++), zi);
        } while (cols.next());
    }

    @Override
    public void XpZd(int pos, int var, DataBlock x, double d) {
        x.addAY(d, Z.row(var));
    }

}
