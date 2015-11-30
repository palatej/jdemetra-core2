/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.implementations;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.DataBlockIterator;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;

/**
 *
 * @author Jean Palate
 */
public class TimeInvariantMeasurement implements ISsfMeasurement {

    private final DataBlock Z;
    private final double var;

    public static TimeInvariantMeasurement of(int dim, ISsfMeasurement measurement) {
        if (!measurement.isTimeInvariant()) {
            return null;
        }
        DataBlock Z = new DataBlock(dim);
        measurement.Z(0, Z);
        if (!measurement.hasErrors()) {
            return new TimeInvariantMeasurement(Z, 0);
        } else {
            return new TimeInvariantMeasurement(Z, measurement.errorVariance(0));
        }

    }

    public TimeInvariantMeasurement(DataBlock Z, double var) {
        this.Z = Z;
        this.var = var;
    }

    @Override
    public boolean isTimeInvariant() {
        return true;
    }

    @Override
    public void Z(int pos, DataBlock z) {
        z.copy(Z);
    }


    @Override
    public boolean hasErrors() {
        return var != 0;
    }

    @Override
    public boolean hasError(int pos) {
        return var != 0;
    }

    @Override
    public double errorVariance(int pos) {
        return var;
    }

    @Override
    public double ZX(int pos, DataBlock m) {
        return Z.dot(m);
    }

    @Override
    public double ZVZ(int pos, SubMatrix V) {
        DataBlock zv = new DataBlock(V.getColumnsCount());
        zv.product(Z, V.columns());
        return zv.dot(Z);
    }

    @Override
    public void VpZdZ(int pos, SubMatrix V, double d) {
        
        DataBlockIterator cols = V.columns();
        DataBlock col = cols.getData();
        DataBlock zi = Z, zj = Z;
        int i = 0;
        do {
            col.addAY(d * zj.get(i++), zi);
        } while (cols.next());
    }

    @Override
    public void XpZd(int pos, DataBlock x, double d) {
        x.addAY(d, Z);
    }

}
