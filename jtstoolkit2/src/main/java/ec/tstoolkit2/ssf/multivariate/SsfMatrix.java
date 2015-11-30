/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.multivariate;

import ec.tstoolkit.data.DescriptiveStatistics;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;

/**
 *
 * @author Jean Palate
 */
public class SsfMatrix implements IMultivariateSsfData {

    private final Matrix x_;

    public SsfMatrix(Matrix x) {
        x_ = x;
    }

    public SsfMatrix(SubMatrix x) {
        x_ = new Matrix(x);
    }

    @Override
    public double get(int pos, int v) {
        return pos < x_.getRowsCount() ? x_.get(pos, v) : Double.NaN;
    }

    @Override
    public boolean isMissing(int pos, int v) {
        if (pos >= x_.getRowsCount()) {
            return true;
        }
        double y = x_.get(pos, v);
        return !DescriptiveStatistics.isFinite(y);
    }

    @Override
    public boolean hasData() {
        return true;
    }

    @Override
    public int getCount() {
        return x_.getRowsCount();
    }

    @Override
    public boolean isHomogeneous() {
        return true;
    }

    @Override
    public int getVarsCount(int pos) {
        return x_.getColumnsCount();
    }

    @Override
    public int getMaxVarsCount() {
        return x_.getColumnsCount();
    }
}
