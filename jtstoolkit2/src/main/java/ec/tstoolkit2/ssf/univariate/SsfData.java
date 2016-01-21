/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.univariate;

import ec.tstoolkit.data.DescriptiveStatistics;
import ec.tstoolkit.data.IReadDataBlock;
import ec.tstoolkit.data.ReadDataBlock;

/**
 *
 * @author Jean Palate
 */
public class SsfData implements ISsfData, IReadDataBlock {

    private final IReadDataBlock x_;

    public SsfData(IReadDataBlock x) {
        x_ = x;
    }

    public SsfData(double[] x) {
        x_ = new ReadDataBlock(x);
    }

    @Override
    public double get(int pos) {
        return pos < x_.getLength() ? x_.get(pos) : Double.NaN;
    }

    @Override
    public boolean isMissing(int pos) {
        if (pos >= x_.getLength()) {
            return true;
        }
        double y = x_.get(pos);
        return !DescriptiveStatistics.isFinite(y);
    }

    @Override
    public boolean hasData() {
        return true;
    }

    @Override
    public int getLength() {
        return x_.getLength();
    }

    @Override
    public void copyTo(double[] buffer, int start) {
        x_.copyTo(buffer, start);
    }

    @Override
    public IReadDataBlock rextract(int start, int length) {
        return x_.rextract(start, length);
    }
}
