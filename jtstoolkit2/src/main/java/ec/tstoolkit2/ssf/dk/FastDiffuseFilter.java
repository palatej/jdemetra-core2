/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.dk;

import ec.tstoolkit2.ssf.univariate.*;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.DataBlockIterator;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.ResultsRange;

/**
 *
 * @author Jean Palate
 */
public class FastDiffuseFilter {

    private final IBaseDiffuseFilteringResults frslts;
    private final ISsfMeasurement measurement;
    private final ISsfDynamics dynamics;
    private final int start, end, enddiffuse;
    private boolean normalized = false;

    public boolean filter(SubMatrix x) {
        return filter(null, x);
    }

    public boolean filter(SubMatrix x0, SubMatrix x) {
        if (x.getColumnsCount() == 1) {
            if (x0 == null) {
                return new FastDiffuseFilter1().filter(x.column(0));
            } else if (x0.getColumnsCount() == 1) {
                return new FastDiffuseFilter1().filter(x0.column(0), x.column(0));
            }else
                return false;
        } else {
            return new FastDiffuseFilterN().filter(x0, x);
        }
    }

    public boolean filter(DataBlock x) {
        return filter(null, x);
    }

    public boolean filter(DataBlock x0, DataBlock x) {
        return new FastDiffuseFilter1().filter(x0, x);
    }
    
    public FastDiffuseFilter(ISsf ssf, IBaseDiffuseFilteringResults frslts, ResultsRange range) {
        this.frslts = frslts;
        measurement = ssf.getMeasurement();
        dynamics = ssf.getDynamics();
        start = range.getStart();
        end = range.getEnd();
        enddiffuse = frslts.getEndDiffusePosition();
    }

    /**
     * @return the normalized
     */
    public boolean isNormalized() {
        return normalized;
    }

    /**
     * @param normalized the normalized to set
     */
    public void setNormalized(boolean normalized) {
        this.normalized = normalized;
    }

    class FastDiffuseFilterN {

        private SubMatrix states;
        // temporaries
        private DataBlock tmp, scol;
        private DataBlockIterator scols;

        boolean filter(SubMatrix x) {
            return filter(null, x);
        }

        boolean filter(SubMatrix x0, SubMatrix x) {
            if (start >= x.getRowsCount()) {
                return true;
            }
            int dim = dynamics.getStateDim();
            if (x0 == null) {
                states = new Matrix(dim, x.getColumnsCount()).subMatrix();
            } else if (x0.getRowsCount() != dim || x0.getColumnsCount() != x.getColumnsCount()) {
                return false;
            } else if (end - start != x.getRowsCount()) {
                return false;
            } else {
                states = x0;
            }
            prepareTmp();
            DataBlockIterator rows = x.rows();
            DataBlock row = rows.getData();
            int pos = start;
            do {
                iterate(pos, row);
            } while (++pos < end && rows.next());
            return true;
        }

        private void prepareTmp() {
            int nvars = states.getColumnsCount();
            tmp = new DataBlock(nvars);
            scols = states.columns();
            scol = scols.getData();
        }

        private void iterate(int i, DataBlock row) {
            boolean missing = !Double.isFinite(frslts.error(i));
            if (!missing) {
                double f = frslts.errorVariance(i);
                double w;
                DataBlock K;
                if (i < enddiffuse) {
                    double fi = frslts.diffuseNorm2(i);
                    if (fi != 0) {
                        w = fi;
                        K = frslts.Mi(i);
                    } else {
                        w = f;
                        K = frslts.M(i);
                    }
                } else {
                    w = f;
                    K = frslts.M(i);
                }

                measurement.ZM(i, states, tmp);
                row.sub(tmp);
                if (normalized && f != 0) {
                    row.mul(1 / Math.sqrt(f));
                }
                // update the states
                scols.begin();
                int j = 0;
                do {
                    scol.addAY(row.get(j++) / w, K);
                    dynamics.TX(i, scol);
                } while (scols.next());
            } else {
                scols.begin();
                do {
                    dynamics.TX(i, scol);
                } while (scols.next());
            }
            //  
        }

    }

    class FastDiffuseFilter1 {

        private DataBlock state;

        boolean filter(DataBlock x) {
            return filter(null, x);
        }

        boolean filter(DataBlock x0, DataBlock x) {
            if (start >= x.getLength()) {
                return true;
            }
            int dim = dynamics.getStateDim();
            if (x0 == null) {
                state = new DataBlock(dim);
            } else if (x0.getLength() != dim) {
                return false;
            } else if (end - start != x.getLength()) {
                return false;
            } else {
                state = x0;
            }
            int pos = start;
            do {
                x.set(pos, iterate(pos, x.get(pos)));
            } while (++pos < end);
            return true;
        }

        private double iterate(int i, double y) {
            boolean missing = !Double.isFinite(frslts.error(i));
            double e = Double.NaN;
            if (!missing) {
                double f = frslts.errorVariance(i);
                double w;
                DataBlock K;
                if (i < enddiffuse) {
                    double fi = frslts.diffuseNorm2(i);
                    if (fi != 0) {
                        w = fi;
                        K = frslts.Mi(i);
                    } else {
                        w = f;
                        K = frslts.M(i);
                    }
                } else {
                    w = f;
                    K = frslts.M(i);
                }

                e = y - measurement.ZX(i, state);
                if (normalized && f != 0) {
                    e /= Math.sqrt(f);
                }
                // update the states
                state.addAY(e / w, K);
            }
            dynamics.TX(i, state);
            return e;
        }

    }
}
