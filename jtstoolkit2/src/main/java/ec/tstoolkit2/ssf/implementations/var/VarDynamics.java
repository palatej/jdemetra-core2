/*
 * Copyright 2013-2014 National Bank of Belgium
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
package ec.tstoolkit2.ssf.implementations.var;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.HouseholderR;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.StateInfo;

/**
 *
 * @author Jean Palate
 */
public class VarDynamics implements ISsfDynamics {

    private final VarDescriptor desc;
    private final Matrix V0;
    private final int neq, nl, nlx;
    private final DataBlock ttmp, xtmp;
    private final StateInfo info;

    public static VarDynamics from(final VarDescriptor desc) {
        return new VarDynamics(desc, desc.getLagsCount());
    }

    public static VarDynamics from(final VarDescriptor desc, final int nlx) {
        if (nlx < desc.getLagsCount()) {
            return null;
        }
        return new VarDynamics(desc, nlx);
    }

    public static VarDynamics conditionalVar(final VarDescriptor desc, final int nlx) {
        if (nlx < desc.getLagsCount()) {
            return null;
        }
        return new VarDynamics(desc, nlx, null, StateInfo.Concurrent);
    }

    public static VarDynamics from(final VarDescriptor desc, final Matrix V0) {
        int nlx = V0.getColumnsCount();
        if (nlx % desc.getEquationsCount() != 0) {
            return null;
        }
        nlx /= desc.getEquationsCount();
        if (nlx < desc.getLagsCount()) {
            return null;
        }
        return new VarDynamics(desc, nlx, V0, StateInfo.Concurrent);
    }

    private VarDynamics(final VarDescriptor desc, final int nlx) {
        this.desc = desc;
        nl = desc.getLagsCount();
        this.nlx = nlx;
        neq = desc.getEquationsCount();
        ttmp = new DataBlock(neq);
        xtmp = new DataBlock(neq * nlx);
        V0 = initialCovariance();
        info = StateInfo.Undefined;
    }

    private VarDynamics(final VarDescriptor desc, final int nlx, final Matrix V0, final StateInfo info) {
        this.desc = desc;
        nl = desc.getLagsCount();
        this.nlx = nlx;
        neq = desc.getEquationsCount();
        if (V0 != null) {
            this.V0 = V0.clone();
        } else {
            this.V0 = null;
        }
        ttmp = new DataBlock(neq);
        xtmp = new DataBlock(neq * nlx);

        this.info = info;
    }

    public VarDescriptor getDescriptor() {
        return desc;
    }

    public int getLagsCount() {
        return nlx;
    }

    @Override
    public int getStateDim() {
        return neq * nlx;
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
        return neq;
    }

    @Override
    public void V(int pos, SubMatrix qm) {
        Matrix v = desc.getInnovationsVariance();
        for (int c = 0; c < neq; ++c) {
            qm.column(c * nlx).extract(0, neq, nlx).copy(v.column(c));
        }
    }

    @Override
    public boolean hasS() {
        return true;
    }

    @Override
    public boolean hasInnovations(int pos) {
        return true;
    }

    @Override
    public void Q(int pos, SubMatrix qm) {
        qm.copy(desc.getInnovationsVariance().subMatrix());
    }

    @Override
    public void S(int pos, SubMatrix sm) {
        for (int i = 0, r = 0; i < neq; ++i, r += nlx) {
            sm.set(r, i, 1);
        }
    }

//    @Override
//    public void addSX(int pos, DataBlock x, DataBlock y) {
//        for (int i = 0, r = 0; i < neq; ++i, r += nlx) {
//            y.add(r, x.get(i));
//        }
//    }
//    
    @Override
    public void T(int pos, SubMatrix tr) {
        Matrix v = desc.getVarMatrix();
        for (int i = 0, r = 0; i < neq; ++i, r += nlx) {
            for (int j = 0, c = 0; j < neq; ++j, c += nlx) {
                SubMatrix B = tr.extract(r, r + nlx, c, c + nlx);
                if (i == j) {
                    B.subDiagonal(-1).set(1);
                }
                B.row(0).range(0, nl).
                        copy(v.row(i).range(j * nl, (j + 1) * nl));
            }
        }
    }

    @Override
    public boolean isDiffuse() {
        return false;
    }

    @Override
    public int getNonStationaryDim() {
        return 0;
    }

    @Override
    public void diffuseConstraints(SubMatrix b) {
    }

    @Override
    public boolean a0(DataBlock a0, StateInfo info) {
        return true;
    }

    @Override
    public boolean Pf0(SubMatrix pf0, StateInfo info) {
        if (info == StateInfo.Smoothed && this.info == StateInfo.Forecast) {
            return false;
        }
        pf0.copy(V0.subMatrix());
        if (this.info == StateInfo.Undefined || this.info == info) {
            return true;
        } else {
            TVT(0, pf0);
            addV(0, pf0);
            return true;
        }
    }

    @Override
    public void TX(int pos, DataBlock x) {
        Matrix v = desc.getVarMatrix();
        // compute first the next item
        for (int i = 0; i < neq; ++i) {
            double r = 0;
            DataBlock p = v.row(i).range(0, nl);
            DataBlock xb = x.range(0, nl);
            for (int j = 0; j < neq; ++j) {
                if (j != 0) {
                    p.move(nl);
                    xb.move(nlx);
                }
                r += p.dot(xb);
            }
            ttmp.set(i, r);
        }
        x.fshift(DataBlock.ShiftOption.Zero);
        x.extract(0, -1, nlx).copy(ttmp);
    }

    @Override
    public void XT(int pos, DataBlock x) {
        Matrix v = desc.getVarMatrix();
        for (int i = 0, k = 0, l = 0; i < neq; ++i) {
            for (int j = 0; j < nl; ++j, ++k) {
                double r = ((k + 1) % nl != 0) ? x.get(k + 1) : 0;
                r += v.column(l++).dot(x.extract(0, neq, nl));
                xtmp.set(k, r);
            }
            for (int j = nl; j < nl - 1; ++j, ++k) {
                xtmp.set(k, x.get(k + 1));
            }
            if (nlx > nl) {
                xtmp.set(k++, 0);
            }
        }
        x.copy(xtmp);
    }

    @Override
    public void addV(final int pos, final SubMatrix v) {
        for (int i = 0; i < neq; ++i) {
            DataBlock cv = v.column(i * nlx).extract(0, neq, nlx);
            cv.add(desc.getInnovationsVariance().column(i));
        }
    }

    private static int pos(int r, int c, int n) {
        return r + c * (2 * n - c - 1) / 2;
    }

    private Matrix initialCovariance() {
        // We have to solve the steady state equation:
        // V = T V T' + Q
        // We consider the nlag*nb, nlag*nb sub-system
        Matrix var = desc.getInnovationsVariance(), coeff = desc.getVarMatrix();

        int n = neq * nl;
        Matrix cov = new Matrix(n, n);
        int np = (n * (n + 1)) / 2;
        Matrix M = new Matrix(np, np);
        double[] b = new double[np];
        // fill the matrix
        for (int c = 0, i = 0; c < n; ++c) {
            for (int r = c; r < n; ++r, ++i) {
                M.set(i, i, 1);
                if (r % nl == 0 && c % nl == 0) {
                    b[i] = var.get(r / nl, c / nl);
                }
                for (int k = 0; k < n; ++k) {
                    for (int l = 0; l < n; ++l) {
                        double zr = 0, zc = 0;
                        if (r % nl == 0) {
                            zr = coeff.get(r / nl, l);
                        } else if (r == l + 1) {
                            zr = 1;
                        }
                        if (c % nl == 0) {
                            zc = coeff.get(c / nl, k);
                        } else if (c == k + 1) {
                            zc = 1;
                        }
                        double z = zr * zc;
                        if (z != 0) {
                            int p = l <= k ? pos(k, l, n) : pos(l, k, n);
                            M.add(i, p, -z);
                        }
                    }
                }
            }
        }
        HouseholderR hous = new HouseholderR(false);
        hous.decompose(M);
        double[] solve = hous.solve(b);
        for (int i = 0, j = 0; i < n; i++) {
            cov.column(i).drop(i, 0).copyFrom(solve, j);
            j += n - i;
        }
        SymmetricMatrix.fromLower(cov);
        Matrix fullCov = new Matrix(getStateDim(), getStateDim());
        for (int r = 0; r < neq; ++r) {
            for (int c = 0; c < neq; ++c) {
                fullCov.subMatrix(r * nlx, r * nlx + nl, c * nlx, c * nlx + nl).copy(cov.subMatrix(r * nl, (r + 1) * nl, c * nl, (c + 1) * nl));
            }
        }
        for (int i = nl; i < nlx; ++i) {
            TVT(0, fullCov.subMatrix());
            addV(0, fullCov.subMatrix());
        }
        return fullCov;
    }

}
