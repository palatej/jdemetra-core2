/*
 * Copyright 2015 National Bank of Belgium
 *  
 * Licensed under the EUPL, Version 1.1 or â€“ as soon they will be approved 
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
/*
 */
package ec.tstoolkit2.ssf.implementations.structural;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.Householder;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.implementations.Measurement;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import ec.tstoolkit2.ssf.univariate.Ssf;

/**
 *
 * @author Jean Palate
 */
public class SsfBsm extends Ssf {

    private SsfBsm(BsmDynamics dynamics, ISsfMeasurement measurement) {
        super(dynamics, measurement);
    }

    public static int searchPosition(BasicStructuralModel model, Component type) {
        int n = 0;
        if (model.nVar > 0) {
            if (type == Component.Noise) {
                return n;
            }
            ++n;
        }
        if (model.cVar >= 0) {
            if (type == Component.Cycle) {
                return n;
            }
            n += 2;
        }
        if (model.lVar >= 0) {
            if (type == Component.Level) {
                return n;
            }
            ++n;
        }
        if (model.sVar >= 0) {
            if (type == Component.Slope) {
                return n;
            }
            ++n;
        }
        if (model.seasVar >= 0 && type == Component.Seasonal) {
            return n;
        } else {
            return -1;
        }
    }

    public static int calcDim(BasicStructuralModel model) {
        int n = 0;
        if (model.nVar > 0) {
            ++n;
        }
        if (model.cVar >= 0) {
            n += 2;
        }
        if (model.lVar >= 0) {
            ++n;
        }
        if (model.sVar >= 0) {
            ++n;
        }
        if (model.seasVar >= 0) {
            n += model.getFrequency() - 1;
        }
        return n;
    }

    /**
     *
     */
    private static int[] calcCmpsIndexes(BasicStructuralModel model) {
        int n = 0;
        if (model.nVar > 0) {
            ++n;
        }
        if (model.cVar >= 0) {
            ++n;
        }
        if (model.lVar >= 0) {
            ++n;
        }
        if (model.seasVar >= 0) {
            ++n;
        }
        int[] cmps = new int[n];
        int i = 0, j = 0;
        if (model.nVar > 0) {
            cmps[i++] = j++;
        }
        if (model.cVar >= 0) {
            cmps[i++] = j;
            j += 2;
        }
        if (model.lVar >= 0) {
            cmps[i++] = j++;
        }
        if (model.sVar >= 0) {
            ++j;
        }
        if (model.seasVar >= 0) {
            cmps[i] = j;
        }
        return cmps;
    }

    public static SsfBsm create(BasicStructuralModel model) {
        int[] idx = calcCmpsIndexes(model);
        ISsfMeasurement measurement = Measurement.create(idx);
        BsmDynamics dynamics = new BsmDynamics(model);
        if (dynamics.isValid()) {
            return new SsfBsm(dynamics, measurement);
        } else {
            return null;
        }
    }

    private static Matrix tsvar(int freq) {
        int n = freq - 1;
        Matrix M = new Matrix(n, freq);
        M.diagonal().set(1);
        M.column(n).set(-1);
        Matrix O = SymmetricMatrix.XXt(M);
        Householder qr = new Householder(false);
        qr.decompose(O);
        Matrix Q = qr.solve(M);
        Matrix H = new Matrix(freq, n);
        // should be improved
        for (int i = 0; i < freq; ++i) {
            double z = 2 * Math.PI * (i + 1) / freq;
            for (int j = 0; j < n / 2; ++j) {
                H.set(i, 2 * j, Math.cos((j + 1) * z));
                H.set(i, 2 * j + 1, Math.sin((j + 1) * z));
            }
            if (n % 2 == 1) {
                H.set(i, n - 1, Math.cos((freq / 2) * z));
            }
        }
        Matrix Z = SymmetricMatrix.XXt(Q.times(H));
        Z.smooth(1e-12);

        return Z;
    }

    /**
     *
     * @param freq
     * @return
     */
    private static synchronized Matrix tsVar(int freq) {
        if (freq == 12) {
            if (g_VTS12 == null) {
                g_VTS12 = tsvar(12);
            }
            return g_VTS12.clone();
        } else if (freq == 4) {
            if (g_VTS4 == null) {
                g_VTS4 = tsvar(4);
            }
            return g_VTS4.clone();
        } else if (freq == 2) {
            if (g_VTS2 == null) {
                g_VTS2 = tsvar(2);
            }
            return g_VTS2.clone();
        } else if (freq == 3) {
            if (g_VTS3 == null) {
                g_VTS3 = tsvar(3);
            }
            return g_VTS3.clone();
        } else if (freq == 6) {
            if (g_VTS6 == null) {
                g_VTS6 = tsvar(6);
            }
            return g_VTS6.clone();
        } else {
            return tsvar(freq);
        }
    }

    static Matrix tsVar(SeasonalModel seasModel, final int freq) {
        if (seasModel == SeasonalModel.Trigonometric) {
            return tsVar(freq);
        } else {
            int n = freq - 1;
            Matrix Q = Matrix.square(n);
            // Dummy
            if (seasModel == SeasonalModel.Dummy) {
                Q.set(n - 1, n - 1, 1);
            } else if (seasModel == SeasonalModel.Crude) {
                Q.set(1);
                //Q.set(0, 0, freq * var);
            } else if (seasModel == SeasonalModel.HarrisonStevens) {
                double v = 1.0 / freq;
                Q.set(-v);
                Q.diagonal().add(1);
            }
            return Q;
        }
    }

    private static Matrix g_VTS2, g_VTS3, g_VTS4, g_VTS6, g_VTS12;

    public static class BsmDynamics implements ISsfDynamics {

        private final SubMatrix tsvar;
        private final double lVar, sVar, seasVar, cVar, nVar, cDump;
        private final double ccos, csin;
        private final int freq;
        private final SeasonalModel seasModel;

        public BsmDynamics(BasicStructuralModel model) {
            lVar = model.lVar;
            sVar = model.sVar;
            seasVar = model.seasVar;
            cVar = model.cVar;
            nVar = model.nVar;
            cDump = model.cDump;
            ccos = model.ccos;
            csin = model.csin;
            seasModel = model.seasModel;
            freq = model.freq;
            if (seasVar > 0) {
                tsvar = tsVar(seasModel, freq).subMatrix();
                tsvar.mul(seasVar);
            } else {
                tsvar = null;
            }
        }

        @Override
        public int getStateDim() {
            int r = 0;
            if (nVar > 0) {
                ++r;
            }
            if (cVar >= 0) {
                r += 2;
            }
            if (lVar >= 0) {
                ++r;
            }
            if (sVar >= 0) {
                ++r;
            }
            if (seasVar >= 0) {
                r += freq - 1;
            }
            return r;
        }

        @Override
        public boolean isTimeInvariant() {
            return true;
        }

        @Override
        public boolean isValid() {
            if (freq == 1 && seasVar >= 0) {
                return false;
            }
            return lVar >= 0 || sVar >= 0 || cVar >= 0 || nVar >= 0;
        }

        @Override
        public int getInnovationsDim() {
            int nr = 0;
            if (seasVar > 0) {
                if (seasModel == SeasonalModel.Dummy) {
//                    || seasModel == SeasonalModel.Crude) {
                    ++nr;
                } else {
                    nr += freq - 1;
                }
            }
            if (nVar > 0) {
                ++nr;
            }
            if (cVar > 0) {
                nr += 2;
            }
            if (lVar > 0) {
                ++nr;
            }
            if (sVar > 0) {
                ++nr;
            }
            return nr;
        }

        @Override
        public void V(int pos, SubMatrix v) {
            int i = 0;
            if (nVar > 0) {
                v.set(i, i, nVar);
                ++i;
            }
            if (cVar >= 0) {
                v.set(i, i, cVar);
                ++i;
                v.set(i, i, cVar);
                ++i;
            }
            if (lVar >= 0) {
                if (lVar != 0) {
                    v.set(i, i, lVar);
                }
                ++i;
            }
            if (sVar >= 0) {
                if (sVar != 0) {
                    v.set(i, i, sVar);
                }
                ++i;
            }
            if (seasVar > 0) {
                if (seasModel == SeasonalModel.Dummy) {
                    v.set(i, i, seasVar);
                } else {
                    tsvar.copyTo(v, i, i);
                }
            }
        }

        @Override
        public boolean hasS() {
            return false;
        }

        @Override
        public boolean hasInnovations(int pos) {
            return true;
        }

        @Override
        public void Q(int pos, SubMatrix q) {
            V(pos, q);
        }

        @Override
        public void S(int pos, SubMatrix sm) {
        }

        @Override
        public void T(int pos, SubMatrix tr) {
            int i = 0;
            if (nVar > 0) {
                ++i;
            }
            if (cVar >= 0) {
                tr.set(i, i, ccos);
                tr.set(i + 1, i + 1, ccos);
                tr.set(i, i + 1, csin);
                tr.set(i + 1, i, -csin);
                i += 2;
            }
            if (lVar >= 0) {
                tr.set(i, i, 1);
                if (sVar >= 0) {
                    tr.set(i, i + 1, 1);
                    ++i;
                    tr.set(i, i, 1);
                }
                ++i;
            }
            if (seasVar >= 0) {
                SubMatrix seas = tr.extract(i, i + freq - 1, i, i + freq - 1);
                seas.row(freq - 2).set(-1);
                seas.subDiagonal(1).set(1);
            }
        }

        @Override
        public boolean isDiffuse() {
            return lVar >= 0 || seasVar >= 0;
        }

        @Override
        public int getNonStationaryDim() {
            int r = 0;
            if (lVar >= 0) {
                ++r;
            }
            if (sVar >= 0) {
                ++r;
            }
            if (seasVar >= 0) {
                r += freq - 1;
            }
            return r;
        }

        @Override
        public void diffuseConstraints(SubMatrix b) {
            int sdim = getStateDim();
            int istart = nVar > 0 ? 1 : 0;
            if (cVar >= 0) {
                istart += 2;
            }
            int iend = sdim;
            for (int i = istart, j = 0; i < iend; ++i, ++j) {
                b.set(i, j, 1);
            }
        }

        @Override
        public void Pi0(SubMatrix p) {
            int sdim = getStateDim();
            int istart = nVar > 0 ? 1 : 0;
            if (cVar >= 0) {
                istart += 2;
            }
            int iend = sdim;
            for (int i = istart; i < iend; ++i) {
                p.set(i, i, 1);
            }
        }

        @Override
        public boolean a0(DataBlock a0, StateInfo info) {
            return true;
        }

        @Override
        public boolean Pf0(SubMatrix p, StateInfo info) {
            if (info == StateInfo.Forecast) {
                int i = 0;
                if (nVar > 0) {
                    p.set(0, 0, nVar);
                    ++i;
                }
                if (cVar > 0) {
                    double q = cVar / (1 - cDump * cDump);
                    p.set(i, i, q);
                    ++i;
                    p.set(i, i, q);
                }
                return true;
            } else {
                // TODO
                return false;
            }
        }

        @Override
        public void TX(int pos, DataBlock x) {
            int i0 = 0;
            if (nVar > 0) {
                x.set(0, 0);
                ++i0;
            }
            if (cVar >= 0) {
                double a = x.get(i0), b = x.get(i0 + 1);
                x.set(i0, a * ccos + b * csin);
                x.set(i0 + 1, -a * csin + b * ccos);
                i0 += 2;
            }
            if (lVar >= 0) {
                if (sVar >= 0) {
                    x.add(i0, x.get(i0 + 1));
                    i0 += 2;
                } else {
                    ++i0;
                }
            }
            if (seasVar >= 0) {
                DataBlock ex = x.extract(i0, freq - 1, 1);
                ex.bshift(DataBlock.ShiftOption.NegSum);
            }
        }

        @Override
        public void XT(int pos, DataBlock x) {
            int i0 = 0;
            if (nVar > 0) {
                x.set(0, 0);
                ++i0;
            }
            if (cVar >= 0) {
                double a = x.get(i0), b = x.get(i0 + 1);
                x.set(i0, a * ccos - b * csin);
                x.set(i0 + 1, a * csin + b * ccos);
                i0 += 2;

            }
            if (lVar >= 0) {
                if (sVar >= 0) {
                    x.add(i0 + 1, x.get(i0));
                    i0 += 2;
                } else {
                    ++i0;
                }
            }
            if (seasVar >= 0) {
                int imax = i0 + freq - 2;
                double xs = x.get(imax);
                for (int i = imax; i > i0; --i) {
                    x.set(i, x.get(i - 1) - xs);
                }
                x.set(i0, -xs);
            }
        }

        @Override
        public void addV(int pos, SubMatrix p) {
            int i = 0;
            if (nVar > 0) {
                p.add(i, i, nVar);
                ++i;
            }
            if (cVar >= 0) {
                p.add(i, i, cVar);
                ++i;
                p.add(i, i, cVar);
                ++i;
            }
            if (lVar >= 0) {
                if (lVar != 0) {
                    p.add(i, i, lVar);
                }
                ++i;
            }
            if (sVar >= 0) {
                if (sVar != 0) {
                    p.add(i, i, sVar);
                }
                ++i;
            }
            if (seasVar > 0) {
                if (seasModel == SeasonalModel.Dummy) {
                    p.add(i, i, seasVar);
                } else {
                    tsvar.addTo(p, i, i);
                }
            }

        }

    }
}
