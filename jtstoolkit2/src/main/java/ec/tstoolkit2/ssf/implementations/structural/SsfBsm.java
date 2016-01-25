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
        Matrix QH = Q.times(H);
        Matrix Z = SymmetricMatrix.XXt(QH);
        Z.smooth(1e-12);

        return Z;
    }

    /**
     *
     * @param freq
     * @return
     */
    private static synchronized Matrix tsVar(int freq) {
        switch (freq) {
            case 12:
                if (VTS12 == null) {
                    VTS12 = tsvar(12);
                }
                return VTS12.clone();
            case 4:
                if (VTS4 == null) {
                    VTS4 = tsvar(4);
                }
                return VTS4.clone();
            case 2:
                if (VTS2 == null) {
                    VTS2 = tsvar(2);
                }
                return VTS2.clone();
            case 3:
                if (VTS3 == null) {
                    VTS3 = tsvar(3);
                }
                return VTS3.clone();
            case 6:
                if (VTS6 == null) {
                    VTS6 = tsvar(6);
                }
                return VTS6.clone();
            default:
                return tsvar(freq);
        }
    }

    private static Matrix hsvar(int freq) {
        Matrix m = Matrix.square(freq - 1);
        m.set(-1.0 / freq);
        m.diagonal().add(1);
        return m;
    }

    private static synchronized Matrix hslVar(int freq) {
        switch (freq) {
            case 12:
                if (LHS12 == null) {
                    LHS12 = hsvar(12);
                    SymmetricMatrix.lcholesky(LHS12);
                }
                return LHS12.clone();
            case 4:
                if (LHS4 == null) {
                    LHS4 = hsvar(4);
                    SymmetricMatrix.lcholesky(LHS4);
                }
                return LHS4.clone();
            case 2:
                if (LHS2 == null) {
                    LHS2 = hsvar(2);
                    SymmetricMatrix.lcholesky(LHS2);
                }
                return LHS2.clone();
            case 3:
                if (LHS3 == null) {
                    LHS3 = hsvar(3);
                    SymmetricMatrix.lcholesky(LHS3);
                }
                return LHS3.clone();
            case 6:
                if (LHS6 == null) {
                    LHS6 = hsvar(6);
                    SymmetricMatrix.lcholesky(LHS6);
                }
                return LHS6.clone();
            default:
                Matrix lhs = hsvar(freq);
                SymmetricMatrix.lcholesky(lhs);
                return lhs;
        }
    }

    public static Matrix tsVar(SeasonalModel seasModel, final int freq) {
        if (seasModel == SeasonalModel.Trigonometric) {
            return tsVar(freq);
        } else {
            int n = freq - 1;
            Matrix Q = Matrix.square(n);
            if (null != seasModel) // Dummy
            {
                switch (seasModel) {
                    case Dummy:
                        Q.set(n - 1, n - 1, 1);
                        break;
                    case Crude:
                        Q.set(1);
                        //Q.set(0, 0, freq * var);
                        break;
                    case HarrisonStevens:
                        double v = 1.0 / freq;
                        Q.set(-v);
                        Q.diagonal().add(1);
                        break;
                    default:
                        break;
                }
            }
            return Q;
        }
    }

    private static synchronized Matrix tslVar(int freq) {
        switch (freq) {
            case 12:
                if (LVTS12 == null) {
                    LVTS12 = tsvar(12);
                    SymmetricMatrix.lcholesky(LVTS12);
                    LVTS12.smooth(1e-12);
                }
                return LVTS12.clone();
            case 4:
                if (LVTS4 == null) {
                    LVTS4 = tsvar(4);
                    LVTS4.smooth(1e-12);
                    SymmetricMatrix.lcholesky(LVTS4);
                }
                return LVTS4.clone();
            case 2:
                if (LVTS2 == null) {
                    LVTS2 = tsvar(2);
                    LVTS2.smooth(1e-12);
                    SymmetricMatrix.lcholesky(LVTS2);
                }
                return LVTS2.clone();
            case 3:
                if (LVTS3 == null) {
                    LVTS3 = tsvar(3);
                    SymmetricMatrix.lcholesky(LVTS3);
                    LVTS3.smooth(1e-12);
                }
                return LVTS3.clone();
            case 6:
                if (LVTS6 == null) {
                    LVTS6 = tsvar(6);
                    SymmetricMatrix.lcholesky(LVTS6);
                    LVTS6.smooth(1e-12);
                }
                return LVTS6.clone();
            default:
                Matrix var = tsvar(freq);
                SymmetricMatrix.lcholesky(var);
                var.smooth(1e-12);
                return var;
        }
    }

    public static Matrix tslVar(SeasonalModel seasModel, final int freq) {
        if (seasModel == SeasonalModel.Trigonometric) {
            return tslVar(freq);
        } else if (seasModel == SeasonalModel.HarrisonStevens) {
            return hslVar(freq);
        } else {
            int n = freq - 1;
            Matrix Q = Matrix.square(n);
            if (null != seasModel) // Dummy
            {
                switch (seasModel) {
                    case Dummy:
                        Q.set(n - 1, n - 1, 1);
                        break;
                    case Crude:
                        Q.set(1);
                        //Q.set(0, 0, freq * var);
                        break;
                    default:
                        break;
                }
            }
            return Q;
        }
    }

    private static Matrix VTS2, VTS3, VTS4, VTS6, VTS12;
    private static Matrix LVTS2, LVTS3, LVTS4, LVTS6, LVTS12, LHS2, LHS3, LHS4, LHS6, LHS12;

    public static class BsmDynamics implements ISsfDynamics {

        private final SubMatrix tsvar, ltsvar;
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
                if (model.seasModel != SeasonalModel.Crude && model.seasModel != SeasonalModel.Dummy) {
                    ltsvar = SsfBsm.tslVar(seasModel, freq).subMatrix();
                    ltsvar.mul(Math.sqrt(seasVar));
                } else {
                    ltsvar = null;
                }
            } else {
                tsvar = null;
                ltsvar = null;
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
                if (seasModel == SeasonalModel.Dummy || seasModel == SeasonalModel.Crude) {
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
        public boolean hasInnovations(int pos) {
            return true;
        }

        @Override
        public void S(int pos, SubMatrix s) {
            int i = 0, j = 0;
            if (nVar > 0) {
                s.set(i++, j++, Math.sqrt(nVar));
            }
            if (cVar > 0) {
                double ce = Math.sqrt(cVar);
                s.set(i++, j++, ce);
                s.set(i++, j++, ce);
            } else if (cVar == 0) {
                i += 2;
            }
            if (lVar > 0) {
                s.set(i++, j++, Math.sqrt(lVar));
            } else if (lVar == 0) {
                ++i;
            }
            if (sVar > 0) {
                s.set(i++, j++, Math.sqrt(sVar));
            } else if (sVar == 0) {
                ++i;
            }
            if (seasVar > 0) {
                if (seasModel == SeasonalModel.Dummy) {
                    s.set(i, j, Math.sqrt(seasVar));
                } else if (seasModel == SeasonalModel.Crude) {
                    s.extract(i, i + freq - 1, j, j + 1).set(Math.sqrt(seasVar));

                } else {
                    s.extract(i, i + freq - 1, j, j + freq - 1).copy(ltsvar);
                }
            }
        }

        @Override
        public void addSU(int pos, DataBlock x, DataBlock u) {
            int i = 0, j = 0;
            if (nVar > 0) {
                x.add(i++, u.get(j++) * Math.sqrt(nVar));
            }
            if (cVar > 0) {
                double ce = Math.sqrt(cVar);
                x.add(i++, u.get(j++) * ce);
                x.add(i++, u.get(j++) * ce);
            } else if (cVar == 0) {
                i += 2;
            }
            if (lVar > 0) {
                x.add(i++, u.get(j++) * Math.sqrt(lVar));
            } else if (lVar == 0) {
                ++i;
            }
            if (sVar > 0) {
                x.add(i++, u.get(j++) * Math.sqrt(sVar));
            } else if (sVar == 0) {
                ++i;
            }
            if (seasVar > 0) {
                if (seasModel == SeasonalModel.Dummy) {
                    x.add(i, u.get(j) * Math.sqrt(seasVar));
                } else if (seasModel == SeasonalModel.Crude) {
                    x.range(i, i + freq - 1).add(Math.sqrt(seasVar) * u.get(j));
                } else {
                    x.range(i, i + freq - 1).addProduct(ltsvar.rows(), u.range(j, j + freq - 1));
                }
            }
        }

        @Override
        public void XS(int pos, DataBlock x, DataBlock xs) {
            int i = 0, j = 0;
            if (nVar > 0) {
                xs.set(j++, x.get(i++) * Math.sqrt(nVar));
            }
            if (cVar > 0) {
                double ce = Math.sqrt(cVar);
                xs.set(j++, x.get(i++) * ce);
                xs.set(j++, x.get(i++) * ce);
            } else if (cVar == 0) {
                i += 2;
            }
            if (lVar > 0) {
                xs.set(j++, x.get(i++) * Math.sqrt(lVar));
            } else if (lVar == 0) {
                ++i;
            }
            if (sVar > 0) {
                xs.set(j++, x.get(i++) * Math.sqrt(lVar));
            } else if (sVar == 0) {
                ++i;
            }
            if (seasVar > 0) {
                if (seasModel == SeasonalModel.Dummy) {
                    xs.set(j, x.get(i) * Math.sqrt(seasVar));
                } else if (seasModel == SeasonalModel.Crude) {
                    xs.set(j, x.range(i, i + freq - 1).sum() * Math.sqrt(seasVar));
                } else {
                    xs.range(j, j + freq - 1).product(x.range(i, i + freq - 1), ltsvar.columns());
                }
            }
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
