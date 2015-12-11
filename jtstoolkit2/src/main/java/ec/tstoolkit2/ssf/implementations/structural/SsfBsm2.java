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
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.implementations.Measurement;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import ec.tstoolkit2.ssf.univariate.Ssf;

/**
 *
 * @author Jean Palate
 */
public class SsfBsm2 extends Ssf {

    private SsfBsm2(Bsm2Dynamics dynamics, ISsfMeasurement measurement) {
        super(dynamics, measurement);
    }

    /**
     *
     */
    private static int[] calcCmpsIndexes(BasicStructuralModel model) {
        int n = 0;
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

    public static SsfBsm2 create(BasicStructuralModel model) {
        int[] idx = calcCmpsIndexes(model);
        ISsfMeasurement measurement = Measurement.create(idx, model.nVar);
        Bsm2Dynamics dynamics = new Bsm2Dynamics(model);
        if (dynamics.isValid()) {
            return new SsfBsm2(dynamics, measurement);
        } else {
            return null;
        }
    }

    public static class Bsm2Dynamics implements ISsfDynamics {

        private final SubMatrix tsvar;
        private final double lVar, sVar, seasVar, cVar, cDump;
        private final double ccos, csin;
        private final int freq;
        private final SeasonalModel seasModel;

        public Bsm2Dynamics(BasicStructuralModel model) {
            lVar = model.lVar;
            sVar = model.sVar;
            seasVar = model.seasVar;
            cVar = model.cVar;
            cDump = model.cDump;
            ccos = model.ccos;
            csin = model.csin;
            seasModel = model.seasModel;
            freq = model.freq;
            if (seasVar > 0) {
                tsvar = SsfBsm.tsVar(seasModel, freq).subMatrix();
                tsvar.mul(seasVar);
            } else {
                tsvar = null;
            }
        }

        @Override
        public int getStateDim() {
            int r = 0;
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
            return lVar >= 0 || sVar >= 0 || cVar >= 0;
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
            int istart = 0;
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
            int istart = 0;
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
                if (cVar > 0) {
                    double q = cVar / (1 - cDump * cDump);
                    p.set(i, i, q);
                    ++i;
                    p.set(i, i, q);
                }
                return true;
            }else{
                // TODO
                return false;
            }
        }

        @Override
        public void TX(int pos, DataBlock x) {
            int i0 = 0;
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
