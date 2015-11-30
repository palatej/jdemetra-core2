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
package ec.tstoolkit2.ssf.implementations.structural;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.StateInfo;

/**
 *
 * @author Jean Palate
 */
public class SeasonalComponent {

    public static class Dynamics implements ISsfDynamics {

        private final SeasonalModel seasModel;
        private final double seasVar;
        private final int freq;
        private final SubMatrix tsvar;

        public Dynamics(final SeasonalModel model, final double seasVar, final int freq) {
            this.seasVar = seasVar;
            this.seasModel = model;
            this.freq = freq;
            if (seasVar > 0) {
                tsvar = SsfBsm.tsVar(seasModel, freq).subMatrix();
                tsvar.mul(seasVar);
            } else {
                tsvar = null;
            }
        }

        @Override
        public int getStateDim() {
            return freq - 1;
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
            if (seasVar > 0) {
                if (seasModel == SeasonalModel.Dummy) {
//                    || seasModel == SeasonalModel.Crude) {
                    return 1;
                } else {
                    return freq - 1;
                }
            } else {
                return 0;
            }
        }

        @Override
        public void V(int pos, SubMatrix v) {
            if (seasVar > 0) {
                if (seasModel == SeasonalModel.Dummy) {
                    v.set(0, 0, seasVar);
                } else {
                    v.copy(tsvar);
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
            if (seasVar >= 0) {
                tr.row(freq - 2).set(-1);
                tr.subDiagonal(1).set(1);
            }
        }

        @Override
        public boolean isDiffuse() {
            return seasVar >= 0;
        }

        @Override
        public int getNonStationaryDim() {
            return freq - 1;
        }

        @Override
        public void diffuseConstraints(SubMatrix b) {
            b.diagonal().set(1);
        }

        @Override
        public void Pi0(SubMatrix p) {
            p.diagonal().set(1);
        }

        @Override
        public boolean a0(DataBlock a0, StateInfo info) {
            return true;
        }

        @Override
        public boolean Pf0(SubMatrix p, StateInfo info) {
            if (info == StateInfo.Forecast) {
                return true;
            } else {
                // TODO
                return false;
            }
        }

        @Override
        public void TX(int pos, DataBlock x) {
            x.bshift(DataBlock.ShiftOption.NegSum);
        }

        @Override
        public void XT(int pos, DataBlock x) {
            int imax = freq - 2;
            double xs = x.get(imax);
            for (int i = imax; i > 0; --i) {
                x.set(i, x.get(i - 1) - xs);
            }
            x.set(0, -xs);

        }

        @Override
        public void addV(int pos, SubMatrix p) {
            if (seasModel == SeasonalModel.Dummy) {
                p.add(0, 0, seasVar);
            } else {
                p.add(tsvar);
            }
        }

    }
}
