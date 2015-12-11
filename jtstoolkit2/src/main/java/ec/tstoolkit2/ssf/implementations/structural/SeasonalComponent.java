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
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.implementations.Measurement;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import ec.tstoolkit2.ssf.univariate.Ssf;

/**
 *
 * @author Jean Palate
 */
public class SeasonalComponent {
    
    public static ISsf create(final SeasonalModel model, final double seasVar, final int period){
        return new Ssf(new Dynamics(model, seasVar, period), Measurement.create(1));
    }

    public static ISsf harrisonStevens(final int period, final double v) {
        return new Ssf(new HarrisonStevensDynamics(period, v), Measurement.cyclical(period));
    }

    public static ISsf harrisonStevens(final double[] var) {
        return new Ssf(new HarrisonStevensDynamics(var), Measurement.cyclical(var.length));
    }
    

    static class Dynamics implements ISsfDynamics {

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


    public static class HarrisonStevensDynamics implements ISsfDynamics {

        private final int period;
        private final double[] var;
        private final Matrix V;

        public HarrisonStevensDynamics(final int period, final double v) {
            this.period = period;
            var = null;
            V = Matrix.square(period - 1);
            V.set(-1.0 / period);
            V.diagonal().add(1);
            V.mul(v);
        }

        public HarrisonStevensDynamics(final double[] var) {
            period = var.length;
            this.var = var.clone();
            Matrix C = new Matrix(period - 1, period);
            C.set(-1.0 / period);
            C.diagonal().add(1);
            Matrix D = Matrix.diagonal(var);
            V = SymmetricMatrix.quadraticFormT(D, C);
        }

        public double[] getVariances() {
            return var;
        }

        @Override
        public int getStateDim() {
            return period - 1;
        }

        @Override
        public boolean isTimeInvariant() {
            return true;
        }

        @Override
        public boolean isValid() {
            return period > 1;
        }

        @Override
        public int getInnovationsDim() {
            return period - 1;
        }

        @Override
        public void V(int pos, SubMatrix qm) {
            qm.copy(V.subMatrix());
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
        public void Q(int pos, SubMatrix qm) {
            qm.copy(V.subMatrix());
        }

        @Override
        public void S(int pos, SubMatrix sm) {
        }

        @Override
        public void T(int pos, SubMatrix tr) {
            tr.diagonal().set(1);
        }

        @Override
        public boolean isDiffuse() {
            return true;
        }

        @Override
        public int getNonStationaryDim() {
            return period - 1;
        }

        @Override
        public void diffuseConstraints(SubMatrix b) {
            b.diagonal().set(1);
        }

        @Override
        public void Pi0(SubMatrix b) {
            b.diagonal().set(1);
        }

        @Override
        public boolean a0(DataBlock a0, StateInfo info) {
            return true;
        }

        @Override
        public boolean Pf0(SubMatrix pf0, StateInfo info) {
            return true;
        }

        @Override
        public void TX(int pos, DataBlock x) {
        }

        @Override
        public void XT(int pos, DataBlock x) {
        }

        @Override
        public void addV(int pos, SubMatrix p) {
            p.add(V.subMatrix());
        }

    }
}
