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
import ec.tstoolkit2.ssf.implementations.Measurement;
import ec.tstoolkit2.ssf.univariate.Ssf;

/**
 * Usual local linear trend y(t)=l(t)+n(t) l(t+1)=s(t)+l(t)+u(t)
 * s(t+1)=s(t)+v(t)
 *
 * @author Jean Palate
 */
public class LocalLinearTrend extends Ssf {

    private final double lv, sv, nv;

    public LocalLinearTrend(double lvar, double svar, double nvar) {
        super(new Dynamics(lvar, svar), Measurement.create(0, nvar));
        lv = lvar;
        sv = svar;
        nv = nvar;
    }

    public double getVariance() {
        return lv;
    }

    public double getSlopeVariance() {
        return sv;
    }

    public double getNoiseVariance() {
        return nv;
    }

    public static class Dynamics implements ISsfDynamics {

        private final double lvar, svar;

        public Dynamics(double var, double svar) {
            lvar = var;
            this.svar = svar;
        }

        public double getVariance() {
            return lvar;
        }

        public double getSlopeVariance() {
            return svar;
        }

        @Override
        public int getStateDim() {
            return 2;
        }

        @Override
        public boolean isTimeInvariant() {
            return true;
        }

        @Override
        public boolean isValid() {
            return lvar >= 0 && svar >= 0;
        }

        @Override
        public int getInnovationsDim() {
            int n = 0;
            if (lvar > 0) {
                ++n;
            }
            if (svar > 0) {
                ++n;
            }
            return n;
        }

        @Override
        public void V(int pos, SubMatrix qm) {
            if (lvar > 0) {
                qm.set(0, 0, lvar);
            }
            if (svar > 0) {
                qm.set(1, 1, svar);
            }
        }

        @Override
        public boolean hasInnovations(int pos) {
            return lvar != 0 || svar != 0;
        }

//        @Override
//        public void Q(int pos, SubMatrix qm) {
//            int i = 0;
//            if (lvar > 0) {
//                qm.set(0, 0, lvar);
//                i = 1;
//            }
//            if (svar > 0) {
//                qm.set(i, i, svar);
//            }
//        }
//
//        @Override
//        public void S(int pos, SubMatrix sm) {
//            if (svar == 0 && lvar != 0) {
//                sm.set(1, 0, 1);
//            } else if (svar != 0 && lvar == 0) {
//                sm.set(0, 0, 1);
//            }
//        }
        @Override
        public void S(int pos, SubMatrix s) {
            //TODO
        }

        @Override
        public void addSU(int pos, DataBlock x, DataBlock u) {
            //TODO
        }
        
        @Override
        public void XS(int pos, DataBlock x, DataBlock xs) {
            //TODO
        }
//        @Override
//        public void addSX(int pos, DataBlock x, DataBlock y) {
//            if (svar == 0 && lvar != 0) {
//                y.add(1, x.get(0));
//            } else if (svar != 0 && lvar == 0) {
//                y.add(0, x.get(0));
//            } else
//            y.add(x);
//        }
//

        @Override
        public void T(int pos, SubMatrix tr) {
            tr.set(0, 0, 1);
            tr.set(0, 1, 1);
            tr.set(1, 1, 1);
        }

        @Override
        public boolean isDiffuse() {
            return true;
        }

        @Override
        public int getNonStationaryDim() {
            return 2;
        }

        @Override
        public void diffuseConstraints(SubMatrix b) {
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
        public void Pi0(SubMatrix pi0) {
            pi0.diagonal().set(1);
        }

        @Override
        public void TX(int pos, DataBlock x) {
            x.add(0, x.get(1));
        }

        @Override
        public void TVT(int pos, SubMatrix vm) {
            double v01 = vm.get(0, 1), v11 = vm.get(1, 1);
            vm.add(0, 0, 2 * v01 + v11);
            vm.add(0, 1, v11);
            vm.add(1, 0, v11);
        }

        @Override
        public void XT(int pos, DataBlock x) {
            x.add(1, x.get(0));
        }

        @Override
        public void addV(int pos, SubMatrix p) {
            if (lvar > 0) {
                p.add(0, 0, lvar);
            }
            if (svar > 0) {
                p.add(1, 1, svar);
            }
        }

    }
}
