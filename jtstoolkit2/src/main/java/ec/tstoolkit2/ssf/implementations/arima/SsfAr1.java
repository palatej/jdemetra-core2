/*
 * Copyright 2015 National Bank of Belgium
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
/*
 */
package ec.tstoolkit2.ssf.implementations.arima;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.implementations.Measurement;
import ec.tstoolkit2.ssf.univariate.Ssf;

/**
 *
 * @author Jean Palate
 */
public class SsfAr1 extends Ssf{
    
    public SsfAr1(final double rho, final double var, final boolean zeroinit){
        super (new Dynamics(rho, var, zeroinit), Measurement.create(0));
    }
    
    private Dynamics dynamics(){
        return (Dynamics) this.dynamics;
    }
    
    public double getRho(){
        return dynamics().rho;
    }

    public double getInnovationVariance(){
        return dynamics().var;
    }
    
    public boolean isZeroInitialization(){
        return dynamics().zeroinit;
    }

    public static class Dynamics implements ISsfDynamics {

        private final boolean zeroinit;
        private final double rho;
        private final double var;

        public Dynamics(double rho, double var, boolean zeroinit) {
            this.rho = rho;
            this.var = var;
            this.zeroinit = zeroinit;
        }

        public boolean isZeroInit() {
            return zeroinit;
        }

        public double getRho() {
            return rho;
        }

        public double getVar() {
            return var;
        }

        @Override
        public int getStateDim() {
            return 1;
        }

        @Override
        public boolean isTimeInvariant() {
            return true;
        }

        @Override
        public boolean isValid() {
            return var > 0;
        }

        @Override
        public int getInnovationsDim() {
            return 1;
        }

        @Override
        public void V(int pos, SubMatrix qm) {
            qm.set(0, 0, var);
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
            qm.set(0, 0, var);
        }

        @Override
        public void S(int pos, SubMatrix sm) {
        }

//        @Override
//        public void addSX(int pos, DataBlock x, DataBlock y) {
//            y.add(x);
//        }
//
        @Override
        public void T(int pos, SubMatrix tr) {
            tr.set(0, 0, rho);
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
            if (zeroinit) {
                if (info == StateInfo.Forecast)
                    pf0.set(0, 0, var);
            } else {
                pf0.set(0, 0, var / (1 - rho * rho));
            }
            return true;
        }

        @Override
        public void TX(int pos, DataBlock x) {
            x.mul(0,rho);
        }

        @Override
        public void TVT(int pos, SubMatrix v) {
            v.mul(0, 0, rho*rho);
        }

        @Override
        public void XT(int pos, DataBlock x) {
            x.mul(0,rho);
        }

        @Override
        public void addV(int pos, SubMatrix p) {
            p.add(0, 0, var);
        }
    }
}
