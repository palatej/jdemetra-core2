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
package ec.tstoolkit2.ssf.implementations;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.IReadDataBlock;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.StateInfo;

/**
 *
 * @author Jean Palate
 */
public class TimeVaryingDynamics {

    public static ISsfDynamics create(IReadDataBlock dvar) {
        return new TimeVaryingDiag(dvar);
    }

    public static ISsfDynamics create(SubMatrix var) {
        return new TimeVaryingFull(var);
    }

    static class TimeVaryingDiag implements ISsfDynamics {

        private final DataBlock var;

        TimeVaryingDiag(final double[] var) {
            this.var = new DataBlock(var);
        }

        TimeVaryingDiag(final IReadDataBlock var) {
            this.var = new DataBlock(var);
        }

        @Override
        public int getStateDim() {
            return var.getLength();
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
            return var.getLength();
        }

        @Override
        public void V(int pos, SubMatrix qm) {
            qm.diagonal().copy(var);
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
            qm.diagonal().copy(var);
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
            tr.diagonal().set(1);
        }

        @Override
        public boolean isDiffuse() {
            return true;
        }

        @Override
        public int getNonStationaryDim() {
            return var.getLength();
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
        public void Pi0(SubMatrix p) {
            p.diagonal().set(1);
        }

        @Override
        public void TX(int pos, DataBlock x) {
        }

        @Override
        public void XT(int pos, DataBlock x) {
        }

        @Override
        public void TVT(int pos, SubMatrix v) {
        }

        @Override
        public void addV(int pos, SubMatrix p) {
            p.diagonal().add(var);
        }
    }

    static class TimeVaryingFull implements ISsfDynamics {

        private final SubMatrix var;

        TimeVaryingFull(final SubMatrix var) {
            this.var = var;
        }

        @Override
        public int getStateDim() {
            return var.getColumnsCount();
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
            return var.getColumnsCount();
        }

        @Override
        public void V(int pos, SubMatrix qm) {
            qm.copy(var);
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
            qm.copy(var);
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
            tr.diagonal().set(1);
        }

        @Override
        public boolean isDiffuse() {
            return true;
        }

        @Override
        public int getNonStationaryDim() {
            return var.getColumnsCount();
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
        public void Pi0(SubMatrix p) {
            p.diagonal().set(1);
        }

        @Override
        public void TX(int pos, DataBlock x) {
        }

        @Override
        public void XT(int pos, DataBlock x) {
        }

        @Override
        public void TVT(int pos, SubMatrix v) {
        }

        @Override
        public void addV(int pos, SubMatrix p) {
            p.add(var);
        }
    }
}
