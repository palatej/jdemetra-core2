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
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.StateInfo;

/**
 *
 * @author Jean Palate
 */
public class ConstantDynamics implements ISsfDynamics {

    private final int dim;

    public ConstantDynamics(final int dim) {
        this.dim = dim;
    }

    @Override
    public int getStateDim() {
        return dim;
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
        return 0;
    }

    @Override
    public void V(int pos, SubMatrix qm) {
    }

    @Override
    public boolean hasS() {
        return false;
    }

    @Override
    public boolean hasInnovations(int pos) {
        return false;
    }

    @Override
    public void Q(int pos, SubMatrix qm) {
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
        return dim;
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
    }

}
