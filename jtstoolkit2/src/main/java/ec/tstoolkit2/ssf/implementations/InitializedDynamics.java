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
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;

/**
 *
 * @author Jean Palate
 */
public final class InitializedDynamics implements ISsfDynamics {

    private final ISsfDynamics dyn;
    private final State start;
    private final int startpos;

    public InitializedDynamics(ISsfDynamics cur, State start, int startpos) {
        this.dyn = cur;
        this.start = start;
        this.startpos = startpos;
    }

    @Override
    public int getStateDim() {
        return dyn.getStateDim();
    }

    @Override
    public boolean isTimeInvariant() {
        return dyn.isTimeInvariant();
    }

    @Override
    public boolean isValid() {
        return dyn.isValid();
    }

    @Override
    public int getInnovationsDim() {
        return dyn.getInnovationsDim();
    }

    @Override
    public void V(int pos, SubMatrix qm) {
        dyn.V(pos + startpos, qm);
    }

    @Override
    public boolean hasS() {
        return dyn.hasS();
    }

    @Override
    public boolean hasInnovations(int pos) {
        return dyn.hasInnovations(pos + startpos);
    }

    @Override
    public void Q(int pos, SubMatrix qm) {
        dyn.Q(pos + startpos, qm);
    }

    @Override
    public void S(int pos, SubMatrix sm) {
        dyn.S(pos + startpos, sm);
    }

    @Override
    public void T(int pos, SubMatrix tr) {
        dyn.T(pos + startpos, tr);
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
        if (info != start.getInfo()) {
            return false;
        }
        a0.copy(start.a());
        return true;
    }

    @Override
    public boolean Pf0(SubMatrix pf0, StateInfo info) {
        if (info != start.getInfo()) {
            return false;
        }
        pf0.copy(start.P().subMatrix());
        return true;
    }

    @Override
    public void TX(int pos, DataBlock x) {
        dyn.TX(pos + startpos, x);
    }

    @Override
    public void XT(int pos, DataBlock x) {
        dyn.XT(pos + startpos, x);
    }

    @Override
    public void addV(int pos, SubMatrix p) {
        dyn.addV(pos + startpos, p);
    }
}
