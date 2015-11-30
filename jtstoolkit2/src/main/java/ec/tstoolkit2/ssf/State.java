/*
 * Copyright 2013 National Bank of Belgium
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
package ec.tstoolkit2.ssf;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.design.Development;
import ec.tstoolkit.maths.matrices.Matrix;

/**
 * Represents a gaussian vector, with its mean and covariance matrix.
 * The way information must be interpreted is given by the state info.
 * This is similar to the NRV (normal random vector) of Snyder/Forbes (apart from 
 * the additional info)
 * @author Jean Palate
 */
@Development(status = Development.Status.Alpha)
public class State {
    
    public static final double ZERO= 1e-9;
    
    public static State of(ISsfDynamics dyn, StateInfo info){
        State state=new State(dyn.getStateDim());
        if (! dyn.a0(state.a, info))
            return null;
        if (!dyn.Pf0(state.P.subMatrix(), info))
            return null;
        state.setInfo(info);
        return state;
    }


    private StateInfo info;
    /**
     * a is the state vector. Its interpretation depends on the considered step
     */
    private final DataBlock a;

    /**
     * P is the covariance of the state vector. Its interpretation depends on
     * the considered step
     */
    private final Matrix P;


    /**
     *
     *
     * @param dim
     */
    public State(final int dim) {
        a = new DataBlock(dim);
        P = Matrix.square(dim);
        info = StateInfo.Undefined;
    }

    public State(final DataBlock a, final Matrix P, final StateInfo info) {
        this.a = a;
        this.P = P;
        this.info = info;
    }
    /**
     *
     * @param state
     */
    public void copy(final State state) {
        a.copy(state.a);
        P.copy(state.P);
        setInfo(state.getInfo());
    }
    
    @Override
    public String toString(){
        StringBuilder builder=new StringBuilder();
        builder.append("mean:\r\n").append(a).append("\r\n");
        builder.append("covariance:\r\n").append(P);
        return builder.toString();
   }

    /**
     * @return the info
     */
    public final StateInfo getInfo() {
        return info;
    }

    /**
     * @param info the info to set
     */
    public final void setInfo(StateInfo info) {
        this.info = info;
    }
    
    public final int getDim(){
        return a.getLength();
    }

    /**
     * @return the a
     */
    public final DataBlock a() {
        return a;
    }

    /**
     * @return the P
     */
    public final Matrix P() {
        return P;
    }
    
}
