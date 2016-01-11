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
package ec.tstoolkit2.ssf;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.DataBlockIterator;
import ec.tstoolkit.maths.matrices.LowerTriangularMatrix;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;

/**
 *
 * @author Jean Palate
 */
public interface ISsfDynamics {

//<editor-fold defaultstate="collapsed" desc="description">
    /**
     * Dimension of the state vector
     *
     * @return
     */
    int getStateDim();

    /**
     * Is this object time invariant
     *
     * @return
     */
    boolean isTimeInvariant();

    /**
     *
     * @return
     */
    boolean isValid();

    /**
     * Dimension of the innovations. See V for further information
     *
     * @return
     */
    int getInnovationsDim(); // E

    /**
     * Variance matrix of the innovations in the transition equation. V is also
     * modelled as
     *
     * V = S*Q*S' 
     *
     * V = d x d where d = getStateDim() Q = q x q where q = getInnovationsDim()
     * S = d x q
     *
     * When S = I (d = q), it should be considered as missing.
     *
     * Another (non unique) modelling is
     *
     * V = U*U'
     * This modelling is useful in square root algorithms. It is used for instance
     * in De Jong.
     *
     * @param pos
     * @param qm
     */
    void V(int pos, SubMatrix qm);

    /**
     * The default implementation computes S*L where L*L' = Q
     * @param pos
     * @param cm
     */
    default void U(int pos, SubMatrix cm){
        if (hasS()){
            S(pos, cm);
            // C = S * L, where L*L'=Q
            int nres=getInnovationsDim();
            Matrix L=Matrix.square(nres);
            Q(pos, L.subMatrix());
            if (nres == 1){
                cm.mul(Math.sqrt(L.get(0, 0)));
            }else{
                SymmetricMatrix.lcholesky(L);
                LowerTriangularMatrix.lmul(L, cm);
            }
        }else{
            Matrix L=Matrix.square(getStateDim());
            V(pos, L.subMatrix());
            SymmetricMatrix.lcholesky(L, State.ZERO);
            cm.copy(L.subMatrix());
        }
    }

    /**
     * @return
     */
    boolean hasS();

    /**
     *
     * @param pos
     * @return
     */
    boolean hasInnovations(int pos);

    /**
     *
     * @param pos
     * @param qm
     */
    void Q(int pos, SubMatrix qm);

    /**
     *
     * @param pos
     * @param sm
     */
    void S(int pos, SubMatrix sm);

    /**
     * Gets the transition matrix.
     *
     * @param pos The position of the model
     * @param tr The sub-matrix that will receive the transition matrix. It must
     * have the dimensions (getStateDim() x getStateDim()). The caller has the
     * responsibility to provide a clean sub-matrix, so that the callee can
     * safely set only the non zero values.
     */
    void T(int pos, SubMatrix tr);

//</editor-fold>
    
//<editor-fold defaultstate="collapsed" desc="initialisation">
    /**
     *
     * @return
     */
    boolean isDiffuse();

    /**
     * Dimension of the non stationary part of the state vector
     *
     * @return
     */
    int getNonStationaryDim();

    /**
     * B = d x nd, where d = getStateDim(), nd = getNonStationaryDim() P(-1,
     * inf) = B * B'
     *
     * @param b
     */
    void diffuseConstraints(SubMatrix b);

    /**
     * Initial state
     *
     * @param a0 Buffer that will contain the initial state
     * @param info
     * @return 
     */
    boolean a0(DataBlock a0, StateInfo info);

    /**
     * Modelling of the stationary variance of the initial state P(-1, *)
     *
     * @param pf0
     * @param info
     * @return 
     */
    boolean Pf0(SubMatrix pf0, StateInfo info);

    /**
     * Modelling of the non stationary part of the initial state P(-1, inf)
     *
     * @param pi0
     * @return 
     */
    default void Pi0(SubMatrix pi0){
        int nd=this.getNonStationaryDim();
        if (nd ==0 )
            return;
        int n=this.getStateDim();
        SubMatrix B=new Matrix(n, nd).subMatrix();
        this.diffuseConstraints(B);
        SymmetricMatrix.XXt(B, pi0);
    }
    
//</editor-fold>    

//<editor-fold defaultstate="collapsed" desc="forward operations">
    /**
     * Computes T(pos) * x
     *
     * @param pos
     * @param x
     */
    void TX(int pos, DataBlock x);

    /**
     * Computes T(pos) * M
     *
     * @param pos
     * @param M
     */
    default void TM(int pos, SubMatrix M) {
        DataBlockIterator cols = M.columns();
        DataBlock col = cols.getData();
        do {
            TX(pos, col);
        } while (cols.next());
    }
    
    default boolean forecast(int pos, State state){
        if (state.getInfo() == StateInfo.Concurrent){
            TX(pos, state.a());
            SubMatrix P= state.P().subMatrix();
            TVT(pos,P);
            addV(pos, P);
            state.setInfo(StateInfo.Forecast);
            return true;
        }else{
            return false;
        }
    }

    /**
     * Computes T V T'
     *
     * @param pos The position of the model
     * @param vm
     */
    default void TVT(int pos, SubMatrix vm) {
        TM(pos, vm);
        TM(pos, vm.transpose());
    }

//</editor-fold>    
    
//<editor-fold defaultstate="collapsed" desc="backward operations">
    /**
     * Computes x * T(pos)
     *
     * @param pos
     * @param x
     */
    void XT(int pos, DataBlock x);

    /**
     * Computes M * T(pos)
     *
     * @param pos
     * @param M
     */
    default void MT(int pos, SubMatrix M) {
        DataBlockIterator rows = M.rows();
        DataBlock row = rows.getData();
        do {
            XT(pos, row);
        } while (rows.next());
    }
    

    /**
     * Adds the variance of the innovations to a given matrix p = p + V(pos)
     *
     * @param pos
     * @param p
     */
    void addV(int pos, SubMatrix p);
//</editor-fold>
    

}
