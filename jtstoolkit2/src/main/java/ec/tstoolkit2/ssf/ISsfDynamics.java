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
     * Another (non unique) modelling is
     *
     * V = S*S'
     * 
     * Wher the m columns of S are independent. m corresponds to innovationsDim 
     * This modelling is useful in square root algorithms. It is used for instance
     * in De Jong. De dimension of the U matrix
     *
     * @param pos
     * @param qm
     */
    void V(int pos, SubMatrix qm);

    /**
     * @param pos
     * @param cm
     */
    void S(int pos, SubMatrix cm);

    /**
     *
     * @param pos
     * @return
     */
    boolean hasInnovations(int pos);


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
     * Computes xs = x*S(pos)
     *
     * @param pos
     * @param x 
     * @param xs 
     */
    void XS(int pos, DataBlock x, DataBlock xs);

    /**
     * Computes x=x+ S(pos) * u. The dimension of u should correspond to the 
     * number of columns of S(po)
     *
     * @param pos
     * @param x 
     * @param u 
     */
    void addSU(int pos, DataBlock x, DataBlock u);

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
