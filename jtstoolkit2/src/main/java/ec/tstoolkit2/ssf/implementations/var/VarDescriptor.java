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
package ec.tstoolkit2.ssf.implementations.var;

import ec.tstoolkit.maths.matrices.Matrix;

/**
 *
 * @author Jean Palate
 */
public class VarDescriptor {

    public static final double AR_DEF = .6;

    /**
     * Number of lags
     */
    private final int nlags, neq;
    /**
     * Parameters of the VAR equations The row i contains the coefficients
     * c(i,k) of fi(t): fi(t)= c(i,0)f0(t-1)+...+c(i,nlags-1)f0(t-nlags)+...
     * +c(i,k)fn(t-1)...+c(i,l)fn(t-nlags))
     */
    private final Matrix varMatrix;
    /**
     * Covariance matrix of the innovations
     */
    private final Matrix covar;
    /**
     * Creates a new descriptor of the transition equation (VAR).
     *
     * @param neq Number of equations
     * @param nlags Number of lags in the VAR model
     */
    public VarDescriptor(int neq, int nlags) {
        varMatrix = new Matrix(neq, neq * nlags);
        covar = new Matrix(neq, neq);
        this.nlags = nlags;
        this.neq=neq;
        setDefault();
    }

    public final void setDefault() {
        covar.set(0);
        covar.diagonal().set(1);
        varMatrix.set(0);
        for (int i = 0; i < varMatrix.getRowsCount(); ++i) {
            varMatrix.set(i, i * nlags, AR_DEF);
        }
    }
    
    public int getLagsCount(){
        return nlags;
    }
    
    public int getEquationsCount(){
        return neq;
    }
    
    public Matrix getVarMatrix(){
        return varMatrix;
    }
    
    public Matrix getInnovationsVariance(){
        return covar;
    }

    /**
     * Gets the matrix of the var parameters corresponding to a given lag
     *
     * @param lag The lag in the var equation. Should belong to [1, nlags]
     * @return A new matrix is returned
     */
    public Matrix getA(int lag) {
        int n = varMatrix.getRowsCount();
        Matrix a = new Matrix(n, n);
        for (int i = 0, j = lag - 1; i < n; ++i, j += nlags) {
            a.column(i).copy(varMatrix.column(j));
        }
        return a;
    }

    /**
     * Sets the matrix of the var parameters corresponding to a given lag
     *
     * @param lag The lag in the var equation. Should belong to [1, nlags]
     * @param a The matrix
     */
    public void setA(int lag, Matrix a) {
        int n = varMatrix.getRowsCount();
        for (int i = 0, j = lag - 1; i < n; ++i, j += nlags) {
            varMatrix.column(j).copy(a.column(i));
        }
    }
}

