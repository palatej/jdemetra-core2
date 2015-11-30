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
package ec.tstoolkit2.ssf.multivariate;

import ec.tstoolkit2.ssf.multivariate.IPredictionErrors;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.Matrix;

/**
 *
 * @author Jean Palate
 */
public class PredictionErrors implements IPredictionErrors {

    /**
     * E is the prediction error (=Y(t)-Z(t)A(t))
     */
    public final DataBlock E;

    /**
     * =(ZPZ'+H) variance/covariance matrix of the
     * prediction errors 
     */
    public final Matrix F;
    
    /**
     * K = P Z' F^-1
     */
    public final Matrix K;
    

    /**
     *
     * @param dim
     * @param nvars
     */
    public PredictionErrors(final int dim, final int nvars) {
        E = new DataBlock(nvars);
        F = new Matrix(nvars, nvars);
        K = new Matrix(dim, nvars);
    }

    @Override
    public DataBlock getPredictionError() {
        return E;
    }

    @Override
    public Matrix getPredictionErrorCovariance() {
        return F;
    }

}
