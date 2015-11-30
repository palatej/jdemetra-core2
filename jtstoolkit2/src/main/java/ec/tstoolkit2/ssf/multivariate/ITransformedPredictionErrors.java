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
package ec.tstoolkit2.ssf.multivariate;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.Matrix;

/**
 * Cholesky form of the prediction error.
 * Suppose that the errors are 
 * E(t) ~ N(0, F) and that LL'=F
 * U(t) = L(t)^-1* E(t) (or L(t)U(t)=E(t)) is then ~ N(0, I).
 * This class provides U(t) and L(t) 
 * @author Jean Palate
 */
public interface ITransformedPredictionErrors {
    /**
     * Provides U(t). See above for the definition of U(t)
     * @return 
     */
    DataBlock getTransformedPredictionErrors();
    /**
     * Provides L(t). See above for the definition of L(t)
     * @return 
     */
    Matrix getCholeskyFactor();
}
