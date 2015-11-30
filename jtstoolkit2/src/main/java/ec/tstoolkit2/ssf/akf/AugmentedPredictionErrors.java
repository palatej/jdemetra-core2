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
package ec.tstoolkit2.ssf.akf;

import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit2.ssf.multivariate.TransformedPredictionErrors;

/**
 * The augmented state only contains information on non missing values.
 * So, it has to be considered with the corresponding data. 
 * @author Jean Palate
 */

public class AugmentedPredictionErrors extends TransformedPredictionErrors {

    /**
     * E is the "prediction error" on the diffuse constraints (=(0-Z(t)A(t))
     * E ~ ndiffuse x nvars
     */
    private final Matrix E;

     /**
     *
     * @param ndiffuse
     * @param nvars
     * @param dim
     */
    public AugmentedPredictionErrors(final int dim, final int nvars, final int ndiffuse) {
        super(dim, nvars);
        E = new Matrix(ndiffuse, nvars);
    }

    public Matrix E() {
        return E;
    }

    public boolean isDiffuse() {
         return !E.isZero();
    }


}
