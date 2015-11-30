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

import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import ec.tstoolkit2.ssf.univariate.Ssf;

/**
 *
 * @author Jean Palate
 */
public class CompositeSsf extends Ssf{
    
    private CompositeSsf(ISsfDynamics dyn, ISsfMeasurement m){
        super(dyn, m);
    }
    
    public static CompositeSsf of(double var, ISsf... ssf ){
        ISsfMeasurement m=CompositeMeasurement.of(var, ssf);
        if (m == null)
            return null;
        ISsfDynamics dyn=CompositeDynamics.of(ssf);
        if (dyn == null)
            return null;
        return new CompositeSsf(dyn, m);
    }
    
}
