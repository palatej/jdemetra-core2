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
package ec.tstoolkit2.ssf.implementations.structural;

import ec.tstoolkit2.ssf.implementations.Measurement;
import ec.tstoolkit2.ssf.implementations.arima.SsfRandomWalk;
import ec.tstoolkit2.ssf.univariate.Ssf;

/**
 * Usual local linear trend 
 * y(t)=l(t)+n(t)
 * l(t+1)=l(t)+u(t) 
 * @author Jean Palate
 */
public class LocalLevel extends Ssf{
    private final double lv,nv;

    public LocalLevel(double lvar, double nvar) {
        super(new SsfRandomWalk.Dynamics(lvar, false), Measurement.create(0, nvar));
        lv=lvar;
        nv=nvar;        
    }

    public double getVariance() {
        return lv;
    }

    public double getNoiseVariance() {
        return nv;
    }

}
