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
package ec.tstoolkit2.ssf.implementations.arima;

import ec.tstoolkit.arima.ArimaModel;
import ec.tstoolkit.ucarima.UcarimaModel;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.implementations.CompositeDynamics;
import ec.tstoolkit2.ssf.implementations.Measurement;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import ec.tstoolkit2.ssf.univariate.Ssf;

/**
 *
 * @author Jean Palate
 */
public class SsfUcarima extends Ssf{
    
    public static SsfUcarima create(final UcarimaModel ucm){
        ucm.simplify();
        ISsfDynamics[] dyn=new ISsfDynamics[ucm.getComponentsCount()] ;
        int[] pos=new int[dyn.length];
        pos[0]=0;
        for (int i=0; i<dyn.length; ++i){
            ArimaModel cmp = ucm.getComponent(i);
            dyn[i]= cmp.isStationary()? new SsfArima.StDynamics(cmp)
                    :new SsfArima.SsfArimaDynamics(cmp);
            if (i>0){
                pos[i]=pos[i-1]+dyn[i-1].getStateDim();
            }                
        }
        return new SsfUcarima(ucm, new CompositeDynamics(dyn), Measurement.create(pos));
    }

    private final UcarimaModel ucm;
    
    private SsfUcarima(UcarimaModel ucm, ISsfDynamics dyn, ISsfMeasurement m){
        super(dyn, m);
        this.ucm=ucm;
    }
    
    public UcarimaModel getUcarimaModel(){
        return ucm;
    }
}
