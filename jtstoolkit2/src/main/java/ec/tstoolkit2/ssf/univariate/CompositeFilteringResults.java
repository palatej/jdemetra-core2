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
package ec.tstoolkit2.ssf.univariate;

import ec.tstoolkit.sarima.SarimaModel;
import ec.tstoolkit.sarima.SarimaSpecification;
import ec.tstoolkit2.ssf.State;

/**
 *
 * @author Jean Palate
 */
public class CompositeFilteringResults implements IFilteringResults {
    
    
    private final IFilteringResults[] subresults;
    public CompositeFilteringResults(final IFilteringResults... subresults){
        this.subresults=subresults;
    }

    @Override
    public void save(int t, PredictionError pe) {
        for (int i=0; i<subresults.length; ++i){
            subresults[i].save(t, pe);
        }
    }

    @Override
    public void save(int t, State state) {
        for (int i=0; i<subresults.length; ++i){
            subresults[i].save(t, state);
        }
    }

    @Override
    public void clear() {
        for (int i=0; i<subresults.length; ++i){
            subresults[i].clear();
        }
    }
    
}
