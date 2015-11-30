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
package ec.tstoolkit2.ssf.akf;

import ec.tstoolkit.eco.ILikelihood;
import ec.tstoolkit2.ssf.univariate.ILikelihoodComputer;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfData;
import ec.tstoolkit2.ssf.univariate.OrdinaryFilter;
import ec.tstoolkit2.ssf.univariate.PredictionErrorDecomposition;

/**
 *
 * @author Jean Palate
 */
public class AkfToolkit {
    private AkfToolkit(){}

    public static ILikelihoodComputer likelihoodComputer(boolean collapsing) {
        return collapsing ? new LLComputer2() : new LLComputer1();
    }

    private static class LLComputer1 implements ILikelihoodComputer {

        @Override
        public ILikelihood compute(ISsf ssf, ISsfData data) {
            AugmentedFilter akf = new AugmentedFilter();
            AugmentedPredictionErrorDecomposition pe = new AugmentedPredictionErrorDecomposition();
            if (!akf.process(ssf, data, pe)) {
                return null;
            }
            return pe.likelihood();
        }

    }

    private static class LLComputer2 implements ILikelihoodComputer {

        @Override
        public ILikelihood compute(ISsf ssf, ISsfData data) {
            AugmentedPredictionErrorDecomposition ipe = new AugmentedPredictionErrorDecomposition();
            AugmentedFilterInitializer initializer = new AugmentedFilterInitializer(ipe);
            OrdinaryFilter filter = new OrdinaryFilter(initializer);
            PredictionErrorDecomposition pe = new PredictionErrorDecomposition(false);
            filter.process(ssf, data, pe);
            AkfDiffuseLikelihood ll = ipe.likelihood();
            ll.add(pe.likelihood());
            return ll;
        }

    }
}
