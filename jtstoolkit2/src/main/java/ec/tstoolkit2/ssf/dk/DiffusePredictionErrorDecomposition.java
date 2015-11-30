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
package ec.tstoolkit2.ssf.dk;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.eco.Determinant;
import ec.tstoolkit.eco.ILikelihood;
import ec.tstoolkit2.ssf.IPredictionErrorDecomposition;
import ec.tstoolkit2.ssf.ResidualsCumulator;
import ec.tstoolkit2.ssf.akf.AugmentedState;
import ec.tstoolkit2.ssf.dk.sqrt.IDiffuseSquareRootFilteringResults;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfData;
import ec.tstoolkit2.ssf.univariate.PredictionErrorDecomposition;

/**
 *
 * @author Jean Palate
 */
public class DiffusePredictionErrorDecomposition extends PredictionErrorDecomposition implements IDiffuseFilteringResults, IDiffuseSquareRootFilteringResults {

    private final Determinant ddet = new Determinant();
    private int nd, n;

    public DiffusePredictionErrorDecomposition(boolean res) {
        super(res);
    }

    @Override
    public ILikelihood likelihood() {
        DkDiffuseLikelihood ll = new DkDiffuseLikelihood();
        int nobs = nd + cumulator.getObsCount();
        ll.set(cumulator.getSsqErr(), cumulator.getLogDeterminant(), ddet.getLogDeterminant(), nobs, nd);
        if (bres) {
//            if (n < res.getLength()) {
//                double[] tmp = new double[n];
//                System.arraycopy(res, 0, tmp, 0, n);
//                ll.setResiduals(tmp);
//            } else {
                ll.setResiduals(res.getData());
//            }
        }
        return ll;
    }

    @Override
    public void close(int pos) {
    }

    @Override
    public void clear() {
        cumulator.clear();
        ddet.clear();
        nd = 0;
        n = 0;
        res = null;
    }

    @Override
    public void save(int t, DiffusePredictionError pe) {
        if (t + 1 > n) {
            n = t + 1;
        }
        if (pe == null || pe.isMissing()) {
            return;
        }
        double d = pe.getDiffuseNorm2();
        if (d != 0) {
            ++nd;
            ddet.add(d);
        } else {
            double e = pe.get();
            cumulator.add(e, pe.getVariance());
            if (bres) {
                res.set(t, e / pe.getStandardDeviation());
            }
        }
    }

    @Override
    public void save(int t, DiffuseState state) {
    }

    @Override
    public void save(int t, AugmentedState state) {
    }

}
