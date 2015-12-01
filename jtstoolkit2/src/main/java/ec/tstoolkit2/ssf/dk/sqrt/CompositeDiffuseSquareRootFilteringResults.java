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
package ec.tstoolkit2.ssf.dk.sqrt;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.univariate.*;
import ec.tstoolkit2.ssf.akf.AugmentedState;
import ec.tstoolkit2.ssf.dk.DiffusePredictionError;

/**
 *
 * @author Jean Palate
 */
public class CompositeDiffuseSquareRootFilteringResults implements IDiffuseSquareRootFilteringResults {
    
    
    private final IDiffuseSquareRootFilteringResults[] subresults;
    public CompositeDiffuseSquareRootFilteringResults(final IDiffuseSquareRootFilteringResults... subresults){
        this.subresults=subresults;
    }

    @Override
    public void close(int pos) {
        for (int i=0; i<subresults.length; ++i){
            subresults[i].close(pos);
        }
    }

    @Override
    public void clear() {
        for (int i=0; i<subresults.length; ++i){
            subresults[i].clear();
        }
    }

    @Override
    public void save(int t, DiffusePredictionError pe) {
        for (int i=0; i<subresults.length; ++i){
            subresults[i].save(t, pe);
        }
    }

    @Override
    public void save(int t, AugmentedState state) {
        for (int i=0; i<subresults.length; ++i){
            subresults[i].save(t, state);
        }
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
    public int getEndDiffusePosition() {
        for (int i = 0; i < subresults.length; ++i) {
            int epos = subresults[i].getEndDiffusePosition();
            if (epos >= 0) {
                return epos;
            }
        }
        return -1;

    }

    @Override
    public DataBlock a(int pos) {
        for (int i = 0; i < subresults.length; ++i) {
            DataBlock a = subresults[i].a(pos);
            if (a != null) {
                return a;
            }
        }
        return null;
    }

    @Override
    public SubMatrix B(int pos) {
        for (int i = 0; i < subresults.length; ++i) {
            SubMatrix P = subresults[i].B(pos);
            if (P != null) {
                return P;
            }
        }
        return null;
    }

    @Override
    public DataBlock M(int pos) {
        for (int i = 0; i < subresults.length; ++i) {
            DataBlock c = subresults[i].M(pos);
            if (c != null) {
                return c;
            }
        }
        return null;
    }

    @Override
    public DataBlock Mi(int pos) {
        for (int i = 0; i < subresults.length; ++i) {
            DataBlock c = subresults[i].Mi(pos);
            if (c != null) {
                return c;
            }
        }
        return null;
    }

    @Override
    public SubMatrix P(int pos) {
        for (int i = 0; i < subresults.length; ++i) {
            SubMatrix P = subresults[i].P(pos);
            if (P != null) {
                return P;
            }
        }
        return null;
    }

    @Override
    public double error(int pos) {
        for (int i = 0; i < subresults.length; ++i) {
            double d = subresults[i].error(pos);
            if (Double.isFinite(d)) {
                return d;
            }
        }
        return Double.NaN;
    }

    @Override
    public double errorVariance(int pos) {
        for (int i = 0; i < subresults.length; ++i) {
            double d = subresults[i].errorVariance(pos);
            if (Double.isFinite(d)) {
                return d;
            }
        }
        return Double.NaN;
    }

    @Override
    public double diffuseNorm2(int pos) {
        for (int i = 0; i < subresults.length; ++i) {
            double d = subresults[i].diffuseNorm2(pos);
            if (Double.isFinite(d)) {
                return d;
            }
        }
        return Double.NaN;
    }
}
