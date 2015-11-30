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
package ec.tstoolkit2.ssf.univariate;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.IReadDataBlock;
import ec.tstoolkit.data.ReadDataBlock;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit2.ssf.DataBlockResults;
import ec.tstoolkit2.ssf.DataResults;
import ec.tstoolkit2.ssf.IStateResults;
import ec.tstoolkit2.ssf.MatrixResults;
import ec.tstoolkit2.ssf.ResultsRange;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;

/**
 *
 * @author Jean Palate
 */
public class DefaultFilteringResults implements IFilteringResults, IStateResults {

    private final DataBlockResults A; // state vector
    private final MatrixResults P;  // P
    private final DataBlockResults C; // C = P*Z'
    private final DataResults e, f; // errors, variances of the errors
    private final ResultsRange range=new ResultsRange();

    protected DefaultFilteringResults(boolean cov) {
        A = new DataBlockResults();
        C = new DataBlockResults();
        P = cov ? new MatrixResults() : null;
        e = new DataResults();
        f = new DataResults();
    }

    public boolean isInitialized(){
        return A.isInitialized();
    }
    
    public ResultsRange getRange(){
        return range;
    }

    public static DefaultFilteringResults full() {
        return new DefaultFilteringResults(true);
    }

    public static DefaultFilteringResults light() {
        return new DefaultFilteringResults(false);
    }

    public void prepare(ISsf ssf, final int start, final int end) {
        int dim = ssf.getStateDim();

        A.prepare(dim, start, end);
        C.prepare(dim, start, end);
        e.prepare(start, end);
        f.prepare(start, end);
        if (P != null) {
            P.prepare(dim, start, end);
        }
    }

    @Override
    public void save(int t, PredictionError pe) {
        e.save(t, pe.get());
        f.save(t, pe.getVariance());
        C.save(t, pe.C());
    }

    @Override
    public void save(int t, State state) {
        if (state.getInfo() != StateInfo.Forecast) {
            return;
        }
        A.save(t, state.a());
        if (P != null) {
            P.save(t, state.P());
        }
        range.add(t);
    }

    public double error(int pos) {
        return e.get(pos);
    }

    public double errorVariance(int pos) {
        return f.get(pos);
    }

    public IReadDataBlock errors() {
        return e;
    }

    public IReadDataBlock errorVariances() {
        return f;
    }

    public DataBlock a(int pos) {
        return A.datablock(pos);
    }

    public DataBlock c(int pos) {
        return C.datablock(pos);
    }

    public SubMatrix P(int pos) {
        return P.subMatrix(pos);
    }

    @Override
    public void clear(){
        e.clear();
        f.clear();
        A.clear();
        C.clear();
        if (P != null)
            P.clear();
        range.clear();
    }
}
