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

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.IReadDataBlock;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit2.ssf.DataBlockResults;
import ec.tstoolkit2.ssf.DataResults;
import ec.tstoolkit2.ssf.MatrixResults;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;

/**
 *
 * @author Jean Palate
 */
public class DefaultSmoothingResults implements ISmoothingResults {

    private final DataBlockResults A;
    private final MatrixResults P;
    private final DataResults e, f;

    private DefaultSmoothingResults(final boolean cov, final boolean err) {
        A = new DataBlockResults();
        P = cov ? new MatrixResults() : null;
        if (err) {
            e = new DataResults();
            f = new DataResults();
        } else {
            e = null;
            f = null;
        }
    }

    public static DefaultSmoothingResults full() {
        return new DefaultSmoothingResults(true, true);
    }

    public static DefaultSmoothingResults light() {
        return new DefaultSmoothingResults(false, false);
    }

    @Override
    public void save(int t, State state) {
        if (state.getInfo() != StateInfo.Smoothed) {
            return;
        }
        A.save(t, state.a());

        if (P != null) {
            P.save(t, state.P());
        }
    }

    public DataBlock getComponent(int pos) {
        return A.item(pos);
    }

    public DataBlock getComponentVariance(int pos) {
        return P.item(pos, pos);
    }

    public void saveSmoothedError(int t, double err, double v) {
        if (e == null) {
            return;
        }
        e.save(t, err);
        f.save(t, v);
    }

    public IReadDataBlock errors() {
        return e == null ? null : e;
    }

    public IReadDataBlock errorVariances() {
        return f == null ? null : f;
    }

    @Override
    public DataBlock a(int pos) {
        return A.datablock(pos);
    }

    @Override
    public SubMatrix P(int pos) {
        return P == null ? null : P.subMatrix(pos);
    }

    public int getStart() {
        return A.getStartSaving();
    }

    public void prepare(ISsf ssf, int start, int end) {
        int dim = ssf.getStateDim();
        if (e != null && ssf.getMeasurement().hasErrors()) {
            e.prepare(start, end);
            f.prepare(start, end);
        } 
        A.prepare(dim, start, end);
        
        if (P != null) {
            P.prepare(dim, start, end);
        }
    }

    public void rescaleVariances(double factor){
        if (P != null){
            P.rescale(factor);
        }
        if (f != null)
            f.rescale(factor);
    }

}
