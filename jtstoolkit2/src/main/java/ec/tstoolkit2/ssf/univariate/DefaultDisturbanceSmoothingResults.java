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
public class DefaultDisturbanceSmoothingResults implements IDisturbanceSmoothingResults {

    private final DataBlockResults U;
    private final MatrixResults UVar;
    private final DataResults e, evar;

    private DefaultDisturbanceSmoothingResults(final boolean cov, final boolean err) {
        U = new DataBlockResults();
        UVar = cov ? new MatrixResults() : null;
        if (err) {
            e = new DataResults();
            if (cov) {
                evar = new DataResults();
            } else {
                evar = null;
            }
        } else {
            e = null;
            evar = null;
        }
    }

    public static DefaultDisturbanceSmoothingResults full(boolean err) {
        return new DefaultDisturbanceSmoothingResults(true, err);
    }

    public static DefaultDisturbanceSmoothingResults light(boolean err) {
        return new DefaultDisturbanceSmoothingResults(false, err);
    }

    @Override
    public void saveSmoothedTransitionDisturbances(int t, DataBlock u, SubMatrix uvar) {
        U.save(t, u);

        if (UVar != null && uvar != null) {
            UVar.save(t, uvar);
        }
    }

    public void saveSmoothedMeasurementDisturbance(int t, double err, double v) {
        if (e == null) {
            return;
        }
        e.save(t, err);
        if (evar != null) {
            evar.save(t, v);
        }
    }

    @Override
    public DataBlock u(int pos) {
        return U.datablock(pos);
    }

    @Override
    public SubMatrix uVar(int pos) {
        return UVar == null ? null : UVar.subMatrix(pos);
    }

    @Override
    public double e(int pos) {
        return e == null ? 0 : e.get(pos);
    }

    @Override
    public double eVar(int pos) {
        return evar == null ? 0 : evar.get(pos);
    }

    public int getStart() {
        return U.getStartSaving();
    }

    public void prepare(ISsf ssf, int start, int end) {
        int dim = ssf.getStateDim();
        if (e != null && ssf.getMeasurement().hasErrors()) {
            e.prepare(start, end);
            evar.prepare(start, end);
        }
        U.prepare(dim, start, end);

        if (UVar != null) {
            UVar.prepare(dim, start, end);
        }
    }

    public void rescaleVariances(double factor) {
        if (UVar != null) {
            UVar.rescale(factor);
        }
        if (evar != null) {
            evar.rescale(factor);
        }
    }

}
