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
package ec.tstoolkit2.ssf.dk;

import ec.tstoolkit2.ssf.univariate.*;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit2.ssf.DataBlockResults;
import ec.tstoolkit2.ssf.DataResults;
import ec.tstoolkit2.ssf.MatrixResults;
import ec.tstoolkit2.ssf.StateInfo;

/**
 *
 * @author Jean Palate
 */
public class DefaultDiffuseFilteringResults extends DefaultFilteringResults implements IDiffuseFilteringResults {

    private final DataBlockResults Ci;
    private final MatrixResults Pi;
    private final DataResults fi;
    private int enddiffuse;

    private DefaultDiffuseFilteringResults(boolean var) {
        super(var);
        Ci = new DataBlockResults();
        fi=new DataResults();
        Pi = var ? new MatrixResults() : null;
    }

    public static DefaultDiffuseFilteringResults full() {
        return new DefaultDiffuseFilteringResults(true);
    }

    public static DefaultDiffuseFilteringResults light() {
        return new DefaultDiffuseFilteringResults(false);
    }
    
    @Override
    public void prepare(ISsf ssf, final int start, final int end) {
        super.prepare(ssf, start, end);
        int dim = ssf.getStateDim(), n = ssf.getDynamics().getNonStationaryDim();
        fi.prepare(start, n);
        Ci.prepare(dim, start, n);
        if (Pi != null) {
            Pi.prepare(dim, start, n);
        }
    }

    @Override
    public void save(int t, DiffusePredictionError pe) {
        super.save(t, pe);
        fi.save(t, pe.getDiffuseNorm2());
        Ci.save(t, pe.Mi());
    }

    @Override
    public void close(int pos) {
        enddiffuse = pos;
    }

    @Override
    public void save(int t, DiffuseState state) {
        if (state.getInfo() != StateInfo.Forecast) {
            return;
        }
        super.save(t, state);
        if (Pi != null) {
            Pi.save(t, state.Pi());
        }
    }


    @Override
    public double diffuseNorm2(int pos) {
        return fi.get(pos);
    }

    @Override
    public DataBlock Mi(int pos) {
        return Ci.datablock(pos);
    }
 
    @Override
    public SubMatrix Pi(int pos) {
        return Pi.subMatrix(pos);
    }

    @Override
    public void clear() {
        super.clear();
        Ci.clear();
        fi.clear();
        Pi.clear();
        enddiffuse = 0;
    }

    @Override
    public int getEndDiffusePosition() {
        return enddiffuse;
    }
}
