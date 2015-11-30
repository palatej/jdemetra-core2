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
package ec.tstoolkit2.ssf.dk.sqrt;

import ec.tstoolkit2.ssf.dk.*;
import ec.tstoolkit2.ssf.univariate.*;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit2.ssf.DataBlockResults;
import ec.tstoolkit2.ssf.DataBlocksResults;
import ec.tstoolkit2.ssf.DataResults;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.akf.AugmentedState;

/**
 *
 * @author Jean Palate
 */
public class DefaultDiffuseSquareRootFilteringResults extends DefaultFilteringResults implements IDiffuseSquareRootFilteringResults {

    private final DataBlockResults Ci;
    private final DataBlocksResults B;
    private DataResults fi;
    private int enddiffuse;

    private DefaultDiffuseSquareRootFilteringResults(boolean var) {
        super(var);
        Ci = new DataBlockResults();
        fi=new DataResults();
        B = var ? new DataBlocksResults() : null;
    }

    public static DefaultDiffuseSquareRootFilteringResults full() {
        return new DefaultDiffuseSquareRootFilteringResults(true);
    }

    public static DefaultDiffuseSquareRootFilteringResults light() {
        return new DefaultDiffuseSquareRootFilteringResults(false);
    }

    @Override
    public void prepare(ISsf ssf, final int start, final int end) {
        super.prepare(ssf, start, end);
        int dim = ssf.getStateDim(), n = ssf.getDynamics().getNonStationaryDim();
        fi.prepare(start, n);
        Ci.prepare(dim, start, n);
        if (B != null) {
            B.prepare(dim, n, n);
        }
    }

    @Override
    public void save(int t, DiffusePredictionError pe) {
        super.save(t, pe);
        fi.save(t, pe.getDiffuseNorm2());
        Ci.save(t, pe.Ci());
    }

    @Override
    public void close(int pos) {
        enddiffuse = pos;
    }

    @Override
    public void save(int t, AugmentedState state) {
        if (state.getInfo() != StateInfo.Forecast) {
            return;
        }
        super.save(t, state);
        if (B != null) {
            B.save(t, state.B());
        }

    }


    public double diffuseNorm(int pos) {
        return fi.get(pos);
    }

    public DataBlock ci(int pos) {
        return Ci.datablock(pos);
    }
 
    public SubMatrix B(int pos) {
        return B.subMatrix(pos);
    }

    @Override
    public void clear() {
        super.clear();
        enddiffuse = 0;
    }

    public int getEndDiffusePosition() {
        return enddiffuse;
    }
}
