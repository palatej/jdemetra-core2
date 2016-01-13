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
import ec.tstoolkit2.ssf.ISsfDynamics;
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

    @Override
    public void saveSmoothedMeasurementDisturbance(int t, double err, double v) {
        if (e == null) {
            return;
        }
        e.save(t, err);
        if (evar != null) {
            evar.save(t, v);
        }
    }
    
    public DataBlock uComponent(int item){
        return U.item(item);
    }

    public DataBlock uComponentVar(int item){
        return UVar.item(item, item);
    }

    public DataBlock e(){
        return e.all();
    }
    
   public DataBlock evar(){
        return evar.all();
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
        ISsfDynamics dynamics = ssf.getDynamics();
        int dim = dynamics.getStateDim(), edim=dynamics.getInnovationsDim();
        if (e != null && ssf.getMeasurement().hasErrors()) {
            e.prepare(start, end);
            evar.prepare(start, end);
        }
        U.prepare(edim, start, end);

        if (UVar != null) {
            UVar.prepare(edim, start, end);
        }
    }

    public void rescale(double factor) {
        double se=Math.sqrt(factor);
        U.rescale(se);
//        if (UVar != null) {
//            UVar.rescale(factor);
//        }
//        if (evar != null) {
//            evar.rescale(factor);
//        }
    }

}
