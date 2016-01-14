/*
 * Copyright 2013 National Bank of Belgium
 *
 * Licensed under the EUPL, Version 1.1 or - as soon they will be approved 
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
import ec.tstoolkit.design.Development;
import ec.tstoolkit.maths.matrices.ElementaryTransformations;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.SsfException;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.akf.AugmentedState;
import ec.tstoolkit2.ssf.dk.DiffusePredictionError;
import ec.tstoolkit2.ssf.dk.IDiffuseFilteringResults;
import ec.tstoolkit2.ssf.univariate.IFilteringResults;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfData;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import ec.tstoolkit2.ssf.univariate.OrdinaryFilter;

/**
 * Mixed algorithm based on the diffuse initializer of Durbin-Koopman and on the
 * (square root) array filter of Kailath for the diffuse part. That solution
 * provides a much more stable estimate of the diffuse part.
 *
 * @author Jean Palate
 */
@Development(status = Development.Status.Preliminary)
public class DiffuseSquareRootInitializer implements OrdinaryFilter.Initializer {

    private final IDiffuseSquareRootFilteringResults results;
    private AugmentedState astate;
    private DiffusePredictionError pe;
    private ISsfMeasurement measurement;
    private ISsfDynamics dynamics;
    private ISsfData data;
    private int pos;
    private DataBlock Z;

    public DiffuseSquareRootInitializer() {
        this.results = null;
    }

    /**
     *
     * @param results
     */
    public DiffuseSquareRootInitializer(IDiffuseSquareRootFilteringResults results) {
        this.results = results;
    }

    /**
     *
     * @param ssf
     * @param data
     * @param state
     * @return
     */
    @Override
    public int initialize(final State state, final ISsf ssf, final ISsfData data) {
        measurement = ssf.getMeasurement();
        dynamics = ssf.getDynamics();
        this.data = data;
        pos = 0;
        int end = data.getCount();
        if (!initState()) {
            return -1;
        }
        pred();
        while (pos < end) {
            if (!astate.isDiffuse()) {
                break;
            }
            if (results != null) {
                results.save(pos, astate);
            }
            if (error()) {
                if (results != null) {
                    results.save(pos, pe);
                }
                update();
            } else {
                if (results != null) {
                    results.save(pos, pe);
                }
                astate.setInfo(StateInfo.Concurrent);
            }
            if (results != null) {
                results.save(pos, astate);
            }
            pred();
            ++pos;
        }

        if (results != null) {
            results.close(pos);
        }
        state.setInfo(StateInfo.Forecast);
        state.P().copy(this.astate.P());
        state.a().copy(this.astate.a());
        return pos;
    }

    private boolean initState() {
        int r = dynamics.getStateDim();
        astate = AugmentedState.of(dynamics, StateInfo.Forecast);
        if (astate == null) {
            return false;
        }
        pe = new DiffusePredictionError(r);
        Z = new DataBlock(astate.getDiffuseDim());
        dynamics.diffuseConstraints(constraints());
        return true;
    }

    /**
     * Computes a(t+1|t), P(t+1|t) from a(t|t), P(t|t) a(t+1|t) = T(t)a(t|t)
     * P(t+1|t) = T(t)P(t|t)T'(t) + V(t)
     */
    protected void pred() {
        if (astate.getInfo() != StateInfo.Forecast) {
            astate.setInfo(StateInfo.Forecast);
            SubMatrix P = astate.P().subMatrix();
            DataBlock a = astate.a();
            dynamics.TX(pos, a);
            dynamics.TVT(pos, P);
            dynamics.addV(pos, P);
            if (astate.isDiffuse()) {
                dynamics.TM(pos, astate.B());
            }
        }
    }

    private void update() {
        if (pe.isDiffuse()) {
            update1();
        } else {
            update0();
        }
        astate.setInfo(StateInfo.Concurrent);
    }

    private void update0() {
        double f = pe.getVariance(), e = pe.get();
        DataBlock C = pe.M();
        Matrix P = astate.P();
        SymmetricMatrix.addXaXt(P, -1 / f, C);

        // state
        // a0 = a0 + f1*Mi*v0.
        if (data.hasData()) {
            double c = e / f;
            astate.a().addAY(c, C);
        }
    }

    private void update1() {
        double fi = pe.getDiffuseNorm2(), f = pe.getVariance(), e = pe.get();
        DataBlock C = pe.M(), Ci = pe.Mi();
        // P = T P T' - 1/f*(TMf)(TMf)'+RQR'+f*(TMf/f-TMi/fi)(TMf/f-TMi/fi)'
        SymmetricMatrix.addXaXt(astate.P(), -1 / f, C);

        DataBlock tmp = C.deepClone();
        tmp.addAY(-f / fi, Ci);
        SymmetricMatrix.addXaXt(astate.P(), 1 / f, tmp);

        if (data.hasData()) {
            // a0 = a0 + f1*Mi*v0. Reuse Mf as temporary buffer
            astate.a().addAY(e / fi, Ci);
        }
    }

    protected boolean error() {
        if (astate.getInfo() != StateInfo.Forecast) {
            throw new SsfException(SsfException.STATUS);
        }
        // calc f and fi
        // fi = Z Pi Z' , f = Z P Z' + H
        preArray();
        double fi = zconstraints().ssq();
        if (fi < State.ZERO) {
            fi = 0;
        }
        pe.setDiffuseNorm2(fi);

        double f = measurement.ZVZ(pos, astate.P().subMatrix());
        if (measurement.hasErrors()) {
            f += measurement.errorVariance(pos);
        }
        pe.setVariance(f);
        if (data.hasData()) {
            double y = data.get(pos);
            if (Double.isNaN(y)) {
                pe.setMissing();
                return false;
            } else {
                pe.set(y - measurement.ZX(pos, astate.a()));
            }
        }

        DataBlock C = pe.M();
        measurement.ZM(pos, astate.P().subMatrix(), C);
        if (pe.isDiffuse()) {
            DataBlock z = zconstraints();
            SubMatrix B = constraints();
            ElementaryTransformations.fastRowGivens(z, B);
            pe.Mi().setAY(z.get(0), B.column(0));
            // move right
            astate.dropDiffuseConstraint();
        }
        return true;
    }

    // Array routines
    //     |R Z*X|
    // X = |     | 
    //     |0   X|
    // XX' = |RR'+ZXX'Z' ZXX'| = |AA'     AB'|
    //       |XX'Z'      XX' | = |BA' BB'+CC'|
    // A = Fi^1/2
    // B = Ci * Fi^-1/2
    // TC = X(t+1)
    private void preArray() {
        DataBlock zconstraints = zconstraints();
        zconstraints.set(0);
        SubMatrix A = constraints();
        measurement.ZM(pos, A, zconstraints);
        //dynamics.TM(pos, A);
    }

    private SubMatrix constraints() {
        return astate.B();
    }

    private DataBlock zconstraints() {
        return Z.range(0, astate.getDiffuseDim());
    }

//    private void checkDiffuse() {
//        SubMatrix C = constraints();
//        for (int c = ndiffuse_ - 1; c >= 0; --c) {
//            if (C.column(c).nrm2() < State.ZERO) {
//                --ndiffuse_;
//            } else {
//                break;
//            }
//        }
//    }
}
