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
import ec.tstoolkit.dstats.Normal;
import ec.tstoolkit.maths.matrices.LowerTriangularMatrix;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit.random.IRandomNumberGenerator;
import ec.tstoolkit.random.XorshiftRNG;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfData;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;

/**
 *
 * @author Jean Palate
 */
public class DiffuseSimulationSmoother {

    private static final Normal N = new Normal();
    private static final IRandomNumberGenerator RNG = XorshiftRNG.fromSystemNanoTime();

    private static void fillRandoms(DataBlock u) {
        synchronized (N) {
            for (int i = 0; i < u.getLength(); ++i) {
                u.set(i, N.random(RNG));
            }
        }
    }

    private static double random() {
        synchronized (N) {
            return N.random(RNG);
        }
    }

    private static final double EPS = 1e-8;

    private Matrix LA;
    private Matrix[] LQ;
    private final ISsf ssf;
    private final ISsfData data;
    

    public DiffuseSimulationSmoother(ISsf ssf, ISsfData data) {
        this.ssf = ssf;
        this.data = data;
        initSsf();
    }

    private Matrix lq(int pos) {
        if (LQ.length == 1) {
            return LQ[0];
        } else if (pos < LQ.length) {
            return LQ[pos];
        } else {
            ISsfDynamics dynamics = ssf.getDynamics();
            Matrix Q = Matrix.square(dynamics.getInnovationsDim());
            dynamics.Q(pos, Q.subMatrix());
            SymmetricMatrix.lcholesky(Q, EPS);
            return Q;
        }
    }

    private double lh(int pos) {
        return Math.sqrt(ssf.getMeasurement().errorVariance(pos));
    }

    private void initSsf() {
        ISsfDynamics dynamics = ssf.getDynamics();
        int dim = dynamics.getStateDim(), resdim = dynamics.getInnovationsDim();
        LA = Matrix.square(dim);
        dynamics.Pf0(LA.subMatrix(), StateInfo.Forecast);
        SymmetricMatrix.lcholesky(LA, EPS);

        if (dynamics.isTimeInvariant()) {
            LQ = new Matrix[1];
        } else {
            LQ = new Matrix[data.getCount()];
        }
        for (int i = 0; i < LQ.length; ++i) {
            Matrix Q = Matrix.square(resdim);
            dynamics.Q(i, Q.subMatrix());
            SymmetricMatrix.lcholesky(Q, EPS);
            LQ[i] = Q;
        }
    }

    private void generateTransitionRandoms(int pos, DataBlock u) {
        fillRandoms(u);
        LowerTriangularMatrix.rmul(lq(pos), u);
    }

    private void generateMeasurementRandoms(DataBlock e) {
        fillRandoms(e);
        e.mul(lh(0));
    }

    private double generateMeasurementRandom(int pos) {
        double e = random();
        return e * lh(pos);
    }

    private void generateInitialState(DataBlock a) {
        fillRandoms(a);
        LowerTriangularMatrix.rmul(LA, a);
    }

}
