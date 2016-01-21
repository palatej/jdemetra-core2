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
import ec.tstoolkit.data.DataBlockStorage;
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
    private Matrix[] S;
    private final ISsf ssf;
    private final ISsfData data;
    private final ISsfDynamics dynamics;
    private final ISsfMeasurement measurement;
    private final Smoothing smoothing;

    public DiffuseSimulationSmoother(ISsf ssf, ISsfData data) {
        this.ssf = ssf;
        dynamics = ssf.getDynamics();
        measurement = ssf.getMeasurement();
        this.data = data;
        initSsf();
        smoothing=new Smoothing();
    }
    
    public Simulation newSimulation(){
        return new Simulation();
    }

    private Matrix lq(int pos) {
        if (LQ.length == 1) {
            return LQ[0];
        } else if (pos < LQ.length) {
            return LQ[pos];
        } else {
            Matrix Q = Matrix.square(dynamics.getInnovationsDim());
            dynamics.Q(pos, Q.subMatrix());
            SymmetricMatrix.lcholesky(Q, EPS);
            return Q;
        }
    }

    private Matrix S(int pos) {
        if (S.length == 1) {
            return S[0];
        } else if (pos < S.length) {
            return S[pos];
        } else if (!dynamics.hasS()) {
            return null;
        } else {
            Matrix s = new Matrix(dynamics.getStateDim(), dynamics.getInnovationsDim());
            dynamics.S(pos, s.subMatrix());
            return s;
        }
    }

    private double lh(int pos) {
        return Math.sqrt(ssf.getMeasurement().errorVariance(pos));
    }

    private void initSsf() {
        int dim = dynamics.getStateDim(), resdim = dynamics.getInnovationsDim();
        LA = Matrix.square(dim);
        dynamics.Pf0(LA.subMatrix(), StateInfo.Forecast);
        SymmetricMatrix.lcholesky(LA, EPS);

        if (dynamics.isTimeInvariant()) {
            LQ = new Matrix[1];
            S = new Matrix[1];
        } else {
            LQ = new Matrix[data.getLength()];
            S = new Matrix[data.getLength()];
        }
        for (int i = 0; i < LQ.length; ++i) {
            Matrix Q = Matrix.square(resdim);
            dynamics.Q(i, Q.subMatrix());
            SymmetricMatrix.lcholesky(Q, EPS);
            LQ[i] = Q;
            Matrix s = new Matrix(dim, resdim);
            dynamics.S(i, s.subMatrix());
            S[i]=s;
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
    
    public class Smoothing {
        final IBaseDiffuseFilteringResults frslts;
        final DataBlockStorage smoothedStates, smoothedInnovations;
        
        Smoothing(){
            frslts=DkToolkit.sqrtFilter(ssf, data, false);
            int dim=dynamics.getStateDim(), resdim=dynamics.getInnovationsDim();
            smoothedStates=new DataBlockStorage(dim, data.getLength());
            smoothedInnovations=new DataBlockStorage(resdim, data.getLength());
        }
        
    }

    public class Simulation {

        final DataBlockStorage states;
        final DataBlockStorage transitionInnovations;
        final double[] measurementErrors;
        private final double[] simulatedData;
        DataBlockStorage smoothedStates;
        DataBlockStorage smoothedInnovations;
        DataBlockStorage R, Ri;

        Simulation() {
            int dim = dynamics.getStateDim();
            int nres = dynamics.getInnovationsDim();
            boolean err = measurement.hasErrors();
            int n = data.getLength();
            states = new DataBlockStorage(dim, n);
            transitionInnovations = new DataBlockStorage(nres, n);
            if (err) {
                measurementErrors = new double[n];
                generateMeasurementRandoms(new DataBlock(measurementErrors));
            } else {
                measurementErrors = null;
            }
            simulatedData = new double[n];
            generateData();
        }

        private void generateData() {
            DataBlock a0f = new DataBlock(dynamics.getStateDim());
            generateInitialState(a0f);
            DataBlock cur = states.block(0);
            dynamics.a0(cur, StateInfo.Forecast);
            cur.add(a0f);
            simulatedData[0] = measurement.ZX(0, cur);
            if (measurementErrors != null) {
                simulatedData[0] += measurementErrors[0];
            }
            // a0 = a(1|0) -> y[1) = Z*a[1|0) + e(1)
            // a(2|1) = T a(1|0) + S * q(1)...
            for (int i = 1; i < simulatedData.length; ++i) {
                DataBlock q = transitionInnovations.block(i - 1);
                generateTransitionRandoms(i - 1, q);
                DataBlock prev = cur;
                cur = states.block(i);
                cur.copy(prev);
                cur.addProduct(S(i).rows(), q);
                simulatedData[i] = measurement.ZX(i, cur);
                if (measurementErrors != null) {
                    simulatedData[i] += measurementErrors[i];
                }
            }
        }

        /**
         * @return the simulatedData
         */
        public double[] getSimulatedData() {
            return simulatedData;
        }
    }
}
