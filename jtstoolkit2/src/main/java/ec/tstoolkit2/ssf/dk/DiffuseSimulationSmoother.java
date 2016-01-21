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
import ec.tstoolkit.data.DescriptiveStatistics;
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
    private Matrix[] Q, S, SQ;
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
        smoothing = new Smoothing();
    }

    public Smoothing getReferenceSmoothing() {
        return smoothing;
    }

    public Simulation newSimulation() {
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

    private Matrix SQ(int pos) {
        if (SQ.length == 1) {
            return SQ[0];
        } else if (pos < SQ.length) {
            return SQ[pos];
        } else if (!dynamics.hasS()) {
            return Q(pos);
        } else {
            Matrix s = new Matrix(dynamics.getStateDim(), dynamics.getInnovationsDim());
            dynamics.S(pos, s.subMatrix());
            return s.times(Q(pos));
        }
    }

    private Matrix Q(int pos) {
        if (Q.length == 1) {
            return Q[0];
        } else if (pos < Q.length) {
            return Q[pos];
        } else {
            Matrix Q = Matrix.square(dynamics.getInnovationsDim());
            dynamics.Q(pos, Q.subMatrix());
            return Q;
        }
    }

    private double lh(int pos) {
        return Math.sqrt(ssf.getMeasurement().errorVariance(pos));
    }

    private double h(int pos) {
        return ssf.getMeasurement().errorVariance(pos);
    }

    private void initSsf() {
        int dim = dynamics.getStateDim(), resdim = dynamics.getInnovationsDim();
        LA = Matrix.square(dim);
        dynamics.Pf0(LA.subMatrix(), StateInfo.Forecast);
        SymmetricMatrix.lcholesky(LA, EPS);

        if (dynamics.isTimeInvariant()) {
            Q = new Matrix[1];
            LQ = new Matrix[1];
            S = new Matrix[1];
            SQ = new Matrix[1];
        } else {
            Q = new Matrix[data.getLength()];
            LQ = new Matrix[data.getLength()];
            S = new Matrix[data.getLength()];
            SQ = new Matrix[data.getLength()];
        }
        for (int i = 0; i < Q.length; ++i) {
            Matrix q = Matrix.square(resdim);
            dynamics.Q(i, q.subMatrix());
            Q[i] = q.clone();
            Matrix s = new Matrix(dim, resdim);
            dynamics.S(i, s.subMatrix());
            S[i] = s;
            SQ[i] = s.times(q);
            SymmetricMatrix.lcholesky(q, EPS);
            LQ[i] = q;
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

    abstract class BaseSimulation {

        protected final IBaseDiffuseFilteringResults frslts;
        protected final DataBlockStorage smoothedInnovations;
        protected final DataBlock esm;
        protected DataBlockStorage smoothedStates;
        protected DataBlock a0;
        protected final int dim, resdim, n, nd;
        protected DataBlock R, Ri;

        protected abstract double getError(int pos);

        protected BaseSimulation(IBaseDiffuseFilteringResults frslts) {
            this.frslts = frslts;
            dim = dynamics.getStateDim();
            resdim = dynamics.getInnovationsDim();
            n = data.getLength();
            nd = frslts.getEndDiffusePosition();
            smoothedInnovations = new DataBlockStorage(resdim, n);
            if (measurement.hasErrors()) {
                esm = new DataBlock(n);
            } else {
                esm = null;
            }
        }

        protected void smooth() {
            R = new DataBlock(dim);
            Ri = new DataBlock(dim);
            // we reproduce here the usual iterations of the smoother
            doNormalSmoothing();
            doDiffuseSmoohing();
            computeInitialState();
        }

        private void doNormalSmoothing() {
            double e, v;
            DataBlock M = new DataBlock(dim);
            DataBlock U = new DataBlock(resdim);
            boolean missing;
            int pos = n;
            while (--pos >= nd) {
                // Get info
                e = getError(pos);
                v = frslts.errorVariance(pos);
                M.copy(frslts.M(pos));
                missing = !DescriptiveStatistics.isFinite(e);
                // Iterate R
                dynamics.XT(pos, R);
                if (!missing && e != 0) {
                    // RT
                    double c = (e - R.dot(M)) / v;
                    measurement.XpZd(pos, R, c);
                }
                // Computes esm, U
                if (esm != null) {
                    if (!missing) {
                        esm.set(pos, h(pos));
                    } else {
                        esm.set(pos, Double.NaN);
                    }
                }
                U.product(R, SQ(pos).columns());
                smoothedInnovations.save(pos, U);
            }
        }

        private void doDiffuseSmoohing() {
            double e, f, fi, c;
            DataBlock C = new DataBlock(dim), Ci = new DataBlock(dim);
            DataBlock U = new DataBlock(resdim);
            boolean missing;
            int pos = nd;
            while (--pos >= 0) {
                // Get info
                e = frslts.error(pos);
                f = frslts.errorVariance(pos);
                fi = frslts.diffuseNorm2(pos);
                C.copy(frslts.M(pos));
                if (fi != 0) {
                    Ci.copy(frslts.Mi(pos));
                    Ci.mul(1 / fi);
                    C.addAY(-f, Ci);
                    C.mul(1 / fi);
                } else {
                    C.mul(1 / f);
                    Ci.set(0);
                }
                missing = !DescriptiveStatistics.isFinite(e);
                // Iterate R
                if (fi == 0) {
                    dynamics.XT(pos, R);
                    dynamics.XT(pos, Ri);
                    if (!missing && f != 0) {
                        c = e / f - R.dot(C);
                        measurement.XpZd(pos, R, c);
                    }
                } else {
                    dynamics.XT(pos, R);
                    dynamics.XT(pos, Ri);
                    if (!missing && f != 0) {
                        c = -Ri.dot(Ci);
                        double ci = e / fi + c - R.dot(C);
                        measurement.XpZd(pos, Ri, ci);
                        double cf = -R.dot(Ci);
                        measurement.XpZd(pos, R, cf);
                    }
                }
                // Computes esm, U
                if (esm != null) {
                    if (!missing) {
                        esm.set(pos, h(pos));
                    } else {
                        esm.set(pos, Double.NaN);
                    }
                }
                U.product(R, SQ(pos).columns());
                smoothedInnovations.save(pos, U);
            }
        }

        private void computeInitialState() {
            // initial state
            a0 = new DataBlock(dim);
            Matrix Pf0 = Matrix.square(dim);
            dynamics.a0(a0, StateInfo.Forecast);
            dynamics.Pf0(Pf0.subMatrix(), StateInfo.Forecast);
            // stationary initialization
            a0.addProduct(R, Pf0.columns());

            // non stationary initialisation
            Matrix Pi0 = Matrix.square(dim);
            dynamics.Pi0(Pi0.subMatrix());
            a0.addProduct(Ri, Pi0.columns());
        }

        public DataBlock getSmoothedInnovations(int pos) {
            return smoothedInnovations.block(pos);
        }

        public DataBlock getSmoothedState(int pos) {
            if (smoothedStates == null) {
                generatesmoothedStates();
            }
            return smoothedStates.block(pos);
        }

        public DataBlockStorage getSmoothedInnovations() {
            return smoothedInnovations;
        }

        public DataBlockStorage getSmoothedStates() {
            if (smoothedStates == null) {
                generatesmoothedStates();
            }
            return smoothedStates;
        }

        private void generatesmoothedStates() {
            smoothedStates = new DataBlockStorage(dim, n);
            smoothedStates.save(0, a0);
            int cur = 1;
            DataBlock a = a0.deepClone();
            do {
                // next: a(t+1) = T a(t) + S*r(t)
                dynamics.TX(cur, a);
                a.addProduct(smoothedInnovations.block(cur), S(cur).rows());
                smoothedStates.save(cur++, a);
            } while (cur < n);
        }
    }

    public class Smoothing extends BaseSimulation {

        Smoothing() {
            super(DkToolkit.sqrtFilter(ssf, data, false));
            smooth();
        }

        @Override
        protected double getError(int pos) {
            return frslts.error(pos);
        }

    }

    public class Simulation extends BaseSimulation{
        
        public Simulation(){
            super(smoothing.frslts);
            boolean err = measurement.hasErrors();
            states = new DataBlockStorage(dim, n);
            transitionInnovations = new DataBlockStorage(resdim, n);
            if (err) {
                measurementErrors = new double[n];
                generateMeasurementRandoms(new DataBlock(measurementErrors));
            } else {
                measurementErrors = null;
            }
            simulatedData = new double[n];
            generateData();
            smooth();
        }

        final DataBlockStorage states;
        final DataBlockStorage transitionInnovations;
        final double[] measurementErrors;
        private final double[] simulatedData;

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

        @Override
        protected double getError(int pos) {
            return simulatedData[pos]-measurement.ZX(pos, states.block(pos));
        }
    }
}
