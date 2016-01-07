/*
 * Copyright 2015 National Bank of Belgium
 *  
 * Licensed under the EUPL, Version 1.1 or â€“ as soon they will be approved 
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
package ec.tstoolkit2.ssf.akf;

import ec.tstoolkit.arima.ArimaModelBuilder;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.eco.ILikelihood;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.sarima.SarimaModel;
import ec.tstoolkit.sarima.SarimaModelBuilder;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.dk.sqrt.DiffuseSquareRootInitializer;
import ec.tstoolkit2.ssf.dk.DurbinKoopmanInitializer;
import ec.tstoolkit2.ssf.implementations.TimeInvariantSsf;
import ec.tstoolkit2.ssf.implementations.arima.SsfArima;
import ec.tstoolkit2.ssf.univariate.FilteringErrors;
import ec.tstoolkit2.ssf.univariate.SsfData;
import static org.junit.Assert.assertTrue;
import org.junit.Ignore;
import org.junit.Test;

/**
 *
 * @author Jean Palate
 */
public class AugmentedKalmanFilterTest {

    private static final int N = 100000, M = 25000;
    static final SarimaModel model;
    static final double[] data;
    private static final double Err, Det, lDet, Ll, DLl;
    private static final Matrix P;
    private static final DataBlock A;

    static {
        SarimaModelBuilder builder = new SarimaModelBuilder();
        model = builder.createArimaModel(12, 3, 1, 1, 0, 1, 1);
        model.setBTheta(1, -.99);
        //model = builder.createAirlineModel(12, -.6, -.8);
        ArimaModelBuilder gbuilder = new ArimaModelBuilder();
        data = gbuilder.generate(model, 500);
        data[5] = Double.NaN;
        data[12] = Double.NaN;
        data[21] = Double.NaN;
        ec.tstoolkit.ssf.arima.SsfArima ssf = new ec.tstoolkit.ssf.arima.SsfArima(model);
        ec.tstoolkit.ssf.Filter filter = new ec.tstoolkit.ssf.Filter();
        ec.tstoolkit.ssf.DiffusePredictionErrorDecomposition pe = new ec.tstoolkit.ssf.DiffusePredictionErrorDecomposition(false);
        filter.setSsf(ssf);
        filter.process(new ec.tstoolkit.ssf.SsfData(data, null), pe);
        Err = pe.getSsqErr();
        Det = pe.getLogDeterminant();
        lDet = pe.getDiffuseLogDeterminant();
        ec.tstoolkit.eco.DiffuseLikelihood dll = new ec.tstoolkit.eco.DiffuseLikelihood();
        ec.tstoolkit.ssf.LikelihoodEvaluation.evaluate(pe, dll);
        Ll = dll.getLogLikelihood();
        DLl = dll.getUncorrectedLogLikelihood();
        ec.tstoolkit.ssf.SsfData ssfdata = new ec.tstoolkit.ssf.SsfData(data, null);
        ec.tstoolkit.ssf.State state = new ec.tstoolkit.ssf.State(ssf.getStateDim(), true);
        ec.tstoolkit.ssf.DurbinKoopmanInitializer dk = new ec.tstoolkit.ssf.DurbinKoopmanInitializer();
        dk.initialize(ssf, ssfdata, state, pe);
        P = state.P;
        A = state.A;
    }

    public AugmentedKalmanFilterTest() {
    }

    @Test
    public void testAkf() {
        AugmentedFilter akf = new AugmentedFilter(true);
        SsfArima ssf = SsfArima.create(model);
        AugmentedPredictionErrorDecomposition pe = new AugmentedPredictionErrorDecomposition(false);
        pe.prepare(ssf, data.length);
        akf.process(ssf, new SsfData(data), pe);
        AugmentedState state = akf.getState();
        assertTrue(state.P().minus(P).nrm2() < 1e-9);
        assertTrue(state.a().distance(A) < 1e-9);
    }

    @Test
    public void testDK() {
        DurbinKoopmanInitializer dk = new DurbinKoopmanInitializer(null);
        SsfArima ssf = SsfArima.create(model);
        FilteringErrors pe = new FilteringErrors(false);
        SsfData ssfdata = new SsfData(data);
        State state = new State(ssf.getStateDim());
        dk.initialize(state, ssf, ssfdata);
        assertTrue(state.P().minus(P).nrm2() < 1e-9);
        assertTrue(state.a().distance(A) < 1e-9);
    }

    @Test
    public void testSquareRoot() {
        DiffuseSquareRootInitializer dk = new DiffuseSquareRootInitializer(null);
        SsfArima ssf = SsfArima.create(model);
        FilteringErrors pe = new FilteringErrors(false);
        SsfData ssfdata = new SsfData(data);
        State state = new State(ssf.getStateDim());
        dk.initialize(state, ssf, ssfdata);
        assertTrue(state.P().minus(P).nrm2() < 1e-9);
        assertTrue(state.a().distance(A) < 1e-9);
    }

    @Ignore
    @Test
    public void testStressInitialisation() {
        long t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            AugmentedFilter akf = new AugmentedFilter(true);
            SsfArima ssf = SsfArima.create(model);
            AugmentedPredictionErrorDecomposition pe = new AugmentedPredictionErrorDecomposition(false);
            pe.prepare(ssf, data.length);
            akf.process(ssf, new SsfData(data), pe);
        }
        long t1 = System.currentTimeMillis();
        System.out.println("akf");
        System.out.println(t1 - t0);
        System.out.println("dkold");
        t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            ec.tstoolkit.ssf.DurbinKoopmanInitializer dk = new ec.tstoolkit.ssf.DurbinKoopmanInitializer();
            ec.tstoolkit.ssf.arima.SsfArima ssf
                    = new ec.tstoolkit.ssf.arima.SsfArima(model);
            ec.tstoolkit.ssf.DiffusePredictionErrorDecomposition pe = new ec.tstoolkit.ssf.DiffusePredictionErrorDecomposition(false);
            ec.tstoolkit.ssf.SsfData ssfdata = new ec.tstoolkit.ssf.SsfData(data, null);
            ec.tstoolkit.ssf.State state = new ec.tstoolkit.ssf.State(ssf.getStateDim(), true);
            dk.initialize(ssf, ssfdata, state, pe);
        }
        t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
        System.out.println("dknew");
        t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            DurbinKoopmanInitializer dk = new DurbinKoopmanInitializer(null);
            SsfArima ssf = SsfArima.create(model);
            FilteringErrors pe = new FilteringErrors(false);
            SsfData ssfdata = new SsfData(data);
            State state = new State(ssf.getStateDim());
            dk.initialize(state, ssf, ssfdata);
        }
        t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
        System.out.println("square root");
        t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            DiffuseSquareRootInitializer dk = new DiffuseSquareRootInitializer(null);
            SsfArima ssf = SsfArima.create(model);
            FilteringErrors pe = new FilteringErrors(false);
            SsfData ssfdata = new SsfData(data);
            State state = new State(ssf.getStateDim());
            dk.initialize(state, ssf, ssfdata);
        }
        t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
        System.out.println("square root - Matrix");
        t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            DiffuseSquareRootInitializer dk = new DiffuseSquareRootInitializer(null);
            SsfArima ssf = SsfArima.create(model);
            FilteringErrors pe = new FilteringErrors(false);
            SsfData ssfdata = new SsfData(data);
            State state = new State(ssf.getStateDim());
            dk.initialize(state, TimeInvariantSsf.of(ssf), ssfdata);
        }
        t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
    }
}
