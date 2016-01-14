/*
 * Copyright 2015 National Bank of Belgium
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
 /*
 */
package ec.tstoolkit2.ssf.implementations.arima;

import ec.tstoolkit.arima.ArimaModel;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit.sarima.SarimaModel;
import ec.tstoolkit.sarima.SarimaModelBuilder;
import ec.tstoolkit.timeseries.simplets.TsData;
import ec.tstoolkit.ucarima.ModelDecomposer;
import ec.tstoolkit.ucarima.SeasonalSelector;
import ec.tstoolkit.ucarima.TrendCycleSelector;
import ec.tstoolkit.ucarima.UcarimaModel;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.akf.AugmentedFilter;
import ec.tstoolkit2.ssf.akf.AugmentedPredictionErrorDecomposition;
import ec.tstoolkit2.ssf.dk.sqrt.DiffuseSquareRootInitializer;
import ec.tstoolkit2.ssf.dk.DurbinKoopmanInitializer;
import ec.tstoolkit2.ssf.implementations.TimeInvariantSsf;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import ec.tstoolkit2.ssf.univariate.FilteringErrors;
import ec.tstoolkit2.ssf.univariate.SsfData;
import java.util.Random;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author Jean Palate
 */
public class SsfUcarimaTest {

    final static int N = 100000;

    public final static SsfUcarima ssf;
    public final static SsfData ssfData;
    public final static ec.tstoolkit.ssf.ucarima.SsfUcarima ossf;
    public final static ec.tstoolkit.ssf.SsfData ossfData;

    final static Matrix P;
    final static DataBlock A;

    static {
        int M = 50;
        TrendCycleSelector tsel = new TrendCycleSelector(.5);
        tsel.setDefaultLowFreqThreshold(12);
        SeasonalSelector ssel = new SeasonalSelector(12, 3);

        ModelDecomposer decomposer = new ModelDecomposer();
        decomposer.add(tsel);
        decomposer.add(ssel);
        TsData x = data.Data.X.clone();
        int[] missing = new int[M];
        Random rng = new Random();
        for (int i = 0; i < M; ++i) {
            missing[i] = rng.nextInt(x.getLength());
        }
        SarimaModel arima = new SarimaModelBuilder().createAirlineModel(12, -.8, -.9);
        UcarimaModel ucm = decomposer.decompose(ArimaModel.create(arima));
        ucm.setVarianceMax(-1);
        ucm.simplify();

        for (int i = 0; i < M; ++i) {
            x.setMissing(missing[i]);
        }

        ssf = SsfUcarima.create(ucm);
        ssfData = new SsfData(x);

        ossf = new ec.tstoolkit.ssf.ucarima.SsfUcarima(ucm);
        ossfData = new ec.tstoolkit.ssf.SsfData(x, null);
        ec.tstoolkit.ssf.DiffusePredictionErrorDecomposition pe = new ec.tstoolkit.ssf.DiffusePredictionErrorDecomposition(false);
        ec.tstoolkit.ssf.State state = new ec.tstoolkit.ssf.State(ssf.getStateDim(), true);
        ec.tstoolkit.ssf.DurbinKoopmanInitializer dk = new ec.tstoolkit.ssf.DurbinKoopmanInitializer();
        dk.initialize(ossf, ossfData, state, pe);
        P = state.P;
        A = state.A;

    }

    public SsfUcarimaTest() {
    }

    @Test
    public void testDynamics() {
        ISsf ref = TimeInvariantSsf.of(ssf);
        int dim = ssf.getStateDim();
        Matrix M1 = new Matrix(dim, dim);
        M1.randomize();
        M1 = SymmetricMatrix.XXt(M1);
        Matrix M2 = M1.clone();
        Matrix M3 = M1.clone();
        DataBlock x1 = new DataBlock(dim);
        x1.randomize();
        DataBlock x2 = x1.deepClone();
        DataBlock x3 = x1.deepClone();
        ISsfDynamics dynref = ref.getDynamics();
        ISsfDynamics dyn = ssf.getDynamics();
        dynref.TX(0, x1);
        dyn.TX(dim, x2);
        ossf.TX(dim, x3);
        assertTrue(x1.distance(x2) < 1e-9);
        assertTrue(x1.distance(x3) < 1e-9);
        dynref.XT(0, x1);
        dyn.XT(0, x2);
        ossf.XT(0, x3);
        assertTrue(x1.distance(x2) < 1e-9);
        assertTrue(x1.distance(x3) < 1e-9);
        dynref.TVT(0, M1.subMatrix());
        dyn.TVT(0, M2.subMatrix());
        ossf.TVT(0, M3.subMatrix());
        assertTrue(M1.distance(M2) < 1e-9);
        assertTrue(M1.distance(M3) < 1e-9);
        dynref.TM(0, M1.subMatrix());
        dyn.TM(0, M2.subMatrix());
        assertTrue(M1.distance(M2) < 1e-9);
        dynref.MT(0, M1.subMatrix());
        dyn.MT(0, M2.subMatrix());
        assertTrue(M1.distance(M2) < 1e-9);
        Matrix Pi0 = Matrix.square(dim);
        Matrix oPi0 = Matrix.square(dim);
        Matrix B = new Matrix(dim, dyn.getNonStationaryDim());
        dyn.Pi0(Pi0.subMatrix());
        dyn.diffuseConstraints(B.subMatrix());
        ossf.Pi0(oPi0.subMatrix());
        assertTrue(Pi0.minus(SymmetricMatrix.XXt(B)).nrm2() < 1e-9);
        assertTrue(Pi0.minus(oPi0).nrm2() < 1e-9);
        Matrix Pf0 = Matrix.square(dim);
        Matrix oPf0 = Matrix.square(dim);
        dyn.Pf0(Pf0.subMatrix(), StateInfo.Forecast);
        ossf.Pf0(oPf0.subMatrix());
        assertTrue(Pf0.minus(oPf0).nrm2() < 1e-9);
    }

    @Test
    public void testMeasurement() {
        ISsf ref = TimeInvariantSsf.of(ssf);
        int dim = ssf.getStateDim();
        Matrix M1 = new Matrix(dim, dim);
        M1.randomize();
        M1 = SymmetricMatrix.XXt(M1);
        Matrix M2 = M1.clone();
        DataBlock x1 = new DataBlock(dim);
        x1.randomize();
        DataBlock x2 = x1.deepClone();
        ISsfMeasurement mref = ref.getMeasurement();
        ISsfMeasurement m = ssf.getMeasurement();
        assertTrue(Math.abs(mref.ZX(0, x1) - m.ZX(0, x2)) < 1e-9);
        assertTrue(Math.abs(mref.ZVZ(0, M1.subMatrix()) - m.ZVZ(0, M1.subMatrix())) < 1e-9);
        mref.VpZdZ(0, M1.subMatrix(), 5);
        m.VpZdZ(0, M2.subMatrix(), 5);
        assertTrue(M1.distance(M2) < 1e-9);
        mref.XpZd(0, x1, 5);
        m.XpZd(0, x2, 5);
        assertTrue(x1.distance(x2) < 1e-9);
        mref.ZM(dim, M1.subMatrix(), x1);
        m.ZM(dim, M1.subMatrix(), x2);
        assertTrue(x1.distance(x2) < 1e-9);
    }

    @Test
    public void testDK() {
        DurbinKoopmanInitializer dk = new DurbinKoopmanInitializer();
        FilteringErrors pe = new FilteringErrors(false);
        State state = new State(ssf.getStateDim());
        dk.initialize(state, ssf, ssfData);
        assertTrue(state.P().minus(P).nrm2() < 1e-6);
        assertTrue(state.a().distance(A) < 1e-6);
    }

    @Test
    public void testSQRDK() {
        DiffuseSquareRootInitializer dk = new DiffuseSquareRootInitializer();
        FilteringErrors pe = new FilteringErrors(false);
        State state = new State(ssf.getStateDim());
        dk.initialize(state, ssf, ssfData);
        assertTrue(state.P().minus(P).nrm2() < 1e-6);
        assertTrue(state.a().distance(A) < 1e-6);
    }

    @Test
    public void testAKF() {
        AugmentedFilter akf = new AugmentedFilter(true);
        AugmentedPredictionErrorDecomposition pe = new AugmentedPredictionErrorDecomposition(false);
        pe.prepare(ssf, ssfData.getCount());
        akf.process(ssf, ssfData, pe);
        assertTrue(akf.getState().P().minus(P).nrm2() < 1e-6);
        assertTrue(akf.getState().a().distance(A) < 1e-6);
    }

    @Ignore
    @Test
    public void testStressInitialisation() {
        long t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            AugmentedFilter akf = new AugmentedFilter(true);
            AugmentedPredictionErrorDecomposition pe = new AugmentedPredictionErrorDecomposition(false);
            pe.prepare(ssf, ssfData.getCount());
            akf.process(ssf, ssfData, pe);
        }
        long t1 = System.currentTimeMillis();
        System.out.println("akf");
        System.out.println(t1 - t0);
        System.out.println("dkold");
        t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            ec.tstoolkit.ssf.DurbinKoopmanInitializer dk = new ec.tstoolkit.ssf.DurbinKoopmanInitializer();
            ec.tstoolkit.ssf.DiffusePredictionErrorDecomposition pe = new ec.tstoolkit.ssf.DiffusePredictionErrorDecomposition(false);
            ec.tstoolkit.ssf.State state = new ec.tstoolkit.ssf.State(ossf.getStateDim(), true);
            dk.initialize(ossf, ossfData, state, pe);
        }
        t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
        System.out.println("dknew");
        t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            DurbinKoopmanInitializer dk = new DurbinKoopmanInitializer();
            FilteringErrors pe = new FilteringErrors(false);
            State state = new State(ssf.getStateDim());
            dk.initialize(state, ssf, ssfData);
        }
        t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
        System.out.println("square root");
        t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            DiffuseSquareRootInitializer dk = new DiffuseSquareRootInitializer();
            FilteringErrors pe = new FilteringErrors(false);
            State state = new State(ssf.getStateDim());
            dk.initialize(state, ssf, ssfData);
        }
        t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
    }

}
