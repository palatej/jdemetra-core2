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
package ec.tstoolkit2.ssf.implementations.structural;

import ec.tstoolkit.arima.ArimaModelBuilder;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.eco.ILikelihood;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit.sarima.SarimaModel;
import ec.tstoolkit.sarima.SarimaModelBuilder;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.akf.AkfToolkit;
import ec.tstoolkit2.ssf.akf.AugmentedFilter;
import ec.tstoolkit2.ssf.akf.AugmentedPredictionErrorDecomposition;
import ec.tstoolkit2.ssf.akf.AugmentedState;
import ec.tstoolkit2.ssf.akf.AkfDiffuseLikelihood;
import ec.tstoolkit2.ssf.dk.DkDiffuseLikelihood;
import ec.tstoolkit2.ssf.dk.sqrt.DiffuseSquareRootInitializer;
import ec.tstoolkit2.ssf.dk.DurbinKoopmanInitializer;
import ec.tstoolkit2.ssf.dk.DkToolkit;
import ec.tstoolkit2.ssf.implementations.TimeInvariantSsf;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import ec.tstoolkit2.ssf.univariate.FilteringErrors;
import ec.tstoolkit2.ssf.univariate.Ssf;
import ec.tstoolkit2.ssf.univariate.SsfData;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author Jean Palate
 */
public class SsfBsmTest {

    BasicStructuralModel model;
    static final double[] data;
    private static final int N = 100000;

    private static final Matrix P;
    private static final DataBlock A;
    private static final ec.tstoolkit.eco.DiffuseLikelihood OldLL;

    static {
        SarimaModelBuilder builder = new SarimaModelBuilder();
        SarimaModel arima = builder.createAirlineModel(12, -.6, -.8);
        ArimaModelBuilder gbuilder = new ArimaModelBuilder();
        data = gbuilder.generate(arima, 50);
        data[5] = Double.NaN;
        data[12] = Double.NaN;
        data[21] = Double.NaN;

        ec.tstoolkit.structural.ModelSpecification spec = new ec.tstoolkit.structural.ModelSpecification();
        spec.useLevel(ec.tstoolkit.structural.ComponentUse.Free);
        spec.useSlope(ec.tstoolkit.structural.ComponentUse.Free);
        spec.useCycle(ec.tstoolkit.structural.ComponentUse.Free);
        spec.useNoise(ec.tstoolkit.structural.ComponentUse.Free);
        spec.setSeasonalModel(ec.tstoolkit.structural.SeasonalModel.Trigonometric);

        ec.tstoolkit.structural.BasicStructuralModel omodel = new ec.tstoolkit.structural.BasicStructuralModel(spec, 12);
        omodel.setCycle(.9, 8);
        omodel.setVariance(ec.tstoolkit.structural.Component.Level, .1);
        omodel.setVariance(ec.tstoolkit.structural.Component.Slope, .2);
        omodel.setVariance(ec.tstoolkit.structural.Component.Cycle, .5);
        omodel.setVariance(ec.tstoolkit.structural.Component.Seasonal, 2);
        omodel.setVariance(ec.tstoolkit.structural.Component.Noise, 1);
        ec.tstoolkit.ssf.DurbinKoopmanInitializer dk = new ec.tstoolkit.ssf.DurbinKoopmanInitializer();

        ec.tstoolkit.ssf.DiffusePredictionErrorDecomposition pe = new ec.tstoolkit.ssf.DiffusePredictionErrorDecomposition(false);
        ec.tstoolkit.ssf.SsfData ssfdata = new ec.tstoolkit.ssf.SsfData(data, null);
        ec.tstoolkit.ssf.State state = new ec.tstoolkit.ssf.State(omodel.getStateDim(), true);
        dk.initialize(omodel, ssfdata, state, pe);
        A = state.A;
        P = state.P;

        ec.tstoolkit.ssf.Filter filter = new ec.tstoolkit.ssf.Filter();
        filter.setSsf(omodel);
        filter.process(ssfdata, pe);

        OldLL = new ec.tstoolkit.eco.DiffuseLikelihood();
        ec.tstoolkit.ssf.LikelihoodEvaluation.evaluate(pe, OldLL);

    }

    public SsfBsmTest() {
        ModelSpecification spec = new ModelSpecification();
        spec.useLevel(ComponentUse.Free);
        spec.useSlope(ComponentUse.Free);
        spec.useCycle(ComponentUse.Free);
        spec.useNoise(ComponentUse.Free);
        spec.setSeasonalModel(SeasonalModel.Trigonometric);

        model = new BasicStructuralModel(spec, 12);
        model.setCycle(.9, 8);
        model.setVariance(Component.Level, .1);
        model.setVariance(Component.Slope, .2);
        model.setVariance(Component.Cycle, .5);
        model.setVariance(Component.Seasonal, 2);
        model.setVariance(Component.Noise, 1);
    }

    @Test
    public void testDynamics() {
        SsfBsm bsm = SsfBsm.create(model);
        ISsf ref = TimeInvariantSsf.of(bsm);
        int dim = bsm.getStateDim();
        Matrix M1 = new Matrix(dim, dim);
        M1.randomize();
        M1 = SymmetricMatrix.XXt(M1);
        Matrix M2 = M1.clone();
        DataBlock x1 = new DataBlock(dim);
        x1.randomize();
        DataBlock x2 = x1.deepClone();
        ISsfDynamics dynref = ref.getDynamics();
        ISsfDynamics dyn = bsm.getDynamics();
        dynref.TX(0, x1);
        dyn.TX(dim, x2);
        assertTrue(x1.distance(x2) < 1e-9);
        dynref.XT(0, x1);
        dyn.XT(0, x2);
        assertTrue(x1.distance(x2) < 1e-9);
        dynref.TVT(0, M1.subMatrix());
        dyn.TVT(0, M2.subMatrix());
        assertTrue(M1.distance(M2) < 1e-9);
        dynref.TM(0, M1.subMatrix());
        dyn.TM(0, M2.subMatrix());
        assertTrue(M1.distance(M2) < 1e-9);
        dynref.MT(0, M1.subMatrix());
        dyn.MT(0, M2.subMatrix());
        assertTrue(M1.distance(M2) < 1e-9);
        Matrix Pi0 = Matrix.square(dim);
        Matrix B = new Matrix(dim, dyn.getNonStationaryDim());
        dyn.Pi0(Pi0.subMatrix());
        dyn.diffuseConstraints(B.subMatrix());
        assertTrue(Pi0.minus(SymmetricMatrix.XXt(B)).nrm2() < 1e-9);
    }

    @Test
    public void testMeasurement() {
        SsfBsm bsm = SsfBsm.create(model);
        ISsf ref = TimeInvariantSsf.of(bsm);
        int dim = bsm.getStateDim();
        Matrix M1 = new Matrix(dim, dim);
        M1.randomize();
        M1 = SymmetricMatrix.XXt(M1);
        Matrix M2 = M1.clone();
        DataBlock x1 = new DataBlock(dim);
        x1.randomize();
        DataBlock x2 = x1.deepClone();
        ISsfMeasurement mref = ref.getMeasurement();
        ISsfMeasurement m = bsm.getMeasurement();
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
    public void testDkLikelihood() {
        Ssf ssf = SsfBsm.create(model);
        SsfData ssfData = new SsfData(data);
        DkDiffuseLikelihood ll = (DkDiffuseLikelihood) DkToolkit.likelihoodComputer(false, false).compute(ssf, ssfData);
        assertTrue(Math.abs(ll.getSsqErr()-OldLL.getSsqErr())<1e-6);
        Ssf ssf2 = SsfBsm2.create(model);
        DkDiffuseLikelihood ll2 = (DkDiffuseLikelihood) DkToolkit.likelihoodComputer(true, false).compute(ssf2, ssfData);
        assertTrue(Math.abs(ll2.getSsqErr()-OldLL.getSsqErr())<1e-6);
        AkfDiffuseLikelihood ll3 = (AkfDiffuseLikelihood) AkfToolkit.likelihoodComputer(true).compute(ssf, ssfData);
        AkfDiffuseLikelihood ll4 = (AkfDiffuseLikelihood) AkfToolkit.likelihoodComputer(false).compute(ssf, ssfData);
        assertEquals(ll.getLogLikelihood(), ll2.getLogLikelihood(), 1e-8);
        assertEquals(ll.getLogLikelihood(), ll3.getLogLikelihood(), 1e-8);
        assertEquals(ll.getLogLikelihood(), ll4.getLogLikelihood(), 1e-8);
    }

    @Test
    public void testAkfLikelihood() {
        Ssf ssf = SsfBsm.create(model);
        SsfData ssfData = new SsfData(data);
        AkfDiffuseLikelihood ll = (AkfDiffuseLikelihood) AkfToolkit.likelihoodComputer(true).compute(ssf, ssfData);
        assertTrue(Math.abs(ll.getSsqErr()-OldLL.getSsqErr())<1e-6);
        assertTrue(Math.abs(ll.getLogDeterminant()+ll.getDiffuseCorrection()-OldLL.getLogDeterminant()-OldLL.getDiffuseLogDeterminant())<1e-6);
    }

    @Test
    public void testDK() {
        DurbinKoopmanInitializer dk = new DurbinKoopmanInitializer();
        //DiffuseSquareRootInitializer dk = new DiffuseSquareRootInitializer();
        Ssf ssf = SsfBsm.create(model);
        FilteringErrors pe = new FilteringErrors(false);
        SsfData ssfdata = new SsfData(data);
        State state = new State(ssf.getStateDim());
        dk.initialize(state, ssf, ssfdata);
        assertTrue(state.P().minus(P).nrm2() < 1e-9);
        assertTrue(state.a().distance(A) < 1e-9);
    }

    @Test
    public void testSQRDK() {
        //DurbinKoopmanInitializer dk = new DurbinKoopmanInitializer();
        DiffuseSquareRootInitializer dk = new DiffuseSquareRootInitializer();
        Ssf ssf = SsfBsm.create(model);
        FilteringErrors pe = new FilteringErrors(false);
        SsfData ssfdata = new SsfData(data);
        State state = new State(ssf.getStateDim());
        dk.initialize(state, ssf, ssfdata);
        assertTrue(state.P().minus(P).nrm2() < 1e-9);
        assertTrue(state.a().distance(A) < 1e-9);
    }

    @Test
    public void testAKF() {
        //DurbinKoopmanInitializer dk = new DurbinKoopmanInitializer();
        AugmentedFilter akf = new AugmentedFilter(true);
        Ssf ssf = SsfBsm.create(model);
        AugmentedPredictionErrorDecomposition pe = new AugmentedPredictionErrorDecomposition();
        akf.process(ssf, new SsfData(data), pe);
        AugmentedState state = akf.getState();
        assertTrue(state.P().minus(P).nrm2() < 1e-9);
        assertTrue(state.a().distance(A) < 1e-9);
    }

    @Test
    @Ignore
    public void testStressInitialisation() {
        long t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            AugmentedFilter akf = new AugmentedFilter(true);
            Ssf ssf = SsfBsm.create(model);
            AugmentedPredictionErrorDecomposition pe = new AugmentedPredictionErrorDecomposition();
            akf.process(ssf, new SsfData(data), pe);
        }
        long t1 = System.currentTimeMillis();
        System.out.println("akf");
        System.out.println(t1 - t0);
        System.out.println("dkold");
        t0 = System.currentTimeMillis();
        ec.tstoolkit.structural.ModelSpecification spec = new ec.tstoolkit.structural.ModelSpecification();
        spec.useLevel(ec.tstoolkit.structural.ComponentUse.Free);
        spec.useSlope(ec.tstoolkit.structural.ComponentUse.Free);
        spec.useCycle(ec.tstoolkit.structural.ComponentUse.Free);
        spec.useNoise(ec.tstoolkit.structural.ComponentUse.Free);
        spec.setSeasonalModel(ec.tstoolkit.structural.SeasonalModel.Crude);

        ec.tstoolkit.structural.BasicStructuralModel omodel = new ec.tstoolkit.structural.BasicStructuralModel(spec, 12);
        omodel.setCycle(.9, 8);
        omodel.setVariance(ec.tstoolkit.structural.Component.Level, .1);
        omodel.setVariance(ec.tstoolkit.structural.Component.Slope, .2);
        omodel.setVariance(ec.tstoolkit.structural.Component.Cycle, .5);
        omodel.setVariance(ec.tstoolkit.structural.Component.Seasonal, 2);
        omodel.setVariance(ec.tstoolkit.structural.Component.Noise, 1);
        for (int i = 0; i < N; ++i) {
            ec.tstoolkit.ssf.DurbinKoopmanInitializer dk = new ec.tstoolkit.ssf.DurbinKoopmanInitializer();
            ec.tstoolkit.ssf.DiffusePredictionErrorDecomposition pe = new ec.tstoolkit.ssf.DiffusePredictionErrorDecomposition(false);
            ec.tstoolkit.ssf.SsfData ssfdata = new ec.tstoolkit.ssf.SsfData(data, null);
            ec.tstoolkit.ssf.State state = new ec.tstoolkit.ssf.State(omodel.getStateDim(), true);
            dk.initialize(omodel, ssfdata, state, pe);
        }
        t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
        System.out.println("dknew");
        t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            DurbinKoopmanInitializer dk = new DurbinKoopmanInitializer();
            Ssf ssf = SsfBsm.create(model);
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
            DiffuseSquareRootInitializer dk = new DiffuseSquareRootInitializer();
            Ssf ssf = SsfBsm.create(model);
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
            DiffuseSquareRootInitializer dk = new DiffuseSquareRootInitializer();
            Ssf ssf = SsfBsm.create(model);
            FilteringErrors pe = new FilteringErrors(false);
            SsfData ssfdata = new SsfData(data);
            State state = new State(ssf.getStateDim());
            dk.initialize(state, TimeInvariantSsf.of(ssf), ssfdata);
        }
        t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
    }

}
