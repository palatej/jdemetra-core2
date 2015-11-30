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
import ec.tstoolkit.eco.ILikelihood;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.sarima.SarimaModel;
import ec.tstoolkit.sarima.SarimaModelBuilder;
import ec.tstoolkit2.ssf.State;
import ec.tstoolkit2.ssf.dk.sqrt.DiffuseSquareRootInitializer;
import ec.tstoolkit2.ssf.dk.DkToolkit;
import static ec.tstoolkit2.ssf.implementations.structural.LocalLinearTrendTest.data;
import ec.tstoolkit2.ssf.univariate.FilteringErrors;
import ec.tstoolkit2.ssf.univariate.SsfData;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author Jean Palate
 */
public class LocalLevelTest {

    static final double[] data;
    private static final int N = 1000000;

    static {
        SarimaModelBuilder builder = new SarimaModelBuilder();
        SarimaModel arima = builder.createAirlineModel(12, -.6, -.8);
        ArimaModelBuilder gbuilder = new ArimaModelBuilder();
        data = gbuilder.generate(arima, 200);
        data[1] = Double.NaN;
        data[2] = Double.NaN;
        data[3] = Double.NaN;
    }

    public LocalLevelTest() {
    }

    @Test
    public void test_LL_Bsm() {
        ModelSpecification spec = new ModelSpecification();
        spec.useLevel(ComponentUse.Free);
        spec.useSlope(ComponentUse.Unused);
        spec.useCycle(ComponentUse.Unused);
        spec.useNoise(ComponentUse.Free);
        spec.setSeasonalModel(SeasonalModel.Unused);

        BasicStructuralModel bsm = new BasicStructuralModel(spec, 1);
        bsm.setVariance(Component.Level, .1);
        bsm.setVariance(Component.Noise, 1);

        SsfBsm ssfbsm = SsfBsm.create(bsm);

        LocalLevel ll = new LocalLevel(.1, 1);

        DiffuseSquareRootInitializer dk = new DiffuseSquareRootInitializer();
        FilteringErrors pe = new FilteringErrors(false);
        SsfData ssfdata = new SsfData(data);
        State state = new State(ssfbsm.getStateDim());
        dk.initialize(state, ssfbsm, ssfdata);
        FilteringErrors lpe = new FilteringErrors(false);
        State lstate = new State(ll.getStateDim());
        dk.initialize(lstate, ll, ssfdata);
        assertTrue(new Matrix(state.P().subMatrix(1, 2, 1, 2)).distance(lstate.P()) < 1e-9);
        assertTrue(state.a().range(1, 2).distance(lstate.a()) < 1e-9);
    }

    @Ignore
    @Test
    public void stressTest_LL_Bsm() {
        ModelSpecification spec = new ModelSpecification();
        spec.useLevel(ComponentUse.Free);
        spec.useSlope(ComponentUse.Unused);
        spec.useCycle(ComponentUse.Unused);
        spec.useNoise(ComponentUse.Free);
        spec.setSeasonalModel(SeasonalModel.Unused);

        BasicStructuralModel bsm = new BasicStructuralModel(spec, 1);
        bsm.setVariance(Component.Level, .1);
        bsm.setVariance(Component.Noise, 1);

        SsfBsm ssfbsm = SsfBsm.create(bsm);

        LocalLevel ll = new LocalLevel(.1, 1);
        long t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            SsfData ssfdata = new SsfData(data);
            double eval = DkToolkit.likelihoodComputer().compute(ssfbsm, ssfdata).getLogLikelihood();
        }
        long t1 = System.currentTimeMillis();
        System.out.println("bsm");
        System.out.println(t1 - t0);
        t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            SsfData ssfdata = new SsfData(data);
            double eval = DkToolkit.likelihoodComputer().compute(ll, ssfdata).getLogLikelihood();
        }
        t1 = System.currentTimeMillis();
        System.out.println("ll");
        System.out.println(t1 - t0);
    }
}
