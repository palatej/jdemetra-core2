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

import data.Data;
import ec.tstoolkit.arima.ArimaModel;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.DataBlockStorage;
import ec.tstoolkit.sarima.SarimaModel;
import ec.tstoolkit.sarima.SarimaModelBuilder;
import ec.tstoolkit.ssf.SmoothingResults;
import ec.tstoolkit.timeseries.simplets.TsData;
import ec.tstoolkit.ucarima.ModelDecomposer;
import ec.tstoolkit.ucarima.SeasonalSelector;
import ec.tstoolkit.ucarima.TrendCycleSelector;
import ec.tstoolkit.ucarima.UcarimaModel;
import ec.tstoolkit2.ssf.implementations.arima.SsfUcarima;
import ec.tstoolkit2.ssf.univariate.DefaultDisturbanceSmoothingResults;
import ec.tstoolkit2.ssf.univariate.DefaultSmoothingResults;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.SsfData;
import java.util.Random;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author Jean Palate
 */
public class DiffuseDisturbanceSmootherTest {

    static final int N = 50;
    static final UcarimaModel ucm;
    static final SsfData data;
    static final ec.tstoolkit.ssf.SsfData odata;

    static {
        TrendCycleSelector tsel = new TrendCycleSelector(.5);
        tsel.setDefaultLowFreqThreshold(12);
        SeasonalSelector ssel = new SeasonalSelector(12, 3);

        ModelDecomposer decomposer = new ModelDecomposer();
        decomposer.add(tsel);
        decomposer.add(ssel);
        TsData x = Data.X.clone();
        int[] missing = new int[N];
        Random rng = new Random(0);
        for (int i = 0; i < N; ++i) {
            missing[i] = rng.nextInt(x.getLength());
        }
        SarimaModel arima = new SarimaModelBuilder().createAirlineModel(12, -.8, -.9);
        ucm = decomposer.decompose(ArimaModel.create(arima));
        ucm.setVarianceMax(-1);
        ucm.simplify();

        for (int i = 0; i < N; ++i) {
            x.setMissing(missing[i]);
        }
        data = new SsfData(x);
        odata = new ec.tstoolkit.ssf.SsfData(x, null);
    }

    public DiffuseDisturbanceSmootherTest() {
    }

    @Test
    public void testStateSmoothing() {
        ec.tstoolkit.ssf.DisturbanceSmoother osmoother = new ec.tstoolkit.ssf.DisturbanceSmoother();
        osmoother.setSsf(new ec.tstoolkit.ssf.ucarima.SsfUcarima(ucm));
        osmoother.process(odata);
        SmoothingResults ss = osmoother.calcSmoothedStates();
        DataBlock ot = new DataBlock(ss.component(0));
        ISsf ssf = SsfUcarima.create(ucm);
        FastStateSmoother smoother = new FastStateSmoother();
        DataBlockStorage str = smoother.process(ssf, data);
        assertTrue(ot.distance(str.item(0)) < 1e-6);
        DefaultSmoothingResults str2 = DkToolkit.smooth(ssf, data, false);
        assertTrue(str2.getComponent(0).distance(str.item(0)) < 1e-6);
    }

    @Test
    @Ignore
    public void stressTestStateSmoothing() {
        int K = 10000;
        ec.tstoolkit.ssf.ucarima.SsfUcarima ossf = new ec.tstoolkit.ssf.ucarima.SsfUcarima(ucm);
        long t0 = System.currentTimeMillis();
        for (int i = 0; i < K; ++i) {
            ec.tstoolkit.ssf.DisturbanceSmoother osmoother = new ec.tstoolkit.ssf.DisturbanceSmoother();
            osmoother.setSsf(ossf);
            osmoother.process(odata);
            SmoothingResults ss = osmoother.calcSmoothedStates();
            DataBlock ot = new DataBlock(ss.component(0));
        }
        long t1 = System.currentTimeMillis();
        System.out.println("Old smoother");
        System.out.println(t1 - t0);
        ISsf ssf = SsfUcarima.create(ucm);
        t0 = System.currentTimeMillis();
        for (int i = 0; i < K; ++i) {
            FastStateSmoother smoother = new FastStateSmoother();
            DataBlockStorage str = smoother.process(ssf, data);
            DataBlock ot = str.item(0);
        }
        t1 = System.currentTimeMillis();
        System.out.println("New smoother");
        System.out.println(t1 - t0);
        t0 = System.currentTimeMillis();
        for (int i = 0; i < K; ++i) {
            DefaultSmoothingResults sr = DkToolkit.sqrtSmooth(ssf, data, false);
            DataBlock ot=sr.getComponent(0);
        }
        t1 = System.currentTimeMillis();
        System.out.println("SQRT DK smoother");
        System.out.println(t1 - t0);
    }
}
