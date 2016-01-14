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

import ec.tstoolkit.arima.ArimaModel;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.sarima.SarimaModel;
import ec.tstoolkit.sarima.SarimaModelBuilder;
import ec.tstoolkit.timeseries.simplets.TsData;
import ec.tstoolkit.ucarima.ModelDecomposer;
import ec.tstoolkit.ucarima.SeasonalSelector;
import ec.tstoolkit.ucarima.TrendCycleSelector;
import ec.tstoolkit.ucarima.UcarimaModel;
import ec.tstoolkit2.ssf.ResultsRange;
import ec.tstoolkit2.ssf.dk.DefaultDiffuseFilteringResults;
import ec.tstoolkit2.ssf.dk.DkToolkit;
import ec.tstoolkit2.ssf.implementations.arima.SsfUcarima;
import java.util.Random;
import static org.junit.Assert.assertTrue;
import org.junit.Ignore;
import org.junit.Test;

/**
 *
 * @author Jean Palate
 */
public class FastFilterTest {

    final static int N = 100000;

    public final static SsfUcarima ssf;
    public final static SsfData ssfData;
    public final static TsData xts;

    static {
        int M = 0;
        TrendCycleSelector tsel = new TrendCycleSelector(.5);
        tsel.setDefaultLowFreqThreshold(12);
        SeasonalSelector ssel = new SeasonalSelector(12, 3);

        ModelDecomposer decomposer = new ModelDecomposer();
        decomposer.add(tsel);
        decomposer.add(ssel);
        xts = data.Data.X.clone();
        int[] missing = new int[M];
        Random rng = new Random();
        for (int i = 0; i < M; ++i) {
            missing[i] = rng.nextInt(xts.getLength());
        }
        SarimaModel arima = new SarimaModelBuilder().createAirlineModel(12, -.8, -.9);
        UcarimaModel ucm = decomposer.decompose(ArimaModel.create(arima));
        ucm.setVarianceMax(-1);
        ucm.simplify();

        for (int i = 0; i < M; ++i) {
            xts.setMissing(missing[i]);
        }

        ssf = SsfUcarima.create(ucm);
        ssfData = new SsfData(xts);
    }

    public FastFilterTest() {
    }

    @Test
    @Ignore
    public void stressTestUcarima() {
        DefaultDiffuseFilteringResults fresults = DkToolkit.filter(ssf, ssfData, false);
        Matrix x = new Matrix(ssfData.getCount(), 10);
        x.randomize(0);
        long t0 = System.currentTimeMillis();
        for (int i = 0; i < 10000; ++i) {
            FastFilter filter = new FastFilter(ssf, fresults, new ResultsRange(fresults.getEndDiffusePosition(), ssfData.getCount()));
            filter.filter(x.subMatrix());
        }
        long t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
    }

    @Test
    public void testUcarima() {
        DefaultDiffuseFilteringResults fresults = DkToolkit.filter(ssf, ssfData, false);
        Matrix x = new Matrix(ssfData.getCount(), 1);
        x.column(0).copy(xts);
        Matrix x0=new Matrix(ssf.getStateDim(), 1);
        int nd=fresults.getEndDiffusePosition();
        x0.column(0).copy(fresults.a(nd));
        FastFilter filter = new FastFilter(ssf, fresults, new ResultsRange(nd, ssfData.getCount()));
        filter.filter(x0.subMatrix(), x.subMatrix(nd,-1,0,1));
        assertTrue(new DataBlock(fresults.errors()).drop(nd, 0)
                .distance(x.column(0).drop(nd, 0))<1e-8);
    }
}