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
package ec.tstoolkit2.ssf.akf;

import ec.tstoolkit.arima.ArimaModelBuilder;
import ec.tstoolkit.arima.estimation.RegArimaModel;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.eco.ILikelihood;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.sarima.SarimaModel;
import ec.tstoolkit.sarima.SarimaModelBuilder;
import ec.tstoolkit.timeseries.simplets.AverageInterpolator;
import ec.tstoolkit.utilities.IntList;
import ec.tstoolkit2.ssf.implementations.arima.SsfArima;
import ec.tstoolkit2.ssf.univariate.SsfData;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author Jean Palate
 */
public class AkfToolkitTest {
    
    private static final int N = 100000, M = 50000;
    static final SarimaModel model;
    static final double[] data;
    private static final double Err, Det, lDet, Ll, DLl, StLl;
    private static final Matrix P;
    private static final DataBlock A;

    static {
        SarimaModelBuilder builder = new SarimaModelBuilder();
        model = builder.createArimaModel(12, 3, 1, 1, 0, 1, 1);
        model.setBTheta(1, -.99);
        //model = builder.createAirlineModel(12, -.6, -.8);
        ArimaModelBuilder gbuilder = new ArimaModelBuilder();
        data = gbuilder.generate(model, 200);
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

        RegArimaModel regarima = new RegArimaModel(model);
        double[] y = data.clone();
        IntList missing = new IntList();
        AverageInterpolator.cleanMissings(y, missing);
        DataBlock Y = new DataBlock(y);
        regarima.setY(Y);
        regarima.setMissings(missing.toArray());
        StLl=regarima.computeLikelihood().getLogLikelihood();
    }
    
    public AkfToolkitTest() {
    }

    @Test
    public void testLikelihood() {
        SsfArima ssf = SsfArima.create(model);
        SsfData ssfData = new SsfData(data);
        AkfDiffuseLikelihood ll = (AkfDiffuseLikelihood) AkfToolkit.likelihoodComputer(false).compute(ssf, ssfData);
        assertTrue(Math.abs(ll.getSsqErr() - Err) / Err < 1e-6);
        assertTrue(Math.abs(ll.getLogDeterminant() + ll.getDiffuseCorrection() - Det) < 1e-6);
        assertTrue(Math.abs(ll.getLogLikelihood() - DLl) < 1e-6);
        assertTrue(Math.abs(ll.getLogLikelihood() - StLl) < 1e-6);
    }

    @Test
    public void testCollapsing() {
        SsfArima ssf = SsfArima.create(model);
        SsfData ssfData = new SsfData(data);
        AkfDiffuseLikelihood ll = (AkfDiffuseLikelihood) AkfToolkit.likelihoodComputer(true).compute(ssf, ssfData);
        assertTrue(Math.abs(ll.getSsqErr() - Err) / Err < 1e-6);
        assertTrue(Math.abs(ll.getLogDeterminant() + ll.getDiffuseCorrection() - Det) < 1e-6);
        assertTrue(Math.abs(ll.getLogLikelihood() - DLl) < 1e-6);
        assertTrue(Math.abs(ll.getLogLikelihood() - StLl) < 1e-6);
    }
    
    @Ignore
    @Test
    public void testStressLikelihood() {
        long t0 = System.currentTimeMillis();
        for (int i = 0; i < M; ++i) {
            SsfArima ssf = SsfArima.create(model);
            ILikelihood ll = AkfToolkit.likelihoodComputer(false).compute(ssf, new SsfData(data));
        }
        long t1 = System.currentTimeMillis();
        System.out.println("AKF");
        System.out.println(t1 - t0);
        t0 = System.currentTimeMillis();
        for (int i = 0; i < M; ++i) {
            SsfArima ssf = SsfArima.create(model);
            ILikelihood ll = AkfToolkit.likelihoodComputer(true).compute(ssf, new SsfData(data));
        }
        t1 = System.currentTimeMillis();
        System.out.println("AKF with collapsing");
        System.out.println(t1 - t0);
        t0 = System.currentTimeMillis();
        for (int i = 0; i < M; ++i) {
            ec.tstoolkit.ssf.arima.SsfArima ssf = new ec.tstoolkit.ssf.arima.SsfArima(model);
            ec.tstoolkit.ssf.Filter filter = new ec.tstoolkit.ssf.Filter();
            ec.tstoolkit.ssf.DiffusePredictionErrorDecomposition pe = new ec.tstoolkit.ssf.DiffusePredictionErrorDecomposition(false);
            filter.setSsf(ssf);
            filter.process(new ec.tstoolkit.ssf.SsfData(data, null), pe);
        }
        t1 = System.currentTimeMillis();
        System.out.println("Old DK Filter");
        System.out.println(t1 - t0);
    }

}
