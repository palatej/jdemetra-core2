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
import ec.tstoolkit.arima.estimation.RegArimaModel;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.eco.ILikelihood;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.realfunctions.IParametricMapping;
import ec.tstoolkit.maths.realfunctions.ISsqFunction;
import ec.tstoolkit.maths.realfunctions.levmar.LevenbergMarquardtMethod;
import ec.tstoolkit.maths.realfunctions.minpack.LevenbergMarquardtMinimizer;
import ec.tstoolkit.sarima.SarimaModel;
import ec.tstoolkit.sarima.SarimaModelBuilder;
import ec.tstoolkit.timeseries.simplets.AverageInterpolator;
import ec.tstoolkit.timeseries.simplets.TsData;
import ec.tstoolkit.utilities.IntList;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.implementations.TimeInvariantSsf;
import ec.tstoolkit2.ssf.implementations.arima.SsfArima;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.SsfData;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author Jean Palate
 */
public class DkToolkitTest {

    private static final int N = 100000, M = 10000;
    static final SarimaModel model;
    static final double[] data;
    private static final double Err, Det, lDet, Ll, DLl, StLl;
    private static final Matrix P;
    private static final DataBlock A;

    static {
        SarimaModelBuilder builder = new SarimaModelBuilder();
        model = builder.createArimaModel(12, 3, 1, 1, 0, 1, 1);
        double[] p = new double[]{-.3, -.3, -.3, -.5, -.9};
        model.setParameters(new DataBlock(p));
        //model = builder.createAirlineModel(12, -.6, -.8);
        TsData s = Data.X;
        data = new double[s.getLength()];
        s.copyTo(data, 0);
        data[2] = Double.NaN;
        data[11] = Double.NaN;
        data[119] = Double.NaN;
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
        StLl = regarima.computeLikelihood().getLogLikelihood();
    }

    public DkToolkitTest() {
    }

    @Test
    public void testLikelihood() {
        SsfArima ssf = SsfArima.create(model);
        SsfData ssfData = new SsfData(data);
        DkDiffuseLikelihood ll = (DkDiffuseLikelihood) DkToolkit.likelihoodComputer(false, true).compute(ssf, ssfData);
        assertTrue(Math.abs(ll.getSsqErr() - Err) / Err < 1e-6);
        assertTrue(Math.abs(ll.getLogDeterminant() + ll.getDiffuseCorrection() - Det) < 1e-6);
        assertTrue(Math.abs(ll.getLogLikelihood() - DLl) < 1e-6);
        assertTrue(Math.abs(ll.getLogLikelihood() - StLl) < 1e-6);
        assertTrue(Math.abs(Ll - DLl + 0.5 * ll.getD() * Math.log(2 * Math.PI)) < 1e-9);
    }

    @Test
    @Ignore
    public void testEstimation() {
        long t0 = System.currentTimeMillis();
        for (int i = 0; i < 1000; ++i) {
            SsfArima ssf = SsfArima.create(model);
            SsfData ssfData = new SsfData(data);
            IParametricMapping<SsfArima> mapping = SsfArima.mapping(model.getSpecification());
            SsfFunction<SsfArima> fn = DkToolkit.likelihoodFunction(ssf, ssfData, mapping);
            //LevenbergMarquardtMinimizer opt = new LevenbergMarquardtMinimizer();
            LevenbergMarquardtMethod opt = new LevenbergMarquardtMethod();
            opt.minimize(fn, mapping.map(ssf));
            SsfFunctionInstance<SsfArima> result = (SsfFunctionInstance) opt.getResult();
//            System.out.println(result.getLikelihood().getLogLikelihood());
//            System.out.println(result.getLikelihood().getSigma());
        }
        long t1 = System.currentTimeMillis();
        System.out.println("Estimation");
        System.out.println(t1 - t0);
    }

    @Test
    public void testSqr() {
        SsfArima ssf = SsfArima.create(model);
        SsfData ssfData = new SsfData(data);
        DkDiffuseLikelihood ll = (DkDiffuseLikelihood) DkToolkit.likelihoodComputer(true, true).compute(ssf, ssfData);
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
            ILikelihood ll = DkToolkit.likelihoodComputer(false).compute(ssf, new SsfData(data));
        }
        long t1 = System.currentTimeMillis();
        System.out.println("DK (normal)");
        System.out.println(t1 - t0);
        t0 = System.currentTimeMillis();
        for (int i = 0; i < M; ++i) {
            SsfArima ssf = SsfArima.create(model);
            ILikelihood ll = DkToolkit.likelihoodComputer(true).compute(ssf, new SsfData(data));
        }
        t1 = System.currentTimeMillis();
        System.out.println("DK (square root form)");
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
        t0 = System.currentTimeMillis();
        for (int i = 0; i < M; ++i) {
            SsfArima ssf = SsfArima.create(model);
            ISsf tssf = TimeInvariantSsf.of(ssf, StateInfo.Forecast);
            ILikelihood ll = DkToolkit.likelihoodComputer(true).compute(tssf, new SsfData(data));
        }

        t1 = System.currentTimeMillis();
        System.out.println("DK Filter. Matrix");
        System.out.println(t1 - t0);
    }

}
