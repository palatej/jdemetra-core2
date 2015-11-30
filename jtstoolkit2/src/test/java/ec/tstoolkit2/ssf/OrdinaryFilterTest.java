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
package ec.tstoolkit2.ssf;

import ec.tstoolkit2.ssf.multivariate.SsfMatrix;
import ec.tstoolkit.arima.estimation.ArmaKF;
import ec.tstoolkit.eco.ILikelihood;
import ec.tstoolkit.eco.Likelihood;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.sarima.SarimaModel;
import ec.tstoolkit.sarima.SarimaSpecification;
import ec.tstoolkit2.ssf.array.ArrayFilter;
import ec.tstoolkit2.ssf.array.MultivariateArrayFilter;
import ec.tstoolkit2.ssf.ckms.FastFilter;
import ec.tstoolkit2.ssf.implementations.MultivariateTimeInvariantSsf;
import ec.tstoolkit2.ssf.implementations.TimeInvariantSsf;
import ec.tstoolkit2.ssf.implementations.arima.SsfArima;
import ec.tstoolkit2.ssf.multivariate.MultivariateOrdinaryFilter;
import ec.tstoolkit2.ssf.multivariate.PredictionErrorsDecomposition;
import ec.tstoolkit2.ssf.multivariate.MultivariateSsf;
import ec.tstoolkit2.ssf.univariate.OrdinaryFilter;
import ec.tstoolkit2.ssf.univariate.PredictionErrorDecomposition;
import ec.tstoolkit2.ssf.univariate.SsfData;
import static org.junit.Assert.assertTrue;
import org.junit.Ignore;
import org.junit.Test;

/**
 *
 * @author Jean Palate
 */
public class OrdinaryFilterTest {

    private static final SarimaModel model;
    private static final Matrix M = new Matrix(360, 1);
    private static final double Err, Det;

    static {
        SarimaSpecification spec = new SarimaSpecification(12);
        spec.setP(1);
        spec.setQ(3);
        //spec.setBP(1);
        spec.setBQ(1);
        model = new SarimaModel(spec);
        M.randomize(0);
        ec.tstoolkit.ssf.Filter filter = new ec.tstoolkit.ssf.Filter();
        ec.tstoolkit.ssf.arima.SsfArima ssf
                = new ec.tstoolkit.ssf.arima.SsfArima(model);
        ec.tstoolkit.ssf.PredictionErrorDecomposition pe = new ec.tstoolkit.ssf.PredictionErrorDecomposition(false);
        filter.setSsf(ssf);
        filter.process(new ec.tstoolkit.ssf.SsfData(M.column(0), null), pe);
        Err = pe.getSsqErr();
        Det = pe.getLogDeterminant();
    }

    public OrdinaryFilterTest() {
    }

    @Test
    public void testNewMethod() {
        MultivariateOrdinaryFilter filter = new MultivariateOrdinaryFilter();
        SsfArima ssf = SsfArima.create(model);
        PredictionErrorsDecomposition pe = new PredictionErrorsDecomposition(false);
        filter.process(MultivariateSsf.proxy(ssf), new SsfMatrix(M), pe);
        ILikelihood ll = pe.likelihood();
        assertTrue(Math.abs(Err - ll.getSsqErr()) / Err < 1e-9);
        assertTrue(Math.abs(Det - ll.getLogDeterminant()) < 1e-9);
    }

    @Test
    public void testNew1Method() {
        OrdinaryFilter filter = new OrdinaryFilter();
        SsfArima ssf = SsfArima.create(model);
        PredictionErrorDecomposition pe = new PredictionErrorDecomposition(false);
        filter.process(ssf, new SsfData(M.column(0)), pe);
        ILikelihood ll = pe.likelihood();
        assertTrue(Math.abs(Err - ll.getSsqErr()) / Err < 1e-9);
        assertTrue(Math.abs(Det - ll.getLogDeterminant()) < 1e-9);
    }

    @Test
    public void testNew2Method() {
        MultivariateOrdinaryFilter filter = new MultivariateOrdinaryFilter();
        SsfArima ssf = SsfArima.create(model);
        PredictionErrorsDecomposition pe = new PredictionErrorsDecomposition(false);
        filter.process(MultivariateTimeInvariantSsf.of(ssf), new SsfMatrix(M), pe);
        ILikelihood ll = pe.likelihood();
        assertTrue(Math.abs(Err - ll.getSsqErr()) / Err < 1e-9);
        assertTrue(Math.abs(Det - ll.getLogDeterminant()) < 1e-9);
    }

    @Test
    public void testNew3Method() {
        ArrayFilter filter = new ArrayFilter();
        SsfArima ssf = SsfArima.create(model);
        PredictionErrorDecomposition pe = new PredictionErrorDecomposition(false);
        filter.process(ssf, new SsfData(M.column(0)), pe);
        ILikelihood ll = pe.likelihood();
        assertTrue(Math.abs(Err - ll.getSsqErr()) / Err < 1e-9);
        assertTrue(Math.abs(Det - ll.getLogDeterminant()) < 1e-9);
    }

    @Test
    public void testNew4Method() {
        MultivariateArrayFilter filter = new MultivariateArrayFilter();
        SsfArima ssf = SsfArima.create(model);
        PredictionErrorsDecomposition pe = new PredictionErrorsDecomposition(false);
        filter.process(MultivariateTimeInvariantSsf.of(ssf), new SsfMatrix(M), pe);
        ILikelihood ll = pe.likelihood();
        assertTrue(Math.abs(Err - ll.getSsqErr()) / Err < 1e-9);
        assertTrue(Math.abs(Det - ll.getLogDeterminant()) < 1e-9);
    }

    @Test
    public void testNew5Method() {
        SsfArima ssf = SsfArima.create(model);
        FastFilter filter = new FastFilter(SsfArima.fastInitializer(ssf));
        PredictionErrorDecomposition pe = new PredictionErrorDecomposition(false);
        filter.process(ssf, new SsfData(M.column(0)), pe);
        ILikelihood ll = pe.likelihood();
        assertTrue(Math.abs(Err - ll.getSsqErr()) < 1e-9);
        assertTrue(Math.abs(Det - ll.getLogDeterminant()) < 1e-9);
    }

    @Ignore
    @Test
    public void stressTest1() {
        for (int q = 0; q < 2; ++q) {
            int N = q == 0 ? 1000 : 50000;

            long t0 = System.currentTimeMillis();
            for (int i = 0; i < N; ++i) {
                ec.tstoolkit.ssf.Filter filter = new ec.tstoolkit.ssf.Filter();
                ec.tstoolkit.ssf.arima.SsfArma ssf
                        = new ec.tstoolkit.ssf.arima.SsfArma(model);
                ec.tstoolkit.ssf.PredictionErrorDecomposition pe = new ec.tstoolkit.ssf.PredictionErrorDecomposition(false);
                filter.setSsf(ssf);
                filter.process(new ec.tstoolkit.ssf.SsfData(M.column(0), null), pe);
            }
            long t1 = System.currentTimeMillis();
            System.out.println("old ordinary filter");
            System.out.println(t1 - t0);
            t0 = System.currentTimeMillis();
            for (int i = 0; i < N; ++i) {
                OrdinaryFilter filter = new OrdinaryFilter();
                SsfArima ssf = SsfArima.create(model);
                PredictionErrorDecomposition pe = new PredictionErrorDecomposition(false);
//            DefaultFilteringResults pe=DefaultFilteringResults.light();
                filter.process(ssf, new SsfData(M.column(0)), pe);
            }
            t1 = System.currentTimeMillis();
            System.out.println("new ordinary filter");
            System.out.println(t1 - t0);
            t0 = System.currentTimeMillis();
            for (int i = 0; i < N; ++i) {
                ArrayFilter filter = new ArrayFilter();
                SsfArima ssf = SsfArima.create(model);
                PredictionErrorDecomposition pe = new PredictionErrorDecomposition(false);
                filter.process(ssf, new SsfData(M.column(0)), pe);
            }
            t1 = System.currentTimeMillis();
            System.out.println("new array filter");
            System.out.println(t1 - t0);
            t0 = System.currentTimeMillis();
            for (int i = 0; i < N; ++i) {
                MultivariateArrayFilter filter = new MultivariateArrayFilter();
                SsfArima ssf = SsfArima.create(model);
                PredictionErrorsDecomposition pe = new PredictionErrorsDecomposition(false);
                filter.process(MultivariateSsf.proxy(ssf), new SsfMatrix(M), pe);
            }
            t1 = System.currentTimeMillis();
            System.out.println("new multivariate array filter");
            System.out.println(t1 - t0);
            t0 = System.currentTimeMillis();
            for (int i = 0; i < N; ++i) {
                SarimaModel arima = new SarimaModel(model.getSpecification());
                SsfArima ssf = SsfArima.create(arima);
                FastFilter filter = new FastFilter(SsfArima.fastInitializer(ssf));
                PredictionErrorDecomposition pe = new PredictionErrorDecomposition(false);
                filter.process(ssf, new SsfData(M.column(0)), pe);
            }
            t1 = System.currentTimeMillis();
            System.out.println("new fast filter");
            System.out.println(t1 - t0);
            t0 = System.currentTimeMillis();
            for (int i = 0; i < N; ++i) {
                SarimaModel arima = new SarimaModel(model.getSpecification());
                ArmaKF kf = new ArmaKF(arima);
                Likelihood ll = new Likelihood();
                kf.process(M.column(0), ll);
            }
            t1 = System.currentTimeMillis();
            System.out.println("specific fast filter");
            System.out.println(t1 - t0);
            t0 = System.currentTimeMillis();
            for (int i = 0; i < N; ++i) {
                OrdinaryFilter filter = new OrdinaryFilter();
                SsfArima ssf = SsfArima.create(model);
                PredictionErrorDecomposition pe = new PredictionErrorDecomposition(false);
                filter.process(TimeInvariantSsf.of(ssf), new SsfData(M.column(0)), pe);
            }
            t1 = System.currentTimeMillis();
            System.out.println("new time invariant ordinary filter");
            System.out.println(t1 - t0);
            t0 = System.currentTimeMillis();
            for (int i = 0; i < N; ++i) {
                MultivariateOrdinaryFilter filter = new MultivariateOrdinaryFilter();
                SsfArima ssf = SsfArima.create(model);
                PredictionErrorsDecomposition pe = new PredictionErrorsDecomposition(false);
                filter.process(MultivariateTimeInvariantSsf.of(ssf), new SsfMatrix(M), pe);
            }
            t1 = System.currentTimeMillis();
            System.out.println("new time invariant multivariate ordinary filter");
            System.out.println(t1 - t0);
        }
    }

}
