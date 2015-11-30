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
package ec.tstoolkit2.ssf.implementations.dfm;

import ec.tstoolkit.eco.ILikelihood;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit2.ssf.array.MultivariateArrayFilter;
import ec.tstoolkit2.ssf.implementations.var.VarDescriptor;
import ec.tstoolkit2.ssf.multivariate.M2uAdapter;
import ec.tstoolkit2.ssf.multivariate.MultivariateOrdinaryFilter;
import ec.tstoolkit2.ssf.multivariate.PredictionErrorsDecomposition;
import ec.tstoolkit2.ssf.multivariate.SsfMatrix;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfData;
import ec.tstoolkit2.ssf.univariate.OrdinaryFilter;
import ec.tstoolkit2.ssf.univariate.PredictionErrorDecomposition;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author Jean Palate
 */
public class SsfDfmTest {

    private static final SsfDfm dfm;
    private static final Matrix data;

    static {
        VarDescriptor vdesc = new VarDescriptor(3, 1);
        vdesc.setDefault();

        MeasurementDescriptor[] mdesc = new MeasurementDescriptor[40];
        for (int i = 0; i < mdesc.length / 2; ++i) {
            mdesc[i] = new MeasurementDescriptor(LevelMeasurement.ML, new double[3], 1);
            mdesc[i].seDefaultCoefficients();
        }
        for (int i = mdesc.length / 2; i < mdesc.length; ++i) {
            mdesc[i] = new MeasurementDescriptor(CumulMeasurement.MC12, new double[3], 1);
            mdesc[i].seDefaultCoefficients();
        }
        dfm = SsfDfm.from(vdesc, mdesc);

        data = new Matrix(240, mdesc.length);
        data.randomize(1);
    }

    public SsfDfmTest() {
    }

    @Test
    public void testLL() {
        MultivariateOrdinaryFilter mfilter = new MultivariateOrdinaryFilter();
        PredictionErrorsDecomposition decomp = new PredictionErrorsDecomposition(false);
        SsfMatrix ssfdata = new SsfMatrix(data);
        mfilter.process(dfm, ssfdata, decomp);
        MultivariateArrayFilter afilter = new MultivariateArrayFilter();
        PredictionErrorsDecomposition adecomp = new PredictionErrorsDecomposition(false);
        afilter.process(dfm, ssfdata, adecomp);
        assertEquals(decomp.likelihood().getLogLikelihood(), adecomp.likelihood().getLogLikelihood(), 1e-6);
        PredictionErrorDecomposition udecomp = new PredictionErrorDecomposition(false);
        ISsf udfm = M2uAdapter.of(dfm);
        ISsfData udata = M2uAdapter.of(ssfdata);
        OrdinaryFilter filter = new OrdinaryFilter();
        filter.process(udfm, udata, udecomp);
        assertEquals(decomp.likelihood().getLogLikelihood(), udecomp.likelihood().getLogLikelihood(), 1e-6);
    }

    int M = 200;

    @Test
    @Ignore
    public void stressTestLL() {
        VarDescriptor vdesc = new VarDescriptor(3, 1);
        vdesc.setDefault();

        MeasurementDescriptor[] mdesc = new MeasurementDescriptor[40];
        for (int i = 0; i < mdesc.length / 2; ++i) {
            mdesc[i] = new MeasurementDescriptor(LevelMeasurement.ML, new double[3], 1);
            mdesc[i].seDefaultCoefficients();
        }
        for (int i = mdesc.length / 2; i < mdesc.length; ++i) {
            mdesc[i] = new MeasurementDescriptor(CumulMeasurement.MC12, new double[3], 1);
            mdesc[i].seDefaultCoefficients();
        }
        long t0 = System.currentTimeMillis();
        ILikelihood ll =null;
        for (int i = 0; i < M; ++i) {
            SsfDfm ssf = SsfDfm.from(vdesc, mdesc);
            MultivariateOrdinaryFilter mfilter = new MultivariateOrdinaryFilter();
            PredictionErrorsDecomposition decomp = new PredictionErrorsDecomposition(false);
            SsfMatrix ssfdata = new SsfMatrix(data);
            mfilter.process(ssf, ssfdata, decomp);
            ll=decomp.likelihood();
        }
        long t1 = System.currentTimeMillis();
        System.out.println("Multivariate");
        System.out.println(t1 - t0);
        System.out.println(ll.getLogLikelihood());
        t0 = System.currentTimeMillis();
        for (int i = 0; i < M; ++i) {
            SsfDfm ssf = SsfDfm.from(vdesc, mdesc);
            PredictionErrorDecomposition decomp = new PredictionErrorDecomposition(false);
            SsfMatrix ssfdata = new SsfMatrix(data);
            ISsf udfm = M2uAdapter.of(ssf);
            ISsfData udata = M2uAdapter.of(ssfdata);
            OrdinaryFilter filter = new OrdinaryFilter();
            filter.process(udfm, udata, decomp);
            ll = decomp.likelihood();
        }
        t1 = System.currentTimeMillis();
        System.out.println("Univariate");
        System.out.println(t1 - t0);
        System.out.println(ll.getLogLikelihood());
    }

}
