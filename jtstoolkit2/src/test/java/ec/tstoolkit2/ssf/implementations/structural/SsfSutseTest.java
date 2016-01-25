/*
 * Copyright 2015 National Bank of Belgium
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
/*
 */
package ec.tstoolkit2.ssf.implementations.structural;

import ec.tstoolkit.eco.ILikelihood;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit.timeseries.simplets.TsData;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.akf.AkfDiffuseLikelihood;
import ec.tstoolkit2.ssf.akf.AkfToolkit;
import ec.tstoolkit2.ssf.akf.AugmentedPredictionErrorsDecomposition;
import ec.tstoolkit2.ssf.akf.MultivariateAugmentedFilter;
import ec.tstoolkit2.ssf.akf.MultivariateAugmentedFilterInitializer;
import ec.tstoolkit2.ssf.dk.DkToolkit;
import ec.tstoolkit2.ssf.multivariate.IMultivariateSsf;
import ec.tstoolkit2.ssf.multivariate.M2uAdapter;
import ec.tstoolkit2.ssf.multivariate.MultivariateOrdinaryFilter;
import ec.tstoolkit2.ssf.multivariate.PredictionErrorsDecomposition;
import ec.tstoolkit2.ssf.multivariate.SsfMatrix;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfData;
import ec.tstoolkit2.ssf.univariate.SsfData;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author Jean Palate
 */
public class SsfSutseTest {

    public static final BasicStructuralModel m1, m2, m3;
    public static final SutseModel sutse;

    public static final int N = 1000;

    static {
        // 3 usual bsm
        ModelSpecification spec = new ModelSpecification();
        spec.useLevel(ComponentUse.Free);
        spec.useSlope(ComponentUse.Free);
        spec.useCycle(ComponentUse.Free);
        spec.useNoise(ComponentUse.Free);
        spec.setSeasonalModel(SeasonalModel.Dummy);

        m1 = new BasicStructuralModel(spec, 12);
        m1.setCycle(.9, 8);
        m1.setVariance(Component.Level, .1);
        m1.setVariance(Component.Slope, .2);
        m1.setVariance(Component.Cycle, .5);
        m1.setVariance(Component.Seasonal, 2);
        m1.setVariance(Component.Noise, 1);

        //spec.useCycle(ComponentUse.Unused);
        m2 = new BasicStructuralModel(spec, 12);
        m2.setVariance(Component.Level, .15);
        m2.setVariance(Component.Slope, .25);
        m2.setVariance(Component.Cycle, .75);
        m2.setCycle(.5, 6);
        m2.setVariance(Component.Seasonal, 2.5);
        m2.setVariance(Component.Noise, 1.5);

        spec.useCycle(ComponentUse.Unused);
        m3 = new BasicStructuralModel(spec, 12);
        m3.setVariance(Component.Level, .15);
        m3.setVariance(Component.Slope, .25);
        m3.setVariance(Component.Seasonal, 2.5);
        m3.setVariance(Component.Noise, 1.5);

        sutse = new SutseModel(new BasicStructuralModel[]{m1, m2, m3});
        Matrix corr = Matrix.square(3);
        corr.diagonal().set(1);
        corr.set(0, 1, .8);
        corr.set(1, 0, .8);
        corr.set(0, 2, .5);
        corr.set(2, 0, .5);
        corr.set(1, 2, .4);
        corr.set(2, 1, .4);
        sutse.setCorrelations(Component.Noise, corr);
        Matrix ccorr = Matrix.square(2);
        ccorr.diagonal().set(1);
        ccorr.set(0, 1, .8);
        ccorr.set(1, 0, .8);
        sutse.setCorrelations(Component.Cycle, ccorr);
    }

    public SsfSutseTest() {
    }

    @Test
    public void testMultivariate() {
        IMultivariateSsf mssf = SsfSutse.of(sutse);
        assertTrue(mssf != null);
        // gets the innovations...
        int rdim = mssf.getDynamics().getInnovationsDim();
        int dim = mssf.getDynamics().getStateDim();
        Matrix V = Matrix.square(dim);
        mssf.getDynamics().V(0, V.subMatrix());
        Matrix S = new Matrix(dim, rdim);
        mssf.getDynamics().S(0, S.subMatrix());
        Matrix W = SymmetricMatrix.XXt(S);
        assertTrue(W.minus(V).nrm2() < 1e-6);
    }

    @Test
    public void testLikelihood() {

        TsData s1 = data.Data.M1;
        TsData s2 = data.Data.M2;
        TsData s3 = data.Data.M3;
        Matrix m = new Matrix(s1.getLength(), 3);
        m.column(0).copy(s1);
        m.column(1).copy(s2);
        m.column(2).copy(s3);
        m.set(2, 1, Double.NaN);
        m.set(3, 0, Double.NaN);
        m.set(5, 1, Double.NaN);
        m.set(5, 0, Double.NaN);
        //SsfMatrix data = new SsfMatrix(m.subMatrix(1, s1.getLength(), 0, 2));
        SsfMatrix data = new SsfMatrix(m);
        long t0 = System.currentTimeMillis();
        AugmentedPredictionErrorsDecomposition idecomp = new AugmentedPredictionErrorsDecomposition();
        MultivariateOrdinaryFilter filter = new MultivariateOrdinaryFilter(new MultivariateAugmentedFilterInitializer(idecomp));
        PredictionErrorsDecomposition decomp = new PredictionErrorsDecomposition(false);
        filter.process(SsfSutse.of(sutse), data, decomp);
        AkfDiffuseLikelihood likelihood = idecomp.likelihood();
        likelihood.add(decomp.likelihood());
        AugmentedPredictionErrorsDecomposition idecomp2 = new AugmentedPredictionErrorsDecomposition();
        MultivariateAugmentedFilter akf = new MultivariateAugmentedFilter(false);
        akf.process(SsfSutse.of(sutse), data, idecomp2);
        AkfDiffuseLikelihood likelihood1 = idecomp2.likelihood();
        assertEquals(likelihood.getLogLikelihood(), likelihood1.getLogLikelihood(), 1e-8);

        ISsfData udata = M2uAdapter.of(data);
        ISsf ussf = M2uAdapter.of(SsfSutse.of(sutse));
        ILikelihood llm2u = DkToolkit.likelihoodComputer().compute(ussf, udata);
        assertEquals(likelihood.getLogLikelihood(), llm2u.getLogLikelihood(), 1e-8);
    }

    @Test
    @Ignore
    public void stressTestLikelihood() {

        TsData s1 = data.Data.M1;
        TsData s2 = data.Data.M2;
        TsData s3 = data.Data.M3;
        Matrix m = new Matrix(s1.getLength(), 3);
        m.column(0).copy(s1);
        m.column(1).copy(s2);
        m.column(2).copy(s3);
        m.set(2, 1, Double.NaN);
        m.set(3, 0, Double.NaN);
        m.set(5, 1, Double.NaN);
        m.set(5, 0, Double.NaN);
        //SsfMatrix data = new SsfMatrix(m.subMatrix(1, s1.getLength(), 0, 2));
        SsfMatrix data = new SsfMatrix(m);
        long t0 = System.currentTimeMillis();
        AkfDiffuseLikelihood likelihood = null;
        for (int i = 0; i < N; ++i) {
            AugmentedPredictionErrorsDecomposition idecomp = new AugmentedPredictionErrorsDecomposition();
            MultivariateOrdinaryFilter filter = new MultivariateOrdinaryFilter(new MultivariateAugmentedFilterInitializer(idecomp));
            PredictionErrorsDecomposition decomp = new PredictionErrorsDecomposition(false);
            filter.process(SsfSutse.of(sutse), data, decomp);
            likelihood = idecomp.likelihood();
            likelihood.add(decomp.likelihood());
        }
        long t1 = System.currentTimeMillis();
        System.out.println("akf (collapsing)");
        System.out.println(t1 - t0);
        System.out.println(likelihood.getLogLikelihood());
        t0 = System.currentTimeMillis();
        AkfDiffuseLikelihood likelihood1 = null;
        for (int i = 0; i < N; ++i) {
            AugmentedPredictionErrorsDecomposition idecomp2 = new AugmentedPredictionErrorsDecomposition();
            MultivariateAugmentedFilter akf = new MultivariateAugmentedFilter(false);
            akf.process(SsfSutse.of(sutse), data, idecomp2);
            likelihood1 = idecomp2.likelihood();
        }
        t1 = System.currentTimeMillis();
        System.out.println("akf (no collapsing)");
        System.out.println(t1 - t0);
        System.out.println(likelihood1.getLogLikelihood());

//        ILikelihood ll1 = AkfToolkit.likelihoodComputer(true).compute(SsfBsm.create(m1), new SsfData(m.column(0)));
//        System.out.println(ll1.getLogLikelihood());
//        ILikelihood ll2 = AkfToolkit.likelihoodComputer(true).compute(SsfBsm.create(m2), new SsfData(m.column(1)));
//        System.out.println(ll2.getLogLikelihood());
        ILikelihood llm2u = null;
        t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            ISsfData udata = M2uAdapter.of(data);
            ISsf ussf = M2uAdapter.of(SsfSutse.of(sutse));
            llm2u = DkToolkit.likelihoodComputer().compute(ussf, udata);
        }
        t1 = System.currentTimeMillis();
        System.out.println("univariate");
        System.out.println(t1 - t0);
        System.out.println(llm2u.getLogLikelihood());
    }

}
