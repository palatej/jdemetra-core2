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

import data.Models;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit2.ssf.ResultsRange;
import ec.tstoolkit2.ssf.dk.DefaultDiffuseFilteringResults;
import ec.tstoolkit2.ssf.dk.DkToolkit;
import static org.junit.Assert.assertTrue;
import org.junit.Ignore;
import org.junit.Test;

/**
 *
 * @author Jean Palate
 */
public class FastFilterTest {

    final static int N = 100000;

    public FastFilterTest() {
    }

    @Test
    @Ignore
    public void stressTestUcarima() {
        ISsf ssf=Models.ssfUcarima;
        SsfData data=Models.ssfXRandom;
        DefaultDiffuseFilteringResults fresults = DkToolkit.filter(ssf, data, false);
        Matrix x = new Matrix(data.getLength(), 10);
        x.randomize(0);
        long t0 = System.currentTimeMillis();
        for (int i = 0; i < 10000; ++i) {
            FastFilter filter = new FastFilter(ssf, fresults, new ResultsRange(fresults.getEndDiffusePosition(), data.getLength()));
            filter.filter(x.subMatrix());
        }
        long t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
    }

    @Test
    public void testUcarima() {
        ISsf ssf=Models.ssfUcarima;
        SsfData data=Models.ssfXRandom;
        DefaultDiffuseFilteringResults fresults = DkToolkit.filter(ssf, data, false);
        Matrix x = new Matrix(data.getLength(), 1);
        x.column(0).copy(data);
        Matrix x0=new Matrix(ssf.getStateDim(), 1);
        int nd=fresults.getEndDiffusePosition();
        x0.column(0).copy(fresults.a(nd));
        FastFilter filter = new FastFilter(ssf, fresults, new ResultsRange(nd, data.getLength()));
        filter.filter(x0.subMatrix(), x.subMatrix(nd,-1,0,1));
        assertTrue(new DataBlock(fresults.errors()).drop(nd, 0)
                .distance(x.column(0).drop(nd, 0))<1e-8);
    }
}
