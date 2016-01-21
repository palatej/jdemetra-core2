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

import data.Models;
import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit2.ssf.ResultsRange;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.SsfData;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author Jean Palate
 */
public class FastDiffuseFilterTest {
    

    public FastDiffuseFilterTest() {
    }

    @Test
    @Ignore
    public void stressTestUcarima() {
        ISsf ssf=Models.ssfUcarima;
        SsfData data=Models.ssfXRandom;
        DefaultDiffuseFilteringResults fresults = DkToolkit.filter(ssf, data, false);
        Matrix x = new Matrix(data.getLength(), 1);
        x.randomize(0);
        long t0 = System.currentTimeMillis();
        for (int i = 0; i < 100000; ++i) {
            FastDiffuseFilter filter = new FastDiffuseFilter(ssf, fresults, new ResultsRange(fresults.getEndDiffusePosition(), data.getLength()));
            filter.filter(x.column(0));
        }
        long t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
    }

    @Test
    public void testUcarima() {
        ISsf ssf=Models.ssfUcarima;
        SsfData data=Models.ssfXRandom;
        DefaultDiffuseFilteringResults fresults = DkToolkit.filter(ssf, data, false);
        Matrix x = new Matrix(data.getLength(), 10);
        x.randomize();
        x.column(0).copy(data);
        FastDiffuseFilter filter = new FastDiffuseFilter(ssf, fresults, new ResultsRange(0, data.getLength()));
        filter.filter(x.subMatrix());
        assertTrue(new DataBlock(fresults.errors()).distance(x.column(0))<1e-9);
    }
}
