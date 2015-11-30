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
package ec.tstoolkit2.ssf.implementations;

import ec.tstoolkit2.ssf.implementations.arima.SsfArima;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.sarima.SarimaModel;
import ec.tstoolkit.sarima.SarimaModelBuilder;
import ec.tstoolkit2.ssf.StateInfo;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Jean Palate
 */
public class SsfArimaTest {

    static final SarimaModel model;
    static final SsfArima ssfnew;
    static final ec.tstoolkit.ssf.arima.SsfArima ssfold;

    static {
        SarimaModelBuilder builder = new SarimaModelBuilder();
        model = builder.createAirlineModel(12, -.6, -.8);
        ssfnew = SsfArima.create(model);
        ssfold = new ec.tstoolkit.ssf.arima.SsfArima(model);
    }

    public SsfArimaTest() {
    }

    @Test
    public void testP0() {
        int n=ssfnew.getStateDim();
        int nd=ssfnew.getDynamics().getNonStationaryDim();
        Matrix P00=Matrix.square(n);
        Matrix P01=Matrix.square(n);
        ssfnew.getDynamics().Pf0(P00.subMatrix(), StateInfo.Forecast);
        ssfold.Pf0(P01.subMatrix());
        assertTrue(P00.equals(P01, 1e-9));
    }

    @Test
    public void testDiffuse() {
        int n=ssfnew.getStateDim();
        int nd=ssfnew.getDynamics().getNonStationaryDim();
        Matrix B0=new Matrix(n, nd);
        Matrix B1=new Matrix(n, nd);
        ssfnew.getDynamics().diffuseConstraints(B0.subMatrix());
        ssfold.diffuseConstraints(B1.subMatrix());
        assertTrue(B0.equals(B1, 1e-9));
    }
}
