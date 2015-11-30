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
package ec.tstoolkit2.ssf.implementations.arima;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SymmetricMatrix;
import ec.tstoolkit.sarima.SarimaModel;
import ec.tstoolkit.sarima.SarimaSpecification;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.implementations.TimeInvariantSsf;
import ec.tstoolkit2.ssf.univariate.ISsf;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Jean Palate
 */
public class SsfArimaTest {

    SarimaModel model;

    public SsfArimaTest() {
        SarimaSpecification spec = new SarimaSpecification(12);
        spec.airline();
        spec.setP(3);
        spec.setBP(1);
        model=new SarimaModel(spec);
        model.setDefault(-.2, -.4);
    }

    @Test
    public void testDynamics() {
        SsfArima arima = SsfArima.create(model);
        ISsf ref = TimeInvariantSsf.of(arima);
        int dim = arima.getStateDim();
        Matrix M1 = new Matrix(dim, dim);
        M1.randomize();
        M1 = SymmetricMatrix.XXt(M1);
        Matrix M2 = M1.clone();
        DataBlock x1 = new DataBlock(dim);
        x1.randomize();
        DataBlock x2 = x1.deepClone();
        ISsfDynamics dynref = ref.getDynamics();
        ISsfDynamics dyn = arima.getDynamics();
        dynref.TX(0, x1);
        dyn.TX(dim, x2);
        assertTrue(x1.distance(x2) < 1e-9);
        dynref.XT(0, x1);
        dyn.XT(0, x2);
        assertTrue(x1.distance(x2) < 1e-9);
        dynref.TVT(0, M1.subMatrix());
        dyn.TVT(0, M2.subMatrix());
        assertTrue(M1.distance(M2) < 1e-9);
        dynref.TM(0, M1.subMatrix());
        dyn.TM(0, M2.subMatrix());
        assertTrue(M1.distance(M2) < 1e-9);
        dynref.MT(0, M1.subMatrix());
        dyn.MT(0, M2.subMatrix());
        assertTrue(M1.distance(M2) < 1e-9);
    }

    @Test
    public void testMeasurement() {
        SsfArima arima = SsfArima.create(model);
        ISsf ref = TimeInvariantSsf.of(arima);
        int dim = arima.getStateDim();
        Matrix M1 = new Matrix(dim, dim);
        M1.randomize();
        M1 = SymmetricMatrix.XXt(M1);
        Matrix M2 = M1.clone();
        DataBlock x1 = new DataBlock(dim);
        x1.randomize();
        DataBlock x2 = x1.deepClone();
        ISsfMeasurement mref = ref.getMeasurement();
        ISsfMeasurement m = arima.getMeasurement();
        assertTrue(Math.abs(mref.ZX(0, x1) - m.ZX(0, x2)) < 1e-9);
        assertTrue(Math.abs(mref.ZVZ(0, M1.subMatrix()) - m.ZVZ(0, M1.subMatrix())) < 1e-9);
        mref.VpZdZ(0, M1.subMatrix(), 5);
        m.VpZdZ(0, M2.subMatrix(), 5);
        assertTrue(M1.distance(M2) < 1e-9);
        mref.XpZd(0, x1, 5);
        m.XpZd(0, x2, 5);
        assertTrue(x1.distance(x2) < 1e-9);
        mref.ZM(dim, M1.subMatrix(), x1);
        m.ZM(dim, M1.subMatrix(), x2);
        assertTrue(x1.distance(x2) < 1e-9);
    }

}
