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

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit2.ssf.univariate.ISsfMeasurement;

/**
 * Shifted measurement
 *
 * @author Jean Palate
 */
public class ShiftedMeasurement implements ISsfMeasurement {

    private final ISsfMeasurement m;
    private final int shift;
    
    public ShiftedMeasurement(ISsfMeasurement m, int shift) {
        this.m = m;
        this.shift = shift;
    }
    
    @Override
    public boolean isTimeInvariant() {
        return m.isTimeInvariant();
    }
    
    @Override
    public void Z(int pos, DataBlock z) {
        m.Z(pos + shift, z);
    }
    
    @Override
    public boolean hasErrors() {
        return m.hasErrors();
    }
    
    @Override
    public boolean hasError(int pos) {
        return m.hasError(pos + shift);
    }
    
    @Override
    public double errorVariance(int pos) {
        return m.errorVariance(pos + shift);
    }
    
    @Override
    public double ZX(int pos, DataBlock b) {
        return m.ZX(pos + shift, b);
    }
    
    @Override
    public double ZVZ(int pos, SubMatrix V) {
        return m.ZVZ(pos + shift, V);
    }
    
    @Override
    public void VpZdZ(int pos, SubMatrix V, double d) {
        m.VpZdZ(pos + shift, V, d);
    }
    
    @Override
    public void XpZd(int pos, DataBlock x, double d) {
        m.XpZd(pos + shift, x, d);
    }
}
