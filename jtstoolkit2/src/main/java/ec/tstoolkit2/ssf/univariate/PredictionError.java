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

import ec.tstoolkit.data.DataBlock;

/**
 *
 * @author Jean Palate
 */
public class PredictionError {

    /**
     * e is the prediction error (=y(t)-Z(t)A(t))
     */
    private double e, f,stde;

    /**
     * C = P Z'
     */
    private final DataBlock M;
    

    /**
     *
     * @param dim
     */
    public PredictionError(final int dim) {
        M = new DataBlock(dim);
    }

    /**
     * 
     * @return 
     */
    public double get() {
        return e;
    }

    /**
     * 
     * @return 
     */
    public double getVariance() {
        return f;
    }

    /**
     * 
     * @return 
     */
    public double getStandardDeviation() {
        return stde;
    }

    /**
     * 
     * @return 
     */
     public DataBlock M() {
        return M;
    }

    /**
     * =(ZPZ'+v) variance of the
     * prediction error 
     * @param f 
     */
    public void setVariance(final double f){
        this.f=f;
        this.stde=Math.sqrt(f);
    };
    
    public void setStandardDeviation(final double e){
        this.f=e*e;
        this.stde=e;
    };
    
    public void set(final double val){
        this.e=val;
    }
    /**
     * 
     * @return
     */
    public boolean isMissing()
    {
	return Double.isNaN(e);
    }
    
    public void setMissing(){
        e=Double.NaN;
    }
}
