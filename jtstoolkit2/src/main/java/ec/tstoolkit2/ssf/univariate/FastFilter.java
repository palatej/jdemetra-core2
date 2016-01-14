/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.univariate;

import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit2.ssf.ISsfDynamics;

/**
 *
 * @author Jean Palate
 */
public class FastFilter {
    
    private final IFilteringResults frslts;
    private final ISsfMeasurement measurement;
    private final ISsfDynamics dynamics;
    
    public FastFilter(ISsf ssf, IFilteringResults frslts){
        this.frslts=frslts;
        measurement=ssf.getMeasurement();
        dynamics=ssf.getDynamics();
    }
    
    public void process(Matrix x0, Matrix x){
        
    }
    
    
}
