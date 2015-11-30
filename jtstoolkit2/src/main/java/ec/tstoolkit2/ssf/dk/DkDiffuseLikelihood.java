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
package ec.tstoolkit2.ssf.dk;

import ec.tstoolkit.eco.ILikelihood;
import ec.tstoolkit.utilities.Arrays2;

/**
 * The diffuse likelihood follows the definition provided in the paper:
 * "Likelihood functions for state space models with diffuse initial conditions"
 * Francke, Koopman, de Vos Journal of Time Series Analysis, July 2010. This
 * definition is slightly different in comparison with JD+ 2.0
 *
 * @author Jean Palate
 */
public class DkDiffuseLikelihood implements ILikelihood {

    /**
     * Respectively: diffuse log-likelihood sum of the squared residuals log
     * determinant of the cov matrix diffuse correction
     */
    private double ll, ssqerr, ldet, lddet;
    private int nobs, nd;
    private double[] res;
    private boolean legacy;

    /**
     *
     */
    public DkDiffuseLikelihood() {
    }

    private int m(){
        return legacy ? nobs : nobs-nd;
    }
    
    /**
     * 
     * @return false by default
     */
    public boolean isLegacy(){
        return legacy;
    }
    
    /**
     * 
     * @param legacy legacy=true should be used only for testing purposes
     */
    public void setLegacy(boolean legacy){
        if (this.legacy != legacy){
            this.legacy=legacy;
            calcll();
        }
    }
    /**
     *
     * @param nparams
     * @return
     */
    public double AIC(final int nparams) {
        return -2 * ll + 2 * nparams;
    }

    /**
     *
     * @param nparams
     * @return
     */
    public double BIC(final int nparams) {
        return -2 * ll + nparams * Math.log(m());
    }

    /**
     *
     */
    public void clear() {
        ll = 0;
        ssqerr = 0;
        ldet = 0;
        lddet = 0;
        nobs = 0;
        nd = 0;
    }

    /**
     *
     * @return
     */
    public int getD() {
        return nd;
    }

    @Override
    public double getFactor() {
        return Math.exp((ldet + lddet) / (m()));
    }

    @Override
    public double getLogLikelihood() {
        return ll;
    }

    @Override
    public int getN() {
        return nobs;
    }

    @Override
    public double[] getResiduals() {
        return res;
    }

    @Override
    public double getLogDeterminant() {
        return ldet;
    }

    /**
     *
     * @return
     */
    public double getSer() {
        return Math.sqrt(ssqerr / (m()));
    }

    @Override
    public double getSigma() {
        return ssqerr / (m());
    }

    @Override
    public double getSsqErr() {
        return ssqerr;
    }

    /**
     * Adjust the likelihood if the data have been pre-multiplied by a given
     * scaling factor
     *
     * @param factor The scaling factor
     */
    public void rescale(final double factor) {
        if (factor == 1) {
            return;
        }
        ssqerr /= factor * factor;
        ll += (m()) * Math.log(factor);
        if (res != null) {
            for (int i = 0; i < res.length; ++i) {
                if (!Double.isNaN(res[i])) {
                    res[i] /= factor;
                }
            }
        }
    }

    public double getDiffuseCorrection() {
        return lddet;
    }

    /**
     * Initialize the diffuse likelihood. We consider below the GLS problem
     * corresponding to a given state space: y = a * X + e, where X is derived
     * from the initial conditions and e ~ N(0, V)
     *
     * The diffuse likelihood is then:
     *
     * -0.5*(m*log(2*pi)+m*log(ssqerr/m)+m+log|V|+log|X'V^-1*X| where m=n-d
     *
     * It should be noted that the usual definition (implemented in JD+ 2.0) is
     * -0.5*(n*log(2*pi)+n*log(ssqerr/n)+n+log|V|+log|X'V^-1*X| The difference
     * is thus -0.5*(d*log(2*pi)+d*log(ssqerr)-n*log(n)+m*log(m))
     *
     * The new definition is more coherent with the marginal likelihood.
     *
     * @param ssqerr The sum of the squared residuals
     * @param ldet The log of the determinant of V
     * @param lddet Diffuse correction
     * @param n The number of observations
     * @param d The number of diffuse constraints
     * @return
     */
    public boolean set(final double ssqerr, final double ldet, final double lddet,
            final int n, final int d) {
        if (d == 0 && lddet != 0) {
            return false;
        }
        this.ssqerr = ssqerr;
        this.ldet = ldet;
        this.lddet = lddet;
        this.nobs = n;
        nd = d;
        calcll();
        return true;
    }

    private void calcll() {
        int m = m();
        if (m <= 0)
            return;
        ll = -.5
                * (m * Math.log(2 * Math.PI) + m
                * (1 + Math.log(ssqerr / m)) + ldet + lddet);

    }

    public void add(ILikelihood ll) {
        nobs += ll.getN();
        ssqerr += ll.getSsqErr();
        ldet += ll.getLogDeterminant();
        if (res != null) {
            double[] nres = ll.getResiduals();
            if (nres != null) {
                res = Arrays2.concat(res, nres);
            }
        }
        calcll();
    }

    void setResiduals(double[] res) {
        this.res = res;
    }

}
