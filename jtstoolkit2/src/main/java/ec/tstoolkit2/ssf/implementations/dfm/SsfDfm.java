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

import ec.tstoolkit2.ssf.implementations.var.VarDescriptor;
import ec.tstoolkit2.ssf.implementations.var.VarDynamics;
import ec.tstoolkit2.ssf.multivariate.MultivariateSsf;

/**
 *
 * @author Jean Palate
 */
public class SsfDfm extends MultivariateSsf {

    public static SsfDfm from(VarDescriptor vdesc, MeasurementDescriptor[] mdesc) {
        int nlx = vdesc.getLagsCount();
        for (int i=0; i<mdesc.length; ++i){
            int n=mdesc[i].getType().getLength();
            if (nlx<n)
                nlx=n;
        }
        return from(vdesc, mdesc, nlx);
    }

    public static SsfDfm from(VarDescriptor vdesc, MeasurementDescriptor[] mdesc, int nlx) {
        int nf = vdesc.getEquationsCount();
        VarDynamics dyn = VarDynamics.from(vdesc, nlx);
        DfmMeasurements m = DfmMeasurements.from(nf, nlx, mdesc);
        return new SsfDfm(dyn, m);
    }

    private SsfDfm(VarDynamics dyn, DfmMeasurements m) {
        super(dyn, m);
    }

    @Override
    public DfmMeasurements getMeasurements() {
        return (DfmMeasurements) super.getMeasurements();
    }

    @Override
    public VarDynamics getDynamics() {
        return (VarDynamics) super.getDynamics();
    }

}
