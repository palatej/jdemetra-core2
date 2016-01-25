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
package ec.tstoolkit2.ssf.implementations.structural;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.maths.matrices.Matrix;
import ec.tstoolkit.maths.matrices.SubMatrix;
import ec.tstoolkit.utilities.DoubleList;
import ec.tstoolkit.utilities.IntList;
import ec.tstoolkit2.ssf.ISsfDynamics;
import ec.tstoolkit2.ssf.StateInfo;
import ec.tstoolkit2.ssf.implementations.CompositeDynamics;
import ec.tstoolkit2.ssf.implementations.CompositeMeasurements;
import ec.tstoolkit2.ssf.multivariate.ISsfMeasurements;
import ec.tstoolkit2.ssf.multivariate.MultivariateSsf;
import ec.tstoolkit2.ssf.univariate.ISsf;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Jean Palate
 */
public class SsfSutse extends MultivariateSsf {

    public static SsfSutse of(SutseModel model) {

        // generate the ssf...
        ISsfDynamics[] dyns = new ISsfDynamics[model.getModelCount()];
        ISsf[] ssfs = new ISsf[model.getModelCount()];
        int nres = 0;
        int nst0 = 0;
        for (int i = 0; i < dyns.length; ++i) {
            BasicStructuralModel m = model.getModel(i);
            SsfBsm ssf = SsfBsm.create(m);
            ModelSpecification spec = m.getSpecification();
            dyns[i] = ssf.getDynamics();
            ssfs[i] = ssf;
            nres += ssf.getDynamics().getInnovationsDim();
            if (spec.hasNoise()) {
                ++nst0;
            }
            if (spec.hasCycle()) {
                nst0 += 2;
            }
        }
        Matrix var = Matrix.square(nres);
        int[] ivar = new int[nres];
        int icur = 0;
        icur += fill(model, Component.Noise, icur, ivar, var);
        icur += fill(model, Component.Cycle, icur, ivar, var);
        icur += fill(model, Component.Level, icur, ivar, var);
        icur += fill(model, Component.Slope, icur, ivar, var);
        icur += fill(model, Component.Seasonal, icur, ivar, var);
        Matrix st0 = Matrix.square(nst0);
        int[] ist0 = new int[nst0];
        icur = 0;
        icur += fill0(model, Component.Noise, icur, ist0, st0);
        icur += fill0(model, Component.Cycle, icur, ist0, st0);
        return new SsfSutse(model, new Dynamics(dyns, var, ivar, st0, ist0),
                CompositeMeasurements.of(ssfs));
    }

    private static int fill(final SutseModel model, final Component cmp, final int start, final int[] ivar, Matrix var) {
        // fill the innovations matrix. Start with noise, cycle, level, slope, seas
        int[] pos = positions(cmp, model);
        double[] v = variances(cmp, model);
        Matrix corr = model.getCorrelations(cmp);
        int icur = start;
        if (corr == null || corr.getColumnsCount() != pos.length) {
            for (int i = 0; i < pos.length; ++i, ++icur) {
                ivar[icur] = pos[i];
                var.set(icur, icur, v[i]);
            }
        } else {
            for (int i = 0; i < pos.length; ++i, ++icur) {
                ivar[icur] = pos[i];
                var.set(icur, icur, v[i]);
                for (int j = 0; j < i; ++j) {
                    double z = corr.get(i, j) * Math.sqrt(v[i] * v[j]);
                    var.set(start + i, start + j, z);
                    var.set(start + j, start + i, z);
                }
            }
        }
        return pos.length;
    }

    private static int fill0(final SutseModel model, final Component cmp, final int start, final int[] ivar, Matrix var) {
        // fill the innovations matrix. Start with noise, cycle, level, slope, seas
        if (cmp == Component.Noise) {
            return fill(model, cmp, start, ivar, var);
        } else if (cmp == Component.Cycle) {
            int[] pos = positions(cmp, model);
            double[] v = stcycles(model);
            Matrix corr = model.getCorrelations(cmp);
            int icur = start;
            if (corr == null || corr.getColumnsCount() != pos.length) {
                for (int i = 0, j=0; i < pos.length; i+=2, icur+=2, ++j) {
                    ivar[icur] = pos[i];
                    ivar[icur + 1] = pos[i + 1];
                    var.set(icur, icur, v[j]);
                    var.set(icur+1, icur+1, v[j]);
                }
            } else {
                for (int i = 0, k=0; i < pos.length; i+=2, icur+=2, ++k) {
                    ivar[icur] = pos[i];
                    ivar[icur + 1] = pos[i + 1];
                    var.set(icur, icur, v[k]);
                    var.set(icur+1, icur+1, v[k]);
                    for (int j = 0, l=0; j < i; j+=2, ++l) {
                        double z = corr.get(i, j) * Math.sqrt(v[k] * v[l]);
                        var.set(start + i, start + j, z);
                        var.set(start + j, start + i, z);
                        var.set(start + i+1, start + j+1, z);
                        var.set(start + j+1, start + i+1, z);
                        var.set(start + i, start + j+1, z);
                        var.set(start + j+1, start + i, z);
                        var.set(start + i+1, start + j, z);
                        var.set(start + j, start + i+1, z);
                    }
                }
            }
            return 2*pos.length;

        } else {
            return 0;
        }
    }

    private static int[] positions(Component cmp, SutseModel sutse) {
        IntList pos = new IntList();
        for (int i = 0, j = 0; i < sutse.getModelCount(); ++i) {
            BasicStructuralModel m = sutse.getModel(i);
            int cpos = SsfBsm.searchPosition(m, cmp);
            if (cpos >= 0) {
                pos.add(j + cpos);
                if (cmp == Component.Cycle) {
                    pos.add(j + cpos + 1);
                }
            }
            j += SsfBsm.calcDim(m);
        }
        return pos.toArray();
    }

    private static double[] variances(Component cmp, SutseModel sutse) {
        DoubleList vars = new DoubleList();
        for (int i = 0, j = 0; i < sutse.getModelCount(); ++i) {
            BasicStructuralModel m = sutse.getModel(i);
            if (m.getSpecification().hasComponent(cmp)) {
                vars.add(m.getVariance(cmp));
                if (cmp == Component.Cycle) {
                    vars.add(m.getVariance(cmp));
                }
            }
        }
        return vars.toArray();
    }

    private static double[] stcycles(SutseModel sutse) {
        DoubleList vars = new DoubleList();
        for (int i = 0, j = 0; i < sutse.getModelCount(); ++i) {
            BasicStructuralModel m = sutse.getModel(i);
            if (m.getSpecification().hasCycle()) {
                double q = m.cVar / (1 - m.cDump * m.cDump);
                vars.add(q);
            }
        }
        return vars.toArray();
    }

    private final SutseModel model;

    private SsfSutse(SutseModel model, ISsfDynamics dyn, ISsfMeasurements m) {
        super(dyn, m);
        this.model = model;
    }

    public SutseModel getModel() {
        return model;
    }

    static class Dynamics extends CompositeDynamics {

        private final Matrix var, st0;
        private final int[] ivar, ist0;

        Dynamics(ISsfDynamics[] dyn, Matrix var, int[] ivar, Matrix st0, int[] ist0) {
            super(dyn);
            this.var = var;
            this.ivar = ivar;
            this.st0 = st0;
            this.ist0 = ist0;
        }

        @Override
        public boolean Pf0(SubMatrix v, StateInfo info) {
            if (info == StateInfo.Forecast) {
                for (int r = 0; r < ist0.length; ++r) {
                    for (int c = 0; c < ist0.length; ++c) {
                        v.set(ist0[r], ist0[c], st0.get(r, c));
                    }
                }
                return true;
            } else {
                return false;
            }
        }

        @Override
        public void V(int pos, SubMatrix v) {
            v.copy(var.subMatrix(), ivar, ivar);
        }

        @Override
        public void addV(int pos, SubMatrix v) {
            for (int c = 0; c < ivar.length; ++c) {
                for (int r = 0; r < ivar.length; ++r) {
                    v.add(ivar[r], ivar[c], var.get(r, c));
                }
            }
        }
    }

}
