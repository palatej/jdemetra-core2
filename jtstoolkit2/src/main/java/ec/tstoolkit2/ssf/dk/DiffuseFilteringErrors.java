/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.dk;

import ec.tstoolkit.data.DataBlockStorage;
import ec.tstoolkit2.ssf.DataResults;
import ec.tstoolkit2.ssf.univariate.*;
import ec.tstoolkit2.ssf.State;

/**
 *
 * @author Jean Palate
 */
public class DiffuseFilteringErrors extends FilteringErrors implements IDiffuseFilteringResults {

    private int enddiffuse;

    public DiffuseFilteringErrors(boolean normalized) {
        super(normalized);
    }

    @Override
    public void close(int pos) {
        enddiffuse = pos;
    }

    @Override
    public void save(int t, DiffuseState state) {
    }

    @Override
    public void save(int t, DiffusePredictionError pe) {
        super.save(t, pe);
    }

    @Override
    public void clear() {
        super.clear();
        enddiffuse = 0;
    }

    @Override
    public int getEndDiffusePosition() {
        return enddiffuse;
    }
}
