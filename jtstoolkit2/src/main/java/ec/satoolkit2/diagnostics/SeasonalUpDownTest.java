/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.satoolkit2.diagnostics;

import ec.tstoolkit.data.DataBlock;
import ec.tstoolkit.data.IReadDataBlock;
import ec.tstoolkit.stats.TestofUpDownRuns;
import ec.tstoolkit.timeseries.simplets.PeriodIterator;
import ec.tstoolkit.timeseries.simplets.TsData;
import ec.tstoolkit.timeseries.simplets.TsDataBlock;
import ec.tstoolkit.timeseries.simplets.TsPeriod;

/**
 *
 * @author Jean Palate
 */
public class SeasonalUpDownTest {

    private final int[] ud, udd;
    private final TsData s, ds;

    public SeasonalUpDownTest(TsData s) {
        this.s = s;
        int freq = s.getFrequency().intValue();
        ds = s.delta(freq);
        ud = new int[freq];
        udd = new int[freq];
        for (int i = 0; i < freq; ++i) {
            IReadDataBlock sdata = seasonalData(s, i);
            IReadDataBlock dsdata = seasonalData(ds, i);
            ud[i]=dcount(sdata);
            udd[i]=dcount(dsdata);
        }
    }

    public TestofUpDownRuns runTest(int period, boolean diff) {
        TsData cur = diff ? ds : s;
        TestofUpDownRuns test = new TestofUpDownRuns();
        test.test(seasonalData(cur, period));
        return test;
    }

    private static IReadDataBlock seasonalData(TsData s, int period) {
        TsPeriod start = s.getStart();
        int i0 = period - start.getPosition();
        int freq=start.getFrequency().intValue();
        if (i0 < 0) {
            i0 += freq;
        }
        DataBlock all=TsDataBlock.all(s).data;
        return all.extract(i0, -1, freq);
    }
    
    private int dcount(IReadDataBlock data){
        int n=data.getLength();
        if (n<3)
            return 0;
        int c=0;
        double prev=data.get(1);
        boolean prevpos = prev>data.get(0);
        for (int i=2; i<n; ++i){
            double cur=data.get(i);
            boolean pos=cur>prev;
            if (pos != prevpos){
                ++c;
            }
            prev=cur;
            prevpos=pos;
        }
        return c;
    }

    public int getUpDownCount(int period) {
        return ud[period];
    }

    public int getUpDownCount() {
        int all = 0;
        for (int i = 0; i < ud.length; ++i) {
            all += ud[i];
        }
        return all;
    }

    public int getUpDownDifferenceCount(int period) {
        return udd[period];
    }

    public int getUpDownDifferenceCount() {
        int all = 0;
        for (int i = 0; i < udd.length; ++i) {
            all += udd[i];
        }
        return all;
    }
}
