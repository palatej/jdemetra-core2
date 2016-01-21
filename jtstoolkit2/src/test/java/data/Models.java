/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package data;

import ec.tstoolkit.arima.ArimaModel;
import ec.tstoolkit.arima.ArimaModelBuilder;
import ec.tstoolkit.sarima.SarimaModel;
import ec.tstoolkit.sarima.SarimaModelBuilder;
import ec.tstoolkit.timeseries.simplets.TsData;
import ec.tstoolkit.ucarima.ModelDecomposer;
import ec.tstoolkit.ucarima.SeasonalSelector;
import ec.tstoolkit.ucarima.TrendCycleSelector;
import ec.tstoolkit.ucarima.UcarimaModel;
import ec.tstoolkit2.ssf.implementations.arima.SsfUcarima;
import ec.tstoolkit2.ssf.univariate.SsfData;
import java.util.Random;

/**
 *
 * @author Admin
 */
public class Models {
    
    public final static SsfUcarima ssfUcarima;
    public final static SsfData ssfRandom, ssfRandomMissing, ssfProd, ssfX, ssfXRandom;
    public final static TsData xts;

    static {
        int M = 50;
        TrendCycleSelector tsel = new TrendCycleSelector(.5);
        tsel.setDefaultLowFreqThreshold(12);
        SeasonalSelector ssel = new SeasonalSelector(12, 3);

        ModelDecomposer decomposer = new ModelDecomposer();
        decomposer.add(tsel);
        decomposer.add(ssel);
        xts = data.Data.X.clone();
        int[] missing = new int[M];
        Random rng = new Random();
        for (int i = 0; i < M; ++i) {
            missing[i] = rng.nextInt(xts.getLength());
        }
        SarimaModel arima = new SarimaModelBuilder().createAirlineModel(12, -.8, -.9);
        UcarimaModel ucm = decomposer.decompose(ArimaModel.create(arima));
        ucm.setVarianceMax(-1);
        ucm.simplify();

        for (int i = 0; i < M; ++i) {
            xts.setMissing(missing[i]);
        }

        ssfUcarima = SsfUcarima.create(ucm);
        ssfXRandom = new SsfData(xts);
        
        ssfProd=new SsfData(Data.P);
        ssfX=new SsfData(Data.X);
        
        ArimaModelBuilder builder=new ArimaModelBuilder();
        double[] rnds = builder.generate(arima,360);
        ssfRandom=new SsfData(rnds.clone());
         for (int i = 0; i < rnds.length; ++i) {
            rnds[rng.nextInt(rnds.length)]=Double.NaN;
        }
        ssfRandomMissing=new SsfData(rnds);
       
    }
 }
