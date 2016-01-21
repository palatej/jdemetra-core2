/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.dk;

import data.Models;
import ec.tstoolkit.data.DataBlock;
import org.junit.Test;

/**
 *
 * @author Admin
 */
public class DiffuseSimulationSmootherTest {

    final static int N = 5000;

    public DiffuseSimulationSmootherTest() {
    }

    @Test
    public void testGenerate() {
        DiffuseSimulationSmoother smoother = new DiffuseSimulationSmoother(Models.ssfUcarima, Models.ssfXRandom);
        DiffuseSimulationSmoother.Simulation simul = smoother.newSimulation();
        simul.getSimulatedData();
        DataBlock item = smoother.getReferenceSmoothing().getSmoothedStates().item(0);
        System.out.println(item);
        System.out.println(DkToolkit.smooth(Models.ssfUcarima, Models.ssfXRandom, false).getComponent(0));
    }

    @Test
    public void stressTestGenerate() {
        long t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            DiffuseSimulationSmoother smoother = new DiffuseSimulationSmoother(Models.ssfUcarima, Models.ssfProd);
            DiffuseSimulationSmoother.Simulation simul = smoother.newSimulation();
            simul.getSimulatedData();
        }
        long t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
        t0 = System.currentTimeMillis();
        for (int i = 0; i < N; ++i) {
            DkToolkit.smooth(Models.ssfUcarima, Models.ssfXRandom, false).getComponent(0);
        }
        t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
    }
}
