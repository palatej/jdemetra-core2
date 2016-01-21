/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf.dk;

import data.Models;
import org.junit.Test;

/**
 *
 * @author Admin
 */
public class DiffuseSimulationSmootherTest {

    final static int N = 10000;

    public DiffuseSimulationSmootherTest() {
    }

    @Test
    public void testGenerate() {
        DiffuseSimulationSmoother smoother = new DiffuseSimulationSmoother(Models.ssfUcarima, Models.ssfProd);
        DiffuseSimulationSmoother.Simulation simul = smoother.newSimulation();
        simul.getSimulatedData();
    }

    @Test
    public void stressTestGenerate() {
        long t0 = System.currentTimeMillis();
        DiffuseSimulationSmoother smoother = new DiffuseSimulationSmoother(Models.ssfUcarima, Models.ssfProd);
        for (int i = 0; i < N; ++i) {
            DiffuseSimulationSmoother.Simulation simul = smoother.newSimulation();
            simul.getSimulatedData();
        }
        long t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
    }
}
