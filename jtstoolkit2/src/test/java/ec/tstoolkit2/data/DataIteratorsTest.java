/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.data;

import ec.tstoolkit.data.DataBlock;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author Jean Palate
 */
public class DataIteratorsTest {
    
    public DataIteratorsTest() {
    }

    //@Ignore
    @Test
    public void stressTestDot() {
        int L=1000, K=1000000;
        DataBlock a=new DataBlock(L);
        a.randomize();
        DataBlock b=new DataBlock(L);
        b.randomize();
        double s=0;
        long t0=System.currentTimeMillis();
        for (int i=0; i<K; ++i){
            s=a.dot(b);
        }
        long t1=System.currentTimeMillis();
        System.out.println(t1-t0);
        System.out.println(s);
        t0=System.currentTimeMillis();
        for (int i=0; i<K; ++i){
            s=dot(a,b);
        }
        t1=System.currentTimeMillis();
        System.out.println(t1-t0);
        System.out.println(s);
    }
    
    public double dot(DataBlock a, DataBlock b){
        IDataIterator ia = DataIterators.iterator(a);
        IDataIterator ib = DataIterators.iterator(b);
        double s=0;
        while (ia.hasNext()){
            s+=ia.next()*ib.next();
        }
        return s;
    }
}
