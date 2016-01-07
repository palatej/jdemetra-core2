/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.data;

/**
 *
 * @author Admin
 */
public interface IDataReader {
    boolean hasNext();
    double next();
    double get();
    int nextCount();
    boolean advance(int n);
}
