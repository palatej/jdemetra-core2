/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.ssf;

import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Jean Palate
 */
public class ResultsRange implements Cloneable{
    
    private int start, end;
    
    public ResultsRange(){}
    
    public ResultsRange(int start, int end){
        this.start=start;
        this.end=end;
    }
    
    @Override
    public ResultsRange clone(){
        try {
            return (ResultsRange) super.clone();
        } catch (CloneNotSupportedException ex) {
            return null;
        }
    }

    public boolean isEmpty(){
        return start>=end;
    }
    
    public int getStart(){
        return start;
    }
    
    public int getEnd(){
        return end;
    }
    
    public void add(int pos){
        if (pos<start)
            start=pos;
        else if (pos>=end)
            end=pos+1;
    }

    public void setRange(int start, int end){
        this.start=start;
        this.end=end;
    }
    
    public void clear(){
        start=end=0;
    }
}
