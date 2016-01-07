/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ec.tstoolkit2.data;

import ec.tstoolkit.data.DataBlock;

/**
 *
 * @author Admin
 */
public class DataIterators {

    public static IDataIterator iterator(DataBlock data) {
        switch (data.getIncrement()) {
            case 1:
                return new Iteratorp1(data);
            case -1:
                return new Iteratorm1(data);
            default:
                return new IteratorN(data);
        }
    }

    static class IteratorN implements IDataIterator {

        private final int inc, end;
        private int cur;
        private final double[] data;

        IteratorN(DataBlock b) {
            inc = b.getIncrement();
            cur = b.getStartPosition();
            end = b.getEndPosition();
            data = b.getData();
        }

        @Override
        public boolean hasNext() {
            return cur != end;
        }

        @Override
        public double next() {
            double val = data[cur];
            cur += inc;
            return val;
        }

        @Override
        public double get() {
            return data[cur];
        }

        @Override
        public void set(double val) {
            data[cur] = val;
        }

        @Override
        public void setAndNext(double val) {
            data[cur] = val;
            cur += inc;
        }

        @Override
        public boolean advance(int n) {
            cur += inc * n;
            if (inc > 0) {
                if (cur > end) {
                    cur = end;
                }
            } else if (cur < end) {
                cur = end;
            }
            return cur != end;
        }

        @Override
        public int nextCount() {
            return (end - cur) / inc;
        }

    }

    static class Iteratorp1 implements IDataIterator {

        private final int end;
        private int cur;
        private final double[] data;

        Iteratorp1(DataBlock b) {
            cur = b.getStartPosition();
            end = b.getEndPosition();
            data = b.getData();
        }

        @Override
        public boolean hasNext() {
            return cur != end;
        }

        @Override
        public double next() {
            return data[cur++];
        }

        @Override
        public double get() {
            return data[cur];
        }

        @Override
        public void set(double val) {
            data[cur] = val;
        }

        @Override
        public void setAndNext(double val) {
            data[cur] = val;
            cur++;
        }

        @Override
        public boolean advance(int n) {
            cur += n;
            if (cur > end) {
                cur = end;
            }
            return cur != end;
        }

        @Override
        public int nextCount() {
            return end - cur;
        }

    }

    static class Iteratorm1 implements IDataIterator {

        private final int end;
        private int cur;
        private final double[] data;

        Iteratorm1(DataBlock b) {
            cur = b.getStartPosition();
            end = b.getEndPosition();
            data = b.getData();
        }

        @Override
        public boolean hasNext() {
            return cur != end;
        }

        @Override
        public double next() {
            return data[cur--];
        }

        @Override
        public double get() {
            return data[cur];
        }

        @Override
        public void set(double val) {
            data[cur] = val;
        }

        @Override
        public void setAndNext(double val) {
            data[cur] = val;
            cur--;
        }

        @Override
        public boolean advance(int n) {
            cur -= n;
            if (cur < end) {
                cur = end;
            }
            return cur != end;

        }

        @Override
        public int nextCount() {
            return cur - end;
        }
    }
}
