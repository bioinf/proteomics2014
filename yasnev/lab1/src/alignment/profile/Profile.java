package alignment.profile;

import common.interfaces.ISequence;
import edu.princeton.cs.introcs.Out;
import edu.princeton.cs.introcs.StdOut;

/**
 * Author: Oleg Yasnev (oyasnev@gmail.com)
 * Date: 12.10.14
 */
public class Profile {

    public char[][] sequences;
    public String[] seqNames;
    public int seqN;
    public int seqL;

    public Profile(ISequence seq) {
        sequences = new char[1][];
        sequences[0] = seq.getSequence().toCharArray();
        seqN = 1;
        seqL = seq.getSequence().length();
        seqNames = new String[1];
        seqNames[0] = seq.getDescription();
    }

    public Profile(StringBuilder[] sbArr, Profile first, Profile second) {
        seqN = sbArr.length;
        sequences = new char[seqN][];
        for (int i = 0; i < sbArr.length; i++) {
            sequences[i] = sbArr[i].toString().toCharArray();
        }
        seqL = sequences[0].length;
        seqNames = new String[seqN];
        for (int i = 0; i < first.seqN; i++) {
            seqNames[i] = first.seqNames[i];
        }
        for (int i = 0; i < second.seqN; i++) {
            seqNames[i + first.seqN] = second.seqNames[i];
        }
    }

    public void print(Out out) {
        for (String name : seqNames) {
            out.println(name);
        }
        for (char[] seq : sequences) {
            out.println(String.valueOf(seq));
        }
    }

    public void print() {
        for (String name : seqNames) {
            StdOut.println(name);
        }
        for (char[] seq : sequences) {
            StdOut.println(String.valueOf(seq));
        }
    }
}
