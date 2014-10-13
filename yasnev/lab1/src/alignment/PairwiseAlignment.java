package alignment;

import common.interfaces.IMatchScoringMatrix;
import common.interfaces.ISequence;
import edu.princeton.cs.introcs.StdOut;

/**
 * Author: Oleg Yasnev (oyasnev@gmail.com)
 * Date: 19.03.14
 */
public class PairwiseAlignment {
    protected IMatchScoringMatrix scoringMatrix;
    protected int[][] matrix;
    protected ISequence[] seqArr;

    public PairwiseAlignment(ISequence[] seqArr, IMatchScoringMatrix scoringMatrix) {
        this.seqArr = seqArr;
        this.scoringMatrix = scoringMatrix;
        align();
    }

    public int[][] getMatrix() {
        return matrix;
    }

    protected void align() {
        int n = seqArr.length;
        matrix = new int[n][n];
        for (int i = 0; i < n; i++) {
            StdOut.printf("%d of %d\n", i + 1, n);
            matrix[i][i] = 0;
            String first = seqArr[i].getSequence();
            for (int j = i + 1; j < n; j++) {
                String second = seqArr[j].getSequence();
                Alignment ed = new Alignment(first, second, scoringMatrix);
                matrix[i][j] = ed.getScore();
                matrix[j][i] = ed.getScore();
            }
        }
    }
}
