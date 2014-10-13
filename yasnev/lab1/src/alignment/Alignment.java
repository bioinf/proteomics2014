package alignment;

import common.interfaces.IMatchScoringMatrix;
import edu.princeton.cs.introcs.StdOut;

enum Direction {
    NULL,
    AB,
    _B, A_
}

/**
 * Author: Oleg Yasnev (oyasnev@gmail.com)
 * Date: 16.02.14
 */
public class Alignment implements IAlignment {
    protected int score;

    protected char[] first;
    protected char[] second;

    protected IMatchScoringMatrix scoringMatrix;
    protected char gap;

    protected int[][] matrix;
    protected Direction[][] dirMatrix;
    protected int m, n;

    public Alignment(String strFirst, String strSecond, IMatchScoringMatrix scoringMatrix) {
        this.scoringMatrix = scoringMatrix;
        gap = scoringMatrix.getGapSymbol();
        // add empty string to the beginning of every string
        first  = (gap + strFirst ).toCharArray();
        second = (gap + strSecond).toCharArray();
        m = first .length;
        n = second.length;
        align();
        //backtracking();
    }

    @Override
    public int getScore() {
        return score;
    }

    @Override
    public String[] getAlignment() {
        return new String[2];
    }

    protected void align() {
        matrix = new int[m][n];
        dirMatrix = new Direction[m][n];

        // Align empty string
        matrix[0][0] = scoringMatrix.getScore(gap, gap);
        dirMatrix[0][0] = Direction.NULL;

        for (int i = 1; i < m; i++) {
            matrix[i][0] = matrix[i-1][0] + scoringMatrix.getScore(first[i], gap);
            dirMatrix[i][0] = Direction.A_;
        }

        for (int i = 1; i < n; i++) {
            matrix[0][i] = matrix[0][i-1] + scoringMatrix.getScore(gap, second[i]);
            dirMatrix[0][i] = Direction._B;
        }

        // Align dinamycally
        for (int i = 1; i < m; i++) {
            for (int j = 1; j < n; j++) {
                int AB = matrix[i-1][j-1] + scoringMatrix.getScore(first[i], second[j]);
                int A_ = matrix[i-1][ j ] + scoringMatrix.getScore(first[i], gap);
                int _B = matrix[ i ][j-1] + scoringMatrix.getScore(gap, second[j]);

                int score = Math.max(AB, Math.max(_B, A_));
                if (score == AB) {
                    dirMatrix[i][j] = Direction.AB;
                } else if (score == _B) {
                    dirMatrix[i][j] = Direction._B;
                } else {
                    dirMatrix[i][j] = Direction.A_;
                }
                matrix[i][j] = score;
            }
        }

        score = matrix[m-1][n-1];
    }

    protected void backtracking() {
        // init
        int i = m - 1;
        int j = n - 1;
        StringBuilder sbFirst  = new StringBuilder();
        StringBuilder sbSecond = new StringBuilder();

        // while not in a cell (0, 0)
        while (dirMatrix[i][j] != Direction.NULL) {
            switch (dirMatrix[i][j]) {
                case AB:
                    sbFirst .append(first [i--]);
                    sbSecond.append(second[j--]);
                    break;
                case A_:
                    sbFirst .append(first[i--] );
                    sbSecond.append(gap        );
                    break;
                case _B:
                    sbFirst .append(gap        );
                    sbSecond.append(second[j--]);
                    break;
                default:
                    StdOut.printf("ERROR while backtracking at (%d, %d)", i, j);
                    return;
            }
        }

        // reverse as we went from the end to beginning
        sbFirst .reverse();
        sbSecond.reverse();
        StdOut.println(sbFirst);
        StdOut.println(sbSecond);
    }
}
