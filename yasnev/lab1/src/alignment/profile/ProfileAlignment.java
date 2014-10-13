package alignment.profile;

import common.interfaces.IMatchScoringMatrix;
import edu.princeton.cs.introcs.StdOut;

enum Direction {
    NULL,
    AB,
    _B, A_
}

/**
 * Author: Oleg Yasnev (oyasnev@gmail.com)
 * Date: 12.10.14
 */
public class ProfileAlignment {
    protected long score;
    protected Profile profile;

    protected Profile first;
    protected Profile second;

    protected IMatchScoringMatrix scoringMatrix;
    protected char gap;

    protected long[][] matrix;
    protected Direction[][] dirMatrix;
    protected int m, n;

    public ProfileAlignment(Profile first, Profile second, IMatchScoringMatrix scoringMatrix) {
        this.scoringMatrix = scoringMatrix;
        gap = scoringMatrix.getGapSymbol();
        this.first = first;
        this.second = second;
        m = first.seqL + 1;
        n = second.seqL + 1;
        align();
        backtracking();
    }

    public long getScore() {
        return score;
    }

    public Profile getProfile() {
        return profile;
    }

    protected void align() {
        matrix = new long[m][n];
        dirMatrix = new Direction[m][n];

        // Align empty string
        matrix[0][0] = first.seqN * second.seqN * scoringMatrix.getScore(gap, gap);
        dirMatrix[0][0] = Direction.NULL;

        for (int i = 1; i < m; i++) {
            matrix[i][0] = matrix[i-1][0] + psp(i, -1);
            dirMatrix[i][0] = Direction.A_;
        }

        for (int i = 1; i < n; i++) {
            matrix[0][i] = matrix[0][i-1] + psp(-1, i);
            dirMatrix[0][i] = Direction._B;
        }

        // Align dinamycally
        for (int i = 1; i < m; i++) {
            for (int j = 1; j < n; j++) {
                long AB = matrix[i-1][j-1] + psp(i, j);
                long A_ = matrix[i-1][ j ] + psp(i, -1);
                long _B = matrix[ i ][j-1] + psp(-1, j);

                long score = Math.max(AB, Math.max(_B, A_));
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
        StringBuilder[] sbArr = new StringBuilder[first.seqN + second.seqN];
        for (int k = 0; k < sbArr.length; k++) {
            sbArr[k] = new StringBuilder();
        }

        // while not in a cell (0, 0)
        while (dirMatrix[i][j] != Direction.NULL) {
            switch (dirMatrix[i][j]) {
                case AB:
                    for (int k = 0; k < first.seqN; k++) {
                        sbArr[k].append(first.sequences[k][i-1]);
                    }
                    for (int k = 0; k < second.seqN; k++) {
                        sbArr[k + first.seqN].append(second.sequences[k][j-1]);
                    }
                    i--;
                    j--;
                    break;
                case A_:
                    for (int k = 0; k < first.seqN; k++) {
                        sbArr[k].append(first.sequences[k][i-1]);
                    }
                    for (int k = 0; k < second.seqN; k++) {
                        sbArr[k + first.seqN].append(gap);
                    }
                    i--;
                    break;
                case _B:

                    for (int k = 0; k < first.seqN; k++) {
                        sbArr[k].append(gap);
                    }
                    for (int k = 0; k < second.seqN; k++) {
                        sbArr[k + first.seqN].append(second.sequences[k][j-1]);
                    }
                    j--;
                    break;
                default:
                    StdOut.printf("ERROR while backtracking at (%d, %d)", i, j);
                    return;
            }
        }

        // reverse as we went from the end to beginning
        for (int k = 0; k < sbArr.length; k++) {
            sbArr[k].reverse();
        }

        profile = new Profile(sbArr, first, second);
    }

    protected long psp(int firstCol, int secondCol) {
        long score = 0;
        if (firstCol == -1) {
            for (int i = 0; i < second.seqN; i++) {
                score += first.seqN * scoringMatrix.getScore(second.sequences[i][secondCol - 1], gap);
            }
        } else if (secondCol == -1) {
            for (int i = 0; i < first.seqN; i++) {
                score += second.seqN * scoringMatrix.getScore(first.sequences[i][firstCol - 1], gap);
            }
        } else {
            for (int i = 0; i < first.seqN; i++) {
                for (int j = 0; j < second.seqN; j++) {
                    score += scoringMatrix.getScore(first.sequences[i][firstCol - 1], second.sequences[j][secondCol - 1]);
                }
            }
        }
        return score;
    }
}
