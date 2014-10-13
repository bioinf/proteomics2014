package guidetree;

import common.Sequence;
import common.interfaces.ISequence;
import edu.princeton.cs.introcs.StdOut;
import settings.Config;

import java.util.ArrayList;

/**
 * Author: Oleg Yasnev (oyasnev@gmail.com)
 * Date: 19.03.14
 */
public class NeighborJoining implements IHierarchicalClustering {
    protected ArrayList<GuideTreeNode> nodeList;
    protected double[][] matrix;
    protected double[] rowSumMatrix;

    public NeighborJoining() { }

    @Override
    public void build(ISequence[] seqArr, int[][] pairwiseAlignmentMatrix) {
        int n = seqArr.length;
        nodeList = new ArrayList<GuideTreeNode>(n);
        // convert each sequence to a single cluster (this is more convenient than star tree)
        for (ISequence seq : seqArr) {
            nodeList.add(new GuideTreeNode(seq));
        }

        int maxScore = Integer.MIN_VALUE;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (pairwiseAlignmentMatrix[i][j] > maxScore) {
                    maxScore = pairwiseAlignmentMatrix[i][j];
                }
            }
        }

        // init score matrix with pairwise distances
        matrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                // convert score to distance (max_score + some_value - score)
                if (i == j) {
                    matrix[i][j] = 0;
                } else {
                    matrix[i][j] = maxScore + 10 - pairwiseAlignmentMatrix[i][j];
                }
            }
        }

        rowSumMatrix = new double[n];

        buildGuideTree();
    }

    @Override
    public GuideTreeNode getRoot() {
        return nodeList.get(0);
    }

    protected void buildGuideTree() {
        // repeat until we have only one cluster
        while (nodeList.size() > 1) {
            int n = nodeList.size();

            calcRowSumMatrix(n);

            if (Config.DEBUG) {
                debugPrintMatrix(n);
            }

            // find min Q
            PairInd minQ = findMinQ(n);

            // join clusters
            GuideTreeNode joinNode = new GuideTreeNode(nodeList.get(minQ.i), nodeList.get(minQ.j), matrix[minQ.i][minQ.j]);
            if (Config.DEBUG) {
                joinNode.seq = new Sequence(joinNode.left.seq.getDescription() + joinNode.right.seq.getDescription(), "");
            }

            // calc scores for the cluster
            double[] newScores = recalcScore(minQ.i, minQ.j, n);

            // remove old nodes and insert the new one
            removeOldInsertNew(minQ.i, minQ.j, joinNode, n, newScores);
        }
    }

    protected double[] recalcScore(int i, int j, int n) {
        double[] newScore = new double[n];
        for (int k = 0; k < n; k++) {
            if (k == i || k == j) {
                continue;
            }
            newScore[k] = (matrix[i][k] + matrix[j][k] - matrix[i][j]) / 2;
        }
        return newScore;
    }

    protected void removeOldInsertNew(int i, int j, GuideTreeNode newNode, int n, double[] newDist) {
        // remove old clusters
        if (i > j) {
            nodeList.remove(i);
            nodeList.remove(j);
        } else {
            nodeList.remove(j);
            nodeList.remove(i);
        }
        // insert the new one to the first position
        nodeList.add(0, newNode);

        // create new matrix
        double[][] newMatrix = new double[n-1][n-1];
        int row = 1, col = 0;
        for (int k = 0; k < n; k++) {
            if (k == i || k == j) {
                // skip old clusters
                continue;
            }
            newMatrix[row][col] = newDist[k];
            newMatrix[col][row] = newDist[k];
            col++;
            for (int l = 0; l < n; l++) {
                if (l == i || l == j) {
                    // skip old clusters
                    continue;
                }
                newMatrix[row][col] = matrix[k][l];
                col++;
            }
            row++;
            col = 0;
        }
        matrix = newMatrix;
    }

    protected PairInd findMinQ(int n) {
        if (Config.DEBUG) {
            StdOut.println("Q: ");
        }
        int minI = -1;
        int minJ = -1;
        double minQ = Double.MAX_VALUE;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double Q = (n - 2) * matrix[i][j] - rowSumMatrix[i] - rowSumMatrix[j];
                if (Config.DEBUG) {
                    StdOut.print(Q + "\t\t");
                }
                if (Q < minQ) {
                    minQ = Q;
                    minI = i;
                    minJ = j;
                }
            }
            if (Config.DEBUG) {
                StdOut.println();
            }
        }
        return new PairInd(minI, minJ);
    }

    protected void calcRowSumMatrix(int n) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += matrix[i][j];
            }
            rowSumMatrix[i] = sum;
        }
    }

    protected void debugPrintMatrix(int n) {
        StdOut.print('\t');
        for (int i = 0; i < n; i++) {
            StdOut.print(nodeList.get(i).seq.getDescription() + "\t\t");
        }
        StdOut.println();
        for (int i = 0; i < n; i++) {
            StdOut.print(nodeList.get(i).seq.getDescription() + "\t");
            for (int j = 0; j < n; j++) {
                StdOut.print(matrix[i][j] + "\t\t");
            }
            StdOut.println();
        }
        StdOut.println();
        StdOut.println();
    }

    private class PairInd {
        public int i, j;
        public PairInd(int i, int j) {
            this.i = i;
            this.j = j;
        }
    }
}
