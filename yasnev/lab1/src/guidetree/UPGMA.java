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
public class UPGMA extends WPGMA implements IHierarchicalClustering {
    public UPGMA() { }

    @Override
    protected double[] recalcScore(int i, int j, int n) {
        GuideTreeNode first = nodeList.get(i);
        GuideTreeNode second = nodeList.get(j);
        double[] newScore = new double[n];
        for (int k = 0; k < n; k++) {
            if (k == i || k == j) {
                continue;
            }
            newScore[k] = (first.size * matrix[i][k] + second.size * matrix[j][k]) / (first.size + second.size);
        }
        return newScore;
    }
}
