package guidetree;

import common.interfaces.ISequence;

/**
 * Author: Oleg Yasnev (oyasnev@gmail.com)
 * Date: 13.10.14
 */
public interface IHierarchicalClustering {
    public void build(ISequence[] seqArr, int[][] pairwiseAlignmentMatrix);
    public GuideTreeNode getRoot();
}
