package guidetree;

import common.interfaces.ISequence;

import java.util.Locale;

/**
 * Author: Oleg Yasnev (oyasnev@gmail.com)
 * Date: 19.03.14
 */
public class GuideTreeNode {
    public boolean isLeaf;
    public double dist;
    public ISequence seq;
    public GuideTreeNode left;
    public GuideTreeNode right;
    public int size;

    public GuideTreeNode(ISequence seq) {
        isLeaf = true;
        this.seq = seq;
        this.size = 1;
    }

    public GuideTreeNode(GuideTreeNode left, GuideTreeNode right, double dist) {
        isLeaf = false;
        this.dist = dist;
        this.left = left;
        this.right = right;
        this.size = left.size + right.size;
    }

    public String toNewick() {
        if (isLeaf) {
            return seq.getDescription();
        } else {
            return String.format(Locale.US, "(%s : %.1f : %s)", left.toNewick(), dist, right.toNewick());
        }
    }
}
