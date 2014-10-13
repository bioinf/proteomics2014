package alignment.progressive;

import alignment.PairwiseAlignment;
import alignment.profile.Profile;
import alignment.profile.ProfileAlignment;
import common.interfaces.IMatchScoringMatrix;
import common.interfaces.ISequence;
import edu.princeton.cs.introcs.StdOut;
import guidetree.GuideTreeNode;
import guidetree.IHierarchicalClustering;
import guidetree.UPGMA;

/**
 * Author: Oleg Yasnev (oyasnev@gmail.com)
 * Date: 13.10.14
 */
public class ProgressiveAlignment {

    protected IMatchScoringMatrix scoringMatrix;
    protected Profile profile;
    protected long score;

    public ProgressiveAlignment(ISequence[] seqArr, IMatchScoringMatrix scoringMatrix, IHierarchicalClustering clustering) {
        this.scoringMatrix = scoringMatrix;

        StdOut.println("Pairwise alignment...");
        PairwiseAlignment alignment = new PairwiseAlignment(seqArr, scoringMatrix);

        StdOut.println("Build guide tree...");
        clustering.build(seqArr, alignment.getMatrix());

        String newick = clustering.getRoot().toNewick();
        StdOut.println(newick);

        profile = alignTree(clustering.getRoot());

        // calc score
        for (int i = 0; i < profile.seqL; i++) {
            for (int j = 0; j < profile.seqN; j++) {
                for (int k = j + 1; k < profile.seqN; k++) {
                    score += scoringMatrix.getScore(profile.sequences[j][i], profile.sequences[k][i]);
                }
            }
        }
    }

    public Profile getProfile() {
        return profile;
    }

    public long getScore() {
        return score;
    }

    protected Profile alignTree(GuideTreeNode node) {
        if (node.isLeaf) {
            return new Profile(node.seq);
        } else {
            return (new ProfileAlignment(alignTree(node.left), alignTree(node.right), scoringMatrix).getProfile());
        }
    }
}
