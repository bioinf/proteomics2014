

import alignment.progressive.ProgressiveAlignment;
import common.ScoringMatrix;
import common.Sequence;
import common.SequenceFactory;
import edu.princeton.cs.introcs.Out;
import edu.princeton.cs.introcs.StdOut;
import guidetree.IHierarchicalClustering;
import guidetree.NeighborJoining;
import guidetree.UPGMA;
import guidetree.WPGMA;
import tools.Fasta;

import java.util.ArrayList;

public class Main {

    @SuppressWarnings("unchecked")
    public static void main(String[] args) {
        ParamParser params = new ParamParser(args);
        if (!params.state) {
            return;
        }

        ArrayList<Sequence> seqList = (ArrayList<Sequence>) Fasta.readSequences(params.inputSeqFile, new SequenceFactory());
        Sequence[] seqArr = new Sequence[seqList.size()];
        seqList.toArray(seqArr);
        ScoringMatrix scoringMatrix = new ScoringMatrix(params.inputMatrixFile, "#");

        IHierarchicalClustering clustering;
        if (params.clusteringMethod.equals("upgma")) {
            clustering = new UPGMA();
            StdOut.println("Progressive alignment with UPGMA");
        } else if (params.clusteringMethod.equals("wpgma")) {
            clustering = new WPGMA();
            StdOut.println("Progressive alignment with WPGMA");
        } else {
            clustering = new NeighborJoining();
            StdOut.println("Progressive alignment with Neighbor-Joining");
        }

        StdOut.println("-------------------------------------------------");
        ProgressiveAlignment progressiveAlignment = new ProgressiveAlignment(seqArr, scoringMatrix, clustering);
        StdOut.println("-------------------------------------------------");

        StdOut.println("Results:");
        StdOut.println("Score: " + progressiveAlignment.getScore());
        progressiveAlignment.getProfile().print();

        Out out = new Out(params.outputFilename);
        out.println("Score: " + progressiveAlignment.getScore());
        progressiveAlignment.getProfile().print(out);

        StdOut.println("\nJob done! The results are available in the output file");
    }
}
