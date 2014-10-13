import edu.princeton.cs.introcs.In;
import edu.princeton.cs.introcs.StdOut;

/**
 * Author: Oleg Yasnev (oyasnev@gmail.com)
 * Date: 25.11.13
 */

public class ParamParser {
    public static final String OUTPUT_FILE_KEY       = "-o";

    public static final String README = "readme.txt";

    public boolean state = false;

    public String inputSeqFile     = "";
    public String inputMatrixFile  = "";
    public String clusteringMethod = "";
    public String outputFilename   = "output.txt";

    public ParamParser(String[] args) {
        if (args.length == 0 || args[0].equals("--help") || args[0].equals("?") || args[0].equals("-h")) {
            printHelp();
            return;
        }

        if (args.length < 3) {
            StdOut.println("Error: three input parameters expected");
            return;
        }

        inputSeqFile    = args[0];
        inputMatrixFile = args[1];
        clusteringMethod = args[2];

        for (int i = 3; i < args.length; i++) {
            String key = args[i];
            if (key.equals(OUTPUT_FILE_KEY)) {
                i++;
                outputFilename = args[i];
            }
        }

        // check state
        state = true;
        if (inputSeqFile.isEmpty()) {
            state = false;
            StdOut.println("Input sequence file must be specified");
        }
        if (inputMatrixFile.isEmpty()) {
            state = false;
            StdOut.println("Input scoring matrix file must be specified");
        }
        if (clusteringMethod.isEmpty()) {
            state = false;
            StdOut.println("Clustering method must be specified");
        }

        if (!state) {
            StdOut.println("For more information run the program with the key '--help'");
        }
    }

    public static void printHelp() {
        In in = new In(README);
        while (!in.isEmpty()) {
            StdOut.println(in.readLine());
        }
    }
}
