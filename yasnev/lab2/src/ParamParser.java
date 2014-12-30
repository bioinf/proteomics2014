import java.io.BufferedReader;
import java.io.FileReader;

/**
 * Author: Oleg Yasnev (oyasnev@gmail.com)
 * Date: 25.11.13
 */

public class ParamParser {
    public static final String OUTPUT_FILE_KEY       = "-o";

    public static final String README = "readme.txt";

    public boolean state = false;

    public String inputPdbFile     = "";
    public double x, y, z;
    public String outputFilename   = "out.pdb";

    public ParamParser(String[] args) {
        if (args.length == 0 || args[0].equals("--help") || args[0].equals("?") || args[0].equals("-h")) {
            printHelp();
            return;
        }

        if (args.length < 4) {
            System.out.println("Error: four (4) input parameters expected");
            return;
        }

        inputPdbFile = args[0];
        x = Double.valueOf(args[1]);
        y = Double.valueOf(args[2]);
        z = Double.valueOf(args[3]);

        for (int i = 4; i < args.length; i++) {
            String key = args[i];
            if (key.equals(OUTPUT_FILE_KEY)) {
                i++;
                outputFilename = args[i];
            }
        }

        // check state
        state = true;
        if (inputPdbFile.isEmpty()) {
            state = false;
            System.out.println("Input PDB file must be specified");
        }

        if (!state) {
            System.out.println("For more information run the program with the key '--help'");
        }
    }

    public static void printHelp() {
        try {
            BufferedReader br = new BufferedReader(new FileReader(README));
            for (String line; (line = br.readLine()) != null; ) {
                System.out.println(line);
            }
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
    }
}
