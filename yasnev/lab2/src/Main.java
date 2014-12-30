import common.Pdb;
import common.Vector;

public class Main {

    public static void main(String[] args) {
        ParamParser params = new ParamParser(args);
        if (!params.state) {
            return;
        }

        System.out.printf("Reading data from %s...\n", params.inputPdbFile);
        Pdb pdb = new Pdb(params.inputPdbFile);

        System.out.println("Processing...");
        new CCD(pdb.atoms, new Vector(params.x, params.y, params.z), 0.001, 1e-6, 100);

        System.out.printf("Writing results to %s...\n", params.outputFilename);
        pdb.writeToFile(params.outputFilename);

        System.out.println("\nJob done! The results are available in the output file");
    }
}
