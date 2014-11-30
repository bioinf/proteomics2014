import org.biojava.bio.structure.*;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;

import java.io.IOException;

public class Main
{
    public static void main(String[] args)
    {
        String filename = args[0];
        double energy_threshold = Double.parseDouble(args[1]);
        PDBFileReader pdb_reader = new PDBFileReader();
        FileParsingParameters params = new FileParsingParameters();
        params.setAlignSeqRes(true);
        pdb_reader.setFileParsingParameters(params);

        try
        {
            Structure structure = pdb_reader.getStructure(filename);
            int helices_assigned = PdbUtils.AssignHelices(structure, 4, energy_threshold, false);
            int beta_sheets_assigned = PdbUtils.AssignBetaSheetsToStructure(structure, energy_threshold);

            System.out.println("Helices assigned: " + helices_assigned);
            System.out.println("Beta sheets assigned: " + beta_sheets_assigned);
        }
        catch (IOException | StructureException e)
        {
            e.printStackTrace();
        }
    }
}
