import org.biojava.bio.structure.*;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;

import java.io.IOException;

public class Main
{
    public static void main(String[] args)
    {
        String filename = "E:\\Dropbox\\3_semester\\structural_biology\\3CSY.pdb";
        PDBFileReader pdb_reader = new PDBFileReader();
        FileParsingParameters params = new FileParsingParameters();
        params.setAlignSeqRes(true);
//        params.setParseSecStruc(true);
//        params.setCreateAtomBonds(true);
        pdb_reader.setFileParsingParameters(params);

        try
        {
            Structure structure = pdb_reader.getStructure(filename);
//            SecStruc sec_struct = new SecStruc();
//            sec_struct.assign(structure);
            int helices_assigned = PdbUtils.AssignHelices(structure, 4, -0.3, false);
            int beta_sheets_assigned = PdbUtils.AssignBetaSheetsToStructure(structure, -0.3);
//            int res = PdbUtils.RemoveExistingBonds(structure);
//            int res = PdbUtils.CalculateHydrogenBonds(structure);
            int a = 5;
        }
        catch (IOException | StructureException e)
        {
            e.printStackTrace();
        }
    }

    public static void TestAssignHelicesPerformance()
    {

    }
}
