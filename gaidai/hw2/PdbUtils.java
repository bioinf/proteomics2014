import org.biojava.bio.structure.*;

import java.util.Collections;
import java.util.List;
import java.util.Map;

public class PdbUtils
{
    /**
     * Calculates approximate position of N-hydrogen
     *
     * @param c  C atom from previous aminoacid
     * @param n  N atom from current aminoacid
     * @param ca CA atom from current aminoacid
     * @return Approximate position of N-hydrogen
     */
    public static Atom CalculateHydrogenPosition(Atom c, Atom n, Atom ca)
    {
        Atom nc = Calc.subtract(n, c);
        Atom nca = Calc.subtract(n, ca);

        Atom u_nc = Calc.unitVector(nc);
        Atom u_nca = Calc.unitVector(nca);

        Atom added = Calc.subtract(u_nc, u_nca);
        //(According to Creighton the distance N-H is 1.03 +/- 0.02 Å.)
        Atom u = Calc.unitVector(added);

        return Calc.add(n, u);
    }

    /**
     * Clears bond structure
     * @param structure structure to clear
     * @return amount of deleted bonds
     */
    public static int RemoveExistingBonds(Structure structure)
    {
        int bonds_removed = 0;
        AtomIterator iterator = new AtomIterator(structure);
        while (iterator.hasNext())
        {
            Atom next = iterator.next();
            List<Bond> bonds = next.getBonds();
            bonds_removed += bonds.size();
            bonds.clear();
        }

        // because bonds store in both atoms so we double-counted them
        return bonds_removed / 2;
    }

    /**
     * Calculates interaction energy between specified amino acids
     * @param amino_acids list of all amino acids in chain
     * @param i first amino acid
     * @param j second amino acid
     * @return interaction energy kcal/mol
     * @throws StructureException
     */
    public static double InteractionEnergy(final List<AminoAcid> amino_acids, int i, int j) throws StructureException
    {
        AminoAcid first = amino_acids.get(i);
        AminoAcid second = amino_acids.get(j);
        AminoAcid previous_second = amino_acids.get(j - 1);
        Atom hydrogen = CalculateHydrogenPosition(previous_second.getC(), second.getN(), second.getCA());

        double ron = Calc.getDistance(first.getO(), second.getN());
        double rch = Calc.getDistance(first.getC(), hydrogen);
        double roh = Calc.getDistance(first.getO(), hydrogen);
        double rcn = Calc.getDistance(first.getC(), second.getN());

        // magic consts - http://en.wikipedia.org/wiki/DSSP_(protein)
        return 0.084 * 332 * (1 / ron + 1 / rch - 1 / roh - 1 / rcn);
    }

    /**
     * Assigns helices as secondary structure to appropriate amino acids in whole structure
     * @param structure structure for assignment
     * @param n type of helix. n = {3, 4, 5}
     * @param energy_threshold any energy below threshold supposed to form a hydrogen bond
     * @param change_direction true if you want to walk through chain in reverse direction
     * @return amount of helices assigned
     * @throws StructureException
     */
    @SuppressWarnings ("unchecked")
    public static int AssignHelices(Structure structure, int n, double energy_threshold, boolean change_direction) throws StructureException
    {
        int structures_assigned = 0;
        for (Chain chain : structure.getChains())
        {
            List<Group> groups = chain.getAtomGroups("amino");
            List<AminoAcid> amino_acids = (List<AminoAcid>) (List<?>) groups;
            if (change_direction)
            {
                Collections.reverse(amino_acids);
            }
            structures_assigned += HelixAssignment(amino_acids, n, energy_threshold);
        }

        return structures_assigned;
    }

    /**
     * Assigns helices as secondary structure to appropriate amino acids in specified list of amino acids
     * @param amino_acids list of all amino acids in chain
     * @param n type of helix. n = {3, 4, 5}
     * @param energy_threshold any energy below threshold supposed to form a hydrogen bond
     * @return amount of helices assigned
     * @throws StructureException
     */
    public static int HelixAssignment(List<AminoAcid> amino_acids, int n, double energy_threshold) throws StructureException
    {
        int structures_assigned = 0;
        for (int i = 2; i < amino_acids.size() - n; ++i)
        {
            AminoAcid amino_acid = amino_acids.get(i);

            double four_turn = InteractionEnergy(amino_acids, i - 1, i - 1 + n);
            double next_four_turn = InteractionEnergy(amino_acids, i, i + n);

            if (four_turn < energy_threshold && next_four_turn < energy_threshold)
            {
                AssignSecondaryStructure(amino_acid, "AlphaHelix");
                structures_assigned += 1;
            }
        }

        return structures_assigned;
    }

    /**
     * Assigns beta sheets as secondary structure to appropriate amino acids in specified list of amino acids
     * @param amino_acids list of all amino acids in chain
     * @param energy_threshold any energy below threshold supposed to form a hydrogen bond
     * @return amount of helices assigned
     * @throws StructureException
     */
    private static int AssignBetaSheets(List<AminoAcid> amino_acids, double energy_threshold) throws StructureException
    {
        int structures_assigned = 0;
        for (int i = 1; i < amino_acids.size() - 1; ++i)
        {
            for (int j = i + 2; j < amino_acids.size() - 1; ++j)
            {
                if (TryAssignParallelBetaSheet(amino_acids, i, j, energy_threshold))
                {
                    structures_assigned += 2;
                }

                if (TryAssignAntiparallelBetaSheet(amino_acids, i, j, energy_threshold))
                {
                    structures_assigned += 2;
                }
            }
        }

        return structures_assigned;
    }

    /**
     * Assigns parallel beta sheet secondary structure between specified amino acids if they conform dssp criteria for it
     * @param amino_acids list of all amino acids in chain
     * @param i first amino acid
     * @param j second amino acid
     * @param energy_threshold any energy below threshold supposed to form a hydrogen bond
     * @return true if specified amino acids form parallel beta sheet
     * @throws StructureException
     */
    public static boolean TryAssignParallelBetaSheet(List<AminoAcid> amino_acids, int i, int j, double energy_threshold) throws StructureException
    {
        double interaction11 = InteractionEnergy(amino_acids, i - 1, j);
        double interaction12 = InteractionEnergy(amino_acids, j, i + 1);
        double interaction21 = InteractionEnergy(amino_acids, j - 1, i);
        double interaction22 = InteractionEnergy(amino_acids, i, j + 1);

        if (interaction11 < energy_threshold && interaction12 < energy_threshold ||
                interaction21 < energy_threshold && interaction22 < energy_threshold)
        {
            String structure_name = "ParallelBetaSheet";
            AminoAcid amino_i = amino_acids.get(i);
            AssignSecondaryStructure(amino_i, structure_name);
            AminoAcid amino_j = amino_acids.get(j);
            AssignSecondaryStructure(amino_j, structure_name);
            return true;
        }

        return false;
    }

    /**
     * Assigns antiparallel beta sheet secondary structure between specified amino acids if they conform dssp criteria for it
     * @param amino_acids list of all amino acids in chain
     * @param i first amino acid
     * @param j second amino acid
     * @param energy_threshold any energy below threshold supposed to form a hydrogen bond
     * @return true if specified amino acids form parallel beta sheet
     * @throws StructureException
     */
    public static boolean TryAssignAntiparallelBetaSheet(List<AminoAcid> amino_acids, int i, int j, double energy_threshold) throws
            StructureException
    {
        double interaction11 = InteractionEnergy(amino_acids, i, j);
        double interaction12 = InteractionEnergy(amino_acids, j, i);
        double interaction21 = InteractionEnergy(amino_acids, i - 1, j + 1);
        double interaction22 = InteractionEnergy(amino_acids, j - 1, i + 1);

        if (interaction11 < energy_threshold && interaction12 < energy_threshold ||
                interaction21 < energy_threshold && interaction22 < energy_threshold)
        {
            String structure_name = "AntiparallelBetaSheet";
            AminoAcid amino_i = amino_acids.get(i);
            AssignSecondaryStructure(amino_i, structure_name);
            AminoAcid amino_j = amino_acids.get(j);
            AssignSecondaryStructure(amino_j, structure_name);
            return true;
        }

        return false;
    }

    /**
     * Assigns specified secondary structure to specified amino acid
     * @param amino_acid amino acid to assign secondary structure
     * @param structure name of structure
     */
    public static void AssignSecondaryStructure(AminoAcid amino_acid, String structure)
    {
        Map<String, String> secondary_structure = amino_acid.getSecStruc();
        secondary_structure.put(structure, "");
        amino_acid.setSecStruc(secondary_structure);
    }

    /**
     * Assigns beta sheets as secondary structure to appropriate amino acids in specified structure
     * @param structure structure for assignment
     * @param energy_threshold any energy below threshold supposed to form a hydrogen bond
     * @return amount of beta sheets assigned
     * @throws StructureException
     */
    @SuppressWarnings ("unchecked")
    public static int AssignBetaSheetsToStructure(Structure structure, double energy_threshold) throws StructureException
    {
        int structures_assigned = 0;
        for (Chain chain : structure.getChains())
        {
            List<Group> groups = chain.getAtomGroups("amino");
            List<AminoAcid> amino_acids = (List<AminoAcid>) (List<?>) groups;
            structures_assigned += AssignBetaSheets(amino_acids, energy_threshold);
        }

        return structures_assigned;
    }
}
