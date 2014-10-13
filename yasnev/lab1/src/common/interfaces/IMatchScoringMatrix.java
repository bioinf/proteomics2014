package common.interfaces;

/**
 * Interface for match scoring matrix (like Blosum62)
 *
 * Author: Oleg Yasnev
 */
public interface IMatchScoringMatrix {
    /**
     * Get a score of substitution for two amino acids (aa) or empty string (-)
     * @param chFirst aa or -
     * @param chSecond aa or -
     * @return a score
     */
    public int getScore(char chFirst, char chSecond);

    /**
     * Get a character symbol used as gap symbol
     * @return a gap symbol
     */
    public char getGapSymbol();

    /**
     * Get a gap score for insertion or deletion
     * @param ch aa
     * @return a gap score
     */
    public int getGapScore(char ch);
}
