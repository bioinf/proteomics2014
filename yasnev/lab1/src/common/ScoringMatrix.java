package common;

import common.interfaces.IMatchScoringMatrix;
import edu.princeton.cs.introcs.In;

import java.util.HashMap;
import java.util.Map;

/**
 * Operates with scoring matrix
 *
 * Author: Oleg Yasnev
 */
public class ScoringMatrix implements IMatchScoringMatrix {
    /**
     * Scoring map. Key is two concatenated amino acids or empty string
     */
    protected Map<String, Integer> map;

    /**
     * A character symbol used as gap symbol
     */
    protected char gapSymbol;

    /**
     * Constructor.
     * Reads scoring matrix from a given file
     * @param filename file with scoring matrix
     * @param commentPrefix prefix string opening line-comment
     */
    public ScoringMatrix(String filename, String commentPrefix) {
        readFromFile(filename, commentPrefix);
    }

    public int getScore(char chFirst, char chSecond) {
        return map.get(String.valueOf(chFirst) + String.valueOf(chSecond));
    }

    public char getGapSymbol() {
        return gapSymbol;
    }

    public int getGapScore(char ch) {
        return map.get(String.valueOf(ch) + String.valueOf(gapSymbol));
    }

    /**
     * Reads scoring matrix from a given file
     * @param filename file with scoring matrix
     * @param commentPrefix prefix string opening line-comment
     */
    protected void readFromFile(String filename, String commentPrefix) {
        In in = new In(filename);

        String str;
        // Skip comments
        do {
            str = in.readLine();
        } while (str.startsWith(commentPrefix));

        // Parse header
        String[] symbols = str.trim().split("\\s+");

        map = new HashMap<String, Integer>(symbols.length * 2);
        // length is twice because AB and BA are two records in map

        // Read score and fill map
        while (!in.isEmpty()) {
            // read header of current line
            char curSymbol = in.readString().charAt(0);
            if (curSymbol < 'A' || curSymbol > 'Z') {
                gapSymbol = curSymbol;
            }
            // add scores for combination row + count
            for (String symbol : symbols) {
                int score = in.readInt();
                map.put(symbol + curSymbol, score);
            }
        }

        in.close();
    }
}
