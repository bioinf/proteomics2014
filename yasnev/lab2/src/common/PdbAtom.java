package common;

import java.util.Locale;

/**
 * Author: Oleg Yasnev (oyasnev@gmail.com)
 * Date: 27.12.2014
 */
public class PdbAtom {

    public String start;
    public String type;
    public String middle;
    public double x, y, z;
    public String end;

    public PdbAtom(String line) {
        start = line.substring(0, 13);
        type = line.substring(13, 16).trim();
        middle = line.substring(16, 26);
        String[] strs = line.substring(26, 54).trim().split("\\s+");
        x = Double.parseDouble(strs[0]);
        y = Double.parseDouble(strs[1]);
        z = Double.parseDouble(strs[2]);
        end = line.substring(55);
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(start);
        sb.append(String.format("%-3s", type));
        sb.append(middle);
        sb.append(String.format(Locale.US, "%12.3f", x));
        sb.append(String.format(Locale.US, "%8.3f", y));
        sb.append(String.format(Locale.US, "%8.3f ", z));
        sb.append(end);
        return sb.toString();
    }
}
