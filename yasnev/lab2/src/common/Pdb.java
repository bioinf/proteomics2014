package common;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;

/**
 * Author: Oleg Yasnev (oyasnev@gmail.com)
 * Date: 27.12.2014
 */
public class Pdb {

    public PdbAtom[] atoms;

    protected String inputFilename;


    public Pdb(String filename) {
        inputFilename = filename;
        readFromFile();
    }

    public void writeToFile(String filename) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(inputFilename));
            PrintWriter pw = new PrintWriter(filename);
            int atomCnt = 0;
            for (String line; (line = br.readLine()) != null; ) {
                if (line.startsWith("ATOM")) {
                    pw.println(atoms[atomCnt].toString());
                    atomCnt++;
                } else {
                    pw.println(line);
                }
            }
            pw.close();
            br.close();
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
    }

    protected void readFromFile() {
        ArrayList<PdbAtom> atoms = new ArrayList<PdbAtom>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(inputFilename));
            for (String line; (line = br.readLine()) != null; ) {
                if (line.startsWith("ATOM")) {
                   atoms.add(new PdbAtom(line));
                }
            }
            br.close();
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
        this.atoms = new PdbAtom[atoms.size()];
        atoms.toArray(this.atoms);
    }
}
