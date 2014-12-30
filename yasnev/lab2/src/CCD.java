import common.PdbAtom;
import common.Vector;

import java.lang.invoke.VolatileCallSite;
import java.util.ArrayList;

/**
 * Author: Oleg Yasnev (oyasnev@gmail.com)
 * Date: 29.12.2014
 */
public class CCD {
    double distEps, equalEps;
    Integer[] boneAtomsIndices;

    public CCD(PdbAtom[] atoms, Vector target, double distEps, double equalEps, int maxIterationCnt) {
        this.distEps = distEps;
        this.equalEps = equalEps;
        getBoneAtomsIndices(atoms);
        int iter = 1;
        while (true) {
            System.out.printf("Iteration %d\n", iter);

            PdbAtom lastAtom = atoms[boneAtomsIndices[boneAtomsIndices.length - 1]];
            for (int i = 0; i < boneAtomsIndices.length - 1; i++) {
                Vector cur = new Vector(lastAtom);
                if (Vector.equal(cur, target, equalEps)) {
                    break;
                }
                PdbAtom iAtom = atoms[boneAtomsIndices[i]];
                Vector theta = getRotationAxis(new Vector(iAtom),
                        new Vector(atoms[boneAtomsIndices[i + 1]]), cur, target);
                if (theta.length() < equalEps) {
                    continue;
                }

                Vector s = computeSCoeffs(theta, new Vector(iAtom), cur, target);
                if (Math.abs(s.y * s.y + s.z * s.z) < equalEps) {
                    continue;
                }
                double cosT = s.y / Math.sqrt(s.y * s.y + s.z * s.z);
                double sinT = s.z / Math.sqrt(s.y * s.y + s.z * s.z);
                for (int j = boneAtomsIndices[i] + 1; j < atoms.length; j++) {
                    PdbAtom jAtom = atoms[j];
                    Vector v = rotate(theta, new Vector(iAtom),
                            new Vector(jAtom), cosT, sinT);
                    jAtom.x = v.x;
                    jAtom.y = v.y;
                    jAtom.z = v.z;
                }
            }

            double dist = (new Vector(new Vector(lastAtom), target)).length();
            System.out.printf("After iteration %d, distance to target: %f\n", iter, dist);
            if (dist < distEps) {
                System.out.println("Distanse < epsilon => finish");
                break;
            }

            if (iter >= maxIterationCnt) {
                System.out.println("Number of iteration > maximum iteration count => finish");
                break;
            }
            iter++;
        }
    }

    public Vector getRotationAxis(Vector point1, Vector point2, Vector cur, Vector target) {
        Vector v1 = new Vector(point1, cur);
        Vector v2 = new Vector(point1, point2);
        Vector theta = v2.unit();
        Vector o = Vector.sum(point1, v2.multByNum(Vector.dotProduct(v1, v2) / Vector.dotProduct(v2, v2)));
        Vector ortho = new Vector(o, cur);
        double ri = ortho.length();
        if (ri < equalEps) {
            Vector axis = Vector.crossProduct(new Vector(o, target), theta);
            if (axis.length() < equalEps || Vector.crossProduct(new Vector(o, cur), theta).length() < equalEps) {
                return new Vector(0, 0, 0);
            }
            return axis.unit();
        }
        Vector r = ortho.unit();
        Vector s = Vector.crossProduct(r, theta);
        Vector fi = new Vector(o, target);
        if (Math.abs(Vector.dotProduct(fi, s)) > equalEps) {
            return theta;
        }
        return s;
    }

    public Vector rotate(Vector theta, Vector point1, Vector c, double cosT, double sinT) {
        Vector v1 = new Vector(point1, c);
        Vector o = Vector.sum(point1, theta.multByNum(Vector.dotProduct(v1, theta)));
        Vector ortho = new Vector(o, c);
        double ri = ortho.length();
        if (ri < equalEps) {
            return c;
        }
        Vector r = ortho.unit();
        Vector s = Vector.crossProduct(r, theta);
        return Vector.sum(o, Vector.sum(s.multByNum(sinT * ri), ortho.multByNum(cosT)));
    }

    public Vector computeSCoeffs(Vector theta, Vector point1, Vector cur, Vector target) {
        Vector v1 = new Vector(point1, cur);
        Vector o = Vector.sum(point1, theta.multByNum(Vector.dotProduct(v1, theta)));
        Vector ri = new Vector(o, cur);
        Vector fi = new Vector(o, target);
        Vector si = Vector.crossProduct(ri.unit(), theta);

        if (ri.length() < equalEps) {
            return new Vector(0, 0, 0);
        }

        Vector s = new Vector();
        s.x = Vector.dotProduct(ri, ri) + Vector.dotProduct(fi, fi);
        s.y = Vector.dotProduct(ri, fi) * 2;
        s.z = Vector.dotProduct(si, fi) * 2 * Math.sqrt(Vector.dotProduct(ri, ri));
        return s;
    }

    public void getBoneAtomsIndices(PdbAtom[] atoms) {
        ArrayList<Integer> arrayList = new ArrayList<Integer>();
        for (int i = 0; i < atoms.length; i++) {
            String type = atoms[i].type;
            if (type.equals("N") || type.equals("CA") || type.equals("C")) {
                arrayList.add(i);
            }
        }
        boneAtomsIndices = new Integer[arrayList.size()];
        arrayList.toArray(boneAtomsIndices);
    }

}
