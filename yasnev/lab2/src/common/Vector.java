package common;

/**
 * Author: Oleg Yasnev (oyasnev@gmail.com)
 * Date: 30.12.2014
 */
public class Vector {
    public double x, y, z;

    public Vector() {}

    public Vector(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public Vector(Vector v1, Vector v2) {
        this.x = v2.x - v1.x;
        this.y = v2.y - v1.y;
        this.z = v2.z - v1.z;
    }

    public Vector(PdbAtom atom) {
        x = atom.x;
        y = atom.y;
        z = atom.z;
    }

    public static Vector sum(Vector v1, Vector v2) {
        return new Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
    }

        public static double dotProduct(Vector v1, Vector v2) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    public static Vector crossProduct(Vector v1, Vector v2) {
        return new Vector(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
    }

    public static boolean equal(Vector v1, Vector v2, double eps) {
        return Math.abs(v1.x - v2.x) < eps && Math.abs(v1.y - v2.y) < eps && Math.abs(v1.z - v2.z) < eps;
    }

    public double length() {
        return Math.sqrt(dotProduct(this, this));
    }

    public Vector multByNum(double num) {
        return new Vector(this.x * num, this.y * num, this.z * num);
    }

    public Vector unit() {
        return multByNum(1.0 / length());
    }
}
