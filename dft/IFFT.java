package dft;

import java.util.Arrays;

/**
 * Schnelle inverse Fourier-Transformation
 *
 * @author Sebastian Rettenberger
 */
public class IFFT {
    /**
     * Schnelle inverse Fourier-Transformation (IFFT).
     * <p>
     * Die Funktion nimmt an, dass die Laenge des Arrays c immer eine
     * Zweierpotenz ist. Es gilt also: c.length == 2^m fuer ein beliebiges m.
     */
    public static Complex[] ifft(Complex[] c) {
        int n = c.length;
        Complex[] v = new Complex[n];

        if (n == 1) {
            v[0] = c[0];
        } else {
            int m = n / 2;
            Complex[] z1 = new Complex[n / 2];
            for (int i = 0; i < z1.length; i++) {
                z1[i] = c[i * 2];
            }
            z1 = ifft(z1);
            Complex[] z2 = new Complex[n / 2];
            for (int i = 0; i < z2.length; i++) {
                z2[i] = c[i * 2 + 1];
            }
            z2 = ifft(z2);
            Complex omega = Complex.fromPolar(1, 2 * Math.PI / n);
            for (int j = 0; j < m; j++) {
                v[j] = omega.power(j).mul(z2[j]).add(z1[j]);
                v[m + j] = z1[j].sub(omega.power(j).mul(z2[j]));
            }
        }
        return v;
    }
}
