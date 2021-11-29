import java.util.Arrays;

/**
 * Die Klasse CubicSpline bietet eine Implementierung der kubischen Splines. Sie
 * dient uns zur effizienten Interpolation von aequidistanten Stuetzpunkten.
 *
 * @author braeckle
 */
public class CubicSpline implements InterpolationMethod {

    /**
     * linke und rechte Intervallgrenze x[0] bzw. x[n]
     */
    double a, b;

    /**
     * Anzahl an Intervallen
     */
    int n;

    /**
     * Intervallbreite
     */
    double h;

    /**
     * Stuetzwerte an den aequidistanten Stuetzstellen
     */
    double[] y;

    /**
     * zu berechnende Ableitunge an den Stuetzstellen
     */
    double yprime[];

    /**
     * {@inheritDoc} Zusaetzlich werden die Ableitungen der stueckweisen
     * Polynome an den Stuetzstellen berechnet. Als Randbedingungen setzten wir
     * die Ableitungen an den Stellen x[0] und x[n] = 0.
     */
    @Override
    public void init(double a, double b, int n, double[] y) {
        this.a = a;
        this.b = b;
        this.n = n;
        h = ((double) b - a) / (n);

        this.y = Arrays.copyOf(y, n + 1);

        /* Randbedingungen setzten */
        yprime = new double[n + 1];
        yprime[0] = 0;
        yprime[n] = 0;

        /* Ableitungen berechnen. Nur noetig, wenn n > 1 */
        if (n > 1) {
            computeDerivatives();
        }
    }

    /**
     * getDerivatives gibt die Ableitungen yprime zurueck
     */
    public double[] getDerivatives() {
        return yprime;
    }

    /**
     * Setzt die Ableitungen an den Raendern x[0] und x[n] neu auf yprime0 bzw.
     * yprimen. Anschliessend werden alle Ableitungen aktualisiert.
     */
    public void setBoundaryConditions(double yprime0, double yprimen) {
        yprime[0] = yprime0;
        yprime[n] = yprimen;
        if (n > 1) {
            computeDerivatives();
        }
    }

    /**
     * Berechnet die Ableitungen der stueckweisen kubischen Polynome an den
     * einzelnen Stuetzstellen. Dazu wird ein lineares System Ax=c mit einer
     * Tridiagonalen Matrix A und der rechten Seite c aufgebaut und geloest.
     * Anschliessend sind die berechneten Ableitungen y1' bis yn-1' in der
     * Membervariable yprime gespeichert.
     * <p>
     * Zum Zeitpunkt des Aufrufs stehen die Randbedingungen in yprime[0] und yprime[n].
     * Speziell bei den "kleinen" Faellen mit Intervallzahlen n = 2
     * oder 3 muss auf die Struktur des Gleichungssystems geachtet werden. Der
     * Fall n = 1 wird hier nicht beachtet, da dann keine weiteren Ableitungen
     * berechnet werden muessen.
     */
    public void computeDerivatives() {
        double[] lower = new double[n - 2];
        double[] diag = new double[n - 1];
        double[] upper = new double[n - 2];
        double[] c = new double[n - 1];

        if (n == 2) {
            c[0] = (3 / h) * (y[2] - y[0]) - yprime[0];
            c[1] = (3 / h) * (y[n] - y[n-2]) - yprime[n];
        } else if (n > 2) {
            c[0] = (3 / h) * (y[2] - y[0] - ((h / 3) * yprime[0]));
            c[n - 2] = (3 / h) * (y[n] - y[n - 2] - ((h / 3) * yprime[n]));
            for (int i = 0; i < n - 2; i++) {
                lower[i] = 1;
                diag[i] = 4;
                upper[i] = 1;
            }
            for (int i = 3; i < n; i++) {
                c[i - 2] = 3 / h * (y[i] - y[i - 2]);
            }
        }
        diag[n - 2] = 4;
        TridiagonalMatrix matrix = new TridiagonalMatrix(lower, diag, upper);
        double[] helper = matrix.solveLinearSystem(c);
        for (int i = 0; i < n - 1; i++) {
            yprime[i + 1] = helper[i];
        }
    }

    /**
     * {@inheritDoc} Liegt z ausserhalb der Stuetzgrenzen, werden die
     * aeussersten Werte y[0] bzw. y[n] zurueckgegeben. Liegt z zwischen den
     * Stuetzstellen x_i und x_i+1, wird z in das Intervall [0,1] transformiert
     * und das entsprechende kubische Hermite-Polynom ausgewertet.
     */
    @Override
    public double evaluate(double z) {
        double result = 0;
        if (z > y[0])
            return y[0];
        else if (z < y[y.length - 1])
            return y[y.length - 1];
        else {
            for (int i = 0; i < y.length - 1; i++) {
                if (y[i] < z && y[i + 1] > z) {
                    result = (z - y[i]) / (y[i + 1] - y[i]);
                    result = hermite(result, i);
                    break;
                }
            }
        }
        return result;
    }

    private double hermite(double t, int i) {
        return y[i] * h0(t) + y[i + 1] * h1(t) + h * yprime[i] * h2(t) + h * yprime[i + 1] * h3(t);
    }

    private double h3(double t) {
        return -Math.pow(t, 2) + Math.pow(t, 3);
    }

    private double h2(double t) {
        return t - 2 * Math.pow(t, 2) + Math.pow(t, 3);
    }

    private double h1(double t) {
        return 3 * Math.pow(t, 2) + 2 * Math.pow(t, 3);
    }

    private double h0(double t) {
        return 1 - 3 * Math.pow(t, 2) + 2 * Math.pow(t, 3);
    }
}
