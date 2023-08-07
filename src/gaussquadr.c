#include <float.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

/* The code below is an adaptation to C from gaussq.f in Netlib:
 * https://www.netlib.org/go/gaussq.f */

#define EV_MAX_ITERS 30 /* used in `gausq2` */

static double class (int kind, int n, double alpha, double beta,
    double *b, double *a) {
  int i;
  double abi, muzero;

  switch (kind) {
    /* Legendre polynomials p(x) on (-1, +1), w(x) = 1 */
    case 1: {
      muzero = 2.;
      for (i = 1; i < n; i++) {
        a[i - 1] = 0.;
        abi = i;
        b[i - 1] = abi / sqrt(4 * abi * abi - 1.);
      }
      a[n - 1] = 0.;
      break;
    }
    /* Chebyshev polynomials of the first kind t(x) on (-1, +1),
     * w(x) = 1 / sqrt(1 - x * x) */
    case 2: {
      muzero = M_PI;
      for (i = 1; i < n; i++) {
        a[i - 1] = 0.;
        b[i - 1] = .5;
      }
      b[0] = M_SQRT1_2;
      a[n - 1] = 0.;
      break;
    }
    /* Chebyshev polynomials of the second kind u(x) on (-1, +1),
     * w(x) = sqrt(1 - x * x) */
    case 3: {
      muzero = M_PI_2;
      for (i = 1; i < n; i++) {
        a[i - 1] = 0.;
        b[i - 1] = .5;
      }
      a[n - 1] = 0.;
      break;
    }
    /* Hermite polynomials h(x) on (-inf, +inf), w(x) = exp(-x * x) */
    case 4: {
      muzero = sqrt(M_PI);
      for (i = 1; i < n; i++) {
        a[i - 1] = 0.;
        b[i - 1] = sqrt(i / 2.);
      }
      a[n - 1] = 0.;
      break;
    }
    /* Jacobi polynomials p(alpha, beta)(x) on (-1, +1),
     * w(x) = (1 - x) ^ alpha + (1 + x) ^ beta, alpha, beta > -1 */
    case 5: {
      double a2b2, ab;
      ab = alpha + beta;
      abi = 2. + ab;
      muzero = exp((ab + 1.) * M_LN2 + lbeta(alpha + 1., beta + 1.));
      /* muzero = pow(2., ab + 1.) * tgamma(alpha + 1.) * tgamma(beta + 1.) /
        tgamma(abi); */
      a[0] = (beta - alpha) / abi;
      b[0] = sqrt(4. * (1. + alpha) * (1. + beta) / ((abi + 1.) * abi * abi));
      a2b2 = beta * beta - alpha * alpha;
      for (i = 2; i < n; i++) {
        abi = 2. * i + ab;
        a[i - 1] = a2b2 / ((abi - 2.) * abi);
        b[i - 1] = sqrt(4. * i * (i + alpha) * (i + beta) * (i + ab) /
            ((abi * abi - 1.) * abi * abi));
      }
      abi = 2. * n + ab;
      a[n - 1] = a2b2 / ((abi - 2.) * abi);
      break;
    }
    /* Laguerre polynomials l(alpha)(x) on (0, +inf),
     * w(x) = x ^ alpha * exp(-x), alpha > -1 */
    case 6: {
      muzero = tgamma(alpha + 1.);
      for (i = 1; i < n; i++) {
        a[i - 1] = 2. * i - 1. + alpha;
        b[i - 1] = sqrt(i * (i + alpha));
      }
      a[n - 1] = 2. * n - 1 + alpha;
      break;
    }
    /* Gaussian density: scaled "probabilist" Hermite on (-inf, inf),
     * w(x) = exp(-x^2 / 2) / sqrt(2 * pi) */
    case 7: {
      muzero = 1.;
      for (i = 1; i < n; i++) {
        a[i - 1] = 0.;
        b[i - 1] = sqrt(i);
      }
      a[n - 1] = 0.;
      break;
    }
  }
  return muzero; /* suppress warnings */
}

static double solve (double shift, int n, double *a, double *b) {
  int i;
  double alpha = a[0] - shift;
  for (i = 2; i < n; i++)
    alpha = a[i - 1] - shift - b[i - 2] * b[i - 2] / alpha;
  return 1. / alpha;
}

static int gausq2 (int n, double *d, double *e, double *z) {
  int i, j, k, l, m, ii;
  double b, c, f, g, p, r, s;

  if (n == 1) return 0;
  e[n - 1] = 0.;
  for (l = 1; l <= n; l++) {
    j = 0;
    for (;;) {
      /* look for small sub-diagonal element */
      for (m = l; m <= n; m++) {
        if (m == n) break;
        if (fabs(e[m - 1]) <= DBL_EPSILON * (fabs(d[m - 1]) + fabs(d[m])))
          break;
      }

      p = d[l - 1];
      if (m == l) break;
      if (j == EV_MAX_ITERS) return l; /* l-th eigenvalue not converged */
      j++;
      /* form shift */
      g = (d[l] - p) / (2. * e[l - 1]);
      r = sqrt(g * g + 1.);
      g = d[m - 1] - p + e[l - 1] / (g + copysign(r, g));
      s = 1.;
      c = 1.;
      p = 0.;

      for (ii = 1; ii <= m - l; ii++) {
        i = m - ii;
        f = s * e[i - 1];
        b = c * e[i - 1];
        if (fabs(f) < fabs(g)) {
          s = f / g;
          r = sqrt(s * s + 1.);
          e[i] = g * r;
          c = 1. / r;
          s *= c;
        } else {
          c = g / f;
          r = sqrt(c * c + 1.);
          e[i] = f * r;
          s = 1. / r;
          c *= s;
        }
        g = d[i] - p;
        r = (d[i - 1] - g) * s + 2. * c * b;
        p = s * r;
        d[i] = g + p;
        g = c * r - b;
        /* form first component of vector */
        f = z[i];
        z[i] = s * z[i - 1] + c * f;
        z[i - 1] = c * z[i - 1] - s * f;
      }
      d[l - 1] -= p;
      e[l - 1] = g;
      e[m - 1] = 0.;
    }
  }

  /* order eigenvalues and eigenvectors */
  for (ii = 2; ii <= n; ii++) {
    i = ii - 1;
    k = i;
    p = d[i - 1];

    for (j = ii; j <= n; j++) {
      if (d[j - 1] >= p) continue;
      k = j;
      p = d[j - 1];
    }

    if (k == i) continue;
    d[k - 1] = d[i - 1];
    d[i - 1] = p;
    p = z[i - 1];
    z[i - 1] = z[k - 1];
    z[k - 1] = p;
  }

  return 0;
}

static int gaussq (int kind, int n, double alpha, double beta, int kpts,
    double *endpts, double *b, double *t, double *w) {
  int i, ierr;
  double muzero = class(kind, n, alpha, beta, b, t);

  if (kpts == 1) /* only t(n) must be changed */
    t[n - 1] = solve(endpts[0], n, t, b) * b[n - 2] * b[n - 2] + endpts[0];
  if (kpts == 2) { /* t(n) and b(n - 1) must be recomputed */
    double t1, gam;
    gam = solve(endpts[0], n, t, b);
    t1 = ((endpts[0] - endpts[1]) / (solve(endpts[1], n, t, b) - gam));
    b[n - 2] = sqrt(t1);
    t[n - 1] = endpts[0] + gam * t1;
  }

  w[0] = 1.;
  for (i = 2; i <= n; i++) w[i - 1] = 0.;
  ierr = gausq2(n, t, b, w);
  for (i = 0; i < n; i++) w[i] = muzero * w[i] * w[i];
  return ierr;
}

SEXP gauss_quad (SEXP skind, SEXP salpha, SEXP sbeta, SEXP skpts,
    SEXP sendpts, SEXP sb, SEXP st, SEXP sw) {
  return ScalarInteger(gaussq(asInteger(skind), length(sb),
        asReal(salpha), asReal(sbeta), asInteger(skpts), REAL(sendpts),
        REAL(sb), REAL(st), REAL(sw)));
}

static const R_CallMethodDef callMethods[] = {
  {"gauss_quad", (DL_FUNC) &gauss_quad, 8},
  {NULL, NULL, 0} /* guard */
};

void R_init_gaussquadr (DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

