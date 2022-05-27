/*
 * Constants copied from DGAM1 in the NSWC library.
 */
/** The constant {@code A0} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_A0 = 0.611609510448141581788e-8;
/** The constant {@code A1} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_A1 = 0.62473083011646551621e-8;
/** The constant {@code B1} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_B1 = 0.2036104140668069873;
/** The constant {@code B2} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_B2 = 0.266205348428949217746e-1;
/** The constant {@code B3} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_B3 = 0.493944979382446875238e-3;
/** The constant {@code B4} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_B4 = -0.851419432440314906588e-5;
/** The constant {@code B5} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_B5 = -0.643045481779353022248e-5;
/** The constant {@code B6} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_B6 = 0.992641840672773722196e-6;
/** The constant {@code B7} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_B7 = -0.607761895722825260739e-7;
/** The constant {@code B8} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_B8 = 0.195755836614639731882e-9;
/** The constant {@code P0} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_P0 = 0.6116095104481415817861e-8;
/** The constant {@code P1} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_P1 = 0.6871674113067198736152e-8;
/** The constant {@code P2} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_P2 = 0.6820161668496170657918e-9;
/** The constant {@code P3} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_P3 = 0.468684332294884803108e-10;
/** The constant {@code P4} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_P4 = 0.1572833027710446286995e-11;
/** The constant {@code P5} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_P5 = -0.1249441572276366213222e-12;
/** The constant {@code P6} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_P6 = 0.4343529937408594255178e-14;
/** The constant {@code Q1} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_Q1 = 0.3056961078365221025009;
/** The constant {@code Q2} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_Q2 = 0.5464213086042296536016e-1;
/** The constant {@code Q3} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_Q3 = 0.495683009382588731202e-2;
/** The constant {@code Q4} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_Q4 = 0.2692369466186361192876e-3;
/** The constant {@code C} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C = -0.422784335098467139393487909917598;
/** The constant {@code C0} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C0 = 0.577215664901532860606512090082402;
/** The constant {@code C1} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C1 = -0.65587807152025388107701951514539;
/** The constant {@code C2} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C2 = -0.420026350340952355290039348754298e-1;
/** The constant {@code C3} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C3 = 0.166538611382291489501700795102105;
/** The constant {@code C4} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C4 = -0.421977345555443367482083012891874e-1;
/** The constant {@code C5} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C5 = -0.96219715278769735621149216723482e-2;
/** The constant {@code C6} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C6 = 0.721894324666309954239501034044657e-2;
/** The constant {@code C7} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C7 = -0.116516759185906511211397108401839e-2;
/** The constant {@code C8} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C8 = -0.215241674114950972815729963053648e-3;
/** The constant {@code C9} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C9 = 0.128050282388116186153198626328164e-3;
/** The constant {@code C10} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C10 = -0.201348547807882386556893914210218e-4;
/** The constant {@code C11} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C11 = -0.125049348214267065734535947383309e-5;
/** The constant {@code C12} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C12 = 0.113302723198169588237412962033074e-5;
/** The constant {@code C13} defined in {@code DGAM1}. */
const INV_GAMMA1P_M1_C13 = -0.205633841697760710345015413002057e-6;
const HALF_LOG_2_PI = 0.5 * Math.log(2.0 * Math.PI);
const LANCZOS_G = 607 / 128;
const LANCZOS = [
  0.99999999999999709182, 57.156235665862923517, -59.597960355475491248,
  14.136097974741747174, -0.49191381609762019978, 0.33994649984811888699e-4,
  0.46523628927048575665e-4, -0.98374475304879564677e-4,
  0.15808870322491248884e-3, -0.21026444172410488319e-3,
  0.2174396181152126432e-3, -0.16431810653676389022e-3,
  0.84418223983852743293e-4, -0.2619083840158140867e-4,
  0.36899182659531622704e-5,
];

const SMALL = 1e-50;

export function regularizedGammaQValue(
  a,
  x,
  epsilon = 1e-15,
  maxIterations = Number.MAX_VALUE
) {
  if (Number.isNaN(a) || Number.isNaN(x) || a <= 0 || x < 0) {
    return Number.NaN;
  } else if (x == 0) {
    return 1;
  } else if (x < a + 1) {
    // P should converge faster in this case.
    return 1 - regularizedGammaPValue(a, x, epsilon, maxIterations);
  } else {
    const getA = function (n, x) {
      return n * (a - n);
    };
    const getB = function (n, x) {
      return 2 * n + 1 - a + x;
    };
    return (
      Math.exp(-x + a * Math.log(x) - logGamma(a)) /
      continuedFraction(x, epsilon, maxIterations, getA, getB)
    );
  }
}

function invGamma1pm1(x) {
  if (x < -0.5 || x > 1.5) {
    throw new Error('Out of range');
  }

  const t = x <= 0.5 ? x : x - 0.5 - 0.5;
  if (t < 0) {
    const a = INV_GAMMA1P_M1_A0 + t * INV_GAMMA1P_M1_A1;
    let b = INV_GAMMA1P_M1_B8;
    b = INV_GAMMA1P_M1_B7 + t * b;
    b = INV_GAMMA1P_M1_B6 + t * b;
    b = INV_GAMMA1P_M1_B5 + t * b;
    b = INV_GAMMA1P_M1_B4 + t * b;
    b = INV_GAMMA1P_M1_B3 + t * b;
    b = INV_GAMMA1P_M1_B2 + t * b;
    b = INV_GAMMA1P_M1_B1 + t * b;
    b = 1.0 + t * b;

    let c = INV_GAMMA1P_M1_C13 + t * (a / b);
    c = INV_GAMMA1P_M1_C12 + t * c;
    c = INV_GAMMA1P_M1_C11 + t * c;
    c = INV_GAMMA1P_M1_C10 + t * c;
    c = INV_GAMMA1P_M1_C9 + t * c;
    c = INV_GAMMA1P_M1_C8 + t * c;
    c = INV_GAMMA1P_M1_C7 + t * c;
    c = INV_GAMMA1P_M1_C6 + t * c;
    c = INV_GAMMA1P_M1_C5 + t * c;
    c = INV_GAMMA1P_M1_C4 + t * c;
    c = INV_GAMMA1P_M1_C3 + t * c;
    c = INV_GAMMA1P_M1_C2 + t * c;
    c = INV_GAMMA1P_M1_C1 + t * c;
    c = INV_GAMMA1P_M1_C + t * c;
    if (x > 0.5) {
      return (t * c) / x;
    } else {
      return x * (c + 0.5 + 0.5);
    }
  } else {
    let p = INV_GAMMA1P_M1_P6;
    p = INV_GAMMA1P_M1_P5 + t * p;
    p = INV_GAMMA1P_M1_P4 + t * p;
    p = INV_GAMMA1P_M1_P3 + t * p;
    p = INV_GAMMA1P_M1_P2 + t * p;
    p = INV_GAMMA1P_M1_P1 + t * p;
    p = INV_GAMMA1P_M1_P0 + t * p;

    let q = INV_GAMMA1P_M1_Q4;
    q = INV_GAMMA1P_M1_Q3 + t * q;
    q = INV_GAMMA1P_M1_Q2 + t * q;
    q = INV_GAMMA1P_M1_Q1 + t * q;
    q = 1.0 + t * q;

    let c = INV_GAMMA1P_M1_C13 + (p / q) * t;
    c = INV_GAMMA1P_M1_C12 + t * c;
    c = INV_GAMMA1P_M1_C11 + t * c;
    c = INV_GAMMA1P_M1_C10 + t * c;
    c = INV_GAMMA1P_M1_C9 + t * c;
    c = INV_GAMMA1P_M1_C8 + t * c;
    c = INV_GAMMA1P_M1_C7 + t * c;
    c = INV_GAMMA1P_M1_C6 + t * c;
    c = INV_GAMMA1P_M1_C5 + t * c;
    c = INV_GAMMA1P_M1_C4 + t * c;
    c = INV_GAMMA1P_M1_C3 + t * c;
    c = INV_GAMMA1P_M1_C2 + t * c;
    c = INV_GAMMA1P_M1_C1 + t * c;
    c = INV_GAMMA1P_M1_C0 + t * c;

    if (x > 0.5) {
      return (t / x) * (c - 0.5 - 0.5);
    } else {
      return x * c;
    }
  }
}

function logGamma1p(x) {
  if (x < -0.5 || x > 1.5) {
    throw new Error('Out of range');
  }

  return -Math.log1p(invGamma1pm1(x));
}

function lanczosApproximation(x) {
  let sum = 0;
  for (let i = LANCZOS.length - 1; i > 0; i--) {
    sum += LANCZOS[i] / (x + i);
  }
  return sum + LANCZOS[0];
}

function logGamma(x) {
  if (Number.isNaN(x) || x <= 0.0) {
    return Number.NaN;
  } else if (x < 0.5) {
    return logGamma1p(x) - Math.log(x);
  } else if (x <= 2.5) {
    return logGamma1p(x - 0.5 - 0.5);
  } else if (x <= 8.0) {
    let n = Math.floor(x - 1.5);
    let prod = 1.0;
    for (let i = 1; i <= n; i++) {
      prod *= x - i;
    }
    return logGamma1p(x - (n + 1)) + Math.log(prod);
  } else {
    let sum = lanczosApproximation(x);
    let tmp = x + LANCZOS_G + 0.5;
    return (x + 0.5) * Math.log(tmp) - tmp + HALF_LOG_2_PI + Math.log(sum / x);
  }
}

export function regularizedGammaPValue(
  a,
  x,
  epsilon = 1e-15,
  maxIterations = Number.MAX_VALUE
) {
  if (Number.isNaN(a) || Number.isNaN(x) || a <= 0 || x < 0) {
    return Number.NaN;
  } else if (x == 0) {
    return 0;
  } else if (x >= a + 1) {
    // Q should converge faster in this case.
    return 1 - regularizedGammaQValue(a, x, epsilon, maxIterations);
  } else {
    // Series.
    let n = 0; // current element index
    let an = 1 / a; // n-th element in the series
    let sum = an; // partial sum
    while (
      Math.abs(an / sum) > epsilon &&
      n < maxIterations &&
      sum < Number.POSITIVE_INFINITY
    ) {
      // compute next element in the series
      n += 1;
      an *= x / (a + n);

      // update partial sum
      sum += an;
    }
    if (n >= maxIterations) {
      throw new Error('Failed to converge');
    } else if (!Number.isFinite(sum)) {
      return 1;
    } else {
      // Ensure result is in the range [0, 1]
      let result = Math.exp(-x + a * Math.log(x) - logGamma(a)) * sum;
      return result > 1.0 ? 1.0 : result;
    }
  }
}

function continuedFraction(x, epsilon, maxIterations, getA, getB) {
  let hPrev = updateIfCloseToZero(getB(0, x));

  let n = 1;
  let dPrev = 0.0;
  let cPrev = hPrev;
  let hN;

  while (n <= maxIterations) {
    const a = getA(n, x);
    const b = getB(n, x);

    let dN = updateIfCloseToZero(b + a * dPrev);
    const cN = updateIfCloseToZero(b + a / cPrev);

    dN = 1 / dN;
    const deltaN = cN * dN;
    hN = hPrev * deltaN;

    if (!Number.isFinite(hN)) {
      throw new Error(
        'Continued fraction convergents diverged to +/- infinity'
      );
    }
    if (Number.isNaN(hN)) {
      throw new Error('Continued fraction diverged to NaN');
    }

    if (Math.abs(deltaN - 1) < epsilon) {
      return hN;
    }

    dPrev = dN;
    cPrev = cN;
    hPrev = hN;
    ++n;
  }

  throw new Error('maximal count exceeded');
}

function updateIfCloseToZero(value) {
  return Math.abs(value - 0.0 <= SMALL) ? SMALL : value;
}
