/**
 * Computes the hypergeometric probability.
 */
function phyper(a, b, c, d) {
    return Math
        .exp((logFactorial(a + b)
            + logFactorial(c + d)
            + logFactorial(a + c) + logFactorial(b + d))
            - (logFactorial(a)
                + logFactorial(b)
                + logFactorial(c)
                + logFactorial(d) + logFactorial(a + b + c + d)));

};

const logFactorials = [0.00000000000000000,
    0.00000000000000000, 0.69314718055994531, 1.79175946922805500,
    3.17805383034794562, 4.78749174278204599, 6.57925121201010100,
    8.52516136106541430, 10.60460290274525023, 12.80182748008146961,
    15.10441257307551530, 17.50230784587388584, 19.98721449566188615,
    22.55216385312342289, 25.19122118273868150, 27.89927138384089157,
    30.67186010608067280, 33.50507345013688888, 36.39544520803305358,
    39.33988418719949404, 42.33561646075348503, 45.38013889847690803,
    48.47118135183522388, 51.60667556776437357, 54.78472939811231919,
    58.00360522298051994, 61.26170176100200198, 64.55753862700633106,
    67.88974313718153498, 71.25703896716800901];

function logFactorial(k) {
    if (k >= 30) { // stirlings approximation
        const C0 = 9.18938533204672742e-01;
        const C1 = 8.33333333333333333e-02;
        const C3 = -2.77777777777777778e-03;
        const C5 = 7.93650793650793651e-04;
        const C7 = -5.95238095238095238e-04;
        const r = 1.0 / k;
        const rr = r * r;
        return (k + 0.5) * Math.log(k) - k + C0 + r
            * (C1 + rr * (C3 + rr * (C5 + rr * C7)));
        // log k! = (k + 1/2)log(k) - k + (1/2)log(2Pi) + stirlingCorrection(k)
    }
    return logFactorials[k];
};

export function fisherTest(a, b, c, d) {
    // match R 2-sided fisher.test
    const p = phyper(a, b, c, d);
    let sum = p;
    for (let _a = 0, n = a + b + c + d; _a <= n; _a++) {
        const _b = a + b - _a;
        const _c = a + c - _a;
        const _d = b + d - _b;
        if (_a !== a && _b >= 0 && _c >= 0 && _d >= 0) {
            const _p = phyper(_a, _b, _c, _d);
            if (_p <= p) {
                sum += _p;
            }
        }
    }
    return Math.min(1, sum);
};
