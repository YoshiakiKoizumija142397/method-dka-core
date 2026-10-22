/**
 * Core Mathematical Engine for Durand-Kerner-Aberth Algorithm
 * 
 * This document provides the complete, self-contained implementation of index.js for method-dka-core.
 * It operates in Node.js using decimal.js for arbitrary-precision complex arithmetic and Yun's square-free factor decomposition algorithm via BigInt arithmetic.
 * 
 * High-Precision Polynomial Root Finder using the Durand-Kerner-Aberth (DKA) Method
 * with Yun's Square-Free Decomposition preprocessing.
 * 
 * Requires: decimal.js (npm install decimal.js)
 */

const fs = require('fs');
const path = require('path');
const Decimal = require('decimal.js');

// Configurable Precision Hierarchy Standards
const PRECISION_PREPROCESS = 320;
const PRECISION_DKA = 260;

/**
 * Configure global Decimal settings to default precision
 */
Decimal.set({
    precision: PRECISION_DKA,
    rounding: Decimal.ROUND_HALF_UP
});

// ============================================================================
// 1. COMPLEX ARITHMETIC HELPER FUNCTIONS
// ============================================================================

function c_add(a, b) {
    return { re: a.re.plus(b.re), im: a.im.plus(b.im) };
}

function c_sub(a, b) {
    return { re: a.re.minus(b.re), im: a.im.minus(b.im) };
}

function c_mul(a, b) {
    return {
        re: a.re.times(b.re).minus(a.im.times(b.im)),
        im: a.re.times(b.im).plus(a.im.times(b.re))
    };
}

function c_div(a, b) {
    const denom = b.re.times(b.re).plus(b.im.times(b.im));
    if (denom.isZero()) {
        throw new Error("Complex division by zero.");
    }
    return {
        re: a.re.times(b.re).plus(a.im.times(b.im)).div(denom),
        im: a.im.times(b.re).minus(a.re.times(b.im)).div(denom)
    };
}

function c_abs(a) {
    return a.re.times(a.re).plus(a.im.times(a.im)).sqrt();
}

/**
 * Evaluates polynomial P(z) = sum_{k=0}^n a_k * z^k using Horner's Method
 */
function evalPoly(coeffs, z) {
    let result = { re: new Decimal(0), im: new Decimal(0) };
    for (let i = coeffs.length - 1; i >= 0; i--) {
        const coeffComp = { re: new Decimal(coeffs[i]), im: new Decimal(0) };
        result = c_add(c_mul(result, z), coeffComp);
    }
    return result;
}

/**
 * Calculates normalized residual: |P(z)| / sum(|a_k| * |z|^k)
 */
function normalizedResidual(coeffs, z) {
    const pVal = c_abs(evalPoly(coeffs, z));
    const zAbs = c_abs(z);

    let denom = new Decimal(0);
    let zPow = new Decimal(1);

    for (let i = 0; i < coeffs.length; i++) {
        const coeffAbs = new Decimal(coeffs[i]).abs();
        denom = denom.plus(coeffAbs.times(zPow));
        zPow = zPow.times(zAbs);
    }

    if (denom.isZero()) return new Decimal(0);
    return pVal.div(denom);
}

/**
 * Evaluates max relative discrepancy between expanded roots and original polynomial coefficients
 */
function coefficientReconstructionResidual(rootsWithMult, origCoeffs) {
    const deg = origCoeffs.length - 1;
    if (deg <= 0) return new Decimal(0);
    let poly = [{ re: new Decimal(1), im: new Decimal(0) }];
    for (const r of rootsWithMult) {
        const negR = { re: r.re.negated(), im: r.im.negated() };
        for (let m = 0; m < r.multiplicity; m++) {
            const nextPoly = new Array(poly.length + 1).fill(null).map(() => ({ re: new Decimal(0), im: new Decimal(0) }));
            for (let i = 0; i < poly.length; i++) {
                nextPoly[i + 1] = c_add(nextPoly[i + 1], poly[i]);
                nextPoly[i] = c_add(nextPoly[i], c_mul(poly[i], negR));
            }
            poly = nextPoly;
        }
    }
    const leadOrig = new Decimal(origCoeffs[deg]);
    let maxRelDiff = new Decimal(0);
    for (let i = 0; i <= deg; i++) {
        const scaledRe = poly[i].re.times(leadOrig);
        const origC = new Decimal(origCoeffs[i]);
        const diff = scaledRe.minus(origC).abs();
        const denom = Decimal.max(origC.abs(), new Decimal(1));
        const relDiff = diff.div(denom);
        if (relDiff.gt(maxRelDiff)) {
            maxRelDiff = relDiff;
        }
    }
    return maxRelDiff;
}

/**
 * Estimates upper bound scale of roots via Cauchy's Bound
 */
function estimateRootScale(coeffs) {
    const n = coeffs.length - 1;
    if (n <= 0) return new Decimal(1);

    const lead = new Decimal(coeffs[n]).abs();
    if (lead.isZero()) return new Decimal(1);
    let maxRatio = new Decimal(0);
    for (let i = 0; i < n; i++) {
        const ratio = new Decimal(coeffs[i]).abs().div(lead);
        if (ratio.gt(maxRatio)) {
            maxRatio = ratio;
        }
    }
    return new Decimal(1).plus(maxRatio);
}

// ============================================================================
// 2. DYNAMIC FLOATING-POINT POLYNOMIAL ALGEBRA & GCD FUNCTIONS
// ============================================================================

function trimLeadingZeros(poly) {
    let end = poly.length;
    while (end > 1 && poly[end - 1].isZero()) {
        end--;
    }
    return poly.slice(0, end);
}

function getPolyDegree(poly) {
    const trimmed = trimLeadingZeros(poly);
    return trimmed.length - 1;
}

function derivativePoly(poly) {
    const trimmed = trimLeadingZeros(poly);
    const deg = trimmed.length - 1;
    if (deg <= 0) return [new Decimal(0)];
    const deriv = [];
    for (let i = 1; i <= deg; i++) {
        deriv.push(trimmed[i].times(i));
    }
    return trimLeadingZeros(deriv);
}

function polyDivide(num, den) {
    let N = trimLeadingZeros(num);
    const D = trimLeadingZeros(den);
    const degN = N.length - 1;
    const degD = D.length - 1;
    if (degD < 0 || (degD === 0 && D[0].isZero())) {
        throw new Error("Polynomial division by zero polynomial.");
    }
    if (degN < degD) {
        return {
            quotient: [new Decimal(0)],
            remainder: N
        };
    }
    const qDeg = degN - degD;
    const quotient = new Array(qDeg + 1).fill(null).map(() => new Decimal(0));
    const leadD = D[degD];
    let rem = N.map(c => new Decimal(c));
    for (let i = qDeg; i >= 0; i--) {
        const curDeg = i + degD;
        if (rem.length - 1 >= curDeg) {
            const coeff = rem[curDeg].div(leadD);
            quotient[i] = coeff;
            for (let j = 0; j <= degD; j++) {
                rem[i + j] = rem[i + j].minus(coeff.times(D[j]));
            }
        }
    }
    return {
        quotient: trimLeadingZeros(quotient),
        remainder: trimLeadingZeros(rem)
    };
}

function normalizePoly(poly) {
    const trimmed = trimLeadingZeros(poly);
    const deg = trimmed.length - 1;
    if (deg < 0 || trimmed[deg].isZero()) return trimmed;
    const lead = trimmed[deg];
    return trimmed.map(c => c.div(lead));
}

function polyGCD(poly1, poly2, eps = new Decimal('1e-50')) {
    let a = normalizePoly(poly1);
    let b = normalizePoly(poly2);
    while (getPolyDegree(b) >= 0 && !b.every(c => c.abs().lt(eps))) {
        const divResult = polyDivide(a, b);
        a = b;
        b = normalizePoly(divResult.remainder);
    }
    return normalizePoly(a);
}

function divideExactly(num, den) {
    const { quotient, remainder } = polyDivide(num, den);
    const remDeg = getPolyDegree(remainder);
    if (remDeg >= 0 && !remainder.every(c => c.abs().lt(new Decimal('1e-40')))) {
        throw new Error("Inexact polynomial division detected.");
    }
    return quotient;
}

// ============================================================================
// 3. EXACT BIGINT POLYNOMIAL ALGEBRA (SQUARE-FREE PREPROCESSING)
// ============================================================================

function bigintAbs(x) {
    return x < 0n ? -x : x;
}

function bigintGcd(a, b) {
    let x = bigintAbs(a);
    let y = bigintAbs(b);
    while (y !== 0n) {
        const t = x % y;
        x = y;
        y = t;
    }
    return x;
}

function trimLeadingZerosBigInt(poly) {
    let end = poly.length;
    while (end > 1 && poly[end - 1] === 0n) {
        end--;
    }
    return poly.slice(0, end);
}

function getDegreeBigInt(poly) {
    return trimLeadingZerosBigInt(poly).length - 1;
}

function primitiveIntegerPoly(poly) {
    const trimmed = trimLeadingZerosBigInt(poly);
    if (trimmed.length === 1 && trimmed[0] === 0n) {
        return { poly: [0n], content: 0n };
    }
    let content = 0n;
    for (const coeff of trimmed) {
        if (coeff !== 0n) {
            content = content === 0n ? bigintAbs(coeff) : bigintGcd(content, coeff);
        }
    }
    if (content === 0n) content = 1n;
    const lead = trimmed[trimmed.length - 1];
    const sign = lead < 0n ? -1n : 1n;
    const divisor = content * sign;
    return {
        poly: trimmed.map(c => c / divisor),
        content: content * sign
    };
}

function integerPolyGcd(poly1, poly2) {
    let A = trimLeadingZerosBigInt(poly1);
    let B = trimLeadingZerosBigInt(poly2);
    if (getDegreeBigInt(A) < getDegreeBigInt(B)) {
        const temp = A;
        A = B;
        B = temp;
    }
    if (getDegreeBigInt(B) < 0 || (B.length === 1 && B[0] === 0n)) {
        return primitiveIntegerPoly(A).poly;
    }
    let { poly: F } = primitiveIntegerPoly(A);
    let { poly: G } = primitiveIntegerPoly(B);
    let g = 1n;
    let h = 1n;
    while (getDegreeBigInt(G) > 0) {
        const df = getDegreeBigInt(F);
        const dg = getDegreeBigInt(G);
        const delta = BigInt(df - dg);
        const leadG = G[G.length - 1];
        let rem = F.map(c => c * (leadG ** (delta + 1n)));
        const qDeg = df - dg;
        for (let i = qDeg; i >= 0; i--) {
            const curDeg = i + dg;
            if (rem.length - 1 >= curDeg && rem[curDeg] !== 0n) {
                const factor = rem[curDeg] / leadG;
                for (let j = 0; j <= dg; j++) {
                    rem[i + j] -= factor * G[j];
                }
            }
        }
        const R = trimLeadingZerosBigInt(rem);
        if (R.length === 1 && R[0] === 0n) {
            return primitiveIntegerPoly(G).poly;
        }
        F = G;
        const denominator = g * (h ** delta);
        G = R.map(c => c / denominator);
        g = leadG;
        h = delta > 1n ? (leadG ** delta) / (h ** (delta - 1n)) : leadG;
    }
    if (G.length === 1 && G[0] !== 0n) {
        return [1n];
    }
    return primitiveIntegerPoly(G).poly;
}

function divideIntegerPolyExactly(num, den) {
    const N = trimLeadingZerosBigInt(num);
    const D = trimLeadingZerosBigInt(den);
    const degN = getDegreeBigInt(N);
    const degD = getDegreeBigInt(D);
    if (degD < 0 || (D.length === 1 && D[0] === 0n)) {
        throw new Error("Exact integer division by zero polynomial.");
    }
    if (degN < degD) {
        throw new Error("Inexact integer polynomial division: degree(num) < degree(den).");
    }
    const qDeg = degN - degD;
    const quotient = new Array(qDeg + 1).fill(0n);
    const rem = [...N];
    const leadD = D[degD];
    for (let i = qDeg; i >= 0; i--) {
        const curDeg = i + degD;
        if (rem[curDeg] !== 0n) {
            if (rem[curDeg] % leadD !== 0n) {
                throw new Error("Inexact integer polynomial division quotient coefficient.");
            }
            const coeff = rem[curDeg] / leadD;
            quotient[i] = coeff;
            for (let j = 0; j <= degD; j++) {
                rem[i + j] -= coeff * D[j];
            }
        }
    }
    for (const c of rem) {
        if (c !== 0n) throw new Error("Inexact integer polynomial division remainder non-zero.");
    }
    return trimLeadingZerosBigInt(quotient);
}

function derivativeIntegerPoly(poly) {
    const trimmed = trimLeadingZerosBigInt(poly);
    const deg = getDegreeBigInt(trimmed);
    if (deg <= 0) return [0n];
    const deriv = [];
    for (let i = 1; i <= deg; i++) {
        deriv.push(trimmed[i] * BigInt(i));
    }
    return trimLeadingZerosBigInt(deriv);
}

function yunDecomposeIntegerPoly(poly) {
    const prim = primitiveIntegerPoly(poly).poly;
    const A = prim;
    const Aprime = derivativeIntegerPoly(A);
    const C = integerPolyGcd(A, Aprime);
    let W = divideIntegerPolyExactly(A, C);
    let Y = divideIntegerPolyExactly(Aprime, C);
    let i = 1;
    const factors = [];
    while (getDegreeBigInt(W) > 0) {
        const YminusWprime = trimLeadingZerosBigInt(
            Y.map((c, idx) => c - (derivativeIntegerPoly(W)[idx] || 0n))
        );
        const G = integerPolyGcd(W, YminusWprime);

        if (getDegreeBigInt(G) > 0) {
            factors.push({ factor: G, multiplicity: i });
        }
        W = divideIntegerPolyExactly(W, G);
        Y = divideIntegerPolyExactly(YminusWprime, G);
        i++;
    }
    return factors;
}

// ============================================================================
// 4. PREPROCESSING & FALLBACK PIPELINE
// ============================================================================

function preprocessSquareFree(inputCoeffs, debug = false) {
    Decimal.set({ precision: PRECISION_PREPROCESS });
    let allInteger = true;
    const bigintCoeffs = [];
    for (const c of inputCoeffs) {
        try {
            const d = new Decimal(c);
            if (d.isInteger()) {
                bigintCoeffs.push(BigInt(d.toFixed(0)));
            } else {
                allInteger = false;
                break;
            }
        } catch (e) {
            allInteger = false;
            break;
        }
    }
    if (allInteger && bigintCoeffs.length > 1) {
        try {
            if (debug) console.log("Attempting Exact BigInt Yun Square-Free Decomposition...");
            const factorsBigInt = yunDecomposeIntegerPoly(bigintCoeffs);

            if (factorsBigInt.length > 0) {
                if (debug) console.log(`Yun's exact integer decomposition successful. Found ${factorsBigInt.length} factor group(s).`);
                return factorsBigInt.map(f => ({
                    factor: f.factor.map(c => new Decimal(c.toString())),
                    multiplicity: f.multiplicity
                }));
            }
        } catch (err) {
            if (debug) console.log("Yun exact decomposition fallback triggered:", err.message);
        }
    }

    // Fallback: Use full normalized polynomial safely without dropping roots
    if (debug) console.log("Using full normalized polynomial as safe fallback.");
    const decCoeffs = inputCoeffs.map(c => new Decimal(c));
    const normalized = normalizePoly(decCoeffs);
    return [{ factor: normalized, multiplicity: 1 }];
}

// ============================================================================
// 5. SOLVER KERNEL (DKA METHOD WITH ABERTH CORRECTION)
// ============================================================================

function solveFactor(factorCoeffs, maxIter = 500, debug = false) {
    Decimal.set({ precision: PRECISION_DKA });
    const coeffs = factorCoeffs.map(c => new Decimal(c));
    const deg = coeffs.length - 1;
    if (deg <= 0) return [];
    if (deg === 1) {
        const root = coeffs[0].negated().div(coeffs[1]);
        return [{ re: root, im: new Decimal(0), iterations: 0, finalCorrection: new Decimal(0) }];
    }
    const R = estimateRootScale(coeffs);
    let roots = [];
    for (let i = 0; i < deg; i++) {
        const angle = new Decimal(2).times(Decimal.acos(-1)).times(i).plus(new Decimal('0.5')).div(deg);
        const radius = R.times(new Decimal('0.9').plus(new Decimal('0.2').times(i).div(deg)));
        roots.push({
            re: radius.times(Decimal.cos(angle)),
            im: radius.times(Decimal.sin(angle))
        });
    }
    let converged = new Array(deg).fill(false);
    let iterations = 0;
    let finalCorrections = new Array(deg).fill(new Decimal(0));
    const leadCoeff = coeffs[deg];
    while (iterations < maxIter && converged.includes(false)) {
        iterations++;
        let maxChange = new Decimal(0);
        for (let i = 0; i < deg; i++) {
            if (converged[i]) continue;
            const z = roots[i];
            const pVal = evalPoly(coeffs, z);

            let pDeriv = { re: new Decimal(0), im: new Decimal(0) };
            for (let k = coeffs.length - 1; k >= 1; k--) {
                const coeffComp = { re: coeffs[k].times(k), im: new Decimal(0) };
                pDeriv = c_add(c_mul(pDeriv, z), coeffComp);
            }
            if (c_abs(pVal).isZero()) {
                converged[i] = true;
                continue;
            }
            const NewtonRatio = c_div(pVal, pDeriv);
            let sumAberth = { re: new Decimal(0), im: new Decimal(0) };
            for (let j = 0; j < deg; j++) {
                if (i !== j) {
                    const diff = c_sub(z, roots[j]);
                    const invDiff = c_div({ re: new Decimal(1), im: new Decimal(0) }, diff);
                    sumAberth = c_add(sumAberth, invDiff);
                }
            }
            const denomAberth = c_sub({ re: new Decimal(1), im: new Decimal(0) }, c_mul(NewtonRatio, sumAberth));
            const correction = c_div(NewtonRatio, denomAberth);
            roots[i] = c_sub(z, correction);
            const corrAbs = c_abs(correction);
            finalCorrections[i] = corrAbs;
            if (corrAbs.gt(maxChange)) {
                maxChange = corrAbs;
            }
            if (corrAbs.lt(new Decimal('1e-250'))) {
                converged[i] = true;
            }
        }
        if (maxChange.lt(new Decimal('1e-250'))) {
            break;
        }
    }
    return roots.map((r, idx) => ({
        re: r.re,
        im: r.im,
        iterations: iterations,
        finalCorrection: finalCorrections[idx]
    }));
}

// ============================================================================
// 6. MAIN EXPORTED FUNCTION
// ============================================================================

/**
 * Solves polynomial root finding using Durand-Kerner-Aberth method with square-free preprocessing.
 * 
 * @param {Array<string|number|Decimal>} coeffs - Coefficients ordered from highest to lowest degree [a_n, a_{n-1}, ..., a_0] (Descending)
 * @param {number} maxIter - Maximum iterations allowed for the DKA loop per factor
 * @param {boolean} debug - Toggles debug output logs
 * @returns {Object} Solution payload containing extracted roots, residuals, and stats.
 */
function solvePolynomial(coeffs, maxIter = 500, debug = false) {
    if (!Array.isArray(coeffs) || coeffs.length < 2) {
        throw new Error("Input coefficients array must contain at least two entries.");
    }
    // Convert external descending order to internal ascending order [a_0, a_1, ..., a_n]
    const origCoeffs = coeffs.map(c => new Decimal(c)).reverse();

    // Step 1: Preprocess to extract square-free factors
    const factors = preprocessSquareFree(origCoeffs, debug);

    // Step 2: Solve each factor
    const solvedRoots = [];
    for (const f of factors) {
        const rawRoots = solveFactor(f.factor, maxIter, debug);
        for (const r of rawRoots) {
            const residual = normalizedResidual(origCoeffs, r);
            solvedRoots.push({
                re: r.re,
                im: r.im,
                residual: residual,
                iterations: r.iterations,
                multiplicity: f.multiplicity,
                finalCorrection: r.finalCorrection
            });
        }
    }

    // Step 3: Global Residual Evaluation
    const maxCoeffResidual = coefficientReconstructionResidual(solvedRoots, origCoeffs);

    return {
        roots: solvedRoots,
        maxCoefficientResidual: maxCoeffResidual,
        factorCount: factors.length
    };
}

module.exports = {
    solvePolynomial,
    c_add,
    c_sub,
    c_mul,
    c_div,
    c_abs,
    evalPoly,
    normalizedResidual,
    coefficientReconstructionResidual,
    estimateRootScale,
    trimLeadingZeros,
    getPolyDegree,
    derivativePoly,
    polyDivide,
    normalizePoly,
    polyGCD,
    divideExactly,
    primitiveIntegerPoly,
    integerPolyGcd,
    divideIntegerPolyExactly,
    derivativeIntegerPoly,
    yunDecomposeIntegerPoly,
    preprocessSquareFree,
    solveFactor
};

// ============================================================================
// 7. CLI EXECUTION & FILE I/O (v3.4.0 EXTENSION)
// ============================================================================

/**
 * Fallback: Generate Wilkinson polynomial of degree 20 (W_20)
 * Evaluates: (x-1)(x-2)...(x-20) and returns coefficients from highest to lowest.
 */
function buildWilkinson20Fallback() {
    let poly = [1n];
    for (let i = 1n; i <= 20n; i++) {
        let nextPoly = new Array(poly.length + 1).fill(0n);
        for (let j = 0; j < poly.length; j++) {
            nextPoly[j + 1] += poly[j];        // * x
            nextPoly[j] -= poly[j] * i;        // * (-i)
        }
        poly = nextPoly;
    }
    // Convert to strings and reverse to get Highest -> Lowest degree
    return poly.reverse().map(c => c.toString());
}

async function main() {
    const inputPath = path.join(__dirname, 'input.txt');
    const outputPath = path.join(__dirname, 'output.txt');
    let coefficients = [];

    // 1. Check for input.txt (Highest priority)
    if (fs.existsSync(inputPath)) {
        const content = fs.readFileSync(inputPath, 'utf-8');
        // Parse by whitespace or newlines, ignoring empty lines
        coefficients = content.split(/\s+/).map(line => line.trim()).filter(line => line.length > 0);
        if (coefficients.length < 2) {
            console.error("[Error] input.txt には最低2つの係数が必要です。");
            process.exit(1);
        }
    } else {
        // 2. Fallback to built-in Wilkinson polynomial generator
        console.log("[Info] input.txt が見つかりません。デフォルトのウィルキンソン多項式(20次)テストを実行します...");
        coefficients = buildWilkinson20Fallback();
    }

    try {
        const startTime = process.hrtime();

        // Execute DKA Solver
        const results = solvePolynomial(coefficients, 5000, false);

        const diff = process.hrtime(startTime);
        const execTimeMs = (diff[0] * 1e3 + diff[1] * 1e-6).toFixed(3);

        // Format Output (Aligned with standard testing format)
        let outputText = `Execution Time: ${execTimeMs} ms\n`;
        outputText += `Max Coefficient Residual: ${results.maxCoefficientResidual.toExponential(10)}\n\n`;

        results.roots.forEach((root, index) => {
            const reStr = root.re.toFixed(20);
            const imStr = root.im.isNegative() ? root.im.toFixed(20) : `+${root.im.toFixed(20)}`;
            const resStr = root.residual.toExponential(6);
            outputText += `Root #${index + 1}: Re=${reStr}, Im=${root.im.toFixed(20)}, Res=${resStr}\n`;
        });

        // Write to output.txt and print to Console
        fs.writeFileSync(outputPath, outputText, 'utf-8');
        console.log("--------------------------------------------------");
        console.log(outputText.trim());
        console.log("--------------------------------------------------");

    } catch (err) {
        console.error(`[Error] 実行時エラー: ${err.message}`);
    }
}

// Execute main() if script is run directly from the command line
if (require.main === module) {
    main();
}