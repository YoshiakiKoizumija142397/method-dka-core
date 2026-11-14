/**
 * Core Mathematical Engine for Durand-Kerner-Aberth Algorithm (v3.7.7)
 * 
 * High-Precision General Polynomial Root Finder using the DKA Method
 * with Yun's Square-Free Decomposition preprocessing, 400-digit internal precision,
 * 180-digit convergence criteria, 180-digit output formatting, and Individual Root Iteration Tracking.
 * 
 * Requires: decimal.js (npm install decimal.js)
 */

const fs = require('fs');
const path = require('path');
const Decimal = require('decimal.js');

// Configurable Precision Hierarchy Standards (v3.7.7)
const PRECISION_PREPROCESS = 420;
const PRECISION_DKA = 400;
const MAX_SAFETY_INTERVENTIONS = 500; // セイフティロック許容上限回数

/**
 * Configure global Decimal settings to 400-digit precision and expanded exponent bounds
 */
Decimal.set({
    precision: PRECISION_DKA,
    maxE: 1e9,
    minE: -1e9,
    rounding: Decimal.ROUND_HALF_UP
});

// ============================================================================
// 1. COMPLEX ARITHMETIC HELPER FUNCTIONS (Type-Safe Decimal Wrapping)
// ============================================================================

function toDec(val) {
    return val instanceof Decimal ? val : new Decimal(val);
}

function c_wrap(re, im) {
    return {
        re: toDec(re),
        im: toDec(im)
    };
}

function c_add(a, b) {
    return {
        re: toDec(a.re).plus(toDec(b.re)),
        im: toDec(a.im).plus(toDec(b.im))
    };
}

function c_sub(a, b) {
    return {
        re: toDec(a.re).minus(toDec(b.re)),
        im: toDec(a.im).minus(toDec(b.im))
    };
}

function c_mul(a, b) {
    const ar = toDec(a.re), ai = toDec(a.im);
    const br = toDec(b.re), bi = toDec(b.im);
    return {
        re: ar.times(br).minus(ai.times(bi)),
        im: ar.times(bi).plus(ai.times(br))
    };
}

function c_div(a, b) {
    const ar = toDec(a.re), ai = toDec(a.im);
    const br = toDec(b.re), bi = toDec(b.im);
    const denom = br.times(br).plus(bi.times(bi));
    if (denom.isZero()) {
        throw new Error("Complex division by zero.");
    }
    return {
        re: ar.times(br).plus(ai.times(bi)).div(denom),
        im: ai.times(br).minus(ar.times(bi)).div(denom)
    };
}

function c_abs(a) {
    const ar = toDec(a.re), ai = toDec(a.im);
    return ar.times(ar).plus(ai.times(ai)).sqrt();
}

/**
 * Evaluates polynomial P(z) = sum_{k=0}^n a_k * z^k using Horner's Method
 */
function evalPoly(coeffs, z) {
    let result = c_wrap(0, 0);
    const zComp = c_wrap(z.re, z.im);
    for (let i = coeffs.length - 1; i >= 0; i--) {
        const coeffComp = c_wrap(coeffs[i], 0);
        result = c_add(c_mul(result, zComp), coeffComp);
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
        const negR = c_wrap(toDec(r.re).negated(), toDec(r.im).negated());
        for (let m = 0; m < r.multiplicity; m++) {
            const nextPoly = new Array(poly.length + 1).fill(null).map(() => c_wrap(0, 0));
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
        const scaledRe = toDec(poly[i].re).times(leadOrig);
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
    while (end > 1 && toDec(poly[end - 1]).isZero()) {
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
        deriv.push(toDec(trimmed[i]).times(i));
    }
    return trimLeadingZeros(deriv);
}

function polyDivide(num, den) {
    let N = trimLeadingZeros(num);
    const D = trimLeadingZeros(den);
    const degN = N.length - 1;
    const degD = D.length - 1;
    if (degD < 0 || (degD === 0 && toDec(D[0]).isZero())) {
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
    const leadD = toDec(D[degD]);
    let rem = N.map(c => new Decimal(c));
    for (let i = qDeg; i >= 0; i--) {
        const curDeg = i + degD;
        if (rem.length - 1 >= curDeg) {
            const coeff = toDec(rem[curDeg]).div(leadD);
            quotient[i] = coeff;
            for (let j = 0; j <= degD; j++) {
                rem[i + j] = toDec(rem[i + j]).minus(coeff.times(toDec(D[j])));
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
    if (deg < 0 || toDec(trimmed[deg]).isZero()) return trimmed;
    const lead = toDec(trimmed[deg]);
    return trimmed.map(c => toDec(c).div(lead));
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
// 4. PREPROCESSING & FALLBACK PIPELINE (v3.7.7)
// ============================================================================

function preprocessSquareFree(inputCoeffs, debug = false) {
    Decimal.set({ precision: PRECISION_PREPROCESS, maxE: 1e9, minE: -1e9 });
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

    if (debug) console.log("Using full normalized polynomial as safe fallback.");
    const decCoeffs = inputCoeffs.map(c => new Decimal(c));
    const normalized = normalizePoly(decCoeffs);
    return [{ factor: normalized, multiplicity: 1 }];
}

// ============================================================================
// 5. SOLVER KERNEL (DKA METHOD WITH INDIVIDUAL ROOT ITERATION TRACKING)
// ============================================================================

function solveFactor(factorCoeffs, maxIter = 5000, debug = false) {
    Decimal.set({ precision: PRECISION_DKA, maxE: 1e9, minE: -1e9 });
    const coeffs = factorCoeffs.map(c => new Decimal(c));
    const deg = coeffs.length - 1;
    if (deg <= 0) return { rawRoots: [], safetyTriggers: 0 };
    if (deg === 1) {
        const root = coeffs[0].negated().div(coeffs[1]);
        return {
            rawRoots: [{ re: root, im: new Decimal(0), iterations: 1, finalCorrection: new Decimal(0) }],
            safetyTriggers: 0
        };
    }
    const R = estimateRootScale(coeffs);
    let roots = [];
    for (let i = 0; i < deg; i++) {
        const angle = new Decimal(2).times(Decimal.acos(-1)).times(i).plus(new Decimal('0.5')).div(deg);
        const radius = R.times(new Decimal('0.9').plus(new Decimal('0.2').times(i).div(deg)));
        roots.push(c_wrap(
            radius.times(Decimal.cos(angle)),
            radius.times(Decimal.sin(angle))
        ));
    }
    let converged = new Array(deg).fill(false);
    let rootConvergedIters = new Array(deg).fill(0); // 各根ごとの個別収束回数
    let iterations = 0;
    let finalCorrections = new Array(deg).fill(new Decimal(0));
    let safetyLockTriggerCount = 0;
    
    const factorLabel = `[Degree ${deg}]`;
    const tolerance = new Decimal('1e-180');

    while (iterations < maxIter && converged.includes(false)) {
        iterations++;
        let maxChange = new Decimal(0);
        let newRoots = roots.map(r => c_wrap(r.re, r.im));

        for (let i = 0; i < deg; i++) {
            if (converged[i]) continue;
            let z_re = toDec(roots[i].re);
            let z_im = toDec(roots[i].im);

            if (!z_re.isFinite() || !z_im.isFinite()) {
                safetyLockTriggerCount++;
                if (safetyLockTriggerCount > MAX_SAFETY_INTERVENTIONS) {
                    throw new Error(`[Fatal Error] Safety lock limit exceeded (${MAX_SAFETY_INTERVENTIONS}). Calculation aborted.`);
                }
                z_re = toDec(R);
                z_im = new Decimal(0);
            }

            const z = c_wrap(z_re, z_im);
            const pVal = evalPoly(coeffs, z);

            let pDeriv = c_wrap(0, 0);
            for (let k = coeffs.length - 1; k >= 1; k--) {
                const coeffComp = c_wrap(coeffs[k].times(k), 0);
                pDeriv = c_add(c_mul(pDeriv, z), coeffComp);
            }
            if (c_abs(pVal).isZero()) {
                converged[i] = true;
                rootConvergedIters[i] = iterations;
                continue;
            }
            const NewtonRatio = c_div(pVal, pDeriv);
            let sumAberth = c_wrap(0, 0);
            for (let j = 0; j < deg; j++) {
                if (i !== j) {
                    const diff = c_sub(z, roots[j]);
                    const invDiff = c_div(c_wrap(1, 0), diff);
                    sumAberth = c_add(sumAberth, invDiff);
                }
            }
            const denomAberth = c_sub(c_wrap(1, 0), c_mul(NewtonRatio, sumAberth));
            let correction = c_div(NewtonRatio, denomAberth);

            if (!toDec(correction.re).isFinite() || !toDec(correction.im).isFinite()) {
                safetyLockTriggerCount++;
                if (safetyLockTriggerCount > MAX_SAFETY_INTERVENTIONS) {
                    throw new Error(`[Fatal Error] Safety lock limit exceeded (${MAX_SAFETY_INTERVENTIONS}). Calculation aborted.`);
                }
                correction = c_wrap(0, 0);
            }

            newRoots[i] = c_sub(z, correction);
            const corrAbs = c_abs(correction);
            finalCorrections[i] = corrAbs;
            if (corrAbs.gt(maxChange)) {
                maxChange = corrAbs;
            }
            if (corrAbs.lt(tolerance)) {
                converged[i] = true;
                rootConvergedIters[i] = iterations;
            }
        }

        // まだ収束していない根についても、ループが終了した時点で現在のイテレーション数を仮設定
        for (let i = 0; i < deg; i++) {
            if (!converged[i]) {
                rootConvergedIters[i] = iterations;
            }
        }

        roots = newRoots;

        if (!debug && (iterations % 10 === 0 || iterations === maxIter || converged.every(Boolean))) {
            const percent = ((iterations / maxIter) * 100).toFixed(1).padStart(5, ' ');
            const progressMsg = `  -> ${factorLabel} Progress: ${percent}% (Iter: ${iterations}/${maxIter}, MaxChange: ${maxChange.toExponential(2)})`;
            process.stdout.write(progressMsg + '\r');
        }

        if (maxChange.lt(tolerance)) {
            break;
        }
    }

    if (!debug) {
        process.stdout.write(`  -> ${factorLabel} Completed in ${iterations} iterations. (Safety Locks Triggered: ${safetyLockTriggerCount})\n`);
    }

    return {
        rawRoots: roots.map((r, idx) => ({
            re: toDec(r.re),
            im: toDec(r.im),
            iterations: rootConvergedIters[idx], // 各根が収束した個別のイテレーション数
            finalCorrection: finalCorrections[idx]
        })),
        safetyTriggers: safetyLockTriggerCount
    };
}

// ============================================================================
// 6. MAIN EXPORTED FUNCTION
// ============================================================================

function solvePolynomial(coeffs, maxIter = 5000, debug = false) {
    if (!Array.isArray(coeffs) || coeffs.length < 2) {
        throw new Error("Input coefficients array must contain at least two entries.");
    }
    const origCoeffs = coeffs.map(c => new Decimal(c)).reverse();

    const factors = preprocessSquareFree(origCoeffs, debug);

    const solvedRoots = [];
    let totalSafetyTriggers = 0;

    for (const f of factors) {
        const factorResult = solveFactor(f.factor, maxIter, debug);
        totalSafetyTriggers += factorResult.safetyTriggers;
        for (const r of factorResult.rawRoots) {
            const residual = normalizedResidual(origCoeffs, r);
            solvedRoots.push({
                re: toDec(r.re),
                im: toDec(r.im),
                residual: residual,
                iterations: r.iterations,
                multiplicity: f.multiplicity,
                finalCorrection: r.finalCorrection
            });
        }
    }

    const maxCoeffResidual = coefficientReconstructionResidual(solvedRoots, origCoeffs);

    return {
        roots: solvedRoots,
        maxCoefficientResidual: maxCoeffResidual,
        factorCount: factors.length,
        safetyTriggers: totalSafetyTriggers
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
    preprocessSquareFree,
    solveFactor
};

// ============================================================================
// 7. CLI EXECUTION & FILE I/O (v3.7.7: Individual Iterations in Report)
// ============================================================================

async function main() {
    const inputPath = path.join(__dirname, 'input.txt');
    const outputPath = path.join(__dirname, 'output.txt');
    let coefficients = [];

    if (fs.existsSync(inputPath)) {
        const content = fs.readFileSync(inputPath, 'utf-8');
        coefficients = content.split(/[\s,]+/).map(line => line.trim()).filter(line => line.length > 0);
        if (coefficients.length < 2) {
            console.error("[Error] input.txt には最低2つの係数が必要です。");
            process.exit(1);
        }
    } else {
        console.error("[Error] input.txt が見つかりません。入力係数ファイルを配置してください。");
        process.exit(1);
    }

    try {
        const startTime = process.hrtime();
        const timestamp = new Date().toISOString();

        const results = solvePolynomial(coefficients, 5000, false);

        const diff = process.hrtime(startTime);
        const execTimeMs = (diff[0] * 1e3 + diff[1] * 1e-6).toFixed(3);

        let outputText = "==================================================\n";
        outputText += " METHOD-DKA-CORE EXECUTION AUDIT REPORT (v3.7.7)\n";
        outputText += "==================================================\n";
        outputText += `Timestamp: ${timestamp}\n`;
        outputText += `Internal Precision: ${PRECISION_DKA} digits (Preprocess: ${PRECISION_PREPROCESS} digits)\n`;
        outputText += `Execution Time: ${execTimeMs} ms\n`;
        outputText += `Safety Lock Interventions Triggered: ${results.safetyTriggers} times (Limit: ${MAX_SAFETY_INTERVENTIONS})\n`;
        outputText += `Max Coefficient Residual: ${results.maxCoefficientResidual.toExponential(10)}\n`;
        outputText += "--------------------------------------------------\n";
        outputText += "Input Coefficients:\n";
        outputText += coefficients.join(", ") + "\n";
        outputText += "--------------------------------------------------\n\n";

        results.roots.forEach((root, index) => {
            const reStr = toDec(root.re).toFixed(180);
            const imVal = toDec(root.im);
            const imSign = imVal.isNegative() ? "" : "+";
            const imStr = imSign + imVal.toFixed(180);
            const resStr = root.residual.toExponential(6);
            const iterVal = root.iterations;
            outputText += `Root #${index + 1}: Re=${reStr}, Im=${imStr}, Res=${resStr}, Iter=${iterVal}\n`;
        });

        outputText += "\n==================================================\n";

        fs.writeFileSync(outputPath, outputText, 'utf-8');
        console.log("--------------------------------------------------");
        console.log(outputText.trim());
        console.log("--------------------------------------------------");

    } catch (err) {
        console.error(`[Error] 実行時エラー: ${err.message}`);
    }
}

if (require.main === module) {
    main();
}