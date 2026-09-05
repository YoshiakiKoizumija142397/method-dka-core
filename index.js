const Decimal = require('decimal.js');

class MethodDka {
    constructor(options = {}) {
        this.precision = options.precision || 230;
        this.eps = new Decimal(options.eps || "1e-170");
        this.maxIter = options.maxIter || 1000;
    }

    c_add(a, b) { return { re: a.re.add(b.re), im: a.im.add(b.im) }; }
    c_sub(a, b) { return { re: a.re.sub(b.re), im: a.im.sub(b.im) }; }
    c_mul(a, b) {
        return {
            re: a.re.mul(b.re).sub(a.im.mul(b.im)),
            im: a.re.mul(b.im).add(a.im.mul(b.re))
        };
    }
    c_div(a, b) {
        const denom = b.re.mul(b.re).add(b.im.mul(b.im));
        return {
            re: a.re.mul(b.re).add(a.im.mul(b.im)).div(denom),
            im: a.im.mul(b.re).sub(a.re.mul(b.im)).div(denom)
        };
    }
    c_abs(a) { return Decimal.sqrt(a.re.mul(a.re).add(a.im.mul(a.im))); }

    evalPoly(A, z) {
        let val = { re: new Decimal(A[0].re), im: new Decimal(A[0].im) };
        for (let i = 1; i < A.length; i++) {
            val = this.c_add(this.c_mul(val, z), A[i]);
        }
        return val;
    }

    /**
     * 多項式の係数配列を受け取り、DKA法で根を計算する
     * @param {Array<string|number>} rawCoeffs 高次から定数項までの係数 [a_n, a_{n-1}, ..., a_0]
     * @returns {Array} 計算された根のリスト
     */
    solve(rawCoeffs) {
        Decimal.set({ precision: this.precision });
        const n = rawCoeffs.length - 1;

        let A = rawCoeffs.map(c => ({ re: new Decimal(c || 0), im: new Decimal(0) }));
        const lead = A[0];
        for (let i = 0; i <= n; i++) {
            A[i] = this.c_div(A[i], lead);
        }

        let z = [];
        let R = new Decimal("150");
        const pi = Decimal.acos(-1);

        for (let k = 0; k < n; k++) {
            let angle = new Decimal(2).mul(pi).mul(k).div(n).add(new Decimal(3).div(new Decimal(4).mul(n)));
            z.push({
                re: R.mul(Decimal.cos(angle)),
                im: R.mul(Decimal.sin(angle)),
                iter: 0,
                converged: false,
                errRadius: new Decimal(Infinity)
            });
        }

        let iter = 1;
        let allConverged = false;

        while (iter <= this.maxIter && !allConverged) {
            allConverged = true;
            for (let i = 0; i < n; i++) {
                if (z[i].converged) continue;

                let f_val = this.evalPoly(A, z[i]);
                let denom = { re: new Decimal(1), im: new Decimal(0) };

                for (let j = 0; j < n; j++) {
                    if (i === j) continue;
                    let diff = this.c_sub(z[i], z[j]);
                    denom = this.c_mul(denom, diff);
                }

                let delta = this.c_div(f_val, denom);
                let errRad = this.c_abs(delta);

                z[i].errRadius = errRad;
                z[i].iter = iter;

                if (errRad.lte(this.eps)) {
                    z[i].converged = true;
                } else {
                    allConverged = false;
                    z[i] = this.c_sub(z[i], delta);
                    z[i].iter = iter;
                    z[i].errRadius = errRad;
                }
            }
            iter++;
        }

        z.sort((a, b) => (a.re.sub(b.re).isNegative() ? -1 : 1));

        return z.map((ans, idx) => ({
            id: idx + 1,
            re: ans.re,
            im: ans.im,
            reStr: ans.re.toFixed(100),
            imStr: ans.im.toFixed(100),
            errRadiusStr: ans.errRadius.toExponential(4),
            iterations: ans.iter,
            converged: ans.converged
        }));
    }
}

module.exports = MethodDka;
