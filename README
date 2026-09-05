# MethodDka

High-precision polynomial root solver using the Durand-Kerner-Aberth (DKA) method and Decimal.js.

## Installation

```bash
npm install method-dka
Usage
JavaScript
const MethodDka = require('method-dka');

// インスタンス化 (必要に応じてオプション指定)
const solver = new MethodDka({
    precision: 230,
    maxIter: 1000
});

// 例: 2次方程式 x^2 - 3x + 2 = 0 の係数 [1, -3, 2]
const coefficients = [1, -3, 2];

const roots = solver.solve(coefficients);

roots.forEach(root => {
    console.log(`Root ${root.id}: ${root.reStr} + ${root.imStr}i (Iter: ${root.iterations})`);
});
