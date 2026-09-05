# MethodDka

High-precision polynomial root-finding library using the Durand-Kerner (DKA) method and Decimal.js.

## 🇯🇵 日本語版 (Japanese)

### 概要

`method-dka` は、デュランド・カーネル法（Durand-Kerner Method / DKA法）と `decimal.js` を用いて、任意の多項式の根を**超高精度**で並行計算するJavaScriptライブラリです。標準の浮動小数点演算では精度の足りない科学技術計算や、極めて高い桁数が求められる数学的検証に最適です。

### 特徴

* **任意の高精度計算**: `decimal.js` をベースにしているため、小数第百位を超えるような高精度な計算が可能です。
* **頑健な反復処理**: 複素数平面上で複数の根を同時に収束させます。
* **エラー追跡**: 各根の反復回数や誤差半径（errRadius）をトラッキングできます。
* **MITライセンス**: 個人開発から商用利用まで自由にご利用いただけます。

### インストール

```bash
npm install method-dka

```

### 使い方

```javascript
const MethodDka = require('method-dka');

// インスタンス化 (精度や最大反復回数はオプションで指定可能)
const solver = new MethodDka({
    precision: 230,
    maxIter: 1000
});

// 多項式の係数を高次から順に指定 [a_n, ..., a_0]
// 例: x^2 - 3x + 2 = 0 の場合
const coefficients = [1, -3, 2];

const roots = solver.solve(coefficients);

roots.forEach(root => {
    const sign = root.im.isNegative() ? '' : '+';
    console.log(`[Root ${root.id}] ${root.reStr} ${sign} ${root.imStr}i (Iter: ${root.iterations}, Error Radius: ${root.errRadiusStr})`);
});

```

---

## 🇬🇧 英語翻訳版 (English Version)

### Overview

`method-dka` is a high-precision JavaScript library designed to find all roots of arbitrary polynomials concurrently using the Durand-Kerner (DKA) method backed by `decimal.js`. It is ideal for scientific computations and mathematical explorations that require arbitrary precision beyond standard floating-point limitations.

### Features

* **Arbitrary High Precision**: Built on `decimal.js`, supporting deep-precision calculations well beyond standard limits (e.g., 230+ digits).
* **Robust Iterative Solver**: Simultaneously converges all roots on the complex plane.
* **Error Tracking**: Tracks iteration counts and error radius (`errRadius`) for each calculated root.
* **MIT License**: Free for open-source and commercial use.

### Installation

```bash
npm install method-dka

```

### Usage

```javascript
const MethodDka = require('method-dka');

// Instantiate with optional configuration
const solver = new MethodDka({
    precision: 230,
    maxIter: 1000
});

// Specify polynomial coefficients from highest degree to constant [a_n, ..., a_0]
// Example: x^2 - 3x + 2 = 0
const coefficients = [1, -3, 2];

const roots = solver.solve(coefficients);

roots.forEach(root => {
    const sign = root.im.isNegative() ? '' : '+';
    console.log(`[Root ${root.id}] ${root.reStr} ${sign} ${root.imStr}i (Iter: ${root.iterations}, Error Radius: ${root.errRadiusStr})`);
});

```
