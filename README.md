# MethodDka Core (v3.4.0)

High-precision polynomial root-finding library using the Durand-Kerner-Aberth (DKA) method, enhanced with Yun's Square-Free Decomposition and `Decimal.js`.

## 🇯🇵 日本語版 (Japanese)

### 概要

`method-dka-core` は、デュランド・カーネル・アバース法（DKA法）と `decimal.js` を用いた多項式数値解法モジュールです。
内部演算には Decimal.js の **260桁精度** を用い、実部は有効数字200桁、虚部は $10^{-200}$ の単位で計算します。**この出力桁数は正しい有効桁数を保証するものではありません。**

**v3.4.0** ではユンの平方無因子分解を導入しました。整数係数はBigIntによる厳密演算で、それ以外の係数は320桁（$\epsilon = 10^{-220}$）で処理し、異なる重複度を持つ根を重複度ごとの平方無因子多項式 $Q_k(x)$ に分離します。検証済みの各因子を既存のDKAエンジンへ個別に渡し、重複度で重み付けした次数の再構成を検査します。検証に失敗した場合は診断ログを出力して元の多項式へフォールバックします。

### 📋 システム仕様 (Technical Specifications)

| 項目 | 詳細仕様 (v3.4.0 拡張) |
| --- | --- |
| **コア解法** | DKA法 (Durand-Kerner-Aberth Simultaneous Root-Finding Method) |
| **演算精度階層** | **階層型マルチプレシジョン設計**<br>

<br>・前処理フェーズ: 整数係数は **BigInt厳密演算**、その他は **320桁固定精度** ($\epsilon = 10^{-220}$)<br>

<br>・後工程ソルバー: **260桁内部演算**（60桁のガード桁、スケール後の補正量閾値 $\epsilon = 10^{-220}$。出力桁数は正しい有効桁数を保証しません） |
| **精度の解釈** | 出力される正規化残差は後退誤差の数値指標です。例えば残差 $10^{-180}$ は「根が180桁正しい」ことを意味せず、有効桁数を認証するものでもありません。あくまで数値検査であり数学的証明ではありません。 |
| **自動前処理 (GCD)** | **ユンの平方無因子分解 (Yun's Square-free Decomposition)**<br>

<br>・微分多項式 $P'(x)$ とのユークリッド互除法による $\gcd(P, P')$ の算出<br>

<br>・整数係数ではBigInt擬似剰余列を使った厳密GCD、それ以外では320桁精度と動的多項式除算 (`polyDivide`) を用いて重複度別の平方無因子群 $Q_k(x)$ を抽出<br>

<br>・重複度 $m$ の自動計算と結果の展開処理 |
| **セーフティガード機構** | 各多項式除算の余り、因子の平方無因子性、重複度を反映した再構成次数を検証。失敗時は元多項式をDKAへ渡す（フォールバック）。 |
| **数値安定化** | 非整数係数の前処理と各DKA因子をそれぞれモニック化し、Fujiwara根上界に基づく変数スケール $x=S y$ を適用。 |
| **多項式評価** | ホーナー法 (Horner's Method) による高精度多項式直評価。 |
| **対応最高次数** | 最高 200 次 ($n \le 200$)。 |

### 🚀 主な機能と特徴

1. **ユンの平方無因子分解**:
整数係数はBigIntによる厳密なGCD・多項式除算で、それ以外は320桁精度と $\epsilon = 10^{-220}$ で処理し、同じ重複度を持つ根の因子群 $Q_k(x)$ を分離します。各群の平方無因子性と元の次数の再構成を検証します。
2. **動的トラッキング型多項式長除法**:
歯抜けの多項式や次数が激しく変動する演算過程においても、常に最高次数を監視し次数差 (`shift`) を動的に計算して商の正しい位置に係数を書き込む最強の除算アルゴリズムを実装。インデックスのズレを完全に防ぎます。
3. **検証付きフォールバック**:
割り算の余り、平方無因子性、再構成次数のいずれかが検証を通らない場合は、元多項式を既存DKAへ渡します。正常に分解できた混在重根は重複度別に個別計算します。
4. **オートスケーリング付きDKA**:
非整数係数の前処理と各平方無因子因子のDKAで、Fujiwara上界 $S=2\max_i \vert{}a_i\vert{}^{1/i}$ を使って $x=S y$ に変数変換します。因子分解後は因子を元の変数へ戻し、DKAは各因子を再スケールして解き、根と最終補正量を元の尺度へ戻します。

### インストール

```bash
npm install method-dka-core

```

*(※依存パッケージとして `decimal.js` が必要です)*

### 使い方

```javascript
const { solvePolynomial } = require('method-dka-core');

// 係数を高次から順に配列で指定 [a_n, ..., a_0]
// 精度を維持するため、数値ではなく「文字列」での指定を推奨します。
// 例: P(x) = x^4 - x^3 - 18x^2 + 52x - 40
const coefficients = ["1", "-1", "-18", "52", "-40"];

// 実行: solvePolynomial(係数配列, 最大反復回数, デバッグ出力の有無)
const result = solvePolynomial(coefficients, 5000, false);

console.log(`因子グループ数: ${result.factorCount}`);
console.log(`最大係数再構成残差: ${result.maxCoefficientResidual.toExponential(4)}`);

result.roots.forEach((root, index) => {
    const sign = root.im.isNegative() ? '' : '+';
    console.log(`[Root ${index + 1}]`);
    console.log(`  解: ${root.re.toFixed(200)} ${sign} ${root.im.toFixed(200)}i`);
    console.log(`  重複度: ${root.multiplicity}`);
    console.log(`  反復回数: ${root.iterations}`);
    console.log(`  正規化残差: ${root.residual.toExponential(4)}`);
    console.log(`  最終補正量: ${root.finalCorrection.toExponential(4)}\n`);
});

```

---
🌐 公式ページ ＆ リポジトリ
Web アプリ (Live Demo): MethodDka Live Demo
GitHub リポジトリ: MethodDka Repository


📜 ライセンス ＆ 開発者情報
開発者: 小泉嘉章 (Yoshiaki Koizumi)
バージョン: v3.4.0
ライセンス: MIT License



## 🇬🇧 English Version

### Overview

`method-dka-core` is a polynomial numerical solver module operating in Node.js. It uses **260-digit internal arithmetic** powered by `Decimal.js`. **Output length does not guarantee correct significant digits.**

Version 3.4.0 adds **Yun's square-free decomposition**, using exact BigInt arithmetic for integer coefficients and 320-digit preprocessing ($\epsilon = 10^{-220}$) otherwise. It separates roots with different multiplicities into square-free factors $Q_k(x)$ and sends each verified factor to the existing DKA engine. Factor square-freeness and reconstructed degree are validated; failures fall back to the original polynomial.

### 📋 System Specifications

| Item | Detailed Specifications |
| --- | --- |
| **Core Solver** | DKA Method (Durand-Kerner-Aberth Simultaneous Root-Finding Method) |
| **Arithmetic Precision Hierarchy** | **Hierarchical Multi-precision Design**<br>

<br>・ Preprocessing Phase: **Exact BigInt arithmetic** for integer coefficients; otherwise **fixed 320-digit precision** ($\epsilon = 10^{-220}$)<br>

<br>・ Post-process Solver: **260-digit internal precision** (60 guard digits; scaled correction threshold $\epsilon = 10^{-220}$; output length does not guarantee correct digits) |
| **Precision Interpretation** | A normalized residual is a numerical backward-error indicator; for example, a residual of $10^{-180}$ does not mean that 180 digits of the root are correct and does not certify any number of correct digits. It is a numerical check, not a mathematical proof. |
| **Automatic Preprocessing (GCD)** | **Yun's Square-free Decomposition**<br>

<br>・ Calculation of $\gcd(P, P')$ using the Euclidean algorithm with derivative polynomial $P'(x)$<br>

<br>・ Extraction of multiplicity-specific square-free factors $Q_k(x)$ using exact BigInt pseudo-remainder sequences for integer coefficients, or 320-digit arithmetic and dynamic polynomial division (`polyDivide`) otherwise. |
| **Safety Guard Mechanism** | Validates division remainders, square-freeness of each factor, and reconstructed degree. On failure, passes the original polynomial to DKA. |
| **Numerical Stabilization** | Non-integer preprocessing and each DKA factor are scaled independently with $x=S y$, using a Fujiwara root bound. |
| **Polynomial Evaluation** | High-precision direct polynomial evaluation via Horner's Method. |
| **Supported Max Degree** | Up to 200th degree ($n \le 200$). |

### 🚀 Key Features

1. **Yun Square-Free Decomposition**:
Exact BigInt GCD and division are used for integer coefficients; otherwise 320-digit precision and $\epsilon = 10^{-220}$ are used. Root groups $Q_k(x)$ are separated by multiplicity, then square-freeness and degree reconstruction are checked.
2. **Dynamic Tracking Polynomial Long Division**:
Implements a powerhouse division algorithm that constantly monitors the highest degree and dynamically calculates degree differences (`shift`) to write coefficients into the correct index positions of the quotient coefficient array `Q`.
3. **Intelligent Safety Guard**:
If division remainders, factor square-freeness, or degree reconstruction fail validation, the original polynomial is passed to the existing DKA engine.
4. **Auto-scaled DKA**:
Non-integer preprocessing and each square-free factor are made monic and scaled with the Fujiwara bound $S=2\max_i \vert{}a_i\vert{}^{1/i}$ using $x=S y$.

### Installation

```bash
npm install method-dka-core

```

*(Requires `decimal.js` as a dependency)*

### Usage

```javascript
const { solvePolynomial } = require('method-dka-core');

// Specify polynomial coefficients from highest degree to constant [a_n, ..., a_0].
// It is highly recommended to use strings to prevent floating-point precision loss.
const coefficients = ["1", "-1", "-18", "52", "-40"];

// Execute: solvePolynomial(coefficients_array, max_iterations, debug_flag)
const result = solvePolynomial(coefficients, 5000, false);

console.log(`Factor Groups Extracted: ${result.factorCount}`);
console.log(`Max Reconstruction Residual: ${result.maxCoefficientResidual.toExponential(4)}`);

result.roots.forEach((root, index) => {
    const sign = root.im.isNegative() ? '' : '+';
    console.log(`[Root ${index + 1}]`);
    console.log(`  Root:         ${root.re.toFixed(200)} ${sign} ${root.im.toFixed(200)}i`);
    console.log(`  Multiplicity: ${root.multiplicity}`);
    console.log(`  Iterations:   ${root.iterations}`);
    console.log(`  Residual:     ${root.residual.toExponential(4)}`);
    console.log(`  Correction:   ${root.finalCorrection.toExponential(4)}\n`);
});

---

🌐 Official Page & Repository
Web App (Live Demo): MethodDka Live Demo
GitHub Repository: MethodDka Repository

📜 License & Developer Information
Developer: Yoshiaki Koizumi
Version: v3.4.0
License: MIT License
