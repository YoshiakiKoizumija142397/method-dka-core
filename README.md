```markdown
# MethodDka Core (v3.4.0)

High-precision polynomial root-finding library and CUI solver using the Durand-Kerner-Aberth (DKA) method, enhanced with Yun's Square-Free Decomposition, 320-digit preprocessing, and `decimal.js`.

---

## 🇯🇵 日本語版 (Japanese)

### 1. 概要 (Overview)

`method-dka-core` は、デュランド・カーネル・アバース法（DKA法）と `decimal.js` を用いた多項式数値解法モジュールおよび CUI ベンチマークエンジンです。
内部演算には Decimal.js の **260桁精度** を用い、実部は有効数字200桁、虚部は $10^{-200}$ の単位で計算します。**※この出力桁数は正しい有効桁数を保証するものではありません。**

**v3.4.0** では、「ユンの平方無因子分解」によるアルゴリズムの根本的な堅牢化と、外部ベンチマークツール連携用の「`input.txt` 優先 I/O インターフェース」を導入しました。整数係数は BigInt による厳密演算で、それ以外の係数は 320 桁（$\epsilon = 10^{-220}$）で処理し、異なる重複度を持つ根を重複度ごとの平方無因子多項式 $Q_k(x)$ に分離します。検証済みの各因子を既存の DKA エンジンへ個別に渡し、重複度で重み付けした次数の再構成を検査します。検証に失敗した場合は診断ログを出力して元の多項式へフォールバックするセーフティガードを備えています。

---

### 2. 📊 驚異的な検証実績 (Benchmark Evidence)

本ソルバーの実用限界を証明するテストとして、**200次のウィルキンソン多項式を200桁精度で完全分離・算出**することに成功しています。この1,025行に及ぶ計算結果の完全な検証ログとベンチマークデータは、以下のGitHub Gistにて公開されています。

👉 **[200次ウィルキンソン多項式の検証結果 (GitHub Gist)]( https://gist.github.com/YoshiakiKoizumija142397/649c89263a53df38e4f6cfbe8637883a )**

---

### 3. 📋 システム仕様 (Technical Specifications)

| 項目 | 詳細仕様 (v3.4.0) |
| --- | --- |
| **コア解法** | DKA法 (Durand-Kerner-Aberth Simultaneous Root-Finding Method) |
| **演算精度階層** | **階層型マルチプレシジョン設計**<br>・前処理フェーズ: 整数係数は **BigInt厳密演算**、その他は **320桁固定精度** ($\epsilon = 10^{-220}$)<br>・後工程ソルバー: **260桁内部演算**（60桁のガード桁、スケール後の補正量閾値 $\epsilon = 10^{-220}$） |
| **精度の解釈** | 出力される正規化残差は後退誤差の数値指標です。例えば残差 $10^{-180}$ は「根が180桁正しい」ことを意味せず、有効桁数を認証するものではありません。あくまで数値検査であり数学的証明ではありません。 |
| **自動前処理 (GCD)** | **ユンの平方無因子分解 (Yun's Square-free Decomposition)**<br>・微分多項式 $P'(x)$ とのユークリッド互除法による $\gcd(P, P')$ の算出<br>・整数係数では BigInt 擬似剰余列を使った厳密GCD、それ以外では 320桁精度と動的多項式除算 (`polyDivide`) を用いて重複度別の平方無因子群 $Q_k(x)$ を抽出。 |
| **セーフティガード機構** | 各多項式除算の余り、因子の平方無因子性、重複度を反映した再構成次数を検証。失敗時は元多項式をそのままDKAへ渡す（フォールバック機能）。 |
| **数値安定化** | 非整数係数の前処理と各DKA因子をそれぞれモニック化し、Fujiwara根上界に基づく変数スケール $x=S y$ を適用。 |
| **多項式評価** | ホーナー法 (Horner's Method) による高精度多項式直評価。 |
| **対応最高次数** | 最高 200 次 ($n \le 200$) |

---

### 4. 入出力仕様 (Input & Output Specification)

本モジュールは、ライブラリAPI呼び出しおよび CUI/CLI 実行の両方に対応します。CUI実行時（`node index.js`）は、外部ツールとの連携のため以下の優先順位で入出力を行います。

#### ① 入力仕様 (Input Specification)
1. **`input.txt` による外部入力（最優先）**:
   - カレントディレクトリに `input.txt` が存在する場合、外部ジェネレータ（`method-dka-roots-generator-js` 等）からの指示と判断し、最優先で読み込みます。
   - **フォーマット**: 改行区切りのテキストデータ。高次から低次の順（例: $a_n, a_{n-1}, \dots, a_0$）。浮動小数点変換誤差を防ぐため、文字列としてパースされます。
2. **内蔵ジェネレータ（フォールバック）**:
   - `input.txt` が存在しない場合は、単独テスト用として内蔵の係数生成ロジック（ウィルキンソン多項式等）を自動起動して処理を実行します。

#### ② 出力仕様 (Output Specification)
計算結果は標準出力（コンソール）に出力されると同時に、同ディレクトリの `output.txt` に自動保存されます。

**出力フォーマット例 (`output.txt`):**
```text
Execution Time: 450.123 ms
Max Coefficient Residual: 3.8907395275e-247

Root #1: Re=0.50000000000000000000, Im=0.00000000000000000000, Res=1.234567e-260
Root #2: Re=0.50003051757812500000, Im=0.00000000000000000000, Res=1.234567e-260
...

```

---

### 5. 🚀 主な機能と特徴

1. **ユンの平方無因子分解**:
整数係数は BigInt による厳密なGCDと多項式除算で、それ以外は 320桁精度で処理し、同じ重複度を持つ根の因子群 $Q_k(x)$ を分離します。各群の平方無因子性と次数の再構成を厳密に検証します。
2. **動的トラッキング型多項式長除法**:
歯抜けの多項式や、演算過程で次数が激しく変動する場合でも、常に最高次数を監視して次数差 (`shift`) を動的に計算。商の正しい位置に係数を書き込み、インデックスのズレを完全に防ぎます。
3. **検証付きフォールバック**:
割り算の余り、平方無因子性、再構成次数のいずれかが検証を通らない場合、元の多項式を安全に既存DKAへ渡します。正常に分解できた混在重根は重複度別に個別計算します。
4. **オートスケーリング付きDKA**:
非整数係数の前処理と各平方無因子因子の DKA において、Fujiwara 上界 $S=2\max_i \vert a_i \vert^{1/i}$ を使って $x=S y$ に変数変換を自動適用します。

---

### 6. インストールと使用方法

#### インストール

```bash
npm install method-dka-core

```

*(※依存パッケージとして `decimal.js` が必要です)*

#### 使い方 A: CUI ベンチマーク実行 (外部ファイル連携)

`input.txt` に係数（高次から低次の順）を配置し、ターミナルで実行します。

```bash
node index.js

```

#### 使い方 B: Node.js ライブラリとしての呼び出し

```javascript
const { solvePolynomial } = require('method-dka-core');

// 係数を高次から順に配列で指定 [a_n, ..., a_0]
// 精度を維持するため、数値ではなく文字列での指定を推奨します。
const coefficients = ["1", "-1", "-18", "52", "-40"];

// solvePolynomial(係数配列, 最大反復回数, デバッグ出力フラグ)
const result = solvePolynomial(coefficients, 5000, false);

```

---

### 7. 変更履歴 (Changelog)

#### v3.4.0

* **仕様拡張**: 外部ベンチマークジェネレータとの自動連携のため、`input.txt` 優先読み込みおよび `output.txt` 出力仕様を新設。
* **機能追加**: ユンの平方無因子分解 (Yun's Square-free Decomposition) を実装し、重根分離性能を向上。
* **機能追加**: 320桁前処理フェーズおよび BigInt 厳密演算 GCD を導入。
* **機能追加**: 動的トラッキング型多項式除算 (`polyDivide`) と Fujiwara オートスケーリングを搭載。

---

## 🇬🇧 English Version

### Overview

`method-dka-core` is a high-precision polynomial root-finding module and CUI benchmark engine operating in Node.js. It uses **260-digit internal arithmetic** powered by `Decimal.js`. **Note: Output length does not guarantee correct significant digits.**

Version 3.4.0 introduces **Yun's square-free decomposition** for robust multiple-root separation, and a new **`input.txt` prioritized File I/O interface** for automated benchmarking tools.

### 📊 Benchmark Evidence

As definitive proof of the solver's practical limits, the engine has successfully resolved a **200th-degree Wilkinson polynomial to 200-digit precision**. The complete 1,025-line verification logs and benchmark dataset are publicly available on GitHub Gist:

👉 **[200th-degree Wilkinson Benchmark on GitHub Gist](https://www.google.com/search?q=%5B%E3%81%93%E3%81%93%E3%81%ABGist%E3%81%AEURL%E3%82%92%E8%B2%BC%E3%82%8A%E4%BB%98%E3%81%91%E3%81%A6%E3%81%8F%E3%81%A0%E3%81%95%E3%81%84%5D)**

### Technical Specifications

| Item | Detailed Specifications |
| --- | --- |
| **Core Solver** | DKA Method (Durand-Kerner-Aberth Simultaneous Root-Finding Method) |
| **Precision Hierarchy** | **Hierarchical Multi-precision Design**<br>

<br>・ Preprocessing Phase: **Exact BigInt arithmetic** for integer coefficients; otherwise **fixed 320-digit precision** ($\epsilon = 10^{-220}$)<br>

<br>・ Post-process Solver: **260-digit internal precision** (60 guard digits; scaled correction threshold $\epsilon = 10^{-220}$) |
| **Precision Interpretation** | A normalized residual is a numerical backward-error indicator. It is a numerical check, not a mathematical proof. |
| **Input / Output** | Prioritizes `input.txt` (line-delimited coefficients) in the working directory for CLI executions (`node index.js`) and outputs results to `output.txt`. Falls back to built-in generators if missing. |
| **Supported Max Degree** | Up to 200th degree ($n \le 200$). |

---

📜 **ライセンス ＆ 開発者情報 (License & Developer Information)**

* **開発者 (Developer)**: 小泉嘉章 (Yoshiaki Koizumi)
* **バージョン (Version)**: v3.4.0
* **ライセンス (License)**: MIT License
* **公式リポジトリ (Repository)**: [MethodDka Repository](https://github.com/YoshiakiKoizumija142397/method-dka-core)

```