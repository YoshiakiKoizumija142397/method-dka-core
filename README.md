# Method-DKA-Core: High-Precision Polynomial Root-Finding Engine

**Software Development Requirements & Specification Document (SRS / README) — v3.7.7**

---

## 🌐 MethodDka Ecosystem 共通クロスリンクバナー

> **MethodDka Ecosystem**: 本リポジトリは超高精度多項式根ソルバープロジェクトの一部です。
> 全6つのディストリビューションを展開しています：
> [📦 Core (CUI)](https://github.com/YoshiakiKoizumija142397/method-dka-core) |
> [🖥️ Win11 (GUI)](https://www.google.com/search?q=https://github.com/YoshiakiKoizumija142397/MethodDkaWindows11) |
> [🌐 Web (PWA)](https://github.com/YoshiakiKoizumija142397/MethodDka) |
> [📱 Android (Java)](https://github.com/YoshiakiKoizumija142397) |
> [📜 N88-BASIC](https://www.google.com/search?q=https://github.com/YoshiakiKoizumija142397/MethodDkaN88BASIC) |
> [📜 VB Express](https://github.com/YoshiakiKoizumija142397/MethodDkaVBExpressEdition)

---

## 1. プロジェクト概要と基本設計思想 (Overview & Philosophy)

`method-dka-core` は、デュランド・ケルナー・アバート法（Durand-Kerner-Aberth Method: DKA法）およびウンの平方分割アルゴリズム（Yun's Square-Free Decomposition）を中核とした、**任意次数・超高精度複素数多項式根ソルバーコア**です。

* **第一原理主義に基づく堅牢性**: 丸め誤差、桁落ち、オーバーフロー、および非数（NaN/Infinity）の発生を徹底的に排除するため、ソフトウェアベースの多倍長浮動小数点数演算（`decimal.js`）を採用しています。
* **完全な監査証跡（Auditability）**: 科学計算や数値解析における再現性を担保するため、入力された係数データから各根の個別収束イテレーション数、正規化残差、係数再構築残差に至るまでを1つの監査レポートファイル（`output.txt`）に自己文書化（セルフドキュメンテーション）します。

---

## 2. 精度階層基準 (Precision Hierarchy Standards)

多倍長演算における計算コストと丸め誤差のバランスを最適化するため、以下の3段階の精度階層を厳格に定義しています。

1. **前処理・因数分解フェーズ (`PRECISION_PREPROCESS`)**: **420桁**
* Yunの平方分割および高精度多項式GCD演算における中間誤差の伝播を防ぎ、重複度を完全に分離するための十分なバッファ精度。


2. **DKA法ソルバー内部演算 (`PRECISION_DKA`)**: **400桁**
* DKA反復における各根の相互作用およびニュートン・アバート補正項の計算精度。


3. **収束判定・出力フォーマット (`Tolerance` & `Output`)**: **180桁** ($10^{-180}$)
* 収束判定閾値は $10^{-180}$ とし、出力される実部・虚部および残差評価も180桁の有効数字フォーマットを標準とします。



---

## 3. アルゴリズム仕様と数学的基盤 (Mathematical Specifications)

### 3.1 Yunの平方分割による前処理（Square-Free Preprocessing）

* 係数がすべて整数の場合、BigIntを用いた厳格なユークリッドの互除法による多項式GCD（最大公約多項式）を計算し、$A(x) = \prod_{i=1}^k P_i(x)^i$ の形に因数分解します。
* これにより、重根を持つ多項式であっても各因数 $P_i(x)$ ごとに単根として独立して解くことができ、DKA法の収束性が劇的に向上します。
* 整数係数以外のフロート係数の場合は、420桁精度による多項式除算と動的トラッキング処理で平方無因子群を抽出します。

### 3.2 総合的な根のスケーリング（Cauchy's Bound）

* コーシーの界（Cauchy's Bound）に基づき、多項式の係数比率から根が存在する最大領域のスケール $R$ を自動推定し、DKA法の初期複素数配置円（Durand-Kerner Initial Circles）を安定的に構築します。

### 3.3 収束判定と残差評価

* **個別根の収束判定**: 各根の補正量 $\Delta z$ の絶対値が $10^{-180}$ を下回った時点で当該根をフリーズ（アーリーイグジット）し、その時点のイテレーション数を記録。
* **正規化残差（Normalized Residual）**: 各根 $z$ における残差を次式で評価します。

$$\text{Res}(z) = \frac{\vert{}P(z)\vert{}}{\sum_{k=0}^n \vert{}a_k\vert{} \vert{}z\vert{}^k}$$



---

## 4. 安全防壁と監査証跡仕様 (Safety Guards & Audit Trail)

### 4.1 CUIセイフティロック機構

* 毎イテレーションにおいて、根の実部・虚部および補正項が有限値（`isFinite()`）であるかを常時監視します。
* 一時的な数値の揺らぎやオーバーフローの兆候を検知した場合、即座に安全なスケール位置へ引き戻すサイレント補正を行いつつ、介入カウンターをインクリメントします。
* 介入回数が上限（`MAX_SAFETY_INTERVENTIONS = 500`）を超えた場合にのみ、致命的エラーとして処理を安全にアボートします。

### 4.2 監査レポートフォーマット (`output.txt`)

計算完了時に自動生成されるレポートの構造は以下の通りです。

```text
==================================================
 METHOD-DKA-CORE EXECUTION AUDIT REPORT (v3.7.7)
==================================================
Timestamp: 2026-11-14T07:18:02.997Z
Internal Precision: 400 digits (Preprocess: 420 digits)
Execution Time: 142890.304 ms
Safety Lock Interventions Triggered: 0 times (Limit: 500)
Max Coefficient Residual: 1.5203333889e-387
--------------------------------------------------
Input Coefficients:
1, -210, 20615, -1256850, 53327946, ...
--------------------------------------------------

Root #1: Re=20.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000, Im=-0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000, Res=3.424672e-401, Iter=412
...
==================================================

```

---

## 5. モジュール構成およびファイル仕様 (Module Structure)

* **`input.txt`**: 計算対象の係数データ（降べきの順、スペースまたはカンマ区切り）。
* **`index.js`**: 複素数ヘルパー、多項式代数、BigInt GCD、DKAソルバーカーネル、およびCLIファイルI/Oを統合したエントリーポイント。
* **`output.txt`**: 監査証跡を記録する出力ファイル。

---

## 6. 実行手順 (Execution Instructions)

1. 依存パッケージのインストール:
```bash
npm install decimal.js

```


2. 係数データの配置 (`input.txt` に記述)
3. ソルバーの実行:
```bash
node index.js

```



---

## 7. ライセンス (License)

MIT License — 開発者: 小泉嘉章 (Yoshiaki Koizumi)