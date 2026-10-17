/**
 * CUI Coefficient Data Generator for MethodDka
 * Generates polynomial coefficients from specified roots and outputs to a text file.
 * Requires: decimal.js (npm install decimal.js)
 */

const fs = require('fs');
const Decimal = require('decimal.js');

// 精度設定
Decimal.set({ precision: 100, rounding: Decimal.ROUND_HALF_UP });

/**
 * 根のリストから多項式の係数（降べきの順: [a_n, a_{n-1}, ..., a_0]）を展開する
 * P(x) = \prod (x - r_i)
 */
function expandPolynomialFromRoots(roots) {
    // 初期状態: P(x) = 1 (係数配列: [1])
    let poly = [new Decimal(1)];

    for (const r of roots) {
        const nextPoly = new Array(poly.length + 1).fill(null).map(() => new Decimal(0));
        const negR = { re: r.re.negated(), im: r.im.negated() }; // 配列展開用 (-r)

        for (let i = 0; i < poly.length; i++) {
            // 次数が1上がるシフト
            nextPoly[i + 1] = nextPoly[i + 1].plus(poly[i]);
            // 定数項側の掛け算
            const term = poly[i].times(negR.re); // 実数根を基本とする拡張
            nextPoly[i] = nextPoly[i].plus(term);
        }
        poly = nextPoly;
    }
    return poly;
}

// --- CUI インターフェイス処理 ---
const args = process.argv.slice(2);
// 使用例: node coeff-generator.js input.txt wilkinson 20
// または: node coeff-generator.js input.txt custom 2,2,2,-5 (重根・カスタム根)

if (args.length < 2) {
    console.error("使用法:");
    console.error("  1. ウィルキンソン多項式: node coeff-generator.js <出力ファイル> wilkinson <次数N>");
    console.error("  2. カスタム根指定:       node coeff-generator.js <出力ファイル> custom <根1,根2,...>");
    process.exit(1);
}

const outputFile = args[0];
const mode = args[1];

let targetRoots = [];

if (mode === 'wilkinson') {
    const n = parseInt(args[2] || "20", 10);
    console.log(`[生成モード] ウィルキンソン多項式 (W_${n}(x), 根: 1 から ${n})`);
    for (let i = 1; i <= n; i++) {
        targetRoots.push({ re: new Decimal(i), im: new Decimal(0) });
    }
} else if (mode === 'custom') {
    const rawRoots = (args[2] || "").split(',').map(s => s.trim()).filter(s => s !== '');
    if (rawRoots.length === 0) {
        console.error("[エラー] 有効な根が指定されていません。");
        process.exit(1);
    }
    console.log(`[生成モード] カスタム根 (指定数: ${rawRoots.length})`);
    targetRoots = rawRoots.map(val => ({ re: new Decimal(val), im: new Decimal(0) }));
} else {
    console.error(`[エラー] 不明なモードです: ${mode}`);
    process.exit(1);
}

try {
    // 多項式係数の展開実行
    const coeffsAsc = expandPolynomialFromRoots(targetRoots);
    // 降べきの順（最高次 → 定数項）に反転
    const coeffsDesc = [...coeffsAsc].reverse();

    // ファイル出力フォーマット（1行に1係数）
    const fileContent = coeffsDesc.map(c => c.toFixed(10)).join('\n');
    fs.writeFileSync(outputFile, fileContent, 'utf-8');

    console.log(`[成功] 係数データ（次数: ${targetRoots.length}）を ${outputFile} に書き出しました。`);
} catch (err) {
    console.error(`[エラー] 係数生成に失敗しました: ${err.message}`);
    process.exit(1);
}