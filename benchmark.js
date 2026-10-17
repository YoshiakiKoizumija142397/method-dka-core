/**
 * CUI Benchmark Script for MethodDka
 * Reads polynomial coefficients from input.txt, solves them using method-dka-core,
 * and measures execution time and precision metrics.
 */

const fs = require('fs');
const { solvePolynomial } = require('./index.js');

// コマンドライン引数の処理
const args = process.argv.slice(2);
const inputFile = args[0] || 'input.txt';
const outputFile = args[1] || null;

console.log("==================================================");
console.log(` MethodDka Benchmark Engine (v3.4.0)             `);
console.log("==================================================");
console.log(`[入力ファイル] ${inputFile}`);

if (!fs.existsSync(inputFile)) {
    console.error(`[エラー] 入力ファイルが見つかりません: ${inputFile}`);
    console.error(`ヒント: 先に 'npm run gen:wilkinson' 等で係数ファイルを生成してください。`);
    process.exit(1);
}

try {
    // 1. input.txt から係数を読み込み
    const fileContent = fs.readFileSync(inputFile, 'utf-8');
    const coeffs = fileContent
        .split('\n')
        .map(line => line.trim())
        .filter(line => line !== '' && !line.startsWith('//'));

    console.log(`[読み込み完了] 係数数: ${coeffs.length} (次数: ${coeffs.length - 1})`);
    console.log(`[ステータス] DKAコアエンジンによる高精度解法を実行中...`);

    // 2. 実行時間の計測開始
    const startTime = process.hrtime.bigint();

    // ソルバーの実行 (最大反復500、デバッグログなし)
    const result = solvePolynomial(coeffs, 500, false);

    // 3. 実行時間の計測終了 (ミリ秒換算)
    const endTime = process.hrtime.bigint();
    const elapsedMs = Number(endTime - startTime) / 1_000_000;

    console.log(`\n================ Benchmark Results ================\n`);
    console.log(` 処理時間 (Execution Time):  ${elapsedMs.toFixed(3)} ms`);
    console.log(` 抽出された因子グループ数:     ${result.factorCount}`);
    console.log(` 最大係数再構成残差:         ${result.maxCoefficientResidual.toExponential(10)}`);
    console.log(`--------------------------------------------------`);

    // 上位の根の詳細表示
    const displayLimit = Math.min(result.roots.length, 25);
    console.log(` 算出された根 (上位 ${displayLimit} 件表示):`);
    
    result.roots.forEach((root, idx) => {
        const sign = root.im.isNegative() ? '' : '+';
        console.log(`   [#${idx + 1}] Re: ${root.re.toFixed(15)}`);
        console.log(`        Im: ${sign}${root.im.toFixed(15)}i | Res: ${root.residual.toExponential(4)} | Iter: ${root.iterations}`);
    });

    if (result.roots.length > displayLimit) {
        console.log(`   ... 他 ${result.roots.length - displayLimit} 件の根が正常に算出されています。`);
    }
    console.log(`==================================================`);

    // 4. 出力ファイルが指定されている場合は結果を書き出し
    if (outputFile) {
        let outputData = `Execution Time: ${elapsedMs.toFixed(3)} ms\n`;
        outputData += `Max Coefficient Residual: ${result.maxCoefficientResidual.toExponential(10)}\n\n`;
        result.roots.forEach((root, idx) => {
            outputData += `Root #${idx + 1}: Re=${root.re.toFixed(20)}, Im=${root.im.toFixed(20)}, Res=${root.residual.toExponential(6)}\n`;
        });
        fs.writeFileSync(outputFile, outputData, 'utf-8');
        console.log(`[出力完了] ベンチマーク結果を ${outputFile} に保存しました。`);
    }

} catch (err) {
    console.error(`[エラー] ベンチマーク実行中に例外が発生しました:`, err);
    process.exit(1);
}