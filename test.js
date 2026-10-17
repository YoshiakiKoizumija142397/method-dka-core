// test.js
const { solvePolynomial } = require('./index.js');

// テストケース: P(x) = (x - 2)^3 * (x + 5) = x^4 - x^3 - 18x^2 + 52x - 40
// 係数配列 (降べきの順)
const coeffs = ["1", "-1", "-18", "52", "-40"];

console.log("==================================================");
console.log("テスト実行: P(x) = x^4 - x^3 - 18x^2 + 52x - 40 ...");
console.log("期待される解: 2 (重複度3), -5 (重複度1)");
console.log("==================================================\n");

try {
    const result = solvePolynomial(coeffs, 500, true);

    console.log("\n================ Solution Summary ================");
    console.log(`抽出された因子グループ数: ${result.factorCount}`);
    console.log(`最大係数再構成残差: ${result.maxCoefficientResidual.toExponential(10)}`);
    console.log("--------------------------------------------------");

    result.roots.forEach((root, idx) => {
        const sign = root.im.isNegative() ? '' : '+';
        console.log(`Root #${idx + 1}:`);
        console.log(`  Real:         ${root.re.toFixed(30)}`);
        console.log(`  Imaginary:    ${sign}${root.im.toFixed(30)}i`);
        console.log(`  Multiplicity: ${root.multiplicity}`);
        console.log(`  Residual:     ${root.residual.toExponential(6)}`);
        console.log(`  Iterations:   ${root.iterations}`);
        console.log(`  Correction:   ${root.finalCorrection.toExponential(6)}\n`);
    });
} catch (error) {
    console.error("テスト実行中にエラーが発生しました:", error);
}