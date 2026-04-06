const katex = require("katex");

const latex = process.argv[2];
const displayMode = process.argv[3] === "true";

try {
    const html = katex.renderToString(latex, {
        displayMode: displayMode,
        throwOnError: false
    });
    console.log(html);
} catch (err) {
    console.error(err);
    process.exit(1);
}
