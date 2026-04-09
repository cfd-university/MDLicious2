import subprocess

def convert_latex_equation(equation, display_mode):
    katex2htmlJS = f'''
    const katex = require("katex");

    const latex = String.raw`{equation}`;
    const displayMode = "{str(display_mode).lower()}" === "true";

    try {{
        const html = katex.renderToString(latex, {{
            displayMode: displayMode,
            throwOnError: true
        }});
        console.log(html);
    }} catch (err) {{
        console.error(err);
        process.exit(1);
    }}
    '''

    result = subprocess.run(
        ["node", "-e", katex2htmlJS],
        capture_output=True,
        text=True
    )

    return result.stdout.strip()