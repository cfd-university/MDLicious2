# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Markdown → HTML converter that produces **WordPress-flavoured HTML** for articles on cfd.university. It wraps `markdown2` and adds LaTeX-style numbered/cross-referenceable figures, tables, code listings and equations, plus privacy-enhanced YouTube embeds. Output carries WordPress block classes (`wp-block-image`, `wp-block-table`, `wp-element-caption`, `wp-block-quote`) — keep these when touching HTML generation.

The full user-facing syntax reference lives in `README.md`; read it before changing any component's parsing.

## Running

```bash
python MDLicious2.py testScripts/example1.json      # config-file is the only CLI argument
mkdir out                                            # output dir must already exist; out/ is gitignored
```

Config JSON: `inputFile` and `outputDirectory` are required; `replace` (a plain find→replace map applied line-by-line before any parsing) is optional.

There is no test suite and no linter. `testScripts/*.md` + matching `*.json` are the manual regression corpus — `example1..3` are stable samples, `temp` and `trouble` are scratch files. Verify a change by converting a sample and diffing the HTML plus `out/stderr.json`.

`.vscode/launch.json` debugs `testScripts/temp.json`.

Setup: `pip install -r requirements.txt`, plus Node.js with `npm install katex` — equation rendering shells out to `node -e` (see `javascriptRuntime.py`), so `require("katex")` must resolve from the working directory. `Dockerfile` builds an image with both toolchains.

The README version badge is bumped one patch level per functional commit. Commits use Conventional Commits (`fix(scope): ...`).

## Pipeline (MDLicious2.py)

The `main()` function is the whole orchestration, and the order is load-bearing:

1. `CommandLineArguments` → `FileProcessor` — reads the markdown as a **list of lines** and applies `replace`.
2. `CaptionMatcher.setup_equation_tags` — scans `$$…$$` blocks and *mutates the line list*, injecting a synthetic `\tag{eq:equation-N}` where the author gave none.
3. `CaptionMatcher.setup_ref_map` — walks the raw markdown and builds `counter_map: tag → number` for all four component types, before any HTML exists.
4. `ComponentManager` + `Preprocessor` — custom components convert their own regions straight to HTML; everything else passes through untouched.
5. `Mark2HTML.convert` — standard markdown2 conversion of the mixed markdown/HTML string, plus post-processing.
6. `CaptionMatcher.substitute` — replaces `\ref{...}` (and KaTeX-mangled equation tags) with numbers from `counter_map`.
7. `CheckManager` — validators over input markdown and output HTML; writes `stderr.json` into `outputDirectory`.
8. `FileProcessor.output` — writes `<inputbasename>.html` into `outputDirectory`.

### Line-list vs. string

Content is a **list of lines** everywhere except between steps 4 and 5: `Preprocessor.processed_content` is a single joined string, and `Mark2HTML.convert()` splits back to a list before returning. Components, `CaptionMatcher`, and every check index into line lists. Getting this wrong is the most common source of breakage.

### Two independent counters

Numbering exists twice and the two must agree:

- `CaptionMatcher.counter` / `counter_map` — built from raw markdown, resolves `\ref{}` to numbers.
- `CaptionExtractor.counter` — each `Component` owns its own instance and emits the visible "Figure N" / "Table N" / "Listing N" caption text during conversion.

If you change what counts as a figure/table/listing, update the detection logic in **both** `CaptionMatcher.__is_*` and the component's `match()`, or references will point at the wrong numbers. Code listings are the asymmetric case: `CaptionExtractor` only numbers a listing when a caption or `\tag` is present, and `CaptionMatcher` only registers one when the preceding line has a `\tag`.

Tag prefixes determine type and are not interchangeable: `fig:`, `tab:`, `code:`, `eq:` (enforced by `RefCheck`).

## Components

`MDLicious2/components/` — each subclasses `Component` (`base.py`) and implements:

- `match(index)` — does the line at `index` start this construct?
- `convert(index)` — return the HTML for the whole construct.
- `increment` — how many source lines were consumed; `Preprocessor` advances by it. **Must be set inside `convert()`** (or via `_find_start_end_based_on_pattern`, which sets it from the closing delimiter), otherwise the block is re-parsed line by line.

`Preprocessor` takes the **first** matching component, so registration order in `main()` matters — equations before code before figures/tables, since `$$` and ```` ``` ```` fences would otherwise be swallowed by generic handling.

To add one: create the class, export it from `MDLicious2/__init__.py`, and register it in `main()`. If it needs numbering, add a `ComponentType` and wire it into both counters described above.

## Gotchas

- **BeautifulSoup re-serialisation**: `Mark2HTML` uses `soup.decode(formatter="minimal")` — a bare `decode()` re-escapes `&lt;` into `<` inside code blocks and has regressed before (commit feeb9d8). `Table.convert` deliberately uses `formatter=None`. BS4 also strips newlines, which is why `__prettify` re-inserts them after block elements.
- **WordPress shortcodes**: `[custom_category_posts_list category_slug="…"]` lines are stashed as placeholders before markdown2 runs and restored afterwards, otherwise they get wrapped in `<p>`.
- **KaTeX tags survive as markup**: after rendering, an equation tag can appear as `<mtext>(eq:…)`, a literal `\tag{eq:…}`, or `<span class="mord">eq:…</span>`. `CaptionMatcher.substitute` handles all three; adding a KaTeX version bump may add a fourth.
- **Inline equations run before markdown2** so it cannot mangle LaTeX; each equation is a separate `node` subprocess, so equation-heavy documents are slow by design.
- **Code styling**: Pygments `nord` style with `linenos="table"`. `code.py` has a commented-out block that regenerates `code.css` — the CSS is not emitted at runtime and is maintained by hand alongside the site theme.
- **Checks are advisory except by convention**: `CheckManager` never raises; it writes warnings/errors to `out/stderr.json`. `RefCheck` flags any `\ref` or malformed `\tag` still present in the final HTML, which is the signal that numbering broke.
