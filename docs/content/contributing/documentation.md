# Maintaining the docs

Write for someone trying to use or maintain the library. Explain what an
operation does, show it on a concrete problem, and state the conditions that
matter to the result. Keep the mathematics precise and the prose direct.

## Build and review

```sh
pixi run --locked -e docs docs-serve
pixi run --locked -e docs docs-build
pixi run --locked -e docs docs-check
pixi run --locked build
pixi run --locked docs-test --output .cache/docs-examples-001 --jobs 1
```

- **`docs-build`** writes the static site to `.cache/docs-site`.
- **`docs-check`** runs the audit, builds in strict mode, then checks every
  rendered local link, anchor and asset, and duplicate HTML ids.
- **`docs-test`** compiles and runs every included example against the prepared
  library. Reuse the existing build when only examples or prose changed.

## Where things live

| Path | Contents |
|---|---|
| `mkdocs.yml` | site configuration and navigation |
| `docs/content/` | the pages, in Markdown |
| `docs/examples/` | every program the site includes |
| `docs/examples.json` | the expected standard output of each included program |
| `docs/hooks.py` | MkDocs hooks: the version and Mojo series for the theme, the audit, directives, source downloads |
| `docs/theme/` | the theme shared by our projects |
| `scripts/docs_contract.py` | the audit and directive rendering |
| `scripts/mojo_api.py` | the reader of Mojo declarations and docstrings that the `api` directive renders |
| `docs/content/reference/` | topic-based reference pages for the public API |

Write tutorials, guides, and the API reference against `poly`. Examples that
use the core `polynomials` API belong on the Core API pages, which serve
developers extending the library.

## Mathematics

Use LaTeX notation for equations: `$x^2$` for inline math, and `$$` on separate
lines for a displayed equation, with a blank line before and after:

```text
$$
\frac{a}{b} + \frac{c}{d} = \frac{ad + bc}{bd}
$$
```

The Markdown extension preserves the notation, and MathJax renders it in
the browser; `\sem{t}` is defined as a macro for $\sem{t}$. Inside tables,
write `\lvert x\rvert` and `\mid` instead of `|`, which would end the cell.
Math inside code spans and fenced blocks stays literal. MathJax loads from a
CDN, so rendering equations requires an internet connection.

Number results per page, with the page's letter (Definition S-1, Theorem
C-6), and label each as on the [theory overview](../theory/index.md).

## Directives

```text
<!-- example: docs/examples/<file>.mojo -->
<!-- source: docs/examples/<file>.mojo -->
<!-- api: src/<package>/<module>.mojo <Name> -->
```

- **`example`** renders a complete program: its file name and download
  link, the code, the run command and the expected output. The program must
  be listed in `docs/examples.json`.
- **`source`** renders a file from `docs/examples/` that is not run on its
  own.

- **`api`** renders a public declaration of a package from its sources:
  a directive naming `src/polynomials/sparse.mojo` and `SparsePoly` gives the declaration's
  signature and docstring, its sections (`Args:`, `Parameters:`,
  `Returns:`, `Raises:`) as lists, and, for a struct, its fields and
  public methods. The reference pages use these directives to stay aligned
  with source declarations and docstrings.

## Rules the audit enforces

The build fails unless:

- every page starts with a `# ` title;
- no page contains a hand-written ` ```mojo ` block (use `text` fences for
  signatures and sketches);
- the set of `example` directives equals the set of programs in
  `docs/examples.json`, and each expected output ends with a newline;
- every `source` directive names a `.mojo` file in `docs/examples/`;
- every `api` directive names a public declaration with a docstring, every
  member it renders has one, and every name a package's `__init__.mojo`
  exports has an `api` directive.

## Adding an example

1. Write a complete program in `docs/examples/`, with a `main` that asserts
   its results with `std.testing`, and a one-line docstring that starts with
   a capital letter.
2. Run it, check its output independently, and record the path and exact
   output in `docs/examples.json`.
3. Include it on the page it illustrates with an `example` directive.
4. Run `docs-test` for it (`--only <name>`) and `docs-check`. Use a fresh
   output directory; the runner records the shared library build ID and tests
   packaged imports at O0 and O3.

Keep each tutorial example focused on one workflow. Assert mathematical results,
not only successful execution, and derive expected values independently of the
code under test. Explain what readers should notice immediately after the
program. Never present an uncompiled sketch as a runnable example.

## Theme

Navigation supports section groups with page links beneath them; keep it to
those two levels. Add every page to navigation or explicitly exclude it in
`not_in_nav` so the strict build can detect accidental omissions.

The theme is the one our projects share. Only the brand strings differ: the
logo text, the eyebrow, the footer, the sidebar note and the theme-preference
storage key.

## Publishing

`.github/workflows/docs.yml` checks the site on every push. It deploys to
GitHub Pages only from `main` of a public repository.
