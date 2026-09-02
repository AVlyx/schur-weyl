# docs

The documentation site for [`schur-weyl`](https://github.com/AVlyx/schur-weyl) — a single page
with interactive figures, built with Vite + React.

```bash
npm install
npm run dev      # http://localhost:5173/schur-weyl/
npm run build    # -> dist/
```

## Layout

| Path | What |
|---|---|
| `src/App.tsx` | The whole page: prose, code samples and the four playgrounds |
| `src/schur-weyl/` | A TypeScript port of `src/schur_weyl/`, so the figures compute for real in the browser |
| `src/components/` | `YoungTableau`, `YoungDiagram`, `PartitionSelector`, plus `Math` (KaTeX) and `CodeBlock` (Prism) |
| `src/components/playgrounds/` | One component per interactive figure |

`YoungTableau` is the primitive everything else draws with: give it a ragged array of numbers
(`0` renders as an empty box) and optionally a `size`, an `isHighlighted(i, j)` predicate and an
`onCellHover` callback.

## Keeping the port in sync

`src/schur-weyl/` mirrors the Python package module for module and is only used by this site. When
you change a function in `src/schur_weyl/`, mirror it here so the figures keep matching the docs.
