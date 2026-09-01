/**
 * LaTeX rendering for the docs, via `react-katex`. Exported as `Tex`/`TexBlock`
 * rather than `Math` so nothing shadows the global `Math`.
 * `katex/dist/katex.min.css` is loaded once in `main.tsx`.
 */
export { InlineMath as Tex, BlockMath as TexBlock } from "react-katex";
