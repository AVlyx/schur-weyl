import { useMemo, useState } from "react";
import { dimSpecht, dimWeyl, hookLength, partitions, type Tableau } from "@/schur-weyl";
import PartitionSelector from "../PartitionSelector";
import YoungTableau from "../YoungTableau";
import { Tex } from "../Math";
import Slider from "./Slider";
import styles from "./Playground.module.css";

/** The hook of (a, b): the cell itself, its arm to the right and its leg below. */
function inHook(a: number, b: number, i: number, j: number): boolean {
  return (i === a && j >= b) || (j === b && i >= a);
}

function DiagramPlayground() {
  const [n, setN] = useState(5);
  const [d, setD] = useState(3);
  const [index, setIndex] = useState(0);
  const [hover, setHover] = useState<[number, number] | null>(null);

  const parts = useMemo(() => [...partitions(n)], [n]);
  const lam = parts[Math.min(index, parts.length - 1)];

  const hooks: Tableau = lam.map((li, i) =>
    Array.from({ length: li }, (_, j) => hookLength(lam, i, j)),
  );
  const product = hooks.flat().reduce((a, b) => a * b, 1);
  const sumOfSquares = parts.reduce((s, p) => s + dimSpecht(p) ** 2, 0);

  return (
    <div className={styles.card}>
      <div className={styles.controls}>
        <Slider label="n" min={1} max={8} value={n} onChange={setN} />
        <Slider label="d" min={1} max={6} value={d} onChange={setD} />
      </div>

      <div className={styles.body}>
        <div style={{ flex: "1 1 300px" }}>
          <PartitionSelector
            partitions={parts}
            selected={lam}
            onChange={(p) => setIndex(parts.indexOf(p))}
          />
        </div>

        <div className={styles.panel}>
          <YoungTableau
            tableau={hooks}
            isHighlighted={hover ? (i, j) => inHook(hover[0], hover[1], i, j) : undefined}
            onCellHover={setHover}
          />
          <div className={styles.hint}>hook lengths — hover a cell to light up its hook</div>
        </div>

        <div className={styles.panel}>
          <div>
            <Tex math={`\\lambda = (${lam.join(",")}) \\vdash ${n}`} />
          </div>
          <div>
            <Tex math={`f^\\lambda = \\frac{${n}!}{${product}} = ${dimSpecht(lam)}`} />
          </div>
          <div>
            <Tex math={`\\dim V_\\lambda^{${d}} = ${dimWeyl(lam, d)}`} />
          </div>
          <div className={styles.hint}>
            {parts.length} partitions of {n}, and{" "}
            <Tex math={`\\textstyle\\sum_\\lambda (f^\\lambda)^2 = ${sumOfSquares} = ${n}!`} />
          </div>
        </div>
      </div>
    </div>
  );
}

export default DiagramPlayground;
