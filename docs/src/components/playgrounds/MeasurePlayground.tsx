import { useMemo, useState } from "react";
import { schurWeylMeasure } from "@/schur-weyl";
import YoungDiagram from "../YoungDiagram";
import { Tex } from "../Math";
import Slider from "./Slider";
import styles from "./Playground.module.css";

function MeasurePlayground() {
  const [p, setP] = useState(0.7);
  const [k, setK] = useState(5);

  const spectrum = useMemo(() => [p, 1 - p], [p]);
  const entries = useMemo(() => schurWeylMeasure(spectrum, k), [spectrum, k]);
  const total = entries.reduce((s, e) => s + e.probability, 0);

  return (
    <div className={styles.card}>
      <div className={styles.controls}>
        <Slider
          label="p"
          min={0.5}
          max={1}
          step={0.01}
          value={p}
          format={(v) => v.toFixed(2)}
          onChange={setP}
        />
        <Slider label="k" min={1} max={8} value={k} onChange={setK} />
        <span className={styles.hint}>
          <Tex math={`\\mathrm{spec}\\,\\rho = (${p.toFixed(2)},\\, ${(1 - p).toFixed(2)})`} />
        </span>
      </div>

      <div className={styles.bars}>
        {entries.map((e, i) => {
          // Partitions outside the support come out as tiny rounding noise.
          const prob = Math.max(0, e.probability);
          return (
            <div key={i} className={styles.bar}>
              <div className={styles.barShape}>
                <YoungDiagram partition={e.partition} size={12} />
              </div>
              <div className={styles.barTrack}>
                <div className={styles.barFill} style={{ width: `${(prob * 100).toFixed(2)}%` }} />
              </div>
              <div className={styles.barValue}>{prob.toFixed(4)}</div>
            </div>
          );
        })}
      </div>

      <div className={styles.hint}>
        Only partitions with at most 2 rows get any mass (the state is a qubit), and the total is{" "}
        {total.toFixed(6)}. As <Tex math="k" /> grows, <Tex math="\lambda / k" /> concentrates on
        the spectrum <Tex math={`(${p.toFixed(2)}, ${(1 - p).toFixed(2)})`} />.
      </div>
    </div>
  );
}

export default MeasurePlayground;
