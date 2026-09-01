import { useMemo, useState } from "react";
import {
  dimSpecht,
  dimWeyl,
  partitions,
  semiStandardYoungTableau,
  standardYoungTableaux,
  type Tableau,
} from "@/schur-weyl";
import PartitionSelector from "../PartitionSelector";
import YoungTableau from "../YoungTableau";
import { Tex } from "../Math";
import Slider from "./Slider";
import styles from "./Playground.module.css";

const MAX_SHOWN = 60;

function TableauxPlayground() {
  const [n, setN] = useState(4);
  const [d, setD] = useState(3);
  const [index, setIndex] = useState(0);
  const [kind, setKind] = useState<"standard" | "semistandard">("standard");

  const parts = useMemo(() => [...partitions(n)], [n]);
  const lam = parts[Math.min(index, parts.length - 1)];

  const all: Tableau[] =
    kind === "standard" ? [...standardYoungTableaux(lam)] : [...semiStandardYoungTableau(lam, d)];

  return (
    <div className={styles.card}>
      <div className={styles.controls}>
        <Slider label="n" min={1} max={6} value={n} onChange={setN} />
        {kind === "semistandard" && (
          <Slider label="d" min={1} max={4} value={d} onChange={setD} />
        )}
        <div className={styles.tabs}>
          <button
            className={kind === "standard" ? `${styles.tab} ${styles.tabActive}` : styles.tab}
            onClick={() => setKind("standard")}
          >
            standard
          </button>
          <button
            className={kind === "semistandard" ? `${styles.tab} ${styles.tabActive}` : styles.tab}
            onClick={() => setKind("semistandard")}
          >
            semistandard
          </button>
        </div>
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
          {kind === "standard" ? (
            <Tex math={`f^\\lambda = ${dimSpecht(lam)}`} />
          ) : (
            <Tex math={`\\dim V_\\lambda^{${d}} = ${dimWeyl(lam, d)}`} />
          )}
          <div className={styles.hint}>
            {all.length} tableaux of shape <Tex math={`(${lam.join(",")})`} />
          </div>
        </div>
      </div>

      <div className={styles.gallery} style={{ marginTop: 24 }}>
        {all.slice(0, MAX_SHOWN).map((t, i) => (
          <YoungTableau key={i} tableau={t} size={26} />
        ))}
      </div>
      {all.length > MAX_SHOWN && (
        <div className={styles.hint}>… and {all.length - MAX_SHOWN} more.</div>
      )}
    </div>
  );
}

export default TableauxPlayground;
