import { useMemo, useState } from "react";
import { characterTable } from "@/schur-weyl";
import YoungDiagram from "../YoungDiagram";
import { Tex } from "../Math";
import Slider from "./Slider";
import styles from "./Playground.module.css";

function CharacterPlayground() {
  const [n, setN] = useState(5);
  const { rows, cols, table } = useMemo(() => characterTable(n), [n]);

  return (
    <div className={styles.card}>
      <div className={styles.controls}>
        <Slider label="n" min={1} max={7} value={n} onChange={setN} />
        <span className={styles.hint}>
          rows are <Tex math="\lambda" />, columns are the cycle type <Tex math="\mu" />
        </span>
      </div>

      <div className={styles.scroll}>
        <table className={styles.table}>
          <thead>
            <tr>
              <th>
                <Tex math="\chi^\lambda(\mu)" />
              </th>
              {cols.map((mu, j) => (
                <th key={j}>
                  <YoungDiagram partition={mu} size={9} />
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {rows.map((lam, i) => (
              <tr key={i}>
                <th>
                  <YoungDiagram partition={lam} size={9} />
                </th>
                {table[i].map((value, j) => (
                  <td key={j} style={value === 0 ? { color: "var(--text)" } : undefined}>
                    {value}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <div className={styles.hint}>
        The first column is <Tex math="\chi^\lambda(1^n) = f^\lambda" />, the dimension of the
        Specht module.
      </div>
    </div>
  );
}

export default CharacterPlayground;
