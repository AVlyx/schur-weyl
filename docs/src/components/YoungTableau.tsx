import type { Tableau } from "@/schur-weyl";
import styles from "./YoungTableau.module.css";

interface YoungTableauProps {
  tableau: Tableau;
}

function YoungTableau({ tableau }: YoungTableauProps) {
  return (
    <div>
      {tableau.map((row, rowIndex) => {
        return (
          <div key={`$row ${rowIndex}`} className={styles.outer}>
            {row.map((box, boxIndex) => {
              return (
                <div key={`box${rowIndex}-${boxIndex}`} className={styles.box}>
                  {box === 0 ? "" : box}
                </div>
              );
            })}
          </div>
        );
      })}
    </div>
  );
}

export default YoungTableau;
