import type { Tableau } from "@/schur-weyl";
import styles from "./YoungTableau.module.css";

interface YoungTableauProps {
  tableau: Tableau;
  /** Box side length in px. */
  size?: number;
  /** Cells for which this returns true are drawn in the accent colour. */
  isHighlighted?: (i: number, j: number) => boolean;
  /** Called with the hovered cell, or null when the pointer leaves. */
  onCellHover?: (cell: [number, number] | null) => void;
}

function YoungTableau({ tableau, size = 40, isHighlighted, onCellHover }: YoungTableauProps) {
  return (
    <div style={{ "--box": `${size}px` } as React.CSSProperties}>
      {tableau.map((row, rowIndex) => {
        return (
          <div key={`$row ${rowIndex}`} className={styles.outer}>
            {row.map((box, boxIndex) => {
              const highlighted = isHighlighted?.(rowIndex, boxIndex) ?? false;
              return (
                <div
                  key={`box${rowIndex}-${boxIndex}`}
                  className={highlighted ? `${styles.box} ${styles.highlighted}` : styles.box}
                  onMouseEnter={() => onCellHover?.([rowIndex, boxIndex])}
                  onMouseLeave={() => onCellHover?.(null)}
                >
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
