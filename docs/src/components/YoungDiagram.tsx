import type { Partition, Tableau } from "@/schur-weyl";
import YoungTableau from "./YoungTableau";

interface YoungDiagramProps {
  partition: Partition;
  /** Box side length in px. */
  size?: number;
  isHighlighted?: (i: number, j: number) => boolean;
  onCellHover?: (cell: [number, number] | null) => void;
}

function YoungDiagram({ partition, size, isHighlighted, onCellHover }: YoungDiagramProps) {
  const convert: Tableau = partition.map((li) => Array(li).fill(0));
  return (
    <YoungTableau
      tableau={convert}
      size={size}
      isHighlighted={isHighlighted}
      onCellHover={onCellHover}
    />
  );
}

export default YoungDiagram;
