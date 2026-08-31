import type { Partition, Tableau } from "@/schur-weyl";
import YoungTableau from "./YoungTableau";

interface YoungDiagramProps {
  partition: Partition;
}

function YoungDiagram({ partition }: YoungDiagramProps) {
  const convert: Tableau = partition.map((li) => Array(li).fill(0));
  return <YoungTableau tableau={convert} />;
}

export default YoungDiagram;
