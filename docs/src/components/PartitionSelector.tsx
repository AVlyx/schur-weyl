import type { Partition } from "@/schur-weyl";
import styles from "./PartitionSelector.module.css";
import YoungDiagram from "./YoungDiagram";

interface PartitionSelectorProps {
  partitions: Partition[];
  selected: Partition; //ref needs to be in partitions array
  onChange: (partition: Partition) => void;
}

function PartitionSelector({ partitions, selected, onChange }: PartitionSelectorProps) {
  const updatePartition = (increment: number) => {
    const selectedIndex = partitions.indexOf(selected);
    const newSelectedIndex = (selectedIndex + increment + partitions.length) % partitions.length;
    onChange(partitions[newSelectedIndex]);
  };

  return (
    <div className={styles.outer}>
      {partitions.length === 0 ? (
        <></>
      ) : (
        <button className={styles.arrow} onClick={() => updatePartition(-1)}>
          {"←"}
        </button>
      )}
      <YoungDiagram partition={selected} />
      {partitions.length === 0 ? (
        <></>
      ) : (
        <button className={styles.arrow} onClick={() => updatePartition(1)}>
          {"→"}
        </button>
      )}
    </div>
  );
}

export default PartitionSelector;
