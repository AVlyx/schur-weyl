import { useState } from "react";
import "./App.css";
import PartitionSelector from "./components/PartitionSelector";
import { type Partition, partitions } from "./schur-weyl";

function App() {
  const partitionsArr = [...partitions(4)];
  const [selectedPartition, setSelectedPartition] = useState<Partition>(partitionsArr[0]);
  return (
    <PartitionSelector
      partitions={partitionsArr}
      selected={selectedPartition}
      onChange={setSelectedPartition}
    />
  );
}

export default App;
