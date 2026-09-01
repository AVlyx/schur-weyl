import styles from "./Playground.module.css";

interface SliderProps {
  label: string;
  min: number;
  max: number;
  step?: number;
  value: number;
  format?: (value: number) => string;
  onChange: (value: number) => void;
}

/** A labelled range input with the current value shown next to it. */
function Slider({ label, min, max, step = 1, value, format, onChange }: SliderProps) {
  return (
    <label className={styles.label}>
      {label}
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={value}
        onChange={(e) => onChange(Number(e.target.value))}
      />
      <span className={styles.value}>{format ? format(value) : value}</span>
    </label>
  );
}

export default Slider;
