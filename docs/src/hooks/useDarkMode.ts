import { useEffect, useState } from "react";

const QUERY = "(prefers-color-scheme: dark)";

/** True when the OS/browser is in dark mode; updates when that changes. */
function useDarkMode(): boolean {
  const [dark, setDark] = useState(() => window.matchMedia(QUERY).matches);

  useEffect(() => {
    const mq = window.matchMedia(QUERY);
    const onChange = (e: MediaQueryListEvent) => setDark(e.matches);
    mq.addEventListener("change", onChange);
    return () => mq.removeEventListener("change", onChange);
  }, []);

  return dark;
}

export default useDarkMode;
