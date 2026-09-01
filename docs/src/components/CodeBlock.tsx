import { Highlight, themes } from "prism-react-renderer";
import useDarkMode from "../hooks/useDarkMode";
import styles from "./CodeBlock.module.css";

interface CodeBlockProps {
  code: string;
  language?: string;
}

/** A syntax-highlighted, read-only code sample. */
function CodeBlock({ code, language = "python" }: CodeBlockProps) {
  const dark = useDarkMode();

  return (
    <Highlight
      code={code.trim()}
      language={language}
      theme={dark ? themes.vsDark : themes.github}
    >
      {({ tokens, getLineProps, getTokenProps }) => (
        <pre className={styles.pre}>
          {tokens.map((line, i) => (
            <span key={i} {...getLineProps({ line, className: styles.line })}>
              {line.map((token, j) => (
                <span key={j} {...getTokenProps({ token })} />
              ))}
            </span>
          ))}
        </pre>
      )}
    </Highlight>
  );
}

export default CodeBlock;
