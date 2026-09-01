import "./App.css";
import CodeBlock from "./components/CodeBlock";
import { Tex, TexBlock } from "./components/Math";
import CharacterPlayground from "./components/playgrounds/CharacterPlayground";
import DiagramPlayground from "./components/playgrounds/DiagramPlayground";
import MeasurePlayground from "./components/playgrounds/MeasurePlayground";
import TableauxPlayground from "./components/playgrounds/TableauxPlayground";

const REPO = "https://github.com/AVlyx/schur-weyl";

const SECTIONS = [
  ["install", "Install"],
  ["start", "Quick start"],
  ["diagrams", "Young diagrams"],
  ["tableaux", "Tableaux"],
  ["characters", "Characters"],
  ["measure", "Schur–Weyl measure"],
  ["reference", "Reference"],
] as const;

const MODULES: [string, string][] = [
  ["young_diagrams", "Partitions, conjugation, hook lengths, dominance order, corners"],
  ["dimensions", "dim_specht (hook-length formula) and dim_weyl (content formula)"],
  ["tableaux", "Standard and semistandard tableaux, Kostka numbers, reading words"],
  ["symmetric_group", "Permutations as one-line tuples: compose, invert, cycle type"],
  ["character", "Irreducible characters of S_n via Murnaghan–Nakayama; character tables"],
  ["symmetric_functions", "Schur polynomials (Jacobi–Trudi) and power sums"],
  ["isotypic", "Isotypic projectors on tensor space"],
  ["sw_measure", "The Schur–Weyl measure for i.i.d. copies of a state"],
  ["young_orthonormal", "Youngs orthogonal form and the Jucys–Murphy elements"],
];

function App() {
  return (
    <div className="page">
      <header className="hero">
        <h1>schur-weyl</h1>
        <p className="lead">
          Partitions, Young tableaux, symmetric-group characters and Schur–Weyl duality, built from
          scratch on <code>numpy</code> and the standard library.
        </p>
        <p className="links">
          <a href={REPO}>GitHub</a>
          <a href="https://pypi.org/project/schur-weyl/">PyPI</a>
        </p>
      </header>

      <nav className="nav">
        {SECTIONS.map(([id, title]) => (
          <a key={id} href={`#${id}`}>
            {title}
          </a>
        ))}
      </nav>

      <section>
        <p>
          Schur–Weyl duality splits <Tex math="(\mathbb{C}^d)^{\otimes n}" /> under the commuting
          actions of <Tex math="S_n" /> (permuting the tensor factors) and <Tex math="GL_d" />{" "}
          (acting diagonally) into blocks indexed by partitions <Tex math="\lambda \vdash n" /> with{" "}
          <Tex math="\ell(\lambda) \le d" />:
        </p>
        <TexBlock math="(\mathbb{C}^d)^{\otimes n} \;\cong\; \bigoplus_{\lambda} S^\lambda \otimes V_\lambda^d" />
        <p>
          This package gives you both sides of that sum — the combinatorics of{" "}
          <Tex math="\lambda" />, and the linear algebra that projects onto each block.
        </p>
      </section>

      <section id="install">
        <h2>Install</h2>
        <CodeBlock language="bash" code="pip install schur-weyl" />
        <p className="note">
          Requires Python ≥ 3.9. The only runtime dependency is <code>numpy</code>.
        </p>
      </section>

      <section id="start">
        <h2>Quick start</h2>
        <p>
          Partitions are plain tuples and tableaux are tuples of tuples. Nothing is a custom class,
          so everything composes with ordinary Python.
        </p>
        <CodeBlock
          code={`
from schur_weyl import partitions, dim_specht, dim_weyl, isotypic_proj
from schur_weyl.character import character
from schur_weyl.sw_measure import schur_weyl_measure

[*partitions(4)]
# [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)]

dim_specht((2, 1))            # 2   -- the standard rep of S_3
dim_weyl((2, 1), d=3)         # 8

character((2, 1), (3,))       # -1  -- Murnaghan-Nakayama

schur_weyl_measure([0.7, 0.3], k=4)
# {(4,): 0.4141, (3, 1): 0.4977, (2, 2): 0.0882, ...}

P = isotypic_proj((2, 1), d=2)   # np.ndarray, P @ P == P
`}
        />
      </section>

      <section id="diagrams">
        <h2>Young diagrams</h2>
        <p>
          A partition <Tex math="\lambda = (\lambda_1 \ge \lambda_2 \ge \cdots)" /> is drawn as rows
          of boxes. The <em>hook</em> of a box is the box itself plus everything to its right and
          everything below it; the hook-length formula turns those numbers into the dimension of the
          Specht module.
        </p>
        <TexBlock math="f^\lambda = \frac{n!}{\prod_{(i,j) \in \lambda} h(i,j)}" />
        <DiagramPlayground />
        <CodeBlock
          code={`
from schur_weyl import partitions, partition_conjugate, hook_lenght, dim_specht

[*partitions(5, max_part=2)]        # [(2, 2, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)]
partition_conjugate((4, 3, 1))      # (3, 2, 2, 1)   -- transpose the diagram
hook_lenght((4, 3, 1), 0, 0)        # 6
dim_specht((3, 2, 1))               # 16
`}
        />
      </section>

      <section id="tableaux">
        <h2>Tableaux</h2>
        <p>
          Fill the diagram with numbers. A <strong>standard</strong> tableau uses{" "}
          <Tex math="1, \dots, n" /> once each, increasing along rows and down columns — there are{" "}
          <Tex math="f^\lambda" /> of them. A <strong>semistandard</strong> tableau draws entries
          from <Tex math="1, \dots, d" /> with repeats allowed, weakly increasing along rows and
          strictly increasing down columns — there are <Tex math="\dim V_\lambda^d" /> of them.
        </p>
        <TableauxPlayground />
        <CodeBlock
          code={`
from schur_weyl import standard_young_tableaux, semi_standard_young_tableau, kostka

[*standard_young_tableaux((2, 1))]
# [((1, 2), (3,)), ((1, 3), (2,))]

[*semi_standard_young_tableau((2, 1), 2)]
# [((1, 1), (2,)), ((1, 2), (2,))]

kostka((2, 1), (1, 1, 1))           # 2
`}
        />
      </section>

      <section id="characters">
        <h2>Characters</h2>
        <p>
          <Tex math="\chi^\lambda(\mu)" /> is the character of the irreducible representation{" "}
          <Tex math="S^\lambda" /> on the conjugacy class of cycle type <Tex math="\mu" />, computed
          by peeling border strips off <Tex math="\lambda" /> (Murnaghan–Nakayama). Distinct rows of
          the table are orthogonal.
        </p>
        <CharacterPlayground />
        <CodeBlock
          code={`
from schur_weyl.character import character, character_table

character((2, 1), (3,))             # -1
character((3, 1), (2, 1, 1))        #  1

rows, cols, table = character_table(5)
`}
        />
      </section>

      <section id="measure">
        <h2>Schur–Weyl measure</h2>
        <p>
          Measure <Tex math="k" /> i.i.d. copies of a state <Tex math="\rho" /> in the Schur basis
          and you get a random partition, with
        </p>
        <TexBlock math="\Pr[\lambda] \;=\; f^\lambda \, s_\lambda(\mathrm{spec}\,\rho)" />
        <p>
          This is the distribution behind spectrum estimation and quantum tomography. Below,{" "}
          <Tex math="\rho" /> is a qubit with eigenvalues <Tex math="(p,\, 1-p)" />.
        </p>
        <MeasurePlayground />
        <CodeBlock
          code={`
from schur_weyl.sw_measure import schur_weyl_measure

schur_weyl_measure(spectrum=[0.7, 0.3], k=4)
# {(4,): 0.4141, (3, 1): 0.4977, (2, 2): 0.0882,
#  (2, 1, 1): 0.0, (1, 1, 1, 1): 0.0}
`}
        />
      </section>

      <section id="reference">
        <h2>Reference</h2>
        <table className="modules">
          <tbody>
            {MODULES.map(([name, what]) => (
              <tr key={name}>
                <td>
                  <a href={`${REPO}/blob/main/src/schur_weyl/${name}.py`}>
                    <code>{name}</code>
                  </a>
                </td>
                <td>{what}</td>
              </tr>
            ))}
          </tbody>
        </table>
        <p className="note">
          Every function carries a docstring with worked examples, and the test suite in{" "}
          <a href={`${REPO}/tree/main/tests`}>
            <code>tests/</code>
          </a>{" "}
          checks closed-form identities (<Tex math="\sum_\lambda (f^\lambda)^2 = n!" />, character
          orthogonality, <Tex math="P_\lambda^2 = P_\lambda" />) rather than fixed outputs.
        </p>
      </section>

      <footer>
        MIT licensed. The interactive figures on this page run a TypeScript port of the library in
        your browser.
      </footer>
    </div>
  );
}

export default App;
