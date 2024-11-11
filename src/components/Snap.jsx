import React from "react";
import { useEffect, useState } from "react";
import { toast } from "react-toastify";
import api, { baseUrl } from "../api";
import { ErrorCard } from "./Error";

function DomainChunk({ domainChunk }) {
  let struct = null;
  if (domainChunk[0].endsWith("CS"))
    [struct, ...domainChunk] = domainChunk;
  let [targetLine, align, queryLine, ...residues] = domainChunk;
  let pp = residues.pop();

  const leadingWhitespaceSize = pp.length - pp.trimStart().length;
  const lastAminoAcidIdx = targetLine.lastIndexOf(" ") - 1;

  function align_class(x) {
    return x === '+' ? "hmmplus" : x === ' ' ? "hmmminus" : "hmmmatch"
  }

  function heat_class(x) {
    return x === '*' ? "heatstar" : x === ' ' ? 'heatgap' : `head${x}`;
  }

  let target = [<span key="target-prefix">{targetLine.slice(0, leadingWhitespaceSize)}</span>];
  {
    let i = leadingWhitespaceSize;
    while (i <= lastAminoAcidIdx) {
      let sequence = "";
      let j = i;
      while (j <= lastAminoAcidIdx && align_class(align[i]) === align_class(align[j])) {
        sequence += targetLine[j];
        j++;
      }
      target.push(<span key={i} className={align_class(align[i])}>{sequence}</span>);
      i = j;
    }
  }
  target.push(<span key="target-suffix">{targetLine.slice(lastAminoAcidIdx + 1, targetLine.length)}</span>);

  let query = [<span key="query-prefix">{queryLine.slice(0, leadingWhitespaceSize)}</span>];
  {
    let i = leadingWhitespaceSize;
    while (i <= lastAminoAcidIdx) {
      let sequence = "";
      let j = i;
      while (j <= lastAminoAcidIdx && heat_class(pp[i]) === heat_class(pp[j])) {
        sequence += queryLine[j];
        j++;
      }
      query.push(<span key={i} className={heat_class(pp[i])}>{sequence}</span>);
      i = j;
    }
  }
  query.push(<span key="target-suffix">{queryLine.slice(lastAminoAcidIdx + 1, queryLine.length)}</span>);

  return (
    <div>
      {struct && <p key="struct" className="alignment-line">{struct}</p>}
      <p key="target" className="alignment-line">{target}</p>
      <p key="align" className="alignment-line">{align}</p>
      <p key="query" className="alignment-line">{query}</p>
      <p key="res0" className="alignment-line">{residues[0]}</p>
      <p key="res1" className="alignment-line">{residues[1]}</p>
      <p key="res2" className="alignment-line">{residues[2]}</p>
      <p key="res3" className="alignment-line">{residues[3]}</p>
      <p key="res4" className="alignment-line">{residues[4]}</p>
      <p key="pp" className="alignment-line">{pp}</p>
    </div>
  );
}

function partitionAt(arr, fn) {
  let chunks = [];
  let i = 0
  while (i < arr.length) {
    let chunk = [];
    let j = i;
    while (j < arr.length && !fn(arr[j])) {
      chunk.push(arr[j]);
      j++;
    }
    if (j < arr.length) chunk.push(arr[j]);
    chunks.push(chunk);
    i = j + 1;
  }
  return chunks;
}

function Domain({ domain }) {
  let chunks = partitionAt(domain.filter(x => x.length > 0), x => x.endsWith("PP"))
  let lines = []
  chunks.forEach((domainChunk, i) => {
    lines.push(<DomainChunk key={`chunk-${i}`} domainChunk={domainChunk}></DomainChunk>);
    lines.push(<React.Fragment key={`break-${i}`}><br /><br /></React.Fragment>);
  })
  return lines;
}

function Alignment({ alignment }) {
  let [headerLine, statLine, ...domain] = alignment.split("\n");
  return (
    <div>
      <p key="header" className={"alignment-line"}>{headerLine || "Alignments for each domain:"}</p>
      <p key="statistics" className={"alignment-line"}>{statLine}</p>
      <Domain domain={domain}></Domain>
    </div>
  );
}

function SnapView({ view }) {
  const alignments = view.split("\n\nAlignments for each domain:");
  return (
    <div className="alignment-output">
      {alignments.map((x, i) => <Alignment key={i} alignment={x}></Alignment>)}
    </div>
  );
};

function SnapLoading() {
  return (
    <div className="alignment-loading">
      <div>                           </div><br />
      <div>                                                             </div><br />
      <div>                                                                                                                  </div><br />
      <div>                                                                                                               </div><br />
      <div>                                                                                                               </div><br />
      <div>                                                                                                               </div><br />
      <div>                                                                                                                   </div><br />
      <div>                                                                                                               </div><br />
      <div>                                                                                                               </div><br />
      <div>                                                                                                               </div><br />
      <div>                                                                                                               </div><br />
      <div>                                                                                                                  </div><br />
      <div></div><br />
      <div></div><br />
      <div>                                                                                               </div><br />
      <div>                                                                                            </div><br />
      <div>                                                                                            </div><br />
      <div>                                                                                            </div><br />
      <div>                                                                                                </div><br />
      <div>                                                                                            </div><br />
      <div>                                                                                            </div><br />
      <div>                                                                                            </div><br />
      <div>                                                                                            </div><br />
      <div>                                                                                               </div>
    </div>);
}

export default function Snap({ scanId }) {
  const [view, setView] = useState(null);
  const [error, setError] = useState();

  useEffect(() => {
    api
      .get(`/scans/${scanId}/snap.dcs/view`)
      .then(x => setView(x.data))
      .catch(error => setError(error))
  }, [scanId]);

  return (
    <div>
      <article className="vf-card vf-card--brand vf-card--bordered">
        <div className="vf-card__content | vf-stack vf-stack--400">
          <div className="vf-sidebar vf-sidebar--end">
            <div className="vf-sidebar__inner">
              <div>
                <h3 className="vf-card__heading">
                  Alignments
                </h3>
                <p className="vf-card__text vf-text-body--5">See <a href="https://hmmer-web-docs.readthedocs.io/en/latest/result.html#alignments" className="vf-link" target="_newtab">HMMER documentation</a> for format.</p>
              </div>
              <div>
                <SnapCopier
                  title={"Alignments"}
                  name="view"
                  scanId={scanId}
                />
              </div>
            </div>
          </div>
          <hr className="vf-divider" />
          {error && (<ErrorCard error={error}></ErrorCard>)}
          {view == null && (<SnapLoading></SnapLoading>)}
          {view != null && (<SnapView view={view}></SnapView>)}
        </div>
      </article>
      <section className="vf-card-container vf-card-container__col-3 | vf-u-background-color--grey--lightest vf-u-fullbleed">
        <div className="vf-card-container__inner">
          <div className="vf-section-header">
            <h2 className="vf-section-header__heading">Results</h2>
            <p className="vf-section-header__text">
              Results files from your search
            </p>
          </div>

          <SnapCard
            title="Alignments"
            description="Alignment of all matches."
            name="view"
            fileFormat="TXT"
            scanId={scanId}
            view={view}
          />

          <SnapCard
            title="GFF"
            description="GFF (General Feature Format) v3 file listing all found matches."
            name="gff"
            fileFormat="GFF"
            scanId={scanId}
          />

          <SnapCard
            title="Original Query"
            description="FA (FASTA) file of matched query subsequences."
            name="queries"
            fileFormat="FASTA"
            scanId={scanId}
          />

          <SnapCard
            title={"Protein sequence matches"}
            description="FAA (FASTA Amino Acids) file of matched amino acid sequences."
            name="aminos"
            fileFormat="FASTA"
            scanId={scanId}
          />

          <SnapCard
            title="DNA of protein sequences"
            description="FA (FASTA) file of codons."
            name="codons"
            fileFormat="FASTA"
            scanId={scanId}
          />

          <SnapCard
            title="HMM Path"
            description="Text file, showing match/insertion/deletion states of matches."
            name="states"
            fileFormat="TXT"
            scanId={scanId}
          />
        </div>
      </section>
    </div>
  );
}

const SnapCopier = ({ name, title, scanId }) => {
  const toastOptions = {
    position: "bottom-left",
    autoClose: 5000,
    hideProgressBar: false,
    closeOnClick: true,
    pauseOnHover: true,
    draggable: true,
    progress: undefined,
  };
  const [data, setData] = useState(null);
  const url = `/scans/${scanId}/snap.dcs/${name}`;

  useEffect(() => {
    api
      .get(url)
      .then((response) => setData(response.data))
  });

  async function onClick() {
    navigator.clipboard.writeText(data)
      .then(() => toast.success(`👍 Copied ${title} to clipboard!`, toastOptions))
  }

  return (
    <div>
      <button
        className={"vf-button vf-button--secondary vf-button--sm"}
        disabled={data === null}
        onClick={onClick}
      >
        <i className="icon icon-common icon-copy" />
        &nbsp;Copy
      </button>
      <button
        className={"vf-button vf-button--secondary vf-button--sm"}
        onClick={() => window.open(`${baseUrl}${url}`, "_blank")}
      >
        <i className="icon icon-common icon-download" />
        &nbsp;View
      </button>
    </div >
  );
};

const SnapCard = ({
  name,
  title,
  description,
  fileFormat,
  scanId,
}) => {
  return (
    <article
      key={name}
      id={`snap-card-${name}`}
      className="vf-card vf-card--brand vf-card--bordered"
    >
      <div className="vf-card__content | vf-stack vf-stack--400">
        <h3 className="vf-card__heading">
          <i className={`icon icon-fileformats icon-${fileFormat}`} />
          &nbsp; {title}
        </h3>
        <p className="vf-card__text">{description}</p>
        <SnapCopier
          title={title}
          name={name}
          scanId={scanId}
        />
      </div>
    </article>
  );
};
