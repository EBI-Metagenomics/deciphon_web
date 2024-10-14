import React from "react";
import { useParams } from "react-router-dom";
import {useEffect, useMemo, useState} from "react";
import api, { baseUrl, pollingInterval } from "../api";
import Loading from "./Loading";
import { toast } from "react-toastify";
import { useInterval } from "react-use";
import { chunk, find } from "lodash";

const UrlCopier = () => {
  const sayCopied = () =>
    toast.success("👍 Copied URL to clipboard!", {
      position: "bottom-left",
      autoClose: 5000,
      hideProgressBar: false,
      closeOnClick: true,
      pauseOnHover: true,
      draggable: true,
      progress: undefined,
    });
  return (
    <div>
      <button
        className={"vf-button vf-button--link vf-button--sm"}
        style={{ fontFamily: "monospace" }}
        onClick={() => {
          navigator.clipboard.writeText(window.location.href);
          sayCopied();
        }}
      >
        {window.location.href}&nbsp;
        <i className="icon icon-common icon-copy" />
      </button>
    </div>
  );
};
const toastOptions = {
  position: "bottom-left",
  autoClose: 5000,
  hideProgressBar: false,
  closeOnClick: true,
  pauseOnHover: true,
  draggable: true,
  progress: undefined,
};

const ResultCopier = ({ resultSuffix, title, scanId }) => {
  return (
    <div>
      <button
        className={"vf-button vf-button--secondary vf-button--sm"}
        onClick={async () => {
          toast.promise(
            () =>
              api
                .get(`/scans/${scanId}/${resultSuffix}`)
                .then((response) => {
                  navigator.clipboard.writeText(response.data);
                }),
            {
              pending: "⏳ Fetching result...",
              success: `👍 Copied ${title} to clipboard!`,
              error: "🚨 Sorry something went wrong.",
            },
            toastOptions
          );
        }}
      >
        <i className="icon icon-common icon-copy" />
        &nbsp;Copy result
      </button>
      <button
        className={"vf-button vf-button--secondary vf-button--sm"}
        onClick={() => {
          window.open(
            `${baseUrl}/scans/${scanId}/${resultSuffix}`,
            "_blank"
          );
        }}
      >
        <i className="icon icon-common icon-download" />
        &nbsp;View
      </button>
    </div>
  );
};

const ResultCard = ({
  resultSuffix,
  title,
  description,
  fileFormat,
  scanId,
}) => {
  return (
    <article
      key={resultSuffix}
      className="vf-card vf-card--brand vf-card--bordered"
    >
      <div className="vf-card__content | vf-stack vf-stack--400">
        <h3 className="vf-card__heading">
          <i className={`icon icon-fileformats icon-${fileFormat}`} />
          &nbsp; {title}
        </h3>
        <p className="vf-card__text">{description}</p>
        <ResultCopier
          title={title}
          resultSuffix={resultSuffix}
          scanId={scanId}
        />
      </div>
    </article>
  );
};

const Result = () => {
  let { jobid } = useParams();
  const [jobState, setJobState] = useState();
  const [errors, setErrors] = useState();
  const [isPolling, setIsPolling] = useState(false);
  const [jobsAhead, setJobsAhead] = useState(null);
  const [scanId, setScanId] = useState(null);
  const [numResults, setNumResults] = useState(null);
  const [alignmentContent, setAlignmentContent] = useState(null);

  const formattedAlignments = useMemo(() => {
    if (!alignmentContent) return null;
    const matches = alignmentContent.split("\n\nAlignments for each domain:");
    let outputLines = [];

    matches.forEach((matchLines, m) => {
        let [headerLine, statLine, ...alignmentsLines] = matchLines.split("\n");
        outputLines.push(<p key={`matches-${m}-header`} className={"alignment-line"}>{headerLine || "Alignments for each domain:"}</p>)
        outputLines.push(<p key={`matches-${m}-stat`} className={"alignment-line"}>{statLine}</p>);
        chunk(alignmentsLines, 11).forEach((alignmentChunk, a) => {
            const [structureLine, targetLine, alignmentLine, queryLine, ...otherLines] = alignmentChunk;
            const leadingWhitespaceSize = structureLine.length - structureLine.trimStart().length;
            const lastAminoAcidIdx = structureLine.lastIndexOf(" ") - 1;

            let formattedTargetLine = [<span key="prefix">{targetLine.slice(0, leadingWhitespaceSize)}</span>];
            let formattedQueryLine = [<span key="prefix">{queryLine.slice(0, leadingWhitespaceSize)}</span>];
            otherLines.pop();
            let ppLine = otherLines.pop();

            outputLines.push(<p key={`matches-${m}-alignmemt-${a}-structure`} className="alignment-line">{structureLine}</p>);

            for (let idx=leadingWhitespaceSize; idx <= lastAminoAcidIdx; idx ++) {
                const alignmentChar = alignmentLine[idx];
                const ppChar = ppLine[idx];
                let alignmentClass = alignmentChar === '+' ? "hmmplus" : alignmentChar === ' ' ? "hmmminus" : "hmmmatch";
                let heatClass = ppChar === '*' ? "heatstar" : ppChar === ' ' ? 'heatgap' : `head${ppChar}`;
                formattedTargetLine.push(<span key={idx} className={alignmentClass}>{targetLine[idx]}</span>);
                formattedQueryLine.push(<span key={idx} className={heatClass}>{queryLine[idx]}</span> );
            }
            outputLines.push(<p key={`matches-${m}-alignmemt-${a}-target`} className="alignment-line">{ formattedTargetLine }</p>);
            outputLines.push(<p key={`matches-${m}-alignmemt-${a}-alignment`} className="alignment-line">{ alignmentLine }</p>);
            outputLines.push(<p key={`matches-${m}-alignmemt-${a}-query`} className="alignment-line">{ formattedQueryLine }</p>);
            outputLines.push(...otherLines.join("\n"));
            outputLines.push(<p key={`matches-${m}-alignmemt-${a}-pp`} className="alignment-line">{ppLine}</p>);
            outputLines.push(<React.Fragment key={`matches-${m}-alignmemt-${a}-brs`}><br/><br/></React.Fragment>);
        });
        outputLines.push(<hr key={`matches-${m}-hr`} className="vf-divider"/>);
    })

    return outputLines;
  }, [alignmentContent]);


  useEffect(() => {
    api.get(`/scans?job_id=${jobid}`).then((response) => { 
      if (response?.data.length) {
        return response.data[0].id;
      }
    }).then((scan_id) => {
        setScanId(scan_id);
    });
  }, [jobid]);

  useInterval(
    () => {
      if (!jobid) return;
      api
        .get(`/jobs/${jobid}`)
        .then((response) => {
          setJobState(response.data);
          if (response.data?.error?.length) {
            setErrors([response.data.error]);
          }
          if (
            response?.data?.state === "done" ||
            response?.data?.state === "fail"
          ) {
            setIsPolling(false);
          } else {
            api
              .get(`/jobs?limit=30`)
              .then((response) => {
                if (response?.data?.length) {
                  const nextPendJob = find(response.data, job => job.state === 'pend' || job.state === 'run');
                  if (nextPendJob != undefined)
                    setJobsAhead(parseInt(jobid) - parseInt(nextPendJob.id));
                  else
                    setJobsAhead(null);
                } else {
                  setJobsAhead(null);
                }
              })
              .catch((err) => {
                console.error(err);
                setJobsAhead(null);
              });
          }
        })
        .catch((err) => setErrors([err?.response?.status]));
    },
    isPolling ? pollingInterval : null
  );

  useEffect(() => {
    if (!jobid) return;
    setIsPolling(true);
  }, [jobid]);

  useEffect(() => {
    if (jobState?.state === "done") {
      api.get(`/scans/${scanId}/snap.dcs/prods`).then((response) => {
        setNumResults(response?.data?.length);
      });
      api.get(`/scans/${scanId}/snap.dcs/view`).then((response) => {
          setAlignmentContent(response?.data)
      })
    }
  }, [jobState, scanId]);

  const finishedAt = jobState
    ? new Date(jobState?.exec_ended).toLocaleString()
    : null;

  return (
    <div className={"vf-stack vf-stack--400"}>
      <nav className="vf-navigation vf-navigation--main | vf-cluster">
        <ul className="vf-navigation__list | vf-list | vf-cluster__inner">
          <li className="vf-navigation__item">
            <a href="/" className="vf-navigation__link" aria-current="page">Query</a>
          </li>
          <li className="vf-navigation__item">
            <a href="/about" className="vf-navigation__link">About</a>
          </li>
        </ul>
      </nav>
      <nav className="vf-breadcrumbs" aria-label="Breadcrumb">
        <ul className="vf-breadcrumbs__list | vf-list vf-list--inline">
          <li className="vf-breadcrumbs__item">
            <a href="/" className="vf-breadcrumbs__link">
              Query
            </a>
          </li>
          <li className="vf-breadcrumbs__item" aria-current="location">
            Job {jobid}
          </li>
        </ul>
      </nav>
      <h1>Results</h1>

      {jobState?.state === "pend" && (
        <>
          <h3>Job is pending</h3>
          <span className="vf-badge vf-badge--primary">live updating</span>
          <UrlCopier />
          {jobsAhead !== null && (
            <p>There are {jobsAhead} jobs ahead of yours in the queue.</p>
          )}
        </>
      )}
      {jobState?.state === "run" && (
        <>
          <h3>Job is running</h3>
          <span className="vf-badge vf-badge--primary">live updating</span>
          <div
            role="progressbar"
            aria-valuenow={jobState.progress}
            aria-valuemin="0"
            aria-valuemax="100"
            className="vf-progress-indicator"
          >
            <div
              className="vf-progress-indicator__mark"
              style={{
                "--vf-progress-indicator__percent": `${jobState.progress}%`,
              }}
            />
            <p className="vf-progress-indicator__helper-text">
              Progress reported by Deciphon Job Server
            </p>
          </div>
          <UrlCopier />
        </>
      )}
      {jobState?.state === "fail" && (
        <>
          <h3>Job has failed</h3>
          <UrlCopier />
        </>
      )}
      <Loading
        isLoading={jobState?.state === "run" || jobState?.state === "pend"}
      />
      {errors && (
        <div className="vf-grid vf-grid__col-2">
          {errors.map((err, idx) => (
            <article
              key={idx}
              className="vf-card vf-card--brand vf-card--bordered"
              style={{ "--vf-card-border-color": "#d41645" }}
            >
              <div className="vf-card__content | vf-stack vf-stack--400">
                <h3 className="vf-card__heading">Error :(</h3>
                <p className="vf-card__subheading">{err}</p>
              </div>
            </article>
          ))}
        </div>
      )}
      {jobState?.state === "done" && (
        <div className="vf-stack vf-stack--400">
          <h4>Job complete &mdash; {numResults || 0} matches found</h4>
          <div>
            <span className="vf-form__helper">Finished at: {finishedAt}</span>
          </div>
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
                            <ResultCopier
                              title={"Alignments"}
                              resultSuffix={"snap.dcs/view"}
                              scanId={scanId}
                            />
                        </div>
                      </div>
                    </div>
                    <hr className="vf-divider"/>
                    <div className="alignment-output">
                        { formattedAlignments }
                    </div>
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

              <ResultCard
                title={"Alignments"}
                description={
                  "Alignment of all matches."
                }
                resultSuffix={"snap.dcs/view"}
                fileFormat={"TXT"}
                scanId={scanId}
              />

              <ResultCard
                title={"GFF"}
                description={
                  "GFF (General Feature Format) v3 file listing all found matches."
                }
                resultSuffix={"snap.dcs/gff"}
                fileFormat={"GFF"}
                scanId={scanId}
              />

              <ResultCard
                title={"Original Query"}
                description={"FA (FASTA) file of matched query subsequences."}
                resultSuffix={"snap.dcs/queries"}
                fileFormat={"FASTA"}
                scanId={scanId}
              />

              <ResultCard
                title={"Protein sequence matches"}
                description={
                  "FAA (FASTA Amino Acids) file of matched amino acid sequences."
                }
                resultSuffix={"snap.dcs/aminos"}
                fileFormat={"FASTA"}
                scanId={scanId}
              />

              <ResultCard
                title={"DNA of protein sequences"}
                description={"FA (FASTA) file of codons."}
                resultSuffix={"snap.dcs/codons"}
                fileFormat={"FASTA"}
                scanId={scanId}
              />

              <ResultCard
                title={"HMM Path"}
                description={
                  "Text file, showing match/insertion/deletion states of matches."
                }
                resultSuffix={"snap.dcs/states"}
                fileFormat={"TXT"}
                scanId={scanId}
              />
            </div>
          </section>
        </div>
      )}
    </div>
  );
};
export default Result;
