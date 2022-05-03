import { useParams } from "react-router-dom";
import { useEffect, useState } from "react";
import api, { baseUrl, pollingInterval } from "../api";
import Loading from "./Loading";
import { toast } from "react-toastify";
import { useInterval } from "react-use";

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
                .get(`/scans/${scanId}/prods/${resultSuffix}`)
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
            `${baseUrl}/scans/${scanId}/prods/${resultSuffix}`,
            "_blank"
          );
        }}
      >
        <i className="icon icon-common icon-download" />
        &nbsp;Download
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
          &nbsp; Download {title}
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

  useEffect(() => {
    api.get(`/jobs/${jobid}/scan`).then((response) => {
      setScanId(response?.data?.id);
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
              .get("/jobs/next_pend")
              .then((response) => {
                setJobsAhead(parseInt(jobid) - response.data.id);
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
      api.get(`/scans/${scanId}/prods`).then((response) => {
        setNumResults(response?.data?.length);
      });
    }
  }, [jobState]);

  const finishedAt = jobState
    ? new Date(1000 * parseInt(jobState?.exec_ended || "0")).toLocaleString()
    : null;

  return (
    <div className={"vf-stack vf-stack--400"}>
      <nav className="vf-breadcrumbs" aria-label="Breadcrumb">
        <ul className="vf-breadcrumbs__list | vf-list vf-list--inline">
          <li className="vf-breadcrumbs__item">
            <a href="/" className="vf-breadcrumbs__link">
              Home
            </a>
          </li>
          <li className="vf-breadcrumbs__item" aria-current="location">
            Job {jobid}
          </li>
        </ul>
      </nav>
      <h1>Query results</h1>

      {jobState?.state === "pend" && (
        <>
          <h3>Job is pending</h3>
          <span className="vf-badge vf-badge--primary">live updating</span>
          <UrlCopier />
          {jobsAhead !== null && jobsAhead > 0 && (
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
        <div className="vf-stack vf-stack--800">
          <h3>Job complete</h3>
          {numResults !== null && (
            <div className="vf-flag vf-flag--bottom vf-flag--200">
              <div className="vf-flag__media">
                <p className="vf-lede">{numResults}</p>
              </div>
              <div className="vf-flag__body">
                <p className="vf-u-type__text-body--3 vf-u-margin--0">
                  matches found
                </p>
              </div>
            </div>
          )}
          <div>
            <span className="vf-form__helper">Finished at: {finishedAt}</span>
          </div>
          <section className="vf-card-container vf-card-container__col-3 | vf-u-background-color--grey--lightest vf-u-fullbleed">
            <div className="vf-card-container__inner">
              <div className="vf-section-header">
                <h2 className="vf-section-header__heading">Downloads</h2>
                <p className="vf-section-header__text">
                  Results files from your queries
                </p>
              </div>

              <ResultCard
                title={"GFF"}
                description={
                  "GFF (General Feature Format) v3 file listing all found features."
                }
                resultSuffix={"gff"}
                fileFormat={"GFF"}
                scanId={scanId}
              />

              <ResultCard
                title={"Fragments"}
                description={"FA (FASTA) file of matched fragment sequences."}
                resultSuffix={"fragment"}
                fileFormat={"FASTA"}
                scanId={scanId}
              />

              <ResultCard
                title={"Amino acids"}
                description={
                  "FAA (FASTA Amino Acids) file of matched amino acid sequences."
                }
                resultSuffix={"amino"}
                fileFormat={"FASTA"}
                scanId={scanId}
              />

              <ResultCard
                title={"Codons"}
                description={"FA (FASTA) file of codons."}
                resultSuffix={"codon"}
                fileFormat={"FASTA"}
                scanId={scanId}
              />

              <ResultCard
                title={"HMM Path"}
                description={
                  "FA (FASTA) file, showing match/insertion/deletion states of matches."
                }
                resultSuffix={"path"}
                fileFormat={"FASTA"}
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
