import { useParams } from "react-router-dom";
import { useEffect, useState } from "react";
import api from "../api";
import Loading from "./Loading";
import { toast } from "react-toastify";

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

const Result = () => {
  let { jobid } = useParams();
  const [jobState, setJobState] = useState();
  const [errors, setErrors] = useState();

  useEffect(() => {
    if (!jobid) return;
    api
      .get(`/jobs/${jobid}`)
      .then((response) => {
        setJobState(response.data);
      })
      .catch((err) => setErrors([err.response.status]));
  }, [jobid]);
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
        </>
      )}
      {jobState?.state === "run" && (
        <>
          <h3>Job is running</h3>
          <span className="vf-badge vf-badge--primary">live updating</span>
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
    </div>
  );
};
export default Result;
