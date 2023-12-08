import {useEffectOnce, useLocalStorage} from "react-use";
import {useState} from "react";
import api from "../api";
import {find, last} from "lodash";

const PreviousJobs = () => {
  const [previousJobs, _, clearPreviousJobs] = useLocalStorage(
    "submittedJobs",
    []
  );

  const [jobsAhead, setJobsAhead] = useState(null);

  const getJobsAhead = () => {
    api
      .get(`/jobs`)
      .then((response) => {
        if (response?.data?.length) {
          const nextPendJob = find(response.data, job => job.state === 'pend');
          const lastSubmittedJob = last(response.data);
          setJobsAhead(parseInt(lastSubmittedJob.id) - parseInt(nextPendJob.id));
        } else {
          setJobsAhead(null);
        }
      })
      .catch((err) => {
        console.error(err);
        setJobsAhead(null);
      });
  }

  useEffectOnce(() => getJobsAhead());

  return (
    <section className="vf-card-container vf-card-container__col-3 | vf-u-background-color--grey--lightest vf-u-fullbleed">
      <div className="vf-card-container__inner">
        {jobsAhead !== null && (
          <div className="vf-section-header">
            <h2 className="vf-section-header__heading">
              Deciphon status
            </h2>
            <p className="vf-text-body vf-text-body--3">
                There are currently {jobsAhead} jobs pending processing.
              <button
                className="vf-button vf-button--link vf-button--sm"
                onClick={() => getJobsAhead()}
              >
                <i className="icon icon-common icon-redo" /> Refresh
              </button>
            </p>
          </div>
        )}
        {!!previousJobs?.length && (
          <>
            <div className="vf-section-header">
              <h2 className="vf-section-header__heading">
                Previously submitted jobs
              </h2>
              <p className="vf-text-body vf-text-body--3">
                Job submission history is stored in your web browser
                <button
                  className="vf-button vf-button--link vf-button--sm"
                  onClick={() => clearPreviousJobs()}
                >
                  <i className="icon icon-common icon-trash" /> Clear history
                </button>
              </p>
            </div>

            {previousJobs.map((job) => (
              <article
                className="vf-card vf-card--brand vf-card--bordered"
                key={job.jobId}
              >
                <div className="vf-card__content | vf-stack vf-stack--400">
                  <h3 className="vf-card__heading">
                    <a className="vf-card__link" href={`/jobs/${job.jobId}`}>
                      Job {job.jobId}{" "}
                      <svg
                        aria-hidden="true"
                        className="vf-card__heading__icon | vf-icon vf-icon-arrow--inline-end"
                        width="1em"
                        height="1em"
                        xmlns="http://www.w3.org/2000/svg"
                      >
                        <path
                          d="M0 12c0 6.627 5.373 12 12 12s12-5.373 12-12S18.627 0 12 0C5.376.008.008 5.376 0 12zm13.707-5.209l4.5 4.5a1 1 0 010 1.414l-4.5 4.5a1 1 0 01-1.414-1.414l2.366-2.367a.25.25 0 00-.177-.424H6a1 1 0 010-2h8.482a.25.25 0 00.177-.427l-2.366-2.368a1 1 0 011.414-1.414z"
                          fill="currentColor"
                          fillRule="nonzero"
                        />
                      </svg>
                    </a>
                  </h3>
                  <p className="vf-card__text">{job.queryNames.join(", ")}</p>
                </div>
              </article>
            ))}
          </>
      )}
      </div>
    </section>
  );
};
export default PreviousJobs;
