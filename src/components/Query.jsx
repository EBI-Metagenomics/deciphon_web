import QuerySequence from "./QuerySequence";
import DatabaseSelection from "./DatabaseSelection";
import { useState } from "react";
import { useNavigate } from "react-router-dom";
import PreviousJobs from "./PreviousJobs";
import { useLocalStorage } from "react-use";
import api from "../api";

const submitJob = (dbId, queryText, handleJobSubmitted) => {
  const seqs = queryText.split(/(?=>)/g);

  const data = {
    db_id: dbId,
    multi_hits: true,
    hmmer3_compat: false,
    seqs: seqs.map((seq) => ({
      name: seq.split("\n")[0].replace(/>/, "").trim(),
      data: seq.split("\n").slice(1).join("").replaceAll(" ", "").toUpperCase(),
    })),
  };

  api.post(`/scans/`, data).then((res) => {
    const jobId = res.data.job.id;
    handleJobSubmitted(
      jobId,
      data.seqs.map((seq) => seq.name)
    );
  });
};

const Query = () => {
  const [selectedDb, setSelectedDb] = useState();
  const [queryText, setQueryText] = useState();
  const nav = useNavigate();
  const [previousJobs, setPreviousJobs] = useLocalStorage("submittedJobs", []);
  return (
    <>
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
      <div className={"vf-stack vf-stack--400"}>
        <div>
          <h1> Query Deciphon </h1>
        </div>
        <div className="vf-grid vf-grid__col-3">
          <div className={"vf-grid__col--span-2"}>
            <QuerySequence onStageSequence={setQueryText} />
          </div>
          <DatabaseSelection
            selectedDb={selectedDb}
            onSelectDb={setSelectedDb}
          />
        </div>
        <div style={{ textAlign: "center" }}>
          <button
            className="vf-button vf-button--primary"
            onClick={() => {
              submitJob(selectedDb, queryText, (jobId, queryNames) => {
                setPreviousJobs([{ jobId, queryNames }, ...previousJobs]);
                nav(`/jobs/${jobId}`);
              });
            }}
            id="submit"
            disabled={!selectedDb || !queryText}
          >
            Submit query
          </button>
        </div>
      </div>
      <PreviousJobs />
    </>
  );
};
export default Query;
