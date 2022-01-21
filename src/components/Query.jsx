import QuerySequence from "./QuerySequence";
import DatabaseSelection from "./DatabaseSelection";
import { useState } from "react";
import axios from "axios";
import config from "../config/config.json";
import { useNavigate } from "react-router-dom";

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

  axios.post(`${config.API_BASE}/jobs`, data).then((res) => {
    const jobId = res.data.id;
    handleJobSubmitted(jobId);
  });
};

const Query = () => {
  const [selectedDb, setSelectedDb] = useState();
  const [queryText, setQueryText] = useState();
  const nav = useNavigate();
  return (
    <div className={"vf-stack vf-stack--400"}>
      <div>
        <h1> Query Deciphon </h1>
      </div>
      <div className="vf-grid vf-grid__col-3">
        <div className={"vf-grid__col--span-2"}>
          <QuerySequence onStageSequence={setQueryText} />
        </div>
        <DatabaseSelection selectedDb={selectedDb} onSelectDb={setSelectedDb} />
      </div>
      <div style={{ textAlign: "center" }}>
        <button
          className="vf-button vf-button--primary"
          onClick={() => {
            submitJob(selectedDb, queryText, (jobId) =>
              nav(`/results/${jobId}`)
            );
          }}
          id="submit"
          disabled={!selectedDb || !queryText}
        >
          Submit query
        </button>
      </div>
    </div>
  );
};
export default Query;
