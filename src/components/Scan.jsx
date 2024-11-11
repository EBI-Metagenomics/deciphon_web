import { useEffect, useState } from "react";
import { useParams } from "react-router-dom";
import api from "../api";
import { SnapComponent } from "./Snap";
import { ErrorCard } from "./Error";


export function ScanResult() {
  let { jobid } = useParams();
  const [scanId, setScanId] = useState(null);
  const [error, setError] = useState();

  useEffect(() => {
    api
      .get(`/scans?job_id=${jobid}`)
      .then(x => {
        if (x.data.length != 1) throw new Error("Expected a single scan.");
        return x.data[0].id;
      })
      .then(scanId => setScanId(scanId))
      .catch(error => setError(error))
  }, [jobid]);
  return (
    <div>
      {error && (<ErrorCard error={error}></ErrorCard>)}
      {scanId != null && (<SnapComponent scanId={scanId} />)}
    </div>
  );
};
