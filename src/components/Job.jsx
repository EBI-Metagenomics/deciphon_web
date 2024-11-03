import React from "react";
import { useInterval } from "react-use";
import api, { pollingInterval } from "../api";
import { find } from "lodash";
import { useEffect, useState } from "react";
import { useParams, useLoaderData } from "react-router-dom";
import { ErrorCard, formatError } from "./Error";
import Snap from "./Snap";
import UrlCopier from "./UrlCopier";
import Progress from "./Progress"
import Spinner from "./Spinner";
import Navigation from "./Navigation";

function JobDone({ job, scan, snapNumProds }) {
  const execEnded = new Date(job.exec_ended).toLocaleString();

  return (
    <div className="vf-stack vf-stack--400">
      <h4>Job complete &mdash; {snapNumProds} matches found</h4>
      <div>
        <span className="vf-form__helper">Finished at: {execEnded}</span>
      </div>
      {<Snap scanId={scan.id} />}
    </div>
  );
}

function JobFail({ job }) {
  return (
    <div>
      <h3>Job has failed</h3>
      {job.error}
      <UrlCopier />
    </div>
  );
}

function JobRun({ job }) {
  return (
    <div>
      <h3>Job is running</h3>
      <Progress progress={job.progress}></Progress>
    </div>
  );
}

function JobPend({ jobsAhead }) {
  return (
    <>
      <h3>Job is pending</h3>
      <span className="vf-badge vf-badge--primary">live updating</span>
      <UrlCopier />
      <p>There are {(jobsAhead == null && <Spinner />) || jobsAhead} jobs ahead of yours in the queue.</p>
    </>
  );
}

async function getJob(jobid) {
  try {
    const x = await api.get(`/jobs/${jobid}`);
    if (x.data.type !== "scan") throw Error(`Job ${jobid} does not exist.`);
    return x.data;
  } catch (error) {
    if (error?.response?.status === 404) throw Error(`Job ${jobid} does not exist.`);
    throw Error(formatError(error))
  }
}

async function getScan(jobId) {
  try {
    const x = await api.get(`/scans?job_id=${jobId}`);
    if (x.data.length !== 1) throw new Error("Expected a single scan.");
    return x.data[0];
  } catch (error) {
    if (error?.response?.status === 404) throw Error(`Job ${jobId} does not exist.`);
    throw Error(formatError(error))
  }
}

async function getSnapNumProds(scanId) {
  try {
    const x = await api.get(`/scans/${scanId}/snap.dcs/prods`);
    return x.data.length;
  } catch (error) {
    if (error?.response?.status === 404) return 0;
    throw Error(formatError(error))
  }
}

export async function loader({ params }) {
  const job = await getJob(params.jobid);
  let scan = null;
  let snapNumProds = null;
  if (job.state === "done") {
    scan = await getScan(params.jobid);
    snapNumProds = await getSnapNumProds(scan.id);
  }
  return { job, scan, snapNumProds };
}

export default function Job() {
  const loaded = useLoaderData()
  const { jobid } = useParams();
  const [jobsAhead, setJobsAhead] = useState();
  const [error, setError] = useState();
  const [job, setJob] = useState(loaded.job);
  const [scan, setScan] = useState(loaded.scan);
  const [snapNumProds, setSnapNumProds] = useState(loaded.snapNumProds);
  const [polling, setPolling] = useState(pollingInterval);

  function fetchJobsAhead(job) {
    api
      .get(`/jobs?limit=30`)
      .then((x) => {
        const nextPend = find(x.data, job => job.state === 'pend' || job.state === 'run');
        if (nextPend) setJobsAhead(parseInt(job.id) - parseInt(nextPend.id));
        else setJobsAhead(0);
      })
      .catch(error => setError(formatError(error)))
  }

  if (job?.type !== "scan") setError(`Job ${job.id} does not exist.`);

  useEffect(() => {
    (async () => {
      if (job?.state === 'pend') fetchJobsAhead(job);
      if (job?.state === 'done' && !scan) {
        const scan = await getScan(job.id);
        setScan(scan);
        setSnapNumProds(await getSnapNumProds(scan.id));
        setPolling(null);
      };
    })();
  }, [job, scan])

  useInterval(() => {
    (async () => {
      if (job?.state === 'pend') {
        fetchJobsAhead(job);
        setJob(await getJob(job.id))
      }
      if (job?.state === 'run') setJob(await getJob(job.id));
      if (job?.state === 'done') {
        const scan = await getScan(job.id);
        setScan(scan);
        setSnapNumProds(await getSnapNumProds(scan.id));
        setPolling(null);
      }
    })();
  }, polling)

  if (error) return <ErrorCard message={error}></ErrorCard>;
  return (<div style={{ minHeight: "348px" }}>
    <Navigation page="query"></Navigation>
    <BreadCrumbs jobid={jobid}></BreadCrumbs>
    <h1>Results</h1>
    {job?.state === 'pend' && <JobPend jobsAhead={jobsAhead}></JobPend>}
    {job?.state === 'run' && <JobRun job={job}></JobRun>}
    {scan && <JobDone job={job} scan={scan} snapNumProds={snapNumProds}></JobDone>}
    {job?.state === 'fail' && <JobFail job={job}></JobFail>}
  </div>)
}

function BreadCrumbs({ jobid }) {
  return (
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
  );
}
