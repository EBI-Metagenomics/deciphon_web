import UrlCopier from "./UrlCopier";

export default function Progress({ progress }) {
  return (
    <div>
      <span className="vf-badge vf-badge--primary">live updating</span>
      <div
        role="progressbar"
        aria-valuenow={progress}
        aria-valuemin="0"
        aria-valuemax="100"
        className="vf-progress-indicator"
      >
        <div
          className="vf-progress-indicator__mark"
          style={{
            "--vf-progress-indicator__percent": `${progress}%`,
          }}
        />
        <p className="vf-progress-indicator__helper-text">
          Progress reported by Deciphon Job Server
        </p>
      </div>
      <UrlCopier />
    </div>
  );
}
