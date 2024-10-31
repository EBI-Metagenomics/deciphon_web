export function ErrorCard({ message }) {
  return (
    <div className="vf-grid vf-grid__col-2">
      <article
        className="vf-card vf-card--brand vf-card--bordered"
        style={{ "--vf-card-border-color": "#d41645" }}
      >
        <div className="vf-card__content | vf-stack vf-stack--400">
          <h3 className="vf-card__heading">Error :(</h3>
          <p className="vf-card__subheading">{message}</p>
        </div>
      </article>
    </div>
  )
}

export function formatError(error) {
  if (typeof a === "string") return error;
  if (error?.response?.data?.detail) return error.response?.data?.detail;
  return error.message;
}
