export default function Spinner({ width, height }) {
  const style = { width: width || "1.0rem", height: height || "1.0rem" };
  return <span className="spinner-border spinner-border-sm" style={style} role="status"></span>;
}
