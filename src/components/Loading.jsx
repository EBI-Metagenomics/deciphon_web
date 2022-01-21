import { PropagateLoader } from "react-spinners";

const Loader = ({ isLoading }) => {
  return (
    <div style={{ textAlign: "center", minHeight: "1em" }}>
      <PropagateLoader color={"#18974c"} loading={isLoading} size={10} />
    </div>
  );
};
export default Loader;
