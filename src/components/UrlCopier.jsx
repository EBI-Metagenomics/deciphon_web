import { toast } from "react-toastify";

export default function UrlCopier() {
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
