import TextareaSequence from "textarea-sequence/dist/textarea-sequence";
import loadWebComponent from "../utils/loadWebComponent";
import { useEffect, useRef, useState } from "react";
import exampleQuery from "../utils/exampleQuery";

const errorBadgeStyle = {
  borderColor: "var(--vf-ui-color--red)",
  color: "var(--vf-ui-color--red)",
};

const workaroundMultipleSequences = (textarea) => {
  const { multipleSequences, ...otherErrors } = textarea.errors;
  if (multipleSequences && Object.values(otherErrors).every((e) => !e)) {
    //  The only "error" is multiple sequences present.
    //  This isn't actually an error, but a Nightingale bug means the textarea is still bordered red.
    document.getElementById("sequence-editor").style.border = "1px solid green";
  }
};

const QuerySequence = ({ onStageSequence }) => {
  loadWebComponent("textarea-sequence", TextareaSequence);
  const textAreaSequenceRef = useRef();
  const [errors, setErrors] = useState({});

  useEffect(() => {
    const handleErrorChange = (e) => {
      setErrors(e.detail.errors);
    };
    if (!textAreaSequenceRef?.current) return;
    textAreaSequenceRef.current.addEventListener(
      "error-change",
      handleErrorChange
    );
    return () =>
      textAreaSequenceRef.current.removeEventListener(
        "error-change",
        handleErrorChange
      );
  }, [textAreaSequenceRef]);
  return (
    <div className="vf-stack vf-stack--400">
      <div className="vf-stack vf-stack--200">
        <button
          className="vf-button vf-button--link vf-button--sm"
          onClick={async () => {
            await textAreaSequenceRef.current.quill.setText(exampleQuery);
            await textAreaSequenceRef.current.cleanUp();
            workaroundMultipleSequences(textAreaSequenceRef.current);
          }}
        >
          Load an example query
        </button>
        <textarea-sequence
          ref={textAreaSequenceRef}
          height="10em"
          min-sequence-length="5"
          single={false}
          caseSensitive="true"
          alphabet="ACTGU "
        />
      </div>
      <div className="vf-cluster vf-cluster--200">
        <div className="vf-cluster__inner">
          {errors.hasInvalidCharacters && (
            <span
              className="vf-badge vf-badge--secondary"
              style={errorBadgeStyle}
            >
              invalid alphabet
            </span>
          )}
          {(errors.missingFirstHeader ||
            errors.headerCheckRequiredForMultipleSequences) && (
            <span
              className="vf-badge vf-badge--secondary"
              style={errorBadgeStyle}
            >
              missing headers
            </span>
          )}
          {errors.tooShort && (
            <span
              className="vf-badge vf-badge--secondary"
              style={errorBadgeStyle}
            >
              sequence length
            </span>
          )}
        </div>
      </div>
      <div>
        <button
          className="vf-button vf-button--secondary vf-button--sm"
          onClick={async () => {
            await textAreaSequenceRef.current.cleanUp();
            workaroundMultipleSequences(textAreaSequenceRef.current);
            onStageSequence(textAreaSequenceRef.current.sequence);
          }}
        >
          Check and autofix queries
        </button>
      </div>
    </div>
  );
};
export default QuerySequence;
