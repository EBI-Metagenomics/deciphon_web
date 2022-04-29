import TextareaSequence from "textarea-sequence/dist/textarea-sequence";
import loadWebComponent from "../utils/loadWebComponent";
import { useEffect, useRef, useState } from "react";

const errorBadgeStyle = {
  borderColor: "var(--vf-ui-color--red)",
  color: "var(--vf-ui-color--red)",
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
      <textarea-sequence
        ref={textAreaSequenceRef}
        height="10em"
        min-sequence-length="5"
        single={false}
        caseSensitive="true"
        alphabet="ACTGU "
      />
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
            onStageSequence(textAreaSequenceRef.current.sequence);
            const { multipleSequences, ...otherErrors } =
              textAreaSequenceRef.current.errors;
            if (
              multipleSequences &&
              Object.values(otherErrors).every((e) => !e)
            ) {
              //  The only "error" is multiple sequences present.
              //  This isn't actually an error, but a Nightingale bug means the textarea is still bordered red.
              document.getElementById("sequence-editor").style.border =
                "1px solid green";
            }
          }}
        >
          Check and autofix queries
        </button>
      </div>
    </div>
  );
};
export default QuerySequence;
