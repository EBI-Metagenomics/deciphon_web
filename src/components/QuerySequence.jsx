import TextareaSequence from "textarea-sequence/dist/textarea-sequence";
import loadWebComponent from "../utils/loadWebComponent";
import { useRef } from "react";

const QuerySequence = ({ onStageSequence }) => {
  loadWebComponent("textarea-sequence", TextareaSequence);
  const textAreaSequenceRef = useRef();
  return (
    <div className="vf-stack vf-stack--400">
      <textarea-sequence
        ref={textAreaSequenceRef}
        height="10em"
        min-sequence-length="5"
        single="false"
        caseSensitive="true"
        alphabet="ACTGU "
      />
      <button
        className="vf-button vf-button--secondary vf-button--sm"
        onClick={async () => {
          console.log(textAreaSequenceRef.current);
          await textAreaSequenceRef.current.cleanUp();
          onStageSequence(textAreaSequenceRef.current.sequence);
        }}
      >
        Check and format query
      </button>
    </div>
  );
};
export default QuerySequence;
