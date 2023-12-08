import { useEffect, useState } from "react";
import Loader from "./Loading";
import api from "../api";

const DatabaseSelection = ({ selectedDb, onSelectDb }) => {
  const [dbs, setDbs] = useState();

  useEffect(() => {
    api
      .get(`dbs`)
      .then((result) => {
        setDbs(result.data);
        if (selectedDb === undefined && result.data.length) {
          onSelectDb(parseInt(result.data[0].id));
        }
      })
      .catch((error) => {
        console.error(error);
      });
  }, [selectedDb, onSelectDb]);

  return (
    <fieldset className="vf-form__fieldset | vf-stack vf-stack--400">
      <legend className="vf-form__legend">Target database</legend>
      <p className="vf-form__helper">
        See the <a href={"/about#databases"} className="vf-link">About</a> page for database details.
      </p>
      <Loader isLoading={!dbs} />
      <div className="vf-cluster vf-cluster--400">
        <div className={"vf-cluster__inner"}>
          {!!dbs &&
            dbs.map((db, idx) => (
              <div className="vf-form__item vf-form__item--radio" key={db.id}>
                <input
                  type="radio"
                  name="target"
                  value={db.id}
                  id={`target_${db.id}`}
                  className={"vf-form__radio"}
                  onChange={(e) => onSelectDb(parseInt(e.target.value))}
                  defaultChecked={
                    db.id === selectedDb || (!selectedDb && idx === 0)
                  }
                />
                <label htmlFor={`target_${db.id}`} className="vf-form__label">
                  {db.file.name.replace(/.dcp$/, "").toUpperCase()}
                </label>
              </div>
            ))}
        </div>
      </div>
    </fieldset>
  );
};
export default DatabaseSelection;
