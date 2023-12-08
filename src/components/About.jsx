import ont from "../img/ont.svg"

const About = () => {
    return (
        <>
          <nav className="vf-navigation vf-navigation--main | vf-cluster">
            <ul className="vf-navigation__list | vf-list | vf-cluster__inner">
              <li className="vf-navigation__item">
                <a href="/" className="vf-navigation__link">Query</a>
              </li>
              <li className="vf-navigation__item">
                <a href="/about" className="vf-navigation__link" aria-current="page">About</a>
              </li>
            </ul>
          </nav>
          <div className={"vf-stack vf-stack--400"}>
            <div>
                <h1> About Deciphon </h1>
                <p className="vf-u-type__text-body--2">
                    Deciphon predicts protein annotations on short- to long-read sequence data, in the presence of technical errors.
                    It uses <a className="vf-link" href="https://github.com/EBI-Metagenomics/imm" target="_newtab">Invisible Markov Models (IMMs)</a> to infer proteins using the Viterbi method.

                </p>
                <p className="vf-u-type__text-body--2">
                    The software can be downloaded from the <a className="vf-link" href="https://github.com/EBI-Metagenomics/deciphon" target="_newtab">EBI-Metagenomics GitHub</a>.
                </p>
                <div className="vf-flag vf-flag--top vf-flag--400">

                  <div className="vf-flag__media">
                    <img src={ont} alt="Oxford Nanopore Technologies logo" height="40px"/>
                  </div>

                  <div className="vf-flag__body">
                    <p className="vf-u-type__text-body--2 vf-u-margin--0">
                        Deciphon’s development has been funded by <a className="vf-link" href="https://nanoporetech.com" target="_newtab">Oxford Nanopore Technologies</a> in collaboration with <a className="vf-link" href="https://www.ebi.ac.uk" target="_newtab">EMBL-EBI</a>.
                    </p>
                  </div>
                </div>

                <h2> Databases </h2>
                <a id="databases"/>
                <h3> Pfam </h3>
                <p className="vf-u-type__text-body--2">
                    Deciphon’s "Pfam" database refers to the <a className="vf-link" href="https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam36.0/Pfam-A.hmm.gz">Pfam-A 36.0</a> protein database processed using the NCBI translation table 11 (the Bacterial, Archaeal and Plan Plastic Code).
                </p>
            </div>
          </div>
        </>
    )
}

export default About;