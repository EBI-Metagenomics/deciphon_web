export default function Footer() {
  return (
    <div class="full-width">
      <footer class="vf-footer">
        <div class="vf-footer__inner">
          <p class="vf-footer__notice">Deciphon: from long reads to protein annotations.</p>
          <p class="vf-footer__notice">
            Deciphon’s development has been funded by <a class="vf-list__link" href="https://nanoporetech.com"
              target="_newtab">Oxford Nanopore Technologies</a> in collaboration with <a class="vf-list__link"
                href="https://www.ebi.ac.uk" target="_newtab">EMBL-EBI</a>.
          </p>
          <div class="vf-footer__links-group vf-grid">
            <div class="vf-links">
              <h4 class="vf-links__heading">Deciphon</h4>
              <ul class="vf-links__list vf-list">
                <li class="vf-list__item">
                  <a class="vf-list__link" href="https://github.com/EBI-Metagenomics/deciphon">Software source</a>
                </li>
              </ul>
            </div>
            <div class="vf-links">
              <h4 class="vf-links__heading">Support</h4>
              <ul class="vf-links__list vf-list">
                <li class="vf-list__item">
                  <a class="vf-list__link" href="mailto:metagenomics-help@ebi.ac.uk">Contact</a>
                </li>
                <li class="vf-list__item">
                  <a class="vf-list__link" href="https://www.ebi.ac.uk/research/finn/">Finn Group</a>
                </li>
              </ul>
            </div>
          </div>
          <section class="vf-footer__legal vf-grid vf-grid__col-1">
            <ul class="vf-footer__list vf-footer__list--legal vf-list vf-list--inline">
              <li class="vf-list__item">
                <a href="https://www.ebi.ac.uk/about/terms-of-use" class="vf-list__link">Terms of use</a>
              </li>
            </ul>
            <p class="vf-footer__legal-text">Copyright © EMBL.</p>
          </section>
        </div>
      </footer>
    </div>
  )
}
