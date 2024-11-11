import "@cypress/code-coverage/support";
import "cypress-real-events";
import "./commands";

beforeEach(() => {
  cy.setCookie("cookies-accepted", "true");
  cy.setCookie("embl-ebi-public-website-v1.0-data-protection-accepted", "true");
});
