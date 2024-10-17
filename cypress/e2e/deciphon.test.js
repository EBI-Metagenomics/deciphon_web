Cypress.Commands.add('assertValueCopiedToClipboard', expected => {
  cy.window().then(win => {
    win.navigator.clipboard.readText().then(text => {
      expect(text).to.contains(expected)
    })
  })
})

describe("Deciphon website tests", () => {
  beforeEach(() => {
    cy.intercept("GET", "http://api/dbs", { fixture: "dbs.json" });

    cy.intercept("GET", "http://api/jobs/99", { fixture: "job_done.json" });
    cy.intercept("GET", "http://api/jobs/100", { fixture: "job_pending.json" });
    cy.intercept("GET", "http://api/jobs/101", { fixture: "job_running.json" });
    cy.intercept("GET", "http://api/jobs/103", { fixture: "job_failed.json" });

    cy.intercept("GET", "http://api/scans?job_id=99", { fixture: "scan_done.json" });
    cy.intercept("GET", "http://api/scans?job_id=100", { fixture: "scan_pending.json" });
    cy.intercept("GET", "http://api/scans?job_id=101", { fixture: "scan_running.json" });
    cy.intercept("GET", "http://api/scans?job_id=103", { fixture: "scan_failed.json" });

    cy.intercept("GET", "http://api/scans/*/snap.dcs/view", {
      fixture: "prod_alignments.txt",
    }).as("alignments");
    cy.intercept("GET", "http://api/scans/*/snap.dcs/gff", {
      fixture: "prod_gff.txt" }
    ).as("gff");
    cy.intercept("GET", "http://api/scans/*/snap.dcs/queries", {
      fixture: "prod_queries.txt" }
    ).as("queries");
    cy.intercept("GET", "http://api/scans/*/snap.dcs/aminos", {
      fixture: "prod_amino.txt",
    }).as("aminos");
    cy.intercept("GET", "http://api/scans/*/snap.dcs/codons", {
      fixture: "prod_codon.txt",
    }).as("codons");
    cy.intercept("GET", "http://api/scans/*/snap.dcs/states", {
      fixture: "prod_hmm.txt",
    }).as("states");
    cy.intercept("GET", "http://api/scans/*/snap.dcs/prods", {
      fixture: "prod_all.json",
    });

    cy.intercept("GET", "http://api/jobs?limit=*", { fixture: "jobs_list.json" });

    cy.intercept("POST", "http://api/scans", {
      fixture: "scan_new.json",
    });

    cy.wrap(
      Cypress.automation("remote:debugger:protocol", {
        command: "Browser.grantPermissions",
        params: {
          permissions: ["clipboardReadWrite", "clipboardSanitizedWrite"],
          origin: "http://localhost:3000",
        },
      })
    );
    cy.window().focus();
  });

  it("loads successfully", () => {
    cy.visit("http://localhost:3000");
    cy.get("h1").contains("Query Deciphon");

    cy.contains("PFAM24").parent().find("input").should("be.checked");

    cy.contains("Submit query").should("be.disabled");
  });

  it("reports query errors", () => {
    cy.visit("http://localhost:3000");
    cy.get(".ql-editor[contenteditable]").should("be.visible");
    cy.get(".ql-editor[contenteditable]").type("nnn");

    cy.contains("invalid alphabet").should("be.visible");
    cy.contains("missing headers").should("be.visible");
    cy.contains("sequence length").should("be.visible");
    cy.get(".ql-editor[contenteditable]").type("aaa");
    cy.contains("invalid alphabet").should("be.visible");

    cy.contains("missing headers").should("be.visible");
    cy.get(".vf-badge.vf-badge--secondary").should("have.length", 2);
    cy.contains("Check and autofix queries").click();
    cy.get(".ql-editor").should("contain.text", "Generated Header");
    cy.get(".vf-badge.vf-badge--secondary").should("have.length", 1);
    cy.contains("sequence length").should("be.visible");
  });

  it("loads example query", () => {
    cy.visit("http://localhost:3000");
    cy.contains("Load an example query").click();
    cy.get(".ql-editor").should("contain.text", "example");

    cy.contains("Check and autofix queries").click();

    cy.contains("Submit query").click();
    cy.location().should((loc) => {
      expect(loc.pathname).to.eq("/jobs/100");
    });
  });

  it("shows pending query", () => {
    cy.visit("http://localhost:3000/jobs/100");
    cy.contains("Job is pending").should("be.visible");
    cy.contains("There are 10 jobs ahead").should("be.visible");
    cy.get("a").first().focus();
    cy.get(".icon-copy").click();
    cy.contains("👍").should("be.visible");
    cy.assertValueCopiedToClipboard("http://localhost:3000/jobs/100");
  });

  it("shows running query", () => {
    cy.visit("http://localhost:3000/jobs/101");
    cy.contains("Job is running").should("be.visible");
    cy.get("[role=progressbar]").should("have.attr", "aria-valuenow", 50);
  });

  it("shows failed query", () => {
    cy.visit("http://localhost:3000/jobs/103");
    cy.contains("Job has failed").should("be.visible");
    cy.contains("Something awful happened").should("be.visible");
  });

  it("shows jobs for successful query", () => {
    cy.visit("http://localhost:3000/jobs/99");
    cy.contains("Job complete").should("be.visible");
    cy.contains("3 matches found").should("be.visible");
    cy.contains("Finished at:").should("be.visible");
    cy.contains("2023").should("be.visible");
    cy.contains("Results files from your search").should("be.visible");

    const titleToProds = {
      "GFF": { clip: "##gff-version", alias: "gff" },
      "Original Query": {
        clip: "ATTTCGACGCTCAAGGAGTCGCTGA",
        alias: "queries",
      },
      "Protein sequence matches": {
        clip: "ISTLKESLIGDRITRIEGILNGTMNYILTEMEEEGASFSEALKEAQQLGYAEADPTDDVE",
        alias: "amino",
      },
      "DNA of protein sequences": {
        clip: "ATTTCGACGCTCAAGGAGTCGCTGATAGGTGACCGTATTACTCGAATCGAAGGGATATTA",
        alias: "codon",
      },
      "HMM Path": {
        clip: "SBM3M4M5M6M7M8M9M10M11M12M13M14M15M16M17M18M19M20M21M22M23M2",
        alias: "states",
      },
    };

    for (const dlTitle in titleToProds) {
      cy.get("a").first().focus();
      cy.get("article")
        .filter(`:contains('${dlTitle}')`)
        .find(".icon-copy")
        .click();
      cy.assertValueCopiedToClipboard(titleToProds[dlTitle].clip)

      cy.get("article")
        .filter(`:contains('${dlTitle}')`)
        .find(".icon-download")
        .parent()
        .should("be.visible");
    }
  });

  it("shows previously submitted job history", () => {
    cy.visit("http://localhost:3000");
    cy.contains("Load an example query").click();
    cy.contains("Check and autofix queries").click();
    cy.contains("Submit query").click();
    cy.visit("http://localhost:3000");
    cy.contains("Previously submitted jobs");
    cy.contains("example");
    cy.get("a[href*='jobs']").should("contain.text", "Job");
    cy.contains("Clear history").click();
    cy.contains("Previously submitted jobs").should("not.exist");
  });
});
