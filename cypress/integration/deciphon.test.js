describe("Deciphon website tests", () => {
  beforeEach(() => {
    cy.intercept("GET", "**/dbs", { fixture: "dbs.json" });

    cy.intercept("GET", "**/jobs/99", { fixture: "job_done.json" });
    cy.intercept("GET", "**/jobs/100", { fixture: "job_pending.json" });
    cy.intercept("GET", "**/jobs/101", { fixture: "job_running.json" });
    cy.intercept("GET", "**/jobs/103", { fixture: "job_failed.json" });

    cy.intercept("GET", "**/scans/9", { fixture: "scan_done.json" });
    cy.intercept("GET", "**/scans/10", { fixture: "scan_pending.json" });
    cy.intercept("GET", "**/scans/11", { fixture: "scan_running.json" });
    cy.intercept("GET", "**/scans/13", { fixture: "scan_failed.json" });

    cy.intercept("GET", "**/jobs/99/scan", { fixture: "scan_done.json" });
    cy.intercept("GET", "**/jobs/100/scan", { fixture: "scan_pending.json" });
    cy.intercept("GET", "**/jobs/101/scan", { fixture: "scan_running.json" });
    cy.intercept("GET", "**/jobs/103/scan", { fixture: "scan_failed.json" });

    cy.intercept("GET", "**/scans/*/prods/gff", { fixture: "prod_gff.txt" }).as(
      "gff"
    );
    cy.intercept("GET", "**/scans/*/prods/fragment", {
      fixture: "prod_frag.txt",
    }).as("frag");
    cy.intercept("GET", "**/scans/*/prods/amino", {
      fixture: "prod_amino.txt",
    }).as("amino");
    cy.intercept("GET", "**/scans/*/prods/codon", {
      fixture: "prod_codon.txt",
    }).as("codon");
    cy.intercept("GET", "**/scans/*/prods/path", {
      fixture: "prod_path.txt",
    }).as("path");
    cy.intercept("GET", "**/scans/*/prods", {
      fixture: "prod_all.json",
    });

    cy.intercept("GET", "**/jobs/next_pend", { fixture: "job_next_pend.json" });

    cy.intercept("POST", "**/scans", {
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
    cy.get(".ql-editor").should("contain.text", "Homoserine_dh-consensus");

    cy.contains("Check and autofix queries").click();

    cy.contains("Submit query").click();
    cy.location().should((loc) => {
      expect(loc.pathname).to.eq("/results/100");
    });
  });

  it("shows pending query", () => {
    cy.visit("http://localhost:3000/results/100");
    cy.contains("Job is pending").should("be.visible");
    cy.contains("There are 10 jobs ahead").should("be.visible");

    cy.get(".icon-copy").click();
    cy.contains("👍").should("be.visible");
    cy.window()
      .its("navigator.clipboard")
      .invoke("readText")
      .should("equal", "http://localhost:3000/results/100");
  });

  it("shows running query", () => {
    cy.visit("http://localhost:3000/results/101");
    cy.contains("Job is running").should("be.visible");
    cy.get("[role=progressbar]").should("have.attr", "aria-valuenow", 50);
  });

  it("shows failed query", () => {
    cy.visit("http://localhost:3000/results/103");
    cy.contains("Job has failed").should("be.visible");
    cy.contains("Something awful happened").should("be.visible");
  });

  it("shows results for successful query", () => {
    cy.visit("http://localhost:3000/results/99");
    cy.contains("Job complete").should("be.visible");
    cy.contains("3").should("be.visible");
    cy.contains("Finished at:").should("be.visible");
    cy.contains("1970").should("be.visible");
    cy.contains("Downloads").should("be.visible");

    const titleToProds = {
      "Download GFF": { clip: "##gff-version", alias: "gff" },
      "Download Fragment": {
        clip: "CCTATCATTTCGACGCTCAAGGAGTCGCTGAC",
        alias: "frag",
      },
      "Download Amino acids": {
        clip: "PIISTLKESLTGDRITRIEGILNGTLNYILTEMEEEGASFSEALKEAQELGYAEADPTDD",
        alias: "amino",
      },
      "Download Codons": {
        clip: "CCTATCATTTCGACGCTCAAGGAGTCGCTGAC",
        alias: "codon",
      },
      "Download HMM Path": {
        clip: "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        alias: "path",
      },
    };

    for (const dlTitle in titleToProds) {
      cy.get("article")
        .filter(`:contains('${dlTitle}')`)
        .find(".icon-copy")
        .click();

      cy.window()
        .its("navigator.clipboard")
        .invoke("readText")
        .should("contain", titleToProds[dlTitle].clip);

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
    cy.contains("Homoserine_dh-consensus");
    cy.get("a[href*='results']").should("contain.text", "Job");
    cy.contains("Clear history").click();
    cy.contains("Previously submitted jobs").should("not.exist");
  });
});
