// Source: https://gist.github.com/mbinic/e75a8910ec51a27a041f967e5b3a5345
Cypress.Commands.add('clipboard', () => {
  if (Cypress.browser.name !== 'electron') {
    cy.wrap(
      Cypress.automation('remote:debugger:protocol', {
        command: 'Browser.grantPermissions',
        params: {
          permissions: ['clipboardReadWrite', 'clipboardSanitizedWrite'],
          origin: window.location.origin,
        },
      }).catch((error) =>
        // Electron (v106) will land here, but that's ok, cause the permissions will be granted anyway
        Cypress.log({ message: `Permission request failed: ${error.message}` })
      )
    );
  }

  return cy.window().then(async (win) => {
    win.focus();
    return win.navigator.clipboard;
  });
});
// usage for e.g. writing text
// cy.clipboard().then(clipboard => clipboard.writeText('some text'))
