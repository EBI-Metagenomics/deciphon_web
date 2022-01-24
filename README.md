![unit tests](https://github.com/EBI-Metagenomics/deciphon_web/actions/workflows/test.yaml/badge.svg)
[![codecov](https://codecov.io/gh/EBI-Metagenomics/deciphon_web/branch/master/graph/badge.svg?token=X15S9LH10H)](https://codecov.io/gh/EBI-Metagenomics/deciphon_web)

# Deciphon Web

The web client for submitting queries to [Deciphon](https://github.com/EBI-Metagenomics/deciphon) and viewing results.

# Architecture

This is a React app.
Deciphon Web interacts with Deciphon via [Deciphon's REST client](https://github.com/EBI-Metagenomics/deciphon/tree/main/rest).

# Development

## The basics

- Check out the repository.
- Install `npm` if you haven't got it. e.g. `brew install node`
- You will need either Deciphon running, to use its REST API, or mock the API using e.g. [Postman](https://www.postman.com).
  - If you're using Postman, there is a collection in this repo `mock_api.postman_collection.json`, that you can import and set up a Mock Server with.
    - This will give you a Mock Server URL, running on Postman's servers, to use for testing.
    - Job ID 99 is a completed job with results. 100 is a pending job, 101 is a running job, and 103 is a failed job.
  - Check the config in `src/config/config.json` to set the API URL (either your local Deciphon or the Postman server).
- `npm install` to install the dependencies.

## Style

Use [Prettier](https://prettier.io) to format code before committing.

E.g. `npx prettier --write .` in the repository’s base directory.

# Testing

TODO
`npm test`

# Use

`npm start`

Browse to [the web interface](http://127.0.0.1:3000).

# Deployment

`npm run build`
