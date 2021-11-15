![unit tests](https://github.com/EBI-Metagenomics/deciphon_web/actions/workflows/test.yaml/badge.svg)
[![codecov](https://codecov.io/gh/EBI-Metagenomics/deciphon_web/branch/master/graph/badge.svg?token=X15S9LH10H)](https://codecov.io/gh/EBI-Metagenomics/deciphon_web)

# Deciphon Web
The web server and client for submitting queries to [Deciphon](https://github.com/EBI-Metagenomics/deciphon).

# Architecture
This is a Django app. 
Deciphon interacts with an sqlite3 database to pull in job, run their queries, and store results.
Deciphon Web wraps this database in [unmanaged Django Models](https://docs.djangoproject.com/en/3.2/ref/models/options/#django.db.models.Options.managed).
It uses Django templates to present a web interface to submit jobs (consisting of queries), and poll for the results. 
We use [Django-Unicorn](https://www.django-unicorn.com) to handle frontend/backend model tying, 
and to do simple partial upsets of the page when polling for job completion.

## Databases
### `default`
A standard Django database in sqlite3. 
This is used for boilerplate Django things like app tables and admin users.

### `deciphon`
The sqlite3 database managed by deciphon itself. 
A [Django DB Router](https://docs.djangoproject.com/en/3.2/topics/db/multi-db/#multiple-databases) makes sure that any `deciphon.models` model read/writes from this database rather than `default`.

**Note** when unit testing, a single database is used so that transactional tests can be used.

# Development
## The basics
- Check out the repository.
- Create and activate a virtualenv or conda env for the project.
- `pip install -r requirements.txt`
- `export DJANGO_SECRET_KEY=<anythingyoulike>`
- `export DECIPHON_WEB_CONFIG=configs/local.yaml`
- `python manage.py migrate`
- `python manage.py collectstatic --noinput`
- `python manage.py runserver 8000`

You can create a new `config.yaml` or modify `configs/local.yaml` to point at different sqlite db files.
E.g. set `django_db_location: /path/to/wherever/deciphon/is/running/deciphon.db`


## Admin interface
The app has Django Admin installed, so you can browse the database models.
- `python manage.py createsuperuser`
- [log in to the admin console](http://127.0.0.1:8000/admin)

## Style
Use [Black](https://black.rtfd.io) to format code before committing.

E.g. `black .` in the repository’s base directory. 

# Testing
`pytest` will run the test suite.

### Selenium integrtation test
`pytest` runs unit tests, as well as a [Selenium](https://pypi.org/project/selenium/) based test of the UI using a browser.
For this test to not fail, you need Chrome and the [Chromedriver](https://chromedriver.chromium.org) installed.
E.g. on a Mac use Homebrew
`brew install --cask chromedriver`
(MacOS will probably complain about launching an unsigned app the first time...
`open /usr/local/Caskroom/chromedriver` to find the chromedriver executable in Finder and then open it once from there to accept the warning.)

# Deployment
