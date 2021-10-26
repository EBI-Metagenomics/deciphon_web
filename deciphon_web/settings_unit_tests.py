from .settings import *

UNIT_TESTING = True

# Use a single database for unit testing, rather than the dual db setup used in production
DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.sqlite3",
        "NAME": BASE_DIR / "databases/tests.db",
    },
}

DATABASE_ROUTERS = []
