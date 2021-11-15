from django.contrib import admin

# Register your models here.
from deciphon.models import TargetDb, Job, QuerySequence, Result


@admin.register(TargetDb)
class TargetAdmin(admin.ModelAdmin):
    list_display = ("id", "name")


class QueryInline(admin.TabularInline):
    model = QuerySequence
    extra = 0


@admin.register(Job)
class JobAdmin(admin.ModelAdmin):
    def query_count(self, obj):
        return obj.queries.count()

    list_display = ("id", "state", "query_count")
    list_filter = (
        "state",
        "multi_hits",
        "hmmer3_compat",
        "target_db",
    )
    inlines = [
        QueryInline,
    ]


@admin.register(QuerySequence)
class QueryAdmin(admin.ModelAdmin):
    list_display = ("id", "name", "job")


@admin.register(Result)
class ResultAdmin(admin.ModelAdmin):
    list_display = ("id", "job", "seq", "match_id", "prof_name", "model")
