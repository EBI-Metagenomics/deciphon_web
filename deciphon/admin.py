from django.contrib import admin

# Register your models here.
from deciphon.models import Alphabet, Target, DeciphonUser, Job, Query, Result


@admin.register(Alphabet)
class AlphabetAdmin(admin.ModelAdmin):
    list_display = ("id", "name")


@admin.register(Target)
class TargetAdmin(admin.ModelAdmin):
    list_display = ("id", "name")


@admin.register(DeciphonUser)
class DeciphonUserAdmin(admin.ModelAdmin):
    list_display = ("id", "username", "name")


class QueryInline(admin.TabularInline):
    model = Query
    extra = 0


@admin.register(Job)
class JobAdmin(admin.ModelAdmin):
    def query_count(self, obj):
        return obj.queries.count()

    list_display = ("id", "sid", "status", "query_count")
    list_filter = (
        "status",
        "multiple_hits",
        "hmmer3_compat",
        "abc",
        "target",
        "user",
        "submission",
        "exec_started",
        "exec_ended",
    )
    inlines = [
        QueryInline,
    ]


@admin.register(Query)
class QueryAdmin(admin.ModelAdmin):
    list_display = ("id", "name")


@admin.register(Result)
class ResultAdmin(admin.ModelAdmin):
    pass
