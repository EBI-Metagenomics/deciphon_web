from django.urls import reverse
from rest_framework import serializers

from deciphon.models import Job, QuerySequence, Result, TargetDb


class QuerySerializer(serializers.ModelSerializer):
    class Meta:
        model = QuerySequence
        fields = ["name", "data", "job"]
        read_only_fields = ["job"]


class TargetDbSerializer(serializers.ModelSerializer):
    class Meta:
        model = TargetDb
        fields = ["id", "name"]
        read_only_fields = ["name"]


class DeciphonJobSerializer(serializers.ModelSerializer):
    queries = QuerySerializer(many=True)
    target_db = TargetDbSerializer()

    class Meta:
        model = Job
        fields = ["state", "error", "target_db", "queries"]
        read_only_fields = ["state", "error"]
