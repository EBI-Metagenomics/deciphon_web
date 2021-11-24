from datetime import datetime

from django import template

register = template.Library()


@register.filter(name="to_datetime")
def datetime_from_timestamp(value, fmt):
    """Formats a unix timestamp as a datetime"""
    return datetime.fromtimestamp(value).strftime(fmt)
