"""One decimal year, one convention - because four of them disagreed.

A decimal year is the calendar year plus the fraction of **that** year already elapsed: 2026.0 is 1 January 2026,
and the denominator is the real length of that year (365 or 366). An average year of 365.25 days is wrong inside a
particular year: it puts a mid-year date about half a day off, and these dates are quoted to the day, so two modules
using different denominators report the same genome on two different days.

Going back to a calendar date the fraction is **floored**: a decimal year names the day it falls on, and rounding
would move the last fraction of a day onto the next one - which is how a figure once showed a date the text did not.

This module imports nothing from the package, so anything may use it.
"""
from datetime import date, timedelta


def year_length(year: int) -> int:
    """Days in that calendar year, by the real leap rule (2000 and 2024 are leap years, 1900 and 2100 are not)."""
    return (date(year + 1, 1, 1) - date(year, 1, 1)).days


def decimal(d: date) -> float:
    """The decimal year of a calendar date. 1 January is the year itself; 31 December is just under the next year."""
    return d.year + (d - date(d.year, 1, 1)).days / year_length(d.year)


# How much of a day a decimal year may fall short of a whole one and still count as that whole one. A date's own decimal
# year is a division that does not come out exactly in binary, so `decimal(d)` lands a few 1e-11 of a day below the day it
# names; flooring that without this would hand back the day before, for about half the days of any year. 1e-6 of a day is
# a tenth of a second: far above that error and far below anything these dates are quoted to.
_EPS_DAYS = 1e-6


def calendar(y: float) -> date:
    """The calendar date a decimal year falls on (the fraction is floored, never rounded)."""
    yr = int(y // 1)
    length = year_length(yr)
    day = int((y - yr) * length + _EPS_DAYS)                 # int() truncates a non-negative fraction, i.e. floors it
    return date(yr, 1, 1) + timedelta(days=min(max(day, 0), length - 1))


def midpoint(a: date, b: date) -> float:
    """The decimal year halfway between two dates - what a genome dated only to a month or a year is given."""
    return (decimal(a) + decimal(b)) / 2
