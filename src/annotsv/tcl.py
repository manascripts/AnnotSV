"""Small helpers that reproduce Tcl semantics the original code relies on.

The Tcl implementation treats tab-separated lines as Tcl lists (``lindex``),
compares values with ``expr`` (numeric when both sides look like numbers,
string comparison otherwise) and prints lists in Tcl's canonical form. Later
modules need the same behaviour to keep outputs byte-comparable.
"""

from __future__ import annotations

import functools
import re
from typing import Iterable, List, Sequence

_WS = " \t\n\r\v\f"


def split_list(text: str) -> List[str]:
    """Parse ``text`` as a Tcl list (whitespace, ``"..."`` and ``{...}`` aware).

    Falls back to a plain whitespace split if the string is not a well formed
    Tcl list (Tcl would raise an error in that case).
    """
    items: List[str] = []
    i, n = 0, len(text)
    try:
        while True:
            while i < n and text[i] in _WS:
                i += 1
            if i >= n:
                return items
            c = text[i]
            if c == "{":
                depth, j = 1, i + 1
                while j < n and depth:
                    if text[j] == "\\":
                        j += 2
                        continue
                    if text[j] == "{":
                        depth += 1
                    elif text[j] == "}":
                        depth -= 1
                    j += 1
                if depth:
                    raise ValueError("unmatched open brace")
                items.append(text[i + 1 : j - 1])
                i = j
            elif c == '"':
                j, buf = i + 1, []
                while j < n and text[j] != '"':
                    if text[j] == "\\" and j + 1 < n:
                        buf.append(_unescape(text[j + 1]))
                        j += 2
                        continue
                    buf.append(text[j])
                    j += 1
                if j >= n:
                    raise ValueError("unmatched quote")
                items.append("".join(buf))
                i = j + 1
            else:
                j, buf = i, []
                while j < n and text[j] not in _WS:
                    if text[j] == "\\" and j + 1 < n:
                        buf.append(_unescape(text[j + 1]))
                        j += 2
                        continue
                    buf.append(text[j])
                    j += 1
                items.append("".join(buf))
                i = j
            if i < n and text[i] not in _WS:
                raise ValueError("list element followed by garbage")
    except ValueError:
        return text.split()


def _unescape(ch: str) -> str:
    return {"n": "\n", "t": "\t", "r": "\r"}.get(ch, ch)


def lindex(text: str, index: int) -> str:
    """``lindex $text $index`` on a string (empty string when out of range)."""
    items = split_list(text)
    return items[index] if 0 <= index < len(items) else ""


def split_chars(text: str, chars: str) -> List[str]:
    """Tcl ``split $text $chars``: every character of ``chars`` is a separator.

    Like Tcl, splitting an empty string returns an empty list.
    """
    if text == "":
        return []
    return re.split("[" + re.escape(chars) + "]", text)


def format_list(items: Iterable[str]) -> str:
    """Render a Python sequence the way Tcl prints a list."""
    out = []
    for item in items:
        item = str(item)
        if item == "":
            out.append("{}")
        elif any(c in item for c in _WS + '{}"\\[]$;'):
            out.append("{" + item + "}")
        else:
            out.append(item)
    return " ".join(out)


def to_number(value: str):
    """Return the value as int/float if Tcl's ``expr`` would see a number."""
    s = value.strip() if isinstance(value, str) else value
    if not isinstance(s, str):
        return s
    try:
        return int(s, 0) if not re.fullmatch(r"[+-]?0\d+", s) else int(s, 8)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        return None


def compare(a: str, b: str) -> int:
    """Three-way comparison with Tcl ``expr`` semantics (numeric if possible)."""
    na, nb = to_number(a), to_number(b)
    if na is not None and nb is not None:
        return (na > nb) - (na < nb)
    sa, sb = str(a), str(b)
    return (sa > sb) - (sa < sb)


def sort_on_element(lines: Sequence[str], index: int, descending: bool = False) -> List[str]:
    """Port of ``lsort -command (Ascending|Descending)SortOnElementN``.

    Tcl's merge sort with those comparators is a stable sort on the Nth list
    element, so a Python stable sort with the same comparison is equivalent.
    """

    def cmp(x: str, y: str) -> int:
        c = compare(lindex(x, index), lindex(y, index))
        return -c if descending else c

    return sorted(lines, key=functools.cmp_to_key(cmp))


def string_is_double(value: str) -> bool:
    """Tcl ``string is double`` (non strict: the empty string is accepted)."""
    if value == "":
        return True
    return to_number(value) is not None
