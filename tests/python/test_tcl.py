from annotsv import tcl


def test_split_list_quotes_and_braces():
    assert tcl.split_list('a "b c" {d {e}} ""') == ["a", "b c", "d {e}", ""]
    assert tcl.split_list("") == []


def test_split_chars_like_tcl():
    assert tcl.split_chars("", ",;") == []
    assert tcl.split_chars("a,b;c", ",;") == ["a", "b", "c"]
    assert tcl.split_chars("a,,b", ",") == ["a", "", "b"]


def test_compare_numeric_then_string():
    assert tcl.compare("9", "10") < 0
    assert tcl.compare("b", "a") > 0
    assert tcl.compare("1.5", "1.50") == 0


def test_sort_on_element_is_stable():
    lines = ["1\t200\t300", "1\t100\t500", "1\t100\t400"]
    out = tcl.sort_on_element(tcl.sort_on_element(lines, 2), 1)
    assert out == ["1\t100\t400", "1\t100\t500", "1\t200\t300"]
    assert tcl.sort_on_element(lines, 1, descending=True)[0] == "1\t200\t300"


def test_format_list():
    assert tcl.format_list(["a", "b c", ""]) == "a {b c} {}"
