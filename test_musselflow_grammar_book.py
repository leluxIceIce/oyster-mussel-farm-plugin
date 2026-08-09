"""Tests pinning the grammar book to the code that actually runs.

The grammar book (``musselflow_grammar_book``) is documentation-as-data.  These
tests guarantee it stays complete and truthful: every key the numerical core
defaults is catalogued, no catalogued key is invented, and the documented
default equals the real default.  A new grammar key therefore cannot be added
to the core without also being documented, and the book cannot drift.
"""

import unittest

import musselflow_ecogrammar_core as core
import musselflow_grammar_book as book


class GrammarBookCompletenessTests(unittest.TestCase):

    def test_book_and_defaults_cover_exactly_the_same_keys(self):
        self.assertEqual(set(book.keys()), set(core.DEFAULTS))

    def test_every_default_value_matches_the_core(self):
        for key, value in core.DEFAULTS.items():
            self.assertEqual(
                book.entry(key)["default"], value,
                "grammar book default for %s disagrees with the core" % key)

    def test_defaults_helper_reproduces_core_defaults(self):
        self.assertEqual(book.defaults(), dict(core.DEFAULTS))

    def test_list_and_bool_and_int_kinds_match_core_key_sets(self):
        kinds = {key: book.entry(key)["kind"] for key in book.keys()}
        self.assertEqual(
            {key for key, kind in kinds.items() if kind == "list"},
            core.LIST_KEYS)
        self.assertEqual(
            {key for key, kind in kinds.items() if kind == "bool"},
            core.BOOL_KEYS)
        self.assertEqual(
            {key for key, kind in kinds.items() if kind == "int"},
            core.INT_KEYS)
        self.assertEqual(
            {key for key, kind in kinds.items() if kind == "string"},
            core.STRING_KEYS)


class GrammarBookIntegrityTests(unittest.TestCase):

    def test_every_entry_has_the_required_fields(self):
        required = {
            "group", "default", "unit", "kind", "summary", "cite", "case",
            "bounds", "bounds_note"}
        for key in book.keys():
            self.assertEqual(set(book.entry(key)) >= required, True,
                             "%s is missing catalog fields" % key)

    def test_case_status_values_are_valid(self):
        for key in book.keys():
            self.assertIn(book.entry(key)["case"], book.CASE_STATUSES)

    def test_citations_resolve_and_carry_a_doi_url(self):
        for key in book.keys():
            record = book.provenance(key)
            if record is not None:
                self.assertIn("url", record)
                self.assertTrue(record["url"].startswith("http"))
                self.assertTrue(record["citation"])
                self.assertTrue(record["equation"])

    def test_every_cite_key_exists_in_the_registry(self):
        for key in book.keys():
            cite = book.entry(key)["cite"]
            if cite is not None:
                self.assertIn(cite, book.CITATIONS)


class GrammarBookCaseGapTests(unittest.TestCase):

    def test_case_gaps_are_the_known_missing_points(self):
        # Only the deferred dry-biomass mode and two compatibility labels remain.
        expected = {
            "profile.species",
            "profile.net",
            "stocking.dry_tissue_kg_per_obstacle",
        }
        self.assertEqual(set(book.case_gaps()), expected)

    def test_render_markdown_produces_a_document(self):
        text = book.render_markdown()
        self.assertIn("# MusselFlow ecological grammar book", text)
        # Every key should appear as a heading.
        for key in book.keys():
            self.assertIn("`%s`" % key, text)


if __name__ == "__main__":
    unittest.main()
