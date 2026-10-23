import argparse
import unittest
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch

import esm
import torch

from esm_pll_pytorch27 import load_local_checkpoint


class SafeCheckpointLoadingTest(unittest.TestCase):
    def test_allowlists_only_argparse_namespace_while_loading(self):
        state = {"inside_safe_globals": False, "allowed": None}
        marker = object()

        @contextmanager
        def safe_globals(values):
            state["allowed"] = values
            state["inside_safe_globals"] = True
            try:
                yield
            finally:
                state["inside_safe_globals"] = False

        def loader(path):
            self.assertTrue(state["inside_safe_globals"])
            self.assertEqual(path, "/trusted/checkpoint.pt")
            return marker, "alphabet"

        with patch.object(torch.serialization, "safe_globals", side_effect=safe_globals), \
             patch.object(esm.pretrained, "load_model_and_alphabet_local", side_effect=loader):
            result = load_local_checkpoint(Path("/trusted/checkpoint.pt"))

        self.assertEqual(state["allowed"], [argparse.Namespace])
        self.assertEqual(result, (marker, "alphabet"))
        self.assertFalse(state["inside_safe_globals"])


if __name__ == "__main__":
    unittest.main()
