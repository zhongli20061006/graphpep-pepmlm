import math
import unittest

from esm_pll import masked_sequences, score_sequence, sum_masked_log_probs


class ESMPLLHelpersTest(unittest.TestCase):
    def test_masks_each_residue_once_without_changing_length(self):
        self.assertEqual(masked_sequences("ACD"), ("<mask>CD", "A<mask>D", "AC<mask>"))

    def test_sums_true_token_log_probabilities_across_positions(self):
        logits = [[0.0, 0.0], [math.log(3.0), 0.0]]
        score = sum_masked_log_probs(logits, [0, 0])
        self.assertAlmostEqual(score, math.log(0.5) + math.log(0.75), places=12)

    def test_rejects_nonstandard_residues(self):
        with self.assertRaises(ValueError):
            masked_sequences("ACX")

    def test_rejects_mismatched_logits_and_targets(self):
        with self.assertRaises(ValueError):
            sum_masked_log_probs([[0.0, 1.0]], [0, 1])

    def test_scores_logits_at_the_masked_residue_position(self):
        import torch

        class Alphabet:
            indices = {"A": 1, "C": 2, "<mask>": 3}

            def get_idx(self, token):
                return self.indices[token]

            def get_batch_converter(self):
                def convert(items):
                    rows = []
                    for _, sequence in items:
                        tokens = [0]
                        while sequence:
                            if sequence.startswith("<mask>"):
                                tokens.append(self.indices["<mask>"])
                                sequence = sequence[6:]
                            else:
                                tokens.append(self.indices[sequence[0]])
                                sequence = sequence[1:]
                        rows.append(tokens + [0])
                    return None, None, torch.tensor(rows)
                return convert

        class Model:
            def __call__(self, tokens):
                logits = torch.zeros((*tokens.shape, 4))
                for batch_index in range(tokens.shape[0]):
                    masked_position = int((tokens[batch_index] == 3).nonzero()[0, 0])
                    target = 1 if masked_position == 1 else 2
                    logits[batch_index, masked_position, target] = math.log(3.0)
                return {"logits": logits}

        values = score_sequence(Model(), Alphabet(), "AC", "cpu")
        self.assertEqual(len(values), 2)
        self.assertAlmostEqual(values[0], math.log(0.5), places=6)
        self.assertAlmostEqual(values[1], math.log(0.5), places=6)


if __name__ == "__main__":
    unittest.main()
