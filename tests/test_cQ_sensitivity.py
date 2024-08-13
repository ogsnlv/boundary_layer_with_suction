import numpy as np
import os
import unittest
import data
from src import constants, solution
from matplotlib import pyplot as plt


class TestcQSensitivity(unittest.case.TestCase):
    def test_cQ_sensitivity(self, K):
