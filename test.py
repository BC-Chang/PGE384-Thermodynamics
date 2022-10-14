# Import unit test packages
import unittest
import numpy as np

# Import test functions
import solve
import eos

class TestSolve(unittest.TestCase):

    def setUp(self):
        # Define test problem
        self.a = 1.
        self.b = -1.0333152016600968
        self.c = 0.06892495528864631
        self.d = -0.0022962487848544105


    def test_get_delta(self):

        _, _, delta = solve._get_delta(self.a, self.b, self.c, self.d)
        np.testing.assert_almost_equal(delta, 3.3053953091446053e-05)

        return

    def test_solve_cardanos(self):
        x = solve.solve_cardanos(self.a, self.b, self.c, self.d)
        np.testing.assert_array_almost_equal(x,
                                             np.array([0.964309,
                                                       0.034503 + 0.034507j,
                                                       0.034503 - 0.034507j]))

    def test_negative_delta(self):

        self.b = -6.
        self.c = 11.
        self.d = -6.
        D, E, delta = solve._get_delta(self.a, self.b, self.c, self.d)
        print(delta)
        x = solve._delta_negative(self.b, D, E)

        np.testing.assert_array_almost_equal(x, np.array([3.000000, 1.000000, 2.000000]))

        return

    def test_positive_delta(self):
        self.b = 2.
        self.c = 3.
        self.d = 4.
        D, _, delta = solve._get_delta(self.a, self.b, self.c, self.d)
        x = solve._delta_positive(self.b, D, delta)

        np.testing.assert_array_almost_equal(x,
                                             np.array([-1.6506291914393882,
                                                       -0.17468540428030588+1.5468688872313963j,
                                                       -0.17468540428030588-1.5468688872313963j]))

        return

    def test_delta_0(self):
        pass


if __name__ == '__main__':
    unittest.main()


