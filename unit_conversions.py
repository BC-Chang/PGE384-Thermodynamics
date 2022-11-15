import numpy as np

class Unit_Converter:
    def __init__(self):
        return
    def make_val_list(self, val):
        if not isinstance(val, np.ndarray):
            val = np.asarray(val)

        return val

    def Pa_to_psi(self, val):
        val = self.make_val_list(val)
        return val * 1.450377E-4

    def psi_to_Pa(self, val):
        val = self.make_val_list(val)
        return val / 1.450377E-4

    def Pa_to_bar(self, val):
        val = self.make_val_list(val)
        return val * 1E-5

    def bar_to_Pa(self, val):
        val = self.make_val_list(val)
        return val / 1E-5

    def C_to_K(self, val):
        val = self.make_val_list(val)
        return val + 273.15

    def K_to_C(self, val):
        val = self.make_val_list(val)
        return val - 273.15

    def C_to_F(self, val):
        val = self.make_val_list(val)
        return 1.80 * val + 32.0

    def F_to_C(self, val):
        val = self.make_val_list(val)
        return (val - 32) * 5.0 / 9.0

    def F_to_K(self, val):
        val = self.make_val_list(val)
        return (val - 32) / 1.8 + 273.15

    def K_to_F(self, val):
        val = self.make_val_list(val)
        return (val - 273.15) * (9.0 / 5.0) + 32







