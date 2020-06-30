from numpy import array as vector
import numpy as np

class SEIIIRD_Model:

    def __init__(self, packed_parameters):
        print("Instantiate the SEIIIRD model ...")
        self.K = packed_parameters[0]
        self.N = packed_parameters[1]
        self.beta_asym = packed_parameters[2]
        self.beta_sym = packed_parameters[3]
        self.beta_sev = packed_parameters[4]
        self.epsilon = packed_parameters[5]
        self.eta = packed_parameters[6]
        self.nu = packed_parameters[7]
        self.sigma = packed_parameters[8]
        self.gamma_asym = packed_parameters[9]
        self.gamma_sym = packed_parameters[10]
        self.gamma_sev_d = packed_parameters[11]
        self.gamma_sev_r = packed_parameters[12]
        self.beds = packed_parameters[13]

        self.numerical_tolerance = 1e-12

    def eval_rhs(self, x):
        S = x[0 : self.K]
        E = x[self.K : 2*self.K]
        I_asym = x[2*self.K : 3*self.K]
        I_sym = x[3*self.K : 4*self.K]
        I_sev = x[4*self.K : 5*self.K]
        R = x[5*self.K : 6*self.K]
        D = x[6*self.K : 7*self.K]

        # Sanity checks
        assert(x.shape[0] == 7 * self.K)
        assert(abs(sum(x) - 1) < self.numerical_tolerance)
        assert(len(S) == self.K)
        assert(len(E) == self.K)
        assert(len(I_asym) == self.K)
        assert(len(I_sym) == self.K)
        assert(len(I_sev) == self.K)
        assert(len(R) == self.K)
        assert(len(D) == self.K)

        f_vec = np.empty(7*self.K)
        idx = 0
        overall_rhs = 0

        # ODE rhs for S
        for k in range(self.K):
            factor = 0.0
            for l in range(self.K):
                factor += self.beta_asym[l,k] * I_asym[l] + self.beta_sym[l,k] * I_sym[l] + self.beta_sev[l,k] * I_sev[l]
            assert(factor >= 0.0)
            f_vec[idx] = -factor * S[k]
            overall_rhs += f_vec[idx]
            idx += 1

        # ODE rhs for E
        for k in range(self.K):
            factor = 0.0
            for l in range(self.K):
                factor += self.beta_asym[l,k] * I_asym[l] + self.beta_sym[l,k] * I_sym[l] + self.beta_sev[l,k] * I_sev[l]
            assert(factor >= 0.0)
            f_vec[idx] = factor * S[k] - self.epsilon[k] * E[k]
            overall_rhs += f_vec[idx]
            idx += 1

        # ODE rhs for I_asym
        for k in range (self.K):
            f_vec[idx] = self.eta[k] * self.epsilon[k] * E[k] - self.gamma_asym * I_asym[k]
            overall_rhs += f_vec[idx]
            idx += 1

        # ODE rhs for I_sym
        for k in range(self.K):
            f_vec[idx] = (1 - self.eta[k]) * (1 - self.nu[k]) * self.epsilon[k] * E[k] - self.gamma_sym[k] * I_sym[k]
            overall_rhs += f_vec[idx]
            idx += 1

        # ODE rhs for I_sev
        for k in range(self.K):
            sigma_k = self._get_sigma_k(I_sev[k], k)
            f_vec[idx] = (1 - self.eta[k]) * self.nu[k] * self.epsilon[k] * E[k] - ((1.0 - sigma_k) * self.gamma_sev_r[k] + sigma_k * self.gamma_sev_d[k]) * I_sev[k]
            overall_rhs += f_vec[idx]
            idx += 1

        # ODE rhs for R
        for k in range(self.K):
            sigma_k = self._get_sigma_k(I_sev[k], k)
            f_vec[idx] = self.gamma_asym * I_asym[k] + self.gamma_sym[k] * I_sym[k] + (1.0 - sigma_k) * self.gamma_sev_r[k] * I_sev[k]
            overall_rhs += f_vec[idx]
            idx += 1

        # ODE rhs for D
        for k in range(self.K):
            sigma_k = self._get_sigma_k(I_sev[k], k)
            f_vec[idx] = sigma_k * self.gamma_sev_d[k] * I_sev[k]
            overall_rhs += f_vec[idx]
            idx += 1

        # Sanity checks
        assert(idx == 7 * self.K)
        assert(abs(overall_rhs) < self.numerical_tolerance)

        return f_vec

    def _get_sigma_k(self, I_sev_k, k):
        assert(len(self.N) == self.K)
        N_total = sum(self.N)
        if (I_sev_k * self.N[k] <= (self.N[k] / N_total) * self.beds):
            return self.sigma[k]
        else:
            sigma_k = (self.sigma[k] * (self.N[k] / N_total) * self.beds + I_sev_k * self.N[k] - (self.N[k] / N_total) * self.beds) / (I_sev_k * self.N[k])
            assert(sigma_k >= 0.0)
            assert(sigma_k <= 1.0)
            return sigma_k

# class SEIIIRD_Model
