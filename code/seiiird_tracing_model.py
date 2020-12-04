from numpy import array as vector
import numpy as np

class SEIIIRD_Tracing_Model:

    def __init__(self, packed_parameters):
        print("Instantiate the SEIIIRD tracing model ...")
        self.K = packed_parameters[0]
        self.N = packed_parameters[1]
        self.N_total = sum(self.N)
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
        self.psi = packed_parameters[13]
        self.beds = packed_parameters[14]

        self.numerical_tolerance_fine = 1e-4 # todo
        self.numerical_tolerance_coarse = 1e-4

    def eval_rhs(self, x):
        S = x[0 : self.K]
        E = x[self.K : 2*self.K]
        E_tracked = x[2*self.K : 3*self.K]
        I_asym = x[3*self.K : 4*self.K]
        I_sym = x[4*self.K : 5*self.K]
        I_sev = x[5*self.K : 6*self.K]
        Q_asym = x[6*self.K : 7*self.K]
        Q_sym = x[7*self.K : 8*self.K]
        Q_sev = x[8*self.K : 9*self.K]
        R = x[9*self.K : 10*self.K]
        D = x[10*self.K : 11*self.K]

        # Sanity checks
        assert(x.shape[0] == 11 * self.K)
        assert(abs(sum(x) - self.N_total) < self.numerical_tolerance_coarse)
        assert(len(S) == self.K)
        assert(len(E) == self.K)
        assert(len(E_tracked) == self.K)
        assert(len(I_asym) == self.K)
        assert(len(I_sym) == self.K)
        assert(len(I_sev) == self.K)
        assert(len(Q_asym) == self.K)
        assert(len(Q_sym) == self.K)
        assert(len(Q_sev) == self.K)
        assert(len(R) == self.K)
        assert(len(D) == self.K)

        for k in range(self.K):
            assert(S[k] >= 0.0)
            assert(S[k] <= self.N_total)
            assert(E[k] >= 0.0)
            assert(E[k] <= self.N_total)
            assert(E_tracked[k] >= 0.0)
            assert(E_tracked[k] <= self.N_total)
            assert(I_asym[k] >= 0.0)
            assert(I_asym[k] <= self.N_total)
            assert(I_sym[k] >= 0.0)
            assert(I_sym[k] <= self.N_total)
            assert(I_sev[k] >= 0.0)
            assert(I_sev[k] <= self.N_total)
            assert(Q_asym[k] >= 0.0)
            assert(Q_asym[k] <= self.N_total)
            assert(Q_sym[k] >= 0.0)
            assert(Q_sym[k] <= self.N_total)
            assert(Q_sev[k] >= 0.0)
            assert(Q_sev[k] <= self.N_total)
            assert(R[k] >= 0.0)
            assert(R[k] <= self.N_total)
            assert(D[k] >= 0.0)
            assert(D[k] <= self.N_total)

        f_vec = np.empty(11*self.K)
        idx = 0
        overall_rhs = 0

        # ODE rhs for S
        for k in range(self.K):
            factor = 0.0
            for l in range(self.K):
                factor += self.beta_asym[l,k] * I_asym[l] + self.beta_sym[l,k] * I_sym[l] + self.beta_sev[l,k] * I_sev[l]
            assert(factor >= 0.0)
            f_vec[idx] = -factor * S[k] / self.N_total
            assert(f_vec[idx] <= 0.0)
            overall_rhs += f_vec[idx]
            idx += 1

        # ODE rhs for E
        for k in range(self.K):
            # First term
            factor = 0.0
            for l in range(self.K):
                factor += self.beta_asym[l,k] * I_asym[l]
            assert(factor >= 0.0)
            first_term = factor * S[k] / self.N_total

            # Second term
            factor = 0.0
            for l in range(self.K):
                psi_lk = self.psi[l] * self.psi[k]
                factor += self.beta_sym[l,k] * (1.0 - psi_lk) * I_sym[l] + self.beta_sev[l,k] * (1.0 - psi_lk) * I_sev[l]
            assert(factor >= 0.0)
            second_term = factor * S[k] / self.N_total

            # Third term and final result
            f_vec[idx] = first_term + second_term - self.epsilon[k] * E[k]
            overall_rhs += f_vec[idx]
            idx += 1

        # ODE rhs for E_tracked
        for k in range(self.K):
            factor = 0.0
            for l in range(self.K):
                psi_lk = self.psi[l] * self.psi[k]
                factor += self.beta_sym[l,k] * psi_lk * I_sym[l] + self.beta_sev[l,k] * psi_lk * I_sev[l]
            assert(factor >= 0.0)
            f_vec[idx] = factor * S[k] / self.N_total - self.epsilon[k] * E_tracked[k]
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
            sigma_k = self._get_sigma_k(I_sev[k], Q_sev[k], k)
            f_vec[idx] = (1 - self.eta[k]) * self.nu[k] * self.epsilon[k] * E[k] - ((1.0 - sigma_k) * self.gamma_sev_r[k] + sigma_k * self.gamma_sev_d[k]) * I_sev[k]
            overall_rhs += f_vec[idx]
            idx += 1

        # ODE rhs for Q_asym
        for k in range(self.K):
            f_vec[idx] = self.eta[k] * self.epsilon[k] * E_tracked[k] - self.gamma_asym * Q_asym[k]
            overall_rhs += f_vec[idx]
            idx += 1

        # ODE rhs for Q_sym
        for k in range(self.K):
            f_vec[idx] = (1.0 - self.eta[k]) * (1.0 - self.nu[k]) * self.epsilon[k] * E_tracked[k] - self.gamma_sym[k] * Q_sym[k]
            overall_rhs += f_vec[idx]
            idx += 1

        # ODE rhs for Q_sev
        for k in range(self.K):
            sigma_k = self._get_sigma_k(I_sev[k], Q_sev[k], k)
            f_vec[idx] = (1.0 - self.eta[k]) * self.nu[k] * self.epsilon[k] * E_tracked[k] - ((1.0 - sigma_k) * self.gamma_sev_r[k] + sigma_k * self.gamma_sev_d[k]) * Q_sev[k]
            overall_rhs += f_vec[idx]
            idx += 1

        # ODE rhs for R
        for k in range(self.K):
            sigma_k = self._get_sigma_k(I_sev[k], Q_sev[k], k)
            I_terms = self.gamma_asym * I_asym[k] + self.gamma_sym[k] * I_sym[k] + (1.0 - sigma_k) * self.gamma_sev_r[k] * I_sev[k]
            Q_terms = self.gamma_asym * Q_asym[k] + self.gamma_sym[k] * Q_sym[k] + (1.0 - sigma_k) * self.gamma_sev_r[k] * Q_sev[k]
            f_vec[idx] = I_terms + Q_terms
            overall_rhs += f_vec[idx]
            idx += 1

        # ODE rhs for D
        for k in range(self.K):
            sigma_k = self._get_sigma_k(I_sev[k], Q_sev[k], k)
            f_vec[idx] = sigma_k * self.gamma_sev_d[k] * I_sev[k] + sigma_k * self.gamma_sev_d[k] * Q_sev[k]
            overall_rhs += f_vec[idx]
            idx += 1

        # Sanity checks
        assert(idx == 11 * self.K)
        assert(abs(overall_rhs) < self.numerical_tolerance_fine)

        return f_vec

    def _get_sigma_k(self, I_sev_k, Q_sev_k, k):
        assert(len(self.N) == self.K)
        beds_for_group_k = (self.N[k] / self.N_total) * self.beds
        if I_sev_k + Q_sev_k <= beds_for_group_k:
            return self.sigma[k]
        else:
            sigma_k = (self.sigma[k] * beds_for_group_k + (I_sev_k + Q_sev_k) - beds_for_group_k) / (I_sev_k + Q_sev_k)
            assert(sigma_k >= 0.0)
            assert(sigma_k <= 1.0)
            return sigma_k

# class SEIIIRD_Tracing_Model
