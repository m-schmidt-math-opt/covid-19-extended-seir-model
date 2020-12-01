from numpy import array as vector
from numpy import matrix as matrix
import numpy as np
import math
import csv

class Data_Reader:

    def read_from_csv_file(self, filename, print_data = False):
        print("Parse data from CSV file " + filename + " ...")

        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=";")
            nr_of_given_tracing_data = 0
            beta_asym_rows = []
            beta_sym_rows = []
            beta_sev_rows = []
            for row in csv_reader:
                if row[0] == "\ufeffN_total": # todo What the heck?!
                    N_total = int(float(row[1].replace(",", ".")))
                    assert(N_total > 0)
                elif row[0] == "K":
                    K = int(float(row[1].replace(",", ".")))
                    assert(K > 0)
                elif row[0] == "Group_labels":
                    pass # todo maybe use later for plots
                elif row[0] == "N":
                    N = vector(self._to_float(row[1:1+K]))
                    assert(sum(N) == N_total)
                    assert(len(N) == K)
                elif row[0] == "beta_asym":
                    beta_asym_rows.append(vector(self._to_float(row[1:1+K])))
                elif row[0] == "beta_sym":
                    beta_sym_rows.append(vector(self._to_float(row[1:1+K])))
                elif row[0] == "beta_sev":
                    beta_sev_rows.append(vector(self._to_float(row[1:1+K])))
                elif row[0] == "epsilon":
                    epsilon = vector(self._to_float(row[1:1+K]))
                    assert(len(epsilon) == K)
                elif row[0] == "eta":
                    eta = vector(self._to_float(row[1:1+K]))
                    assert(len(eta) == K)
                elif row[0] == "nu":
                    nu = vector(self._to_float(row[1:1+K]))
                    assert(len(nu) == K)
                elif row[0] == "sigma":
                    sigma = vector(self._to_float(row[1:1+K]))
                    assert(len(sigma) == K)
                elif row[0] == "psi":
                    psi = vector(self._to_float(row[1:1+K]))
                    assert(len(psi) == K)
                    nr_of_given_tracing_data += 1
                elif row[0] == "gamma_asym":
                    gamma_asym = float(row[1].replace(",","."))
                elif row[0] == "gamma_sym":
                    gamma_sym = vector(self._to_float(row[1:1+K]))
                    assert(len(gamma_sym) == K)
                elif row[0] == "gamma_sev_d_hat":
                    gamma_sev_d_hat = vector(self._to_float(row[1:1+K]))
                    assert(len(gamma_sev_d_hat) == K)
                elif row[0] == "gamma_sev_r_hat":
                    gamma_sev_r_hat = vector(self._to_float(row[1:1+K]))
                    assert(len(gamma_sev_r_hat) == K)
                elif row[0] == "beds":
                    beds = float(row[1].replace(",","."))
                    assert(beds > 0)
                elif row[0] == "S":
                    S = vector(self._to_float(row[1:1+K]))
                    assert(len(S) == K)
                elif row[0] == "E":
                    E = vector(self._to_float(row[1:1+K]))
                    assert(len(E) == K)
                elif row[0] == "E_T":
                    E_tracked = vector(self._to_float(row[1:1+K]))
                    assert(len(E_tracked) == K)
                    nr_of_given_tracing_data += 1
                elif row[0] == "I_asym":
                    I_asym = vector(self._to_float(row[1:1+K]))
                    assert(len(I_asym) == K)
                elif row[0] == "I_sym":
                    I_sym = vector(self._to_float(row[1:1+K]))
                    assert(len(I_sym) == K)
                elif row[0] == "I_sev":
                    I_sev = vector(self._to_float(row[1:1+K]))
                    assert(len(I_sev) == K)
                elif row[0] == "Q_asym":
                    Q_asym = vector(self._to_float(row[1:1+K]))
                    assert(len(Q_asym) == K)
                    nr_of_given_tracing_data += 1
                elif row[0] == "Q_sym":
                    Q_sym = vector(self._to_float(row[1:1+K]))
                    assert(len(Q_sym) == K)
                    nr_of_given_tracing_data += 1
                elif row[0] == "Q_sev":
                    Q_sev = vector(self._to_float(row[1:1+K]))
                    assert(len(Q_sev) == K)
                    nr_of_given_tracing_data += 1
                elif row[0] == "R":
                    R = vector(self._to_float(row[1:1+K]))
                    assert(len(R) == K)
                elif row[0] == "D":
                    D = vector(self._to_float(row[1:1+K]))
                    assert(len(D) == K)
                else:
                    print("Unknown data field: " + str(row[0]))
                    assert(False)

        beta_asym = self._parse_symmetric_matrix(beta_asym_rows, K)
        beta_sym = self._parse_symmetric_matrix(beta_sym_rows, K)
        beta_sev = self._parse_symmetric_matrix(beta_sev_rows, K)

        numerical_tolerance = 1e-12
        tracing_data_given = (nr_of_given_tracing_data >= 1)
        if tracing_data_given:
            assert(nr_of_given_tracing_data == 5)
            x0_total = np.concatenate((S, E, E_tracked, I_asym, I_sym, I_sev, Q_asym, Q_sym, Q_sev, R, D))
            assert(abs(sum(x0_total) - N_total) < numerical_tolerance)
            x0_share = vector([x / N_total for x in x0_total])
            assert(x0_total.shape[0] == 11 * K)
            assert(x0_share.shape[0] == 11 * K)
            if print_data:
                self._print_all_data(tracing_data_given,
                                     N_total,
                                     K,
                                     N,
                                     beta_asym,
                                     beta_sym,
                                     beta_sev,
                                     epsilon,
                                     eta,
                                     nu,
                                     sigma,
                                     gamma_asym,
                                     gamma_sym,
                                     gamma_sev_d_hat,
                                     gamma_sev_r_hat,
                                     psi,
                                     beds,
                                     x0_total)
            packed_data = [tracing_data_given,
                           K,
                           N,
                           beta_asym,
                           beta_sym,
                           beta_sev,
                           epsilon,
                           eta,
                           nu,
                           sigma,
                           gamma_asym,
                           gamma_sym,
                           gamma_sev_d_hat,
                           gamma_sev_r_hat,
                           psi,
                           beds,
                           x0_total]
        else:
            assert(nr_of_given_tracing_data == 0)
            x0_total = np.concatenate((S, E, I_asym, I_sym, I_sev, R, D))
            x0_share = vector([x / N_total for x in x0_total])
            assert(abs(sum(x0_total) - N_total) < numerical_tolerance)
            assert(x0_total.shape[0] == 7 * K)
            assert(x0_share.shape[0] == 7 * K)
            packed_data = [tracing_data_given,
                           K,
                           N,
                           beta_asym,
                           beta_sym,
                           beta_sev,
                           epsilon,
                           eta,
                           nu,
                           sigma,
                           gamma_asym,
                           gamma_sym,
                           gamma_sev_d_hat,
                           gamma_sev_r_hat,
                           beds,
                           x0_total]

        return packed_data

    def _to_float(self, my_list):
        return [float(x.replace(",",".")) if x is not '' else float('nan') for x in my_list]

    def _parse_symmetric_matrix(self, matrix_data, n):
        sym_matrix = matrix(matrix_data)
        # replace nan entries at [i,j] with entries at [j,i]
        for i in range(n):
            for j in range(n):
                if math.isnan(sym_matrix[i,j]):
                    # check only lower diagonal part is nan
                    assert(i > j)
                    assert(not math.isnan(sym_matrix[j,i]))
                    sym_matrix[i,j] = sym_matrix[j,i]
        return sym_matrix

    def _print_all_data(self,
                        tracing_data_given,
                        N_total,
                        K,
                        N,
                        beta_asym,
                        beta_sym,
                        beta_sev,
                        epsilon,
                        eta,
                        nu,
                        sigma,
                        gamma_asym,
                        gamma_sym,
                        gamma_sev_d_hat,
                        gamma_sev_r_hat,
                        psi,
                        beds,
                        x0_total):
        print("tracing_data_given:")
        print(tracing_data_given)
        print("N_total:")
        print(N_total)
        print("K:")
        print(K)
        print("N:")
        print(N)
        print("beta_asym:")
        print(beta_asym)
        print("beta_sym:")
        print(beta_sym)
        print("beta_sev:")
        print(beta_sev)
        print("epsilon:")
        print(epsilon)
        print("eta:")
        print(eta)
        print("nu:")
        print(nu)
        print("sigma:")
        print(sigma)
        print("gamma_asym:")
        print(gamma_asym)
        print("gamma_sym:")
        print(gamma_sym)
        print("gamma_sev_d_hat:")
        print(gamma_sev_d_hat)
        print("gamma_sev_r_hat:")
        print(gamma_sev_r_hat)
        print("psi:")
        print(psi)
        print("beds:")
        print(beds)
        print("x0_total:")
        print(x0_total)

# class Data_Reader
