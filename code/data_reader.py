from numpy import array as vector
from numpy import matrix as matrix
import numpy as np
import csv

class Data_Reader:

    def read_from_csv_file(self, filename):
        print("Parse data from CSV file " + filename + " ...")

        with open(filename, encoding='utf-8-sig') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=";")
            nr_of_given_tracing_data = 0
            for row in csv_reader:
                if row[0] == "N_total":
                    N_total = int(float(row[1].replace(",", ".")))
                    assert(N_total > 0)
                elif row[0] == "K":
                    K = int(float(row[1].replace(",", ".")))
                    assert(K > 0)
                elif row[0] == "N":
                    N = vector(self._to_float(row[1:1+K]))
                    assert(sum(N) == N_total)
                    assert(len(N) == K)
                elif row[0] == "beta_asym":
                    beta_asym = matrix(vector(self._to_float(row[1:1+K*K]))).reshape((K, K))
                elif row[0] == "beta_sym":
                    beta_sym = matrix(vector(self._to_float(row[1:1+K*K]))).reshape((K, K))
                elif row[0] == "beta_sev":
                    beta_sev = matrix(vector(self._to_float(row[1:1+K*K]))).reshape((K, K))
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
                elif row[0] == "psi":
                    psi = vector(self._to_float(row[1:1+K]))
                    assert(len(psi) == K)
                    nr_of_given_tracing_data += 1
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

        tracing_data_given = (nr_of_given_tracing_data >= 1)
        if tracing_data_given:
            assert(nr_of_given_tracing_data == 5)
            x0_total = np.concatenate((S, E, E_tracked, I_asym, I_sym, I_sev, Q_asym, Q_sym, Q_sev, R, D))
            assert(sum(S) + sum(E) + sum(E_tracked) + sum(I_asym) + sum(I_sym) + sum(I_sev) + sum(Q_asym) + sum(Q_sym) + sum(Q_sev) + sum(R) + sum(D) == N_total)
            x0_share = vector([x / N_total for x in x0_total])
            assert(x0_total.shape[0] == 11 * K)
            assert(x0_share.shape[0] == 11 * K)
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
            assert(sum(S) + sum(E) + sum(I_asym) + sum(I_sym) + sum(I_sev) + sum(R) + sum(D) == N_total)
            x0_share = vector([x / N_total for x in x0_total])
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
        return [float(x.replace(",",".")) for x in my_list]

# class Data_Reader
