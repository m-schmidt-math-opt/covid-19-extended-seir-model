class Result_Analyzer:

    def __init__(self, tracing_model_solved, results, stepsize, N, K, beds, t_start, t_end, filename):
        # Store given values
        self.tracing_model_solved = tracing_model_solved
        self.stepsize = stepsize
        self.N = N
        self.K = K
        self.beds = beds
        self.t_start = t_start
        self.t_end = t_end
        self.filename = filename

        # Sanity checks
        assert(self.t_start < self.t_end)

        # Preparing lists
        self.t_vals = []
        self.S_vals = [[] for k in range(K)]
        self.E_vals = [[] for k in range(K)]
        self.I_asym_vals = [[] for k in range(K)]
        self.I_sym_vals = [[] for k in range(K)]
        self.I_sev_vals = [[] for k in range(K)]
        self.R_vals = [[] for k in range(K)]
        self.D_vals = [[] for k in range(K)]

        # Preparing lists that are only used if a tracing model was solved
        if self.tracing_model_solved:
            self.E_tracked_vals = [[] for k in range(K)]
            self.Q_asym_vals = [[] for k in range(K)]
            self.Q_sym_vals = [[] for k in range(K)]
            self.Q_sev_vals = [[] for k in range(K)]

        # Extract and separate results
        for entry in results:
            assert(len(entry) == 2) # entry is of type [t, x]
            self.t_vals.append(entry[0])
            x = entry[1]
            if self.tracing_model_solved:
                assert(len(x) == 11 * K)
            else:
                assert(len(x) == 7 * K)
            idx = 0
            for k in range(K):
                assert(x[idx] >= 0.0)
                assert(x[idx] <= 1.0)
                self.S_vals[k].append(x[idx])
                idx += 1
            for k in range(K):
                assert(x[idx] >= 0.0)
                assert(x[idx] <= 1.0)
                self.E_vals[k].append(x[idx])
                idx += 1
            if self.tracing_model_solved:
                for k in range(K):
                    assert(x[idx] >= 0.0)
                    assert(x[idx] <= 1.0)
                    self.E_tracked_vals[k].append(x[idx])
                    idx += 1
            for k in range(K):
                assert(x[idx] >= 0.0)
                assert(x[idx] <= 1.0)
                self.I_asym_vals[k].append(x[idx])
                idx += 1
            for k in range(K):
                assert(x[idx] >= 0.0)
                assert(x[idx] <= 1.0)
                self.I_sym_vals[k].append(x[idx])
                idx += 1
            for k in range(K):
                assert(x[idx] >= 0.0)
                assert(x[idx] <= 1.0)
                self.I_sev_vals[k].append(x[idx])
                idx += 1
            if self.tracing_model_solved:
                for k in range(K):
                    assert(x[idx] >= 0.0)
                    assert(x[idx] <= 1.0)
                    self.Q_asym_vals[k].append(x[idx])
                    idx += 1
                for k in range(K):
                    assert(x[idx] >= 0.0)
                    assert(x[idx] <= 1.0)
                    self.Q_sym_vals[k].append(x[idx])
                    idx += 1
                for k in range(K):
                    assert(x[idx] >= 0.0)
                    assert(x[idx] <= 1.0)
                    self.Q_sev_vals[k].append(x[idx])
                    idx += 1
            for k in range(K):
                self.R_vals[k].append(x[idx])
                idx += 1
            for k in range(K):
                self.D_vals[k].append(x[idx])
                idx += 1

            # Sanity checks
            if self.tracing_model_solved:
                assert(idx == 11 * K)
            else:
                assert(idx == 7 * K)

        # Sanity checks
        for k in range(K):
            assert(len(self.t_vals) == len(self.S_vals[k]))
            assert(len(self.t_vals) == len(self.E_vals[k]))
            assert(len(self.t_vals) == len(self.I_asym_vals[k]))
            assert(len(self.t_vals) == len(self.I_sym_vals[k]))
            assert(len(self.t_vals) == len(self.I_sev_vals[k]))
            assert(len(self.t_vals) == len(self.R_vals[k]))
            assert(len(self.t_vals) == len(self.D_vals[k]))
            if self.tracing_model_solved:
                assert(len(self.t_vals) == len(self.E_tracked_vals[k]))
                assert(len(self.t_vals) == len(self.Q_asym_vals[k]))
                assert(len(self.t_vals) == len(self.Q_sym_vals[k]))
                assert(len(self.t_vals) == len(self.Q_sev_vals[k]))

    def analyze(self):

        # Compute total number of deaths
        nr_of_deaths = 0.0
        for k in range(self.K):
            # Sanity checks for the final death values
            assert(self.D_vals[k][-1] >= 0.0)
            assert(self.D_vals[k][-1] <= 1.0)
            nr_of_deaths += self.D_vals[k][-1]
        N_total = sum(self.N)
        nr_of_deaths *= N_total

        # Compute number of days with exceeded ICU bed capacity
        sev_total = []
        for index, t in enumerate(self.t_vals):
            sev_total_at_t = 0.0
            for k in range(self.K):
                assert(len(self.t_vals) == len(self.I_sev_vals[k]))
                assert(self.I_sev_vals[k][index] >= 0.0)
                assert(self.I_sev_vals[k][index] <= 1.0)
                sev_total_at_t += self.I_sev_vals[k][index]
                if self.tracing_model_solved:
                    assert(len(self.t_vals) == len(self.Q_sev_vals[k]))
                    assert(self.Q_sev_vals[k][index] >= 0.0)
                    assert(self.Q_sev_vals[k][index] <= 1.0)
                    sev_total_at_t += self.Q_sev_vals[k][index]
            sev_total.append(sev_total_at_t * N_total)
        assert(len(sev_total) == len(self.t_vals))
        t_w_exceeded_beds = [val for val in sev_total if val > self.beds]
        days_w_exceeded_beds = len(t_w_exceeded_beds) * self.stepsize

        # Compute duration of pandemic
        threshold = 0.05 * self.beds
        time_points_below_threshold = 0
        sev_below_at_end = (sev_total[-1] < threshold)
        assert(sev_below_at_end)
        for val in reversed(sev_total):
            if val < threshold:
                time_points_below_threshold += 1
            else:
                break
        number_of_days = self.t_end - self.t_start
        days_below_threshold = time_points_below_threshold * self.stepsize
        duration_of_pandemic = number_of_days - days_below_threshold

        result_file = open(self.filename + ".res", "w")
        result_file.write("Total number of deaths:              " + str(nr_of_deaths) + "\n")
        result_file.write("Days with exceeded ICU bed capacity: " + str(days_w_exceeded_beds) + "\n")
        result_file.write("Duration of pandemic (in days):      " + str(duration_of_pandemic) + "\n")
        result_file.close()

# class Result_Analyzer
