class Result_Analyzer:

    def __init__(self, tracing_model_solved, results, stepsize, N, K, beds, t_start, t_end, filename):
        # Store given values
        self.tracing_model_solved = tracing_model_solved
        self.results = results
        self.stepsize = stepsize
        self.N = N
        self.N_total = sum(N)
        self.K = K
        self.beds = beds
        self.t_start = t_start
        self.t_end = t_end
        self.filename = filename

        # Sanity check
        assert(self.t_start < self.t_end)

        # Preparing empty lists (of lists)
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

    def get_extracted_results_as_dict(self):
        result_dict = {}
        result_dict["t_vals"] = self.t_vals
        result_dict["S_vals"] = self.S_vals
        result_dict["E_vals"] = self.E_vals
        result_dict["I_asym_vals"] = self.I_asym_vals
        result_dict["I_sym_vals"] = self.I_sym_vals
        result_dict["I_sev_vals"] = self.I_sev_vals
        result_dict["R_vals"] = self.R_vals
        result_dict["D_vals"] = self.D_vals
        result_dict["sev_total"] = self.sev_total
        if self.tracing_model_solved:
            result_dict["E_tracked_vals"] = self.E_tracked_vals
            result_dict["Q_asym_vals"] = self.Q_asym_vals
            result_dict["Q_sym_vals"] = self.Q_sym_vals
            result_dict["Q_sev_vals"] = self.I_sev_vals
        return result_dict

    def extract_results(self):
        # Extract and separate results
        for entry in self.results:
            assert(len(entry) == 2) # entry is of type [t, x]
            self.t_vals.append(entry[0])
            x = entry[1]
            if self.tracing_model_solved:
                assert(len(x) == 11 * self.K)
            else:
                assert(len(x) == 7 * self.K)
            idx = 0
            for k in range(self.K):
                assert(x[idx] >= 0.0)
                assert(x[idx] <= self.N_total)
                self.S_vals[k].append(x[idx])
                idx += 1
            for k in range(self.K):
                assert(x[idx] >= 0.0)
                assert(x[idx] <= self.N_total)
                self.E_vals[k].append(x[idx])
                idx += 1
            if self.tracing_model_solved:
                for k in range(self.K):
                    assert(x[idx] >= 0.0)
                    assert(x[idx] <= self.N_total)
                    self.E_tracked_vals[k].append(x[idx])
                    idx += 1
            for k in range(self.K):
                assert(x[idx] >= 0.0)
                assert(x[idx] <= self.N_total)
                self.I_asym_vals[k].append(x[idx])
                idx += 1
            for k in range(self.K):
                assert(x[idx] >= 0.0)
                assert(x[idx] <= self.N_total)
                self.I_sym_vals[k].append(x[idx])
                idx += 1
            for k in range(self.K):
                assert(x[idx] >= 0.0)
                assert(x[idx] <= self.N_total)
                self.I_sev_vals[k].append(x[idx])
                idx += 1
            if self.tracing_model_solved:
                for k in range(self.K):
                    assert(x[idx] >= 0.0)
                    assert(x[idx] <= self.N_total)
                    self.Q_asym_vals[k].append(x[idx])
                    idx += 1
                for k in range(self.K):
                    assert(x[idx] >= 0.0)
                    assert(x[idx] <= self.N_total)
                    self.Q_sym_vals[k].append(x[idx])
                    idx += 1
                for k in range(self.K):
                    assert(x[idx] >= 0.0)
                    assert(x[idx] <= self.N_total)
                    self.Q_sev_vals[k].append(x[idx])
                    idx += 1
            for k in range(self.K):
                self.R_vals[k].append(x[idx])
                idx += 1
            for k in range(self.K):
                self.D_vals[k].append(x[idx])
                idx += 1

            # Sanity checks
            if self.tracing_model_solved:
                assert(idx == 11 * self.K)
            else:
                assert(idx == 7 * self.K)

        # Sanity checks
        for k in range(self.K):
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

        self.sev_total = self._compute_total_severe_cases()

    def compute_and_write_results(self):
        self.nr_of_deaths = self._compute_total_nr_deaths()
        self.days_w_exceeded_beds = self._compute_exceeded_icu_bed_days()
        self.duration_of_pandemic = self._compute_duration_of_pandemic()
        self._write_result_file()

    def _compute_total_nr_deaths(self):
        # Compute total number of deaths
        nr_of_deaths = 0.0
        for k in range(self.K):
            # Sanity checks for the final death values
            assert(self.D_vals[k][-1] >= 0.0)
            assert(self.D_vals[k][-1] <= self.N_total)
            nr_of_deaths += self.D_vals[k][-1]
        return(nr_of_deaths)

    def _compute_total_severe_cases(self):
        sev_total = []
        for index, t in enumerate(self.t_vals):
            sev_total_at_t = 0.0
            for k in range(self.K):
                assert(len(self.t_vals) == len(self.I_sev_vals[k]))
                assert(self.I_sev_vals[k][index] >= 0.0)
                assert(self.I_sev_vals[k][index] <= self.N_total)
                sev_total_at_t += self.I_sev_vals[k][index]
                if self.tracing_model_solved:
                    assert(len(self.t_vals) == len(self.Q_sev_vals[k]))
                    assert(self.Q_sev_vals[k][index] >= 0.0)
                    assert(self.Q_sev_vals[k][index] <= self.N_total)
                    sev_total_at_t += self.Q_sev_vals[k][index]
            sev_total.append(sev_total_at_t)
        assert(len(sev_total) == len(self.t_vals))
        return sev_total

    def _compute_exceeded_icu_bed_days(self):
        # Compute number of days with exceeded ICU bed capacity
        t_w_exceeded_beds = [val for val in self.sev_total if val > self.beds]
        days_w_exceeded_beds = len(t_w_exceeded_beds) * self.stepsize
        return days_w_exceeded_beds

    def _compute_duration_of_pandemic(self):
        # Compute duration of pandemic
        threshold = 0.05 * self.beds
        time_points_below_threshold = 0
        sev_below_at_end = (self.sev_total[-1] < threshold)
        assert(sev_below_at_end)
        for val in reversed(self.sev_total):
            if val < threshold:
                time_points_below_threshold += 1
            else:
                break
        number_of_days = self.t_end - self.t_start
        days_below_threshold = time_points_below_threshold * self.stepsize
        duration_of_pandemic = number_of_days - days_below_threshold
        return duration_of_pandemic

    def _write_result_file(self):
        result_file = open(self.filename + ".res", "w")
        result_file.write("Total number of deaths:              " + str(self.nr_of_deaths) + "\n")
        result_file.write("Days with exceeded ICU bed capacity: " + str(self.days_w_exceeded_beds) + "\n")
        result_file.write("Duration of pandemic (in days):      " + str(self.duration_of_pandemic) + "\n")
        result_file.close()

# class Result_Analyzer
