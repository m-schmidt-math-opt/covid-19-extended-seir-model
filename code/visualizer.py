import matplotlib.pyplot as plot
import numpy as np


class Visualizer:

    def __init__(self, tracing_model_solved, results, N, N_total, K, beds, t_start, t_end, filename):
        # Store given values
        self.tracing_model_solved = tracing_model_solved
        self.N = N
        self.N_total = N_total
        assert(self.N_total == sum(self.N))
        self.K = K
        self.beds = beds
        self.t_start = t_start
        self.t_end = t_end
        self.filename = filename

        # Sanity checks
        assert(self.t_start < self.t_end)

        # Prepare lists
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
                assert(x[idx] >= 0.0)
                assert(x[idx] <= 1.0)
                self.R_vals[k].append(x[idx])
                idx += 1
            for k in range(K):
                assert(x[idx] >= 0.0)
                assert(x[idx] <= 1.0)
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

        # Compute total number of severe cases
        self.sev_total_vals = []
        for index, t in enumerate(self.t_vals):
            sev_total_at_t = 0.0
            for k in range(self.K):
                assert(len(self.t_vals) == len(self.I_sev_vals[k]))
                val = self.I_sev_vals[k][index]
                assert(val >= 0.0)
                assert(val <= 1.0)
                sev_total_at_t += val
                if self.tracing_model_solved:
                    assert(len(self.t_vals) == len(self.Q_sev_vals[k]))
                    val = self.Q_sev_vals[k][index]
                    assert(val >= 0.0)
                    assert(val <= 1.0)
                    sev_total_at_t += val
            self.sev_total_vals.append(sev_total_at_t)

        # Sanity checks
        assert(len(self.sev_total_vals) == len(self.t_vals))
        assert(len(self.t_vals) == len(self.I_sev_vals[k]))



    def plot(self):
        figure, axes = plot.subplots()
        figure.subplots_adjust(bottom = 0.15)
        axes.grid(linestyle = ':', linewidth = 0.5, color = "#808080")

        group_linestyles = {0: "solid", 1: "dotted"} # well, yes, that's hard-coded for K = 2

        for k in range(self.K):
            S_plot, = axes.plot(self.t_vals, self.S_vals[k],
                                color = "lightskyblue", linestyle = group_linestyles[k])
            S_plot.set_label("S_" + str(k))

            E_plot, = axes.plot(self.t_vals, self.E_vals[k],
                                color = "teal", linestyle = group_linestyles[k])
            E_plot.set_label("E_" + str(k))

            if self.tracing_model_solved:
                E_tracked_plot, = axes.plot(self.t_vals, self.E_tracked_vals[k],
                                            color = "deeppink", linestyle = group_linestyles[k])
                E_tracked_plot.set_label("E_tracked_" + str(k))

            I_asym_plot, = axes.plot(self.t_vals, self.I_asym_vals[k],
                                     color = "gold", linestyle = group_linestyles[k])
            I_asym_plot.set_label("I_asym_" + str(k))

            I_sym_plot, = axes.plot(self.t_vals, self.I_sym_vals[k],
                                    color = "orange", linestyle = group_linestyles[k])
            I_sym_plot.set_label("I_sym_" + str(k))

            I_sev_plot, = axes.plot(self.t_vals, self.I_sev_vals[k],
                                    color = "chocolate", linestyle = group_linestyles[k])
            I_sev_plot.set_label("I_sev_" + str(k))

            if self.tracing_model_solved:
                Q_asym_plot, = axes.plot(self.t_vals, self.Q_asym_vals[k],
                                         color = "plum", linestyle = group_linestyles[k])
                Q_asym_plot.set_label("Q_asym_" + str(k))

                Q_sym_plot, = axes.plot(self.t_vals, self.Q_sym_vals[k],
                                         color = "blueviolet", linestyle = group_linestyles[k])
                Q_sym_plot.set_label("Q_sym_" + str(k))

                Q_sev_plot, = axes.plot(self.t_vals, self.Q_sev_vals[k],
                                         color = "navy", linestyle = group_linestyles[k])
                Q_sev_plot.set_label("Q_sev_" + str(k))

            R_plot, = axes.plot(self.t_vals, self.R_vals[k],
                                color = "limegreen", linestyle = group_linestyles[k])
            R_plot.set_label("R_" + str(k))

            D_plot, = axes.plot(self.t_vals, self.D_vals[k],
                                color = "firebrick", linestyle = group_linestyles[k])
            D_plot.set_label("D_" + str(k))

        sev_total_plot, = axes.plot(self.t_vals, self.sev_total_vals,
                                    color = "chocolate", linestyle = "dashed")
        sev_total_plot.set_label("sev_total")

        # Plot horizontal lines for beds
        beds_plot = axes.hlines(float(self.beds / self.N_total), self.t_start, self.t_end,
                                color = "cornflowerblue", linestyle = "dashed")
        beds_plot.set_label("beds")

        # Shrink current axis by 20%
        box = axes.get_position()
        axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plot.savefig(self.filename + "_full.pdf")
        ymax = 0.2
        axes.set_ylim([-ymax * 0.05, ymax])
        plot.savefig(self.filename + "_small.pdf")
        ymax = 0.005
        axes.set_ylim([-ymax * 0.05, ymax])
        plot.savefig(self.filename + "_super_small.pdf")
        plot.close("all")

    def paper_plot_figure_1(self):
        figure, axes = plot.subplots()
        figure.subplots_adjust(bottom = 0.15)
        axes.grid(linestyle = ':', linewidth = 0.5, color = "#808080")

        # 0: old; 1: young
        group_linestyles = {0: "solid", 1: "dotted"} # well, yes, that's hard-coded for K = 2

        # hard-coded for K = 2
        I_0_aggr_vals = [self.N_total * (x + y + z) for x, y, z in zip(self.I_asym_vals[0], self.I_sym_vals[0], self.I_sev_vals[0])]
        I_0_aggr_plot, = axes.plot(self.t_vals, I_0_aggr_vals, color = "tab:orange", linewidth=3.5, linestyle = group_linestyles[0])

        # hard-coded for K = 2
        I_1_aggr_vals = [self.N_total * (x + y + z) for x, y, z in zip(self.I_asym_vals[1], self.I_sym_vals[1], self.I_sev_vals[1])]
        I_1_aggr_plot, = axes.plot(self.t_vals, I_1_aggr_vals, color = "tab:orange", linewidth=3.5, linestyle = group_linestyles[1])

        for k in range(self.K):
            S_vals_k_abs = [self.N_total * self.S_vals[k][t] for t in range(len(self.S_vals[k]))]
            S_plot, = axes.plot(self.t_vals, S_vals_k_abs,
                                color = "tab:green", linewidth=3.5, linestyle = group_linestyles[k])
            #S_plot.set_label("S_" + str(k))

            R_vals_k_abs = [self.N_total * self.R_vals[k][t] for t in range(len(self.R_vals[k]))]
            R_plot, = axes.plot(self.t_vals, R_vals_k_abs,
                                color = "tab:blue", linewidth=3.5, linestyle = group_linestyles[k])
            #R_plot.set_label("R_" + str(k))

        plot.yticks(np.arange(0, 8e7, 1e7), ['0',
                                             r'$10\,$M',
                                             r'$20\,$M',
                                             r'$30\,$M',
                                             r'$40\,$M',
                                             r'$50\,$M',
                                             r'$60\,$M',
                                             r'$70\,$M'])
        plot.savefig(self.filename + "_figure_1.pdf")
        plot.close("all")

    def paper_plot_figure_3(self):
        figure, axes = plot.subplots()
        figure.subplots_adjust(bottom = 0.15)
        axes.grid(linestyle = ':', linewidth = 0.5, color = "#808080")

        # 0: old; 1: young
        group_linestyles = {0: "solid", 1: "dotted"} # well, yes, that's hard-coded for K = 2

        for k in range(self.K):

            I_and_Q_sym_total_vals_k = [self.N_total * (x + y) for x, y in zip(self.I_sym_vals[k], self.Q_sym_vals[k])]
            I_and_Q_sym_total_plot, = axes.plot(self.t_vals, I_and_Q_sym_total_vals_k,
                                                color = "tab:blue", linewidth=3.5, linestyle = group_linestyles[k])

            I_and_Q_asym_total_vals_k = [self.N_total * (x + y) for x, y in zip(self.I_asym_vals[k], self.Q_asym_vals[k])]
            I_and_Q_asym_total_plot, = axes.plot(self.t_vals, I_and_Q_asym_total_vals_k,
                                                 color = "tab:green", linewidth=3.5, linestyle = group_linestyles[k])

            I_and_Q_sev_total_vals_k = [self.N_total * (x + y) for x, y in zip(self.I_sev_vals[k], self.Q_sev_vals[k])]
            I_and_Q_sev_total_plot, = axes.plot(self.t_vals, I_and_Q_sev_total_vals_k,
                                                color = "tab:purple", linewidth=3.5, linestyle = group_linestyles[k])

            D_total_vals_k = [self.N_total * x for x in self.D_vals[k]]
            D_plot, = axes.plot(self.t_vals, D_total_vals_k, color = "tab:red", linewidth=3.5, linestyle = group_linestyles[k])

        I_and_Q_sev_total_vals = [self.N_total * (a + b + c + d) for a, b, c, d in zip(self.I_sev_vals[0],
                                                                                       self.I_sev_vals[1],
                                                                                       self.Q_sev_vals[0],
                                                                                       self.Q_sev_vals[1])]
        I_and_Q_sev_total_plot, = axes.plot(self.t_vals, I_and_Q_sev_total_vals, color = "tab:orange", linewidth=3.5, linestyle = "dashed")


        # Plot horizontal line for beds
        beds_plot = axes.hlines(self.beds, self.t_start, self.t_end, color = "black", linewidth=3.5, linestyle = "dashed")


        plot.yticks(np.arange(0, 400000, 50000), ['0',
                                                  r'$50\,$K',
                                                  r'$100\,$K',
                                                  r'$150\,$K',
                                                  r'$200\,$K',
                                                  r'$250\,$K',
                                                  r'$300\,$K',
                                                  r'$350\,$K'])
        ymax = 375000
        axes.set_ylim([-ymax * 0.05, ymax])
        plot.savefig(self.filename + "_figure_3.pdf")
        plot.close("all")

    def paper_plot_figure_4(self):
        figure, axes = plot.subplots()
        figure.subplots_adjust(bottom = 0.15)
        axes.grid(linestyle = ':', linewidth = 0.5, color = "#808080")

        # 0: old, 1: young
        group_linestyles = {0: "solid", 1: "dotted"} # well, yes, that's hard-coded for K = 2

        for k in range(self.K):

            I_and_Q_sev_total_vals_k = [self.N_total * (x + y) for x, y in zip(self.I_sev_vals[k], self.Q_sev_vals[k])]
            I_and_Q_sev_total_plot, = axes.plot(self.t_vals, I_and_Q_sev_total_vals_k,
                                                color = "tab:purple", linewidth=3.5, linestyle = group_linestyles[k])

            D_total_vals_k = [self.N_total * x for x in self.D_vals[k]]
            D_plot, = axes.plot(self.t_vals, D_total_vals_k, color = "tab:red", linewidth=3.5, linestyle = group_linestyles[k])

        I_and_Q_sev_total_vals = [self.N_total * (a + b + c + d) for a, b, c, d in zip(self.I_sev_vals[0],
                                                                                       self.I_sev_vals[1],
                                                                                       self.Q_sev_vals[0],
                                                                                       self.Q_sev_vals[1])]
        I_and_Q_sev_total_plot, = axes.plot(self.t_vals, I_and_Q_sev_total_vals, color = "tab:orange", linewidth=3.5, linestyle = "dashed")


        # Plot horizontal line for beds
        beds_plot = axes.hlines(self.beds, self.t_start, self.t_end, color = "black", linewidth=3.5, linestyle = "dashed")


        plot.yticks(np.arange(0, 120000, 20000), ['0',
                                                  r'$20\,$K',
                                                  r'$40\,$K',
                                                  r'$60\,$K',
                                                  r'$80\,$K',
                                                  r'$100\,$K'])
        ymax = 105000
        axes.set_ylim([-ymax * 0.05, ymax])
        plot.savefig(self.filename + "_figure_4.pdf")
        plot.close("all")

# class Visualizer
