# Global imports
import os

# Local imports
from data_reader import Data_Reader
from seiiird_model import SEIIIRD_Model
from seiiird_tracing_model import SEIIIRD_Tracing_Model
from explicit_euler import Explicit_Euler
from result_analyzer import Result_Analyzer
from visualizer import Visualizer

# Start and end time, stepsize
t_start = 0
t_end = 500
stepsize = 1e-2

# Preparing data parsing
data_set_name = "Figure-4"
data_directory_name = "../data/" + data_set_name + "/"
data_directory = os.fsencode(data_directory_name)

for file in os.listdir(data_directory):
    data_filename = os.fsdecode(file)
    if data_filename.endswith(".csv"):

        # Parsing the data
        data_filename_prefix = data_filename.split(".")[0]
        data_reader = Data_Reader()
        packed_data = data_reader.read_from_csv_file(data_directory_name + data_filename_prefix + ".csv")

        # Unpacking some of the parsed data
        tracing_data_given = packed_data[0]
        assert(isinstance(tracing_data_given, bool))
        K = packed_data[1]
        N = packed_data[2]
        beds = packed_data[-2]
        x0 = packed_data[-1]

        # Instantiate the ODE system class
        if tracing_data_given:
            ode_system = SEIIIRD_Tracing_Model(packed_data[1:-1])
        else:
            ode_system = SEIIIRD_Model(packed_data[1:-1])

        # Solve the ODE system
        explicit_euler = Explicit_Euler(ode_system, stepsize)
        results = explicit_euler.solve(t_start, x0, t_end)

        # Analyze results
        result_analyzer = Result_Analyzer(tracing_data_given, results, stepsize, N, K, beds, t_start, t_end,
                                          "../results/" + data_set_name + "/" + data_filename_prefix)
        result_analyzer.analyze()

        # Visualize results
        N_total = sum(N)
        visualizer = Visualizer(tracing_data_given, results, N, N_total, K, beds, t_start, t_end,
                                "../results/" + data_set_name + "/" + data_filename_prefix)
        visualizer.plot()
        visualizer.paper_plot_figure_1()
        if tracing_data_given:
            visualizer.paper_plot_figure_3()
            visualizer.paper_plot_figure_4()
