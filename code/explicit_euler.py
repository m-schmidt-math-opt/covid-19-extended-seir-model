class Explicit_Euler:

    def __init__(self, ode, stepsize):
        self.ode = ode
        self.stepsize = stepsize

    def solve(self, t_start, x_start, t_end):
        print("Solve SEIR model with explicit Euler scheme ...")
        t = t_start
        x = x_start
        results = [[t, x]]
        nr_of_steps = int((t_end - t_start) / self.stepsize) + 1
        print("Number of steps: " + str(nr_of_steps))
        for k in range(0, nr_of_steps):
            x = x + self.stepsize * self.ode.eval_rhs(x)
            t = t_start + k * self.stepsize
            results.append([t, x])
        return results

# class Explicit_Euler
