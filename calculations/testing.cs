using calculations;

namespace testing
{
    public class Testing
    {
        const double eps = 1e-5;
        // grid step size sensitivity

        public static int Step_size()
        {
            double[] eta, u, error;

            double u_error = 1, previous_u_0 = 0;

            int K = 1000;

            while (u_error > eps)
            {
                (eta, u, error) = Calculations.Calculate_velocity_field(K);
                using StreamWriter writer_u = new(Convert.ToString($"/Users/levonoganesyan/Desktop/Projects/GitHub/CP1/results/grid_step_size_sensitivity_study/data/testing_K={K}.txt"));

                for (int i = 0; i < K; i++)
                {
                    writer_u.WriteLine("{0} {1}", eta[i], u[i]);
                }

                u_error = Math.Abs(u[0] - previous_u_0);

                previous_u_0 = u[0];

                K += 2000;


            }

            return K;

        }

        public static void C_Q(int K)
        {

            double[] eta, u;

            for (double c_Q = -1.5; c_Q <= 1.5; c_Q += 0.5)
            {

                (eta, u) = Calculations.Calculate_velocity_field(c_Q, K);
            }
        }
    }
}
