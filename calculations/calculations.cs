using System.ComponentModel;
using System.Reflection;
using System.Reflection.Metadata.Ecma335;

namespace calculations
{
    public class Calculations
    {
        const double u_inf = 1;
        const double u_0 = 0;
        const int a = 0;
        const int b = 10;
        const double eps = 1e-7;

        static double Trapz(double[] u, double[] eta, double h, int k)
        {
            double[] u_slice = new double[k];
            double[] eta_slice = new double[k];

            Array.Copy(eta, eta_slice, k);
            Array.Copy(u, u_slice, k);

            double integral = 0;

            for (int i = 1; i < eta_slice.Length; i++)
                integral += (u_slice[i] + u_slice[i - 1]) / 2.0;

            return h * integral;
        }

        static double[] Gradient(double[] u, double h)
        {
            double[] gradient = new double[u.Length];

            for (int i = 1; i < gradient.Length - 1; i++)
                gradient[i] = (u[i + 1] - u[i - 1]) / (2.0 * h);

            return gradient;
        }


        static public double Error(double[] u, double[] previous_u){
            
            double[] difference_squared = u.Zip(previous_u, (x, y) => Math.Pow(x - y, 2)).ToArray();
    
            return Math.Sqrt(difference_squared.Average());
        }
        public static (double[], double[], double[]) Calculate_velocity_field(int K)
        {
            double c_Q = 0;

            double[] A = new double[K];
            double[] B = new double[K];
            double[] C = new double[K];
            double[] F = new double[K];
            double[] alpha = new double[K];
            double[] beta = new double[K];

            double[] u = new double[K];
            double[] u_prime = new double[K];
            double[] previous_f = new double[K];
            double[] previous_u = new double[K];
            double[] difference_squared = new double[K];
            double[] eta = new double[K];

            List<double> iteration_error_list = [];

            double h = (b - a) / Convert.ToDouble(K - 1);

            for (int i = 0; i < K; i++)
                eta[i] = a + h * i;

            B[0] = 1;
            C[0] = 0;
            F[0] = u_0 / u_inf;

            previous_f[0] = c_Q;

            A[K - 1] = 0;
            B[K - 1] = 1;
            F[K - 1] = 1;
            alpha[1] = -C[0] / B[0];
            beta[1] = F[0] / B[0];

            iteration_error_list.Add(1);

            while (iteration_error_list.Last() > eps)
            {
                for (int k = 0; k < K - 1; k++)
                {
                    A[k] = 4 - h * previous_f[k];
                    B[k] = -8;
                    C[k] = 4 + h * previous_f[k];
                    alpha[k + 1] = -C[k] / (A[k] * alpha[k] + B[k]);
                    beta[k + 1] = (F[k] - A[k] * beta[k]) / (A[k] * alpha[k] + B[k]);
                    u[K - 1] =
                        (F[K - 1] - A[K - 1] * beta[K - 1]) / (B[K - 1] + A[K - 1] * alpha[K - 1]);
                }

                for (int k = K - 2; k >= 0; k--)
                    u[k] = alpha[k + 1] * u[k + 1] + beta[k + 1];

                for (int k = 0; k < K - 2; k++)
                    previous_f[k] = Trapz(u, eta, h, k);

                iteration_error_list.Add(Error(u,previous_u));

                Array.Copy(u, previous_u, K);

            }

            return (eta, u, iteration_error_list.ToArray());
        }

        public static (double[], double[])  Calculate_velocity_field(double c_Q, int K)
        {
            double[] A = new double[K];
            double[] B = new double[K];
            double[] C = new double[K];
            double[] F = new double[K];
            double[] alpha = new double[K];
            double[] beta = new double[K];

            double[] u = new double[K];
            double[] u_prime = new double[K];
            double[] previous_f = new double[K];
            double[] previous_u = new double[K];
            double[] difference_squared = new double[K];
            double[] eta = new double[K];

            List<double> error = new List<double>();

            double h = (b - a) / Convert.ToDouble(K - 1);

            for (int i = 0; i < K; i++)
                eta[i] = a + h * i;

            B[0] = 1;
            C[0] = 0;
            F[0] = u_0 / u_inf;

            previous_f[0] = c_Q;

            A[K - 1] = 0;
            B[K - 1] = 1;
            F[K - 1] = 1;
            alpha[1] = -C[0] / B[0];
            beta[1] = F[0] / B[0];

            error.Add(1);

            while (error.Last() > eps)
            {
                for (int k = 0; k < K - 1; k++)
                {
                    A[k] = 4 - h * previous_f[k];
                    B[k] = -8;
                    C[k] = 4 + h * previous_f[k];
                    alpha[k + 1] = -C[k] / (A[k] * alpha[k] + B[k]);
                    beta[k + 1] = (F[k] - A[k] * beta[k]) / (A[k] * alpha[k] + B[k]);
                    u[K - 1] =
                        (F[K - 1] - A[K - 1] * beta[K - 1]) / (B[K - 1] + A[K - 1] * alpha[K - 1]);
                }

                for (int k = K - 2; k >= 0; k--)
                    u[k] = alpha[k + 1] * u[k + 1] + beta[k + 1];

                for (int k = 0; k < K - 2; k++)
                    previous_f[k] = Trapz(u, eta, h, k);

                difference_squared = u.Zip(previous_u, (x, y) => Math.Pow(x - y, 2)).ToArray();
                error.Add(Math.Sqrt(difference_squared.Average()));

                Array.Copy(u, previous_u, K);

                u_prime = Gradient(u, h);
            }

            return (eta, u);
        }
    }
}
