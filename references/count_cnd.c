// Sumber; https://labex.io/tutorials/c-compute-a-cumulative-distribution-function-cdf-in-c-435339
// SUMBER INI BERFUNGSI SEBAGAI REFERENSI MENGHITUNG CUMULATIVE NORMAL DISTRIBUTION

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Standard Normal CDF approximation function (Abramowitz and Stegun)
double standard_normal_cdf(double x) {
    const double a1 =  0.254829592;
    const double a2 = -0.284496736;
    const double a3 =  1.421413741;
    const double a4 = -1.453152027;
    const double a5 =  1.061405429;
    const double p  =  0.3275911;

    // Handle negative values
    int sign = (x < 0) ? -1 : 1;
    x = fabs(x);

    // Approximation formula
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

    return 0.5 * (1.0 + sign * y);
}

// Compute CDF for normal distribution
double normal_cdf(double x, double mean, double std_dev) {
    // Z-score calculation
    double z_score = (x - mean) / std_dev;
    return standard_normal_cdf(z_score);
}

int main() {
    // Distribution parameters
    double mean, std_dev;
    double x_value;

    // Prompt user for distribution parameters
    printf("Enter the mean (μ): ");
    scanf("%lf", &mean);

    printf("Enter the standard deviation (σ): ");
    scanf("%lf", &std_dev);

    // Prompt user for x value
    printf("Enter the x value to compute CDF: ");
    scanf("%lf", &x_value);

    // Calculate and print CDF
    double cdf_value = normal_cdf(x_value, mean, std_dev);

    printf("\nCDF Calculation Results:\n");
    printf("Mean (μ): %.2f\n", mean);
    printf("Standard Deviation (σ): %.2f\n", std_dev);
    printf("X Value: %.2f\n", x_value);
    printf("CDF P(X ≤ x): %.4f\n", cdf_value);

    return 0;
}