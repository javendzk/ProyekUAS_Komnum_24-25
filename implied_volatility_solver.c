#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "references/brent.h"

#define MAKS_PNJNG_NPUT 1000
#define MAKS_INPUT 100
#define PHI 3.14159265358979323846

typedef struct {
    char jenisOption[20];
    double stockPrice;
    double strikePrice;
    double waktuSisaExp;
    double riskFreeRate;
    double marketPrice;
    double sigLow;
    double sigHigh;
    double tolerance;
    int maxIterations;
} OptionData;

typedef struct {
    double impliedVolatility;
    int iterations;
    int converged;
    double finalError;
    double calculatedPrice;
    double priceDifference;
    double waktuExec;
} ResultData;

OptionData currentOption;

double standard_normal_cdf(double x);
double hitungD1(double S, double K, double T, double r, double sigma);
double hitungD2(double S, double K, double T, double r, double sigma);
double blackScholesCall(double S, double K, double T, double r, double sigma);
double blackScholesPut(double S, double K, double T, double r, double sigma);
double objectiveFunction(double sigma);
int readCsv(const char* filename, OptionData* options);
void writeHasil(const char* filename, OptionData* options, ResultData* results, int count);
int validateData(OptionData* option);
ResultData hitungIV(OptionData option);

int main() {
    OptionData options[MAKS_INPUT];
    ResultData results[MAKS_INPUT];
    int numberOfOptions;
    int i;
    
    printf("PROYEK UAS KOMNUM - Javana Muhammad Dzaki 2306161826\n");
    
    numberOfOptions = readCsv("data-io/data.csv", options);
    
    if (numberOfOptions <= 0) {
        printf("[!] Data kosong\n");
        return 1;
    }
        
    for (i = 0; i < numberOfOptions; i++) {
        printf("Processing options ke-%d (%s, K=%.2f)", i+1, options[i].jenisOption, options[i].strikePrice);
        
        results[i] = hitungIV(options[i]);
        
        if (results[i].converged) {
            printf("> done hitung iv= %.4f\n", results[i].impliedVolatility);
        } else {
            printf("[!] Fagal converge");
        }
    }
    
    writeHasil("data-io/results.csv", options, results, numberOfOptions);
    
    printf("- total options diproses: %d\n", numberOfOptions);
    
    int successCount = 0;
    for (i = 0; i < numberOfOptions; i++) {
        if (results[i].converged) successCount++;
    }
    
    printf("- berhasil konvergen: %d\n", successCount);
    printf("- fagal konvergen: %d\n", numberOfOptions - successCount);
    
    return 0;
}


double standard_normal_cdf(double x) {
    const double a1 =  0.254829592;
    const double a2 = -0.284496736;
    const double a3 =  1.421413741;
    const double a4 = -1.453152027;
    const double a5 =  1.061405429;
    const double p  =  0.3275911;
    
    int sign;
    double t;
    double y;
    
    sign = (x < 0) ? -1 : 1;
    x = fabs(x);
    
    t = 1.0 / (1.0 + p * x);
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);
    
    return 0.5 * (1.0 + sign * y);
}

double hitungD1(double S, double K, double T, double r, double sigma) {
    return (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
}

double hitungD2(double S, double K, double T, double r, double sigma) {
    return hitungD1(S, K, T, r, sigma) - sigma * sqrt(T);
}

double blackScholesCall(double S, double K, double T, double r, double sigma) {
    double d1;
    double d2;
    
    if (T <= 0.0 || sigma <= 0.0) {
        return 0.0;
    }
    
    d1 = hitungD1(S, K, T, r, sigma);
    d2 = hitungD2(S, K, T, r, sigma);
    
    return S * standard_normal_cdf(d1) - K * exp(-r * T) * standard_normal_cdf(d2);
}

double blackScholesPut(double S, double K, double T, double r, double sigma) {
    double d1;
    double d2;
    
    if (T <= 0.0 || sigma <= 0.0) {
        return 0.0;
    }
    
    d1 = hitungD1(S, K, T, r, sigma);
    d2 = hitungD2(S, K, T, r, sigma);
    
    return K * exp(-r * T) * standard_normal_cdf(-d2) - S * standard_normal_cdf(-d1);
}

double objectiveFunction(double sigma) {
    double calculatedPrice;
    
    if (strcmp(currentOption.jenisOption, "option_call") == 0) {
        calculatedPrice = blackScholesCall(currentOption.stockPrice, currentOption.strikePrice,
                                         currentOption.waktuSisaExp, currentOption.riskFreeRate, sigma);
    } else {
        calculatedPrice = blackScholesPut(currentOption.stockPrice, currentOption.strikePrice,
                                        currentOption.waktuSisaExp, currentOption.riskFreeRate, sigma);
    }
    
    return calculatedPrice - currentOption.marketPrice;
}

int readCsv(const char* filename, OptionData* options) {
    FILE* file;
    char line[MAKS_PNJNG_NPUT];
    char* token;
    int count;
    int i;
    
    count = 0;
    file = fopen(filename, "r");
    
    if (file == NULL) {
        printf("[!] Error file data null]\n");
        return -1;
    }
        
    if (fgets(line, sizeof(line), file) == NULL) {
        printf("[!] File data kosong\n");
        fclose(file);
        return -1;
    }
    
    while (fgets(line, sizeof(line), file) && count < MAKS_INPUT) {
        i = 0;
        token = strtok(line, ",");
        
        while (token != NULL && i < 10) {
            switch (i) {
                case 0:
                    strcpy(options[count].jenisOption, token);
                    break;
                case 1:
                    options[count].stockPrice = atof(token);
                    break;
                case 2:
                    options[count].strikePrice = atof(token);
                    break;
                case 3:
                    options[count].waktuSisaExp = atof(token);
                    break;
                case 4:
                    options[count].riskFreeRate = atof(token);
                    break;
                case 5:
                    options[count].marketPrice = atof(token);
                    break;
                case 6:
                    options[count].sigLow = atof(token);
                    break;
                case 7:
                    options[count].sigHigh = atof(token);
                    break;
                case 8:
                    options[count].tolerance = atof(token);
                    break;
                case 9:
                    options[count].maxIterations = atoi(token);
                    break;
            }
            token = strtok(NULL, ",");
            i++;
        }
        count++;
    }
    
    fclose(file);
    return count;
}

void writeHasil(const char* filename, OptionData* options, ResultData* results, int count) {
    FILE* file;
    int i;
    
    file = fopen(filename, "w");
    
    if (file == NULL) {
        printf("[!] Error file adalah null");
        return;
    }
    
    fprintf(file, "problem_type,S,K,T,r,C_market,implied_volatility,iterations,converged,error,calculated_price,price_difference,time_ms\n");
    
    for (i = 0; i < count; i++) {
        fprintf(file, "%s,%.2f,%.2f,%.6f,%.3f,%.2f,%.6f,%d,%s,%.2e,%.6f,%.6f,%.2f\n",
                options[i].jenisOption,
                options[i].stockPrice,
                options[i].strikePrice,
                options[i].waktuSisaExp,
                options[i].riskFreeRate,
                options[i].marketPrice,
                results[i].impliedVolatility,
                results[i].iterations,
                results[i].converged ? "yes" : "no",
                results[i].finalError,
                results[i].calculatedPrice,
                results[i].priceDifference,
                results[i].waktuExec);
    }
    
    fclose(file);
    printf("[v] Program selesai\n");
}

int validateData(OptionData* option) {
    double intrinsicValue;
    
    if (option->stockPrice <= 0 || option->strikePrice <= 0 || 
        option->waktuSisaExp <= 0 || option->riskFreeRate < 0 || 
        option->marketPrice <= 0) {
        return 0;
    }
    
    if (strcmp(option->jenisOption, "option_call") == 0) {
        intrinsicValue = fmax(option->stockPrice - option->strikePrice * exp(-option->riskFreeRate * option->waktuSisaExp), 0.0);
        if (option->marketPrice < intrinsicValue) {
            return 0;
        }
    } else if (strcmp(option->jenisOption, "option_put") == 0) {
        intrinsicValue = fmax(option->strikePrice * exp(-option->riskFreeRate * option->waktuSisaExp) - option->stockPrice, 0.0);
        if (option->marketPrice < intrinsicValue) {
            return 0;
        }
    }
    
    return 1;
}

ResultData hitungIV(OptionData option) {
    ResultData result;
    double machineEpsilon;
    double fLow;
    double fHigh;
    double impliedVol;
    clock_t startTime;
    clock_t endTime;
    
    result.converged = 0;
    result.iterations = 0;
    result.finalError = 1.0;
    result.impliedVolatility = 0.0;
    result.calculatedPrice = 0.0;
    result.priceDifference = 0.0;
    result.waktuExec = 0.0;
    
    startTime = clock();
    
    if (!validateData(&option)) {
        return result;
    }
    
    currentOption = option;
    
    fLow = objectiveFunction(option.sigLow);
    fHigh = objectiveFunction(option.sigHigh);
    
    if (fLow * fHigh > 0) {
        printf("[!] bracket tidak valid\n");
        return result;
    }
    
    machineEpsilon = 2.22e-16;
    
    impliedVol = zero(option.sigLow, option.sigHigh, machineEpsilon, 
                      option.tolerance, objectiveFunction);
    
    result.impliedVolatility = impliedVol;
    result.finalError = fabs(objectiveFunction(impliedVol));
    result.converged = (result.finalError < option.tolerance) ? 1 : 0;
    
    if (strcmp(option.jenisOption, "option_call") == 0) {
        result.calculatedPrice = blackScholesCall(option.stockPrice, option.strikePrice,
                                                option.waktuSisaExp, option.riskFreeRate, impliedVol);
    } else {
        result.calculatedPrice = blackScholesPut(option.stockPrice, option.strikePrice,
                                               option.waktuSisaExp, option.riskFreeRate, impliedVol);
    }
    
    result.priceDifference = fabs(result.calculatedPrice - option.marketPrice);
    
    endTime = clock();
    result.waktuExec = ((double)(endTime - startTime)) / CLOCKS_PER_SEC * 1000.0;
    
    return result;
}
