#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "references/brent.h"

#define MAX_LINE_LENGTH 1000
#define MAX_OPTIONS 100
#define PI 3.14159265358979323846

typedef struct {
    char problemType[20];
    double stockPrice;
    double strikePrice;
    double timeToExpiry;
    double riskFreeRate;
    double marketPrice;
    double sigmaLow;
    double sigmaHigh;
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
    double executionTime;
} ResultData;

OptionData currentOption;

double cumulativeNormalDistribution(double x) {
    double t;
    double z;
    double ans;
    
    z = fabs(x);
    t = 1.0 / (1.0 + 0.5 * z);
    
    ans = t * exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + 
          t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + 
          t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + 
          t * 0.17087277)))))))));
    
    if (x >= 0.0) {
        return 1.0 - ans;
    } else {
        return ans;
    }
}

double calculateD1(double S, double K, double T, double r, double sigma) {
    return (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
}

double calculateD2(double S, double K, double T, double r, double sigma) {
    return calculateD1(S, K, T, r, sigma) - sigma * sqrt(T);
}

double blackScholesCall(double S, double K, double T, double r, double sigma) {
    double d1;
    double d2;
    
    if (T <= 0.0 || sigma <= 0.0) {
        return 0.0;
    }
    
    d1 = calculateD1(S, K, T, r, sigma);
    d2 = calculateD2(S, K, T, r, sigma);
    
    return S * cumulativeNormalDistribution(d1) - K * exp(-r * T) * cumulativeNormalDistribution(d2);
}

double blackScholesPut(double S, double K, double T, double r, double sigma) {
    double d1;
    double d2;
    
    if (T <= 0.0 || sigma <= 0.0) {
        return 0.0;
    }
    
    d1 = calculateD1(S, K, T, r, sigma);
    d2 = calculateD2(S, K, T, r, sigma);
    
    return K * exp(-r * T) * cumulativeNormalDistribution(-d2) - S * cumulativeNormalDistribution(-d1);
}

double objectiveFunction(double sigma) {
    double calculatedPrice;
    
    if (strcmp(currentOption.problemType, "option_call") == 0) {
        calculatedPrice = blackScholesCall(currentOption.stockPrice, currentOption.strikePrice,
                                         currentOption.timeToExpiry, currentOption.riskFreeRate, sigma);
    } else {
        calculatedPrice = blackScholesPut(currentOption.stockPrice, currentOption.strikePrice,
                                        currentOption.timeToExpiry, currentOption.riskFreeRate, sigma);
    }
    
    return calculatedPrice - currentOption.marketPrice;
}

int readCsvData(const char* filename, OptionData* options) {
    FILE* file;
    char line[MAX_LINE_LENGTH];
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
    
    while (fgets(line, sizeof(line), file) && count < MAX_OPTIONS) {
        i = 0;
        token = strtok(line, ",");
        
        while (token != NULL && i < 10) {
            switch (i) {
                case 0:
                    strcpy(options[count].problemType, token);
                    break;
                case 1:
                    options[count].stockPrice = atof(token);
                    break;
                case 2:
                    options[count].strikePrice = atof(token);
                    break;
                case 3:
                    options[count].timeToExpiry = atof(token);
                    break;
                case 4:
                    options[count].riskFreeRate = atof(token);
                    break;
                case 5:
                    options[count].marketPrice = atof(token);
                    break;
                case 6:
                    options[count].sigmaLow = atof(token);
                    break;
                case 7:
                    options[count].sigmaHigh = atof(token);
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

void writeResults(const char* filename, OptionData* options, ResultData* results, int count) {
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
                options[i].problemType,
                options[i].stockPrice,
                options[i].strikePrice,
                options[i].timeToExpiry,
                options[i].riskFreeRate,
                options[i].marketPrice,
                results[i].impliedVolatility,
                results[i].iterations,
                results[i].converged ? "yes" : "no",
                results[i].finalError,
                results[i].calculatedPrice,
                results[i].priceDifference,
                results[i].executionTime);
    }
    
    fclose(file);
    printf("[v] Program selesai\n");
}

int validateOptionData(OptionData* option) {
    double intrinsicValue;
    
    if (option->stockPrice <= 0 || option->strikePrice <= 0 || 
        option->timeToExpiry <= 0 || option->riskFreeRate < 0 || 
        option->marketPrice <= 0) {
        return 0;
    }
    
    if (strcmp(option->problemType, "option_call") == 0) {
        intrinsicValue = fmax(option->stockPrice - option->strikePrice * exp(-option->riskFreeRate * option->timeToExpiry), 0.0);
        if (option->marketPrice < intrinsicValue) {
            return 0;
        }
    } else if (strcmp(option->problemType, "option_put") == 0) {
        intrinsicValue = fmax(option->strikePrice * exp(-option->riskFreeRate * option->timeToExpiry) - option->stockPrice, 0.0);
        if (option->marketPrice < intrinsicValue) {
            return 0;
        }
    }
    
    return 1;
}

ResultData calculateImpliedVolatility(OptionData option) {
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
    result.executionTime = 0.0;
    
    startTime = clock();
    
    if (!validateOptionData(&option)) {
        return result;
    }
    
    currentOption = option;
    
    fLow = objectiveFunction(option.sigmaLow);
    fHigh = objectiveFunction(option.sigmaHigh);
    
    if (fLow * fHigh > 0) {
        printf("[!] bracket tidak valid\n");
        return result;
    }
    
    machineEpsilon = 2.22e-16;
    
    impliedVol = zero(option.sigmaLow, option.sigmaHigh, machineEpsilon, 
                      option.tolerance, objectiveFunction);
    
    result.impliedVolatility = impliedVol;
    result.finalError = fabs(objectiveFunction(impliedVol));
    result.converged = (result.finalError < option.tolerance) ? 1 : 0;
    
    if (strcmp(option.problemType, "option_call") == 0) {
        result.calculatedPrice = blackScholesCall(option.stockPrice, option.strikePrice,
                                                option.timeToExpiry, option.riskFreeRate, impliedVol);
    } else {
        result.calculatedPrice = blackScholesPut(option.stockPrice, option.strikePrice,
                                               option.timeToExpiry, option.riskFreeRate, impliedVol);
    }
    
    result.priceDifference = fabs(result.calculatedPrice - option.marketPrice);
    
    endTime = clock();
    result.executionTime = ((double)(endTime - startTime)) / CLOCKS_PER_SEC * 1000.0;
    
    return result;
}

int main() {
    OptionData options[MAX_OPTIONS];
    ResultData results[MAX_OPTIONS];
    int numberOfOptions;
    int i;
    
    printf("PROYEK UAS KOMNUM - Javana Muhammad Dzaki 2306161826\n");
    
    numberOfOptions = readCsvData("data-io/data.csv", options);
    
    if (numberOfOptions <= 0) {
        printf("[!] Data kosong\n");
        return 1;
    }
        
    for (i = 0; i < numberOfOptions; i++) {
        printf("Processing options ke-%d (%s, K=%.2f)", i+1, options[i].problemType, options[i].strikePrice);
        
        results[i] = calculateImpliedVolatility(options[i]);
        
        if (results[i].converged) {
            printf("> done hitung iv= %.4f\n", results[i].impliedVolatility);
        } else {
            printf("[!] Fagal converge");
        }
    }
    
    writeResults("data-io/results.csv", options, results, numberOfOptions);
    
    printf("- total options diproses: %d\n", numberOfOptions);
    
    int successCount = 0;
    for (i = 0; i < numberOfOptions; i++) {
        if (results[i].converged) successCount++;
    }
    
    printf("- berhasil konvergen: %d\n", successCount);
    printf("- fagal konvergen: %d\n", numberOfOptions - successCount);
    
    return 0;
}