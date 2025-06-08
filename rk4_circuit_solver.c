#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_LINE_LENGTH 1024
#define MAX_CIRCUITS 50
#define MAX_TIME_STEPS 100000

typedef struct {
    char circuitType[20];
    double R;
    double C;
    double L;
    double V_source;
    double I_s;
    double n;
    double V_t;
    double V0;
    double I0;
    double t_start;
    double t_end;
    double step_size;
    double output_interval;
} CircuitParams;

typedef struct {
    double time;
    double V_capacitor;
    double V_diode;
    double I_circuit;
    double I_inductor;
    int step_count;
    double execution_time_ms;
} SimulationResult;

double calculateDiodeVoltage(double current, double I_s, double n, double V_t) {
    if (current <= 0) {
        return 0.0;
    }
    
    double ratio = current / I_s;
    if (ratio < 1e-15) ratio = 1e-15;
    
    return n * V_t * log(ratio + 1.0);
}

double calculateCircuitCurrent(double V_cap, double V_source, double R, double I_s, double n, double V_t) {
    double I_guess = (V_source - V_cap) / R;
    
    for (int iter = 0; iter < 10; iter++) {
        double V_diode = calculateDiodeVoltage(I_guess, I_s, n, V_t);
        double I_new = (V_source - V_cap - V_diode) / R;
        
        if (fabs(I_new - I_guess) < 1e-9) break;
        I_guess = I_new;
        
        if (I_guess < 0) I_guess = 0;
    }
    
    return I_guess;
}

double rcDiodeDerivative(double t, double V, CircuitParams* params) {
    double I = calculateCircuitCurrent(V, params->V_source, params->R, 
                                     params->I_s, params->n, params->V_t);
    
    return I / params->C;
}

void rlcDiodeDerivatives(double t, double V, double I_L, double* dVdt, double* dIdt, CircuitParams* params) {
    double V_diode = calculateDiodeVoltage(I_L, params->I_s, params->n, params->V_t);
    
    *dVdt = I_L / params->C;
    
    *dIdt = (params->V_source - V - V_diode - I_L * params->R) / params->L;
    
    if (I_L < 0) {
        *dIdt = 0;
    }
}

void rungeKutta4_RC(CircuitParams* params, SimulationResult* results, int* result_count) {
    double t = params->t_start;
    double V = params->V0;
    double h = params->step_size;
    double next_output = params->t_start;
    int step_count = 0;
    *result_count = 0;
    
    clock_t start_time = clock();
    
    printf("R=%.0fΩ, C=%.2e F, V_sumber=%.1fV\n", params->R, params->C, params->V_source);
    printf("Dioda: I_s=%.2e A, n=%.2f, V_t=%.5f V\n", params->I_s, params->n, params->V_t);
    
    while (t <= params->t_end && *result_count < MAX_TIME_STEPS) {
        if (t >= next_output) {
            double I = calculateCircuitCurrent(V, params->V_source, params->R, 
                                             params->I_s, params->n, params->V_t);
            double V_diode = calculateDiodeVoltage(I, params->I_s, params->n, params->V_t);
            
            results[*result_count].time = t;
            results[*result_count].V_capacitor = V;
            results[*result_count].V_diode = V_diode;
            results[*result_count].I_circuit = I;
            results[*result_count].I_inductor = 0;
            results[*result_count].step_count = step_count;
            
            clock_t current_time = clock();
            results[*result_count].execution_time_ms = 
                (double)(current_time - start_time) * 1000.0 / CLOCKS_PER_SEC;
            
            (*result_count)++;
            next_output += params->output_interval;
        }
        
        double k1 = h * rcDiodeDerivative(t, V, params);
        double k2 = h * rcDiodeDerivative(t + h/2, V + k1/2, params);
        double k3 = h * rcDiodeDerivative(t + h/2, V + k2/2, params);
        double k4 = h * rcDiodeDerivative(t + h, V + k3, params);
        
        V = V + (k1 + 2*k2 + 2*k3 + k4) / 6.0;
        t = t + h;
        step_count++;
        
        if (V < 0) V = 0;
        if (V > params->V_source) V = params->V_source;
    }
    
    printf("Simulasi dl: %d langkah, %d data output\n", step_count, *result_count);
}

void rungeKutta4_RLC(CircuitParams* params, SimulationResult* results, int* result_count) {
    double t = params->t_start;
    double V = params->V0;
    double I_L = params->I0;
    double h = params->step_size;
    double next_output = params->t_start;
    int step_count = 0;
    *result_count = 0;
    
    clock_t start_time = clock();
    
    printf("R=%.0fΩ, L=%.2e H, C=%.2e F, V_source=%.1fV\n", 
           params->R, params->L, params->C, params->V_source);
    printf("Dioda: I_s=%.2e A, n=%.2f, V_t=%.5f V\n", params->I_s, params->n, params->V_t);
    
    while (t <= params->t_end && *result_count < MAX_TIME_STEPS) {
        if (t >= next_output) {
            double V_diode = calculateDiodeVoltage(I_L, params->I_s, params->n, params->V_t);
            
            results[*result_count].time = t;
            results[*result_count].V_capacitor = V;
            results[*result_count].V_diode = V_diode;
            results[*result_count].I_circuit = I_L;
            results[*result_count].I_inductor = I_L;
            results[*result_count].step_count = step_count;
            
            clock_t current_time = clock();
            results[*result_count].execution_time_ms = 
                (double)(current_time - start_time) * 1000.0 / CLOCKS_PER_SEC;
            
            (*result_count)++;
            next_output += params->output_interval;
        }
        
        double dVdt1, dIdt1, dVdt2, dIdt2, dVdt3, dIdt3, dVdt4, dIdt4;
        
        rlcDiodeDerivatives(t, V, I_L, &dVdt1, &dIdt1, params);
        double k1_V = h * dVdt1;
        double k1_I = h * dIdt1;
        
        rlcDiodeDerivatives(t + h/2, V + k1_V/2, I_L + k1_I/2, &dVdt2, &dIdt2, params);
        double k2_V = h * dVdt2;
        double k2_I = h * dIdt2;
        
        rlcDiodeDerivatives(t + h/2, V + k2_V/2, I_L + k2_I/2, &dVdt3, &dIdt3, params);
        double k3_V = h * dVdt3;
        double k3_I = h * dIdt3;
        
        rlcDiodeDerivatives(t + h, V + k3_V, I_L + k3_I, &dVdt4, &dIdt4, params);
        double k4_V = h * dVdt4;
        double k4_I = h * dIdt4;
        
        V = V + (k1_V + 2*k2_V + 2*k3_V + k4_V) / 6.0;
        I_L = I_L + (k1_I + 2*k2_I + 2*k3_I + k4_I) / 6.0;
        
        t = t + h;
        step_count++;
        
        if (V < 0) V = 0;
        if (I_L < 0) I_L = 0;
    }
    
    printf("Simulasi done: %d step, %d data output\n", step_count, *result_count);
}

int readTransientCsvFile(const char* filename, CircuitParams* circuits) {
    FILE* file = fopen(filename, "r");
    
    char line[MAX_LINE_LENGTH];
    int circuitCount = 0;
    int lineNum = 0;
    
    while (fgets(line, MAX_LINE_LENGTH, file) != NULL && circuitCount < MAX_CIRCUITS) {
        if (lineNum == 0) {
            lineNum++;
            continue;
        }
        
        CircuitParams* circuit = &circuits[circuitCount];
        char* token;
        int tokenIndex = 0;
        
        token = strtok(line, ",");
        while (token != NULL && tokenIndex < 12) {
            switch (tokenIndex) {
                case 0: strcpy(circuit->circuitType, token); break;
                case 1: circuit->R = atof(token); break;
                case 2: circuit->C = atof(token); break;
                case 3: 
                    if (strcmp(circuit->circuitType, "rlc_diode") == 0) {
                        circuit->L = atof(token);
                    } else {
                        circuit->V_source = atof(token);
                    }
                    break;
                case 4: 
                    if (strcmp(circuit->circuitType, "rlc_diode") == 0) {
                        circuit->V_source = atof(token);
                    } else {
                        circuit->I_s = atof(token);
                    }
                    break;
                case 5: 
                    if (strcmp(circuit->circuitType, "rlc_diode") == 0) {
                        circuit->I_s = atof(token);
                    } else {
                        circuit->n = atof(token);
                    }
                    break;
                case 6: 
                    if (strcmp(circuit->circuitType, "rlc_diode") == 0) {
                        circuit->n = atof(token);
                    } else {
                        circuit->V_t = atof(token);
                    }
                    break;
                case 7: 
                    if (strcmp(circuit->circuitType, "rlc_diode") == 0) {
                        circuit->V_t = atof(token);
                    } else {
                        circuit->V0 = atof(token);
                    }
                    break;
                case 8: 
                    if (strcmp(circuit->circuitType, "rlc_diode") == 0) {
                        circuit->V0 = atof(token);
                    } else {
                        circuit->t_start = atof(token);
                    }
                    break;
                case 9: 
                    if (strcmp(circuit->circuitType, "rlc_diode") == 0) {
                        circuit->t_start = atof(token);
                    } else {
                        circuit->t_end = atof(token);
                    }
                    break;
                case 10: 
                    if (strcmp(circuit->circuitType, "rlc_diode") == 0) {
                        circuit->t_end = atof(token);
                    } else {
                        circuit->step_size = atof(token);
                    }
                    break;
                case 11: 
                    if (strcmp(circuit->circuitType, "rlc_diode") == 0) {
                        circuit->step_size = atof(token);
                    } else {
                        circuit->output_interval = atof(token);
                    }
                    break;
            }
            token = strtok(NULL, ",");
            tokenIndex++;
        }
        
        if (strcmp(circuit->circuitType, "rlc_diode") == 0) {
            circuit->I0 = 0.0;
            if (token != NULL) {
                circuit->output_interval = atof(token);
            }
        }
        
        circuitCount++;
        lineNum++;
    }
    
    fclose(file);
    return circuitCount;
}

void writeTransientResults(const char* filename, CircuitParams* circuits, int circuitCount, 
                          SimulationResult** allResults, int* resultCounts) {
    FILE* file = fopen(filename, "w");
    
    fprintf(file, "circuit_type,R,C,L,V_source,I_s,n,V_t,V0,I0,time,V_capacitor,V_diode,I_circuit,I_inductor,step_count,execution_time_ms\n");
    
    for (int i = 0; i < circuitCount; i++) {
        CircuitParams* circuit = &circuits[i];
        SimulationResult* results = allResults[i];
        int count = resultCounts[i];
        
        for (int j = 0; j < count; j++) {
            fprintf(file, "%s,%.1f,%.2e,%.2e,%.2f,%.2e,%.2f,%.5f,%.3f,%.3f,%.6e,%.6f,%.6f,%.6e,%.6e,%d,%.3f\n",
                   circuit->circuitType,
                   circuit->R,
                   circuit->C,
                   (strcmp(circuit->circuitType, "rlc_diode") == 0) ? circuit->L : 0.0,
                   circuit->V_source,
                   circuit->I_s,
                   circuit->n,
                   circuit->V_t,
                   circuit->V0,
                   (strcmp(circuit->circuitType, "rlc_diode") == 0) ? circuit->I0 : 0.0,
                   results[j].time,
                   results[j].V_capacitor,
                   results[j].V_diode,
                   results[j].I_circuit,
                   results[j].I_inductor,
                   results[j].step_count,
                   results[j].execution_time_ms);
        }
    }
    
    fclose(file);
    printf("Hasil simulasi berhasil disimpan ke %s\n", filename);
}

int main() {
    printf("ANALISIS RANGKAIAN DIODA TRANSIEN - METODE RK4\n");
    printf("Komputasi Numerik (ENCE604015) - UAS 2024/2025\n");
    printf("Javana Muhammad Dzaki - NPM 2306161826\n\n");
    
    CircuitParams circuits[MAX_CIRCUITS];
    SimulationResult* allResults[MAX_CIRCUITS];
    int resultCounts[MAX_CIRCUITS];
    
    int circuitCount = readTransientCsvFile("./data-io/transient_data.csv", circuits);
    if (circuitCount == 0) {
        printf("[!] Data kosong\n");
        return 1;
    }
    
    for (int i = 0; i < circuitCount; i++) {
        allResults[i] = (SimulationResult*)malloc(MAX_TIME_STEPS * sizeof(SimulationResult));
        if (!allResults[i]) {
            printf("[!] Fail allocate memori \n");
            return 1;
        }
    }
    
    for (int i = 0; i < circuitCount; i++) {
        printf("\n> Rangkaian ke %d/%d: %s\n", i+1, circuitCount, circuits[i].circuitType);
        
        if (strcmp(circuits[i].circuitType, "rc_diode") == 0) {
            rungeKutta4_RC(&circuits[i], allResults[i], &resultCounts[i]);
        } else if (strcmp(circuits[i].circuitType, "rlc_diode") == 0) {
            rungeKutta4_RLC(&circuits[i], allResults[i], &resultCounts[i]);
        } else {
            printf("[!] Circuit type '%s' tidak dikenal\n", circuits[i].circuitType);
            resultCounts[i] = 0;
        }
        
        if (resultCounts[i] > 0) {
            SimulationResult* lastResult = &allResults[i][resultCounts[i]-1];
            printf("Hasil akhir: V_cap=%.4fV, I=%.6fA, V_diode=%.4fV\n", 
                   lastResult->V_capacitor, lastResult->I_circuit, lastResult->V_diode);
            printf("Const waktu τ = RC = %.6f s\n", circuits[i].R * circuits[i].C);
        }
    }
    
    writeTransientResults("./data-io/transient_results.csv", circuits, circuitCount, allResults, resultCounts);
    
    for (int i = 0; i < circuitCount; i++) {
        free(allResults[i]);
    }
    
    printf("Total rangkaian dianalisis: %d\n", circuitCount);
    
    return 0;
}
