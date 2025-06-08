/*
Analisis Rangkaian RC Dioda Transiet dengan Metode Runge-Kutta Order 4
Proyek UAS Komnum Javana Muhammad Dzaki NPM 2306161826
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_PANJANG 1024
#define MAX_DATA 50
#define MAX_TIME_STEPS 100000

typedef struct { 
    double R;                   // resistansi (ohm)
    double C;                   // kapasitansi (farad)
    double V_sumber;           
    double I_s;                 // arus saturasi dioda (ampere)
    double n;                   // faktor idealitas diodanya
    double V_t;                 
    double V0;               
    double t_start;             // time mulai simulasi 
    double t_end;               // time akhir simulasi 
    double step_size;           // step size rk4nya
    double output_interval;    
} inputRangkaian;

typedef struct {
    double time;
    double V_kapasitor;
    double V_diode;
    double I_rangkaian;
    int step_count;
    double waktu_ms;
} hasilSimulasi;

double hitungVDioda(double current, double I_s, double n, double V_t);
double hitungArus(double V_cap, double V_sumber, double R, double I_s, double n, double V_t);
double hitungTurunanDioda(double t, double V, inputRangkaian* params);
void rungeKutta4_RC(inputRangkaian* params, hasilSimulasi* results, int* result_count);
int readInput(const char* filename, inputRangkaian* circuits);
void writeOutput(const char* filename, inputRangkaian* circuits, int circuitCount, hasilSimulasi** allResults, int* resultCounts);

int main() {
    printf("ANALISIS RANGKAIAN RC DIODA TRANSIEN - METODE RK4\n");
    printf("Komputasi Numerik (ENCE604015) - Proyek UAS 2024/2025\n");
    printf("Javana Muhammad Dzaki - NPM 2306161826\n\n");
    
    inputRangkaian circuits[MAX_DATA];
    hasilSimulasi* allResults[MAX_DATA];
    int resultCounts[MAX_DATA];
    
    // read params dari input
    int circuitCount = readInput("./data-io/transient_data.csv", circuits);
    if (circuitCount == 0) {
        printf("[!] Data kosong\n");
        return 1;
    }
    
    for (int i = 0; i < circuitCount; i++) {
        allResults[i] = (hasilSimulasi*)malloc(MAX_TIME_STEPS * sizeof(hasilSimulasi));
        if (!allResults[i]) {
            printf("[!] Fail allocate memori \n");
            return 1;
        }
    }
    
    // processing setiap rangkaian RC
    for (int i = 0; i < circuitCount; i++) {
        printf("\n> Rangkaian ke %d/%d: RC-Dioda\n", i+1, circuitCount);
        
        rungeKutta4_RC(&circuits[i], allResults[i], &resultCounts[i]);
        
        // show final result
        if (resultCounts[i] > 0) {
            hasilSimulasi* lastResult = &allResults[i][resultCounts[i]-1];
            printf("Hasil akhir: V_cap=%.4fV, I=%.6fA, V_diode=%.4fV\n", 
                   lastResult->V_kapasitor, lastResult->I_rangkaian, lastResult->V_diode);
            printf("Const waktu τ = RC = %.6f s\n", circuits[i].R * circuits[i].C);
        }
    }
    
    writeOutput("./data-io/transient_results.csv", circuits, circuitCount, allResults, resultCounts);
    
    // clean memori
    for (int i = 0; i < circuitCount; i++) {
        free(allResults[i]);
    }
            
    return 0;
}

// function untuk hitung tegangan dioda menggunakan shockley diode equ
double hitungVDioda(double current, double I_s, double n, double V_t) {
    if (current <= 0) {
        return 0.0;  // dioda reverse bias
    }
    
    // V_dioda = n * V_t * ln((I/I_s) + 1)
    double ratio = current / I_s;
    if (ratio < 1e-15) ratio = 1e-15; 
    
    return n * V_t * log(ratio + 1.0);
}

// fungsi untuk menghitung arus rangkaian dari tegangan kapasitor
double hitungArus(double V_cap, double V_sumber, double R, double I_s, double n, double V_t) {
    // loop untuk mencari arus yang konsisten dengan tegangan dioda
    double I_guess = (V_sumber - V_cap) / R;  // estimasi awal tanpa dioda
    
    for (int iter = 0; iter < 10; iter++) {
        double V_diode = hitungVDioda(I_guess, I_s, n, V_t);
        double I_new = (V_sumber - V_cap - V_diode) / R;
        
        if (fabs(I_new - I_guess) < 1e-9) break;
        I_guess = I_new;
        
        if (I_guess < 0) I_guess = 0;  // dioda tidak bisa reverse current
    }
    
    return I_guess;
}

// fungsi turunan untuk RC diode circuit yaitu  dV/dt
double hitungTurunanDioda(double t, double V, inputRangkaian* params) {
    double I = hitungArus(V, params->V_sumber, params->R, params->I_s, params->n, params->V_t);
    
    return I / params->C;
}

// implementasi rk4 untuk rc
void rungeKutta4_RC(inputRangkaian* params, hasilSimulasi* results, int* result_count) {
    double t = params->t_start;
    double V = params->V0;
    double h = params->step_size;
    double next_output = params->t_start;
    int step_count = 0;
    *result_count = 0;
    
    clock_t start_time = clock();
    
    printf("R=%.0fΩ, C=%.2e F, V_sumber=%.1fV\n", params->R, params->C, params->V_sumber);
    printf("Dioda: I_s=%.2e A, n=%.2f, V_t=%.5f V\n", params->I_s, params->n, params->V_t);
    
    while (t <= params->t_end && *result_count < MAX_TIME_STEPS) {
        if (t >= next_output) {
            double I = hitungArus(V, params->V_sumber, params->R, 
                                             params->I_s, params->n, params->V_t);
            double V_diode = hitungVDioda(I, params->I_s, params->n, params->V_t);
            
            results[*result_count].time = t;
            results[*result_count].V_kapasitor = V;
            results[*result_count].V_diode = V_diode;
            results[*result_count].I_rangkaian = I;
            results[*result_count].step_count = step_count;
            
            clock_t current_time = clock();
            results[*result_count].waktu_ms = 
                (double)(current_time - start_time) * 1000.0 / CLOCKS_PER_SEC;
            
            (*result_count)++;
            next_output += params->output_interval;
        }
        
        // rk4
        double k1 = h * hitungTurunanDioda(t, V, params);
        double k2 = h * hitungTurunanDioda(t + h/2, V + k1/2, params);
        double k3 = h * hitungTurunanDioda(t + h/2, V + k2/2, params);
        double k4 = h * hitungTurunanDioda(t + h, V + k3, params);
        
        V = V + (k1 + 2*k2 + 2*k3 + k4) / 6.0;
        t = t + h;
        step_count++;
        
        // tegangan kapasitor tidak boleh negatif melebihi src nya
        if (V < 0) V = 0;
        if (V > params->V_sumber) V = params->V_sumber;
    }
    
    printf("Simulasi done: %d langkah, %d data output\n", step_count, *result_count);
}

// func untuk read input
int readInput(const char* filename, inputRangkaian* circuits) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("[!] Unable buka file %s\n", filename);
        return 0;
    }
    
    char line[MAX_PANJANG];
    int circuitCount = 0;
    int lineNum = 0;
    
    while (fgets(line, MAX_PANJANG, file) != NULL && circuitCount < MAX_DATA) {
        if (lineNum == 0) {  
            lineNum++;
            continue;
        }
        
        inputRangkaian* circuit = &circuits[circuitCount];
        char* token;
        int tokenIndex = 0;
        
        // parsing input
        token = strtok(line, ",");
        while (token != NULL && tokenIndex < 11) {
            switch (tokenIndex) {
                case 0: circuit->R = atof(token); break;
                case 1: circuit->C = atof(token); break;
                case 2: circuit->V_sumber = atof(token); break;
                case 3: circuit->I_s = atof(token); break;
                case 4: circuit->n = atof(token); break;
                case 5: circuit->V_t = atof(token); break;
                case 6: circuit->V0 = atof(token); break;
                case 7: circuit->t_start = atof(token); break;
                case 8: circuit->t_end = atof(token); break;
                case 9: circuit->step_size = atof(token); break;
                case 10: circuit->output_interval = atof(token); break;
            }
            token = strtok(NULL, ",");
            tokenIndex++;
        }
        
        circuitCount++;
        lineNum++;
    }
    
    fclose(file);
    printf("Berhasil membaca %d rangkaian RC dari %s\n", circuitCount, filename);
    return circuitCount;
}

// func untuk write ke csv
void writeOutput(const char* filename, inputRangkaian* circuits, int circuitCount, hasilSimulasi** allResults, int* resultCounts) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("[!] Unable create file tempat output %s\n", filename);
        return;
    }
    
    fprintf(file, "R,C,V_sumber,I_s,n,V_t,V0,time,V_kapasitor,V_diode,I_rangkaian,step_count,waktu_ms\n");
    
    for (int i = 0; i < circuitCount; i++) {
        inputRangkaian* circuit = &circuits[i];
        hasilSimulasi* results = allResults[i];
        int count = resultCounts[i];
        
        for (int j = 0; j < count; j++) {
            fprintf(file, "%.1f,%.2e,%.2f,%.2e,%.2f,%.5f,%.3f,%.6e,%.6f,%.6f,%.6e,%d,%.3f\n",
                   circuit->R,
                   circuit->C,
                   circuit->V_sumber,
                   circuit->I_s,
                   circuit->n,
                   circuit->V_t,
                   circuit->V0,
                   results[j].time,
                   results[j].V_kapasitor,
                   results[j].V_diode,
                   results[j].I_rangkaian,
                   results[j].step_count,
                   results[j].waktu_ms);
        }
    }
    
    fclose(file);
    printf("Hasil simulasi berhasil disimpan ke %s\n", filename);
}
