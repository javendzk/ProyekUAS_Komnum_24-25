# Analisis Transient Rangkaian RC Dioda dengan Metode Runge-Kutta Order 4

**Proyek UAS Komputasi Numerik 2024/2025**  
**Author: Javana Muhammad Dzaki - NPM 2306161826**

##  Deskripsi Proyek
Metode Runge-Kutta Order 4 (RK4) adalah metode numerik untuk menyelesaikan persamaan diferensial biasa karena memberikan keseimbangan yang baik antara akurasi dan efisiensi komputasi.

Proyek ini mengimplementasikan metode **Runge-Kutta Order 4 (RK4)** pada buku chapra untuk menganalisis respons transient pada rangkaian RC yang mengandung dioda. Program ini dapat mensimulasikan bagaimana tegangan kapasitor dan arus rangkaian berubah seiring waktu ketika ada perubahan mendadak pada sumber tegangan.

##  Latar Belakang 

Dalam rangkaian RC-dioda, kita perlu menyelesaikan persamaan diferensial:

```
dV/dt = (V_sumber - V - V_dioda) / (R × C)
```

Di mana tegangan dioda mengikuti persamaan Shockley:
```
V_dioda = n × V_t × ln((I/I_s) + 1)
```

Metode RK4 digunakan karena memberikan akurasi tinggi (orde ke-4) untuk menyelesaikan persamaan diferensial ini secara numerik.

##  Cara Kerja Program

Program pada proyek UAS Ini menggunakan metode RK4 untuk menyelesaikan persamaan diferensial rangkaian. Proses implementasinya adalah sebagai berikut:

1. **Baca Data**: Program membaca parameter rangkaian dari file CSV yang berisi nilai resistansi, kapasitansi, tegangan sumber, dan parameter dioda.

2. **Inisialisasi**: Set nilai awal tegangan kapasitor = 0V dan mempersiapkan variabel-variabel yang diperlukan untuk simulasi.

3. **Loop Simulasi**: Untuk setiap step waktu, program akan:
   - Menghitung k1, k2, k3, k4 menggunakan rumus RK4
   - Memperbarui tegangan kapasitor dengan rumus: V_baru = V_lama + (k1+2k2+2k3+k4)/6
   - Menghitung arus dan tegangan dioda berdasarkan karakteristik non-linear
   - Menyimpan hasil iterasi ke file output

4. **Output**: Hasil simulasi disimpan ke file CSV untuk analisis dan menbuat grafik


##  Cara Menjalankan

### 1. Compile Program
```bash
gcc -o rk4_circuit_solver rk4_circuit_solver.c -lm
```

### 2. Run Simulasi
```bash
./rk4_circuit_solver
```

##  Input Data

File `data-io/transient_data.csv` berisi parameter rangkaian:
- **R**: Resistansi (Ω)
- **C**: Kapasitansi (F) 
- **V_source**: Tegangan sumber (V)
- **I_s**: Arus saturasi dioda (A)
- **n**: Faktor idealitas dioda
- **V_t**: Tegangan termal (V)
- **V0**: Kondisi awal tegangan kapasitor (V)

**Contoh:**
```
R,C,V_source,I_s,n,V_t,V0,t_start,t_end,step_size,output_interval
1000,1e-6,5.0,1e-12,1.5,0.02585,0.0,0.0,0.001,1e-6,1e-4
2200,2.2e-6,3.3,2e-12,1.2,0.02585,0.5,0.0,0.002,5e-7,1e-4
```

##  Output Data

File `data-io/transient_results.csv` berisi hasil simulasi:
- **time**: Waktu simulasi (s)
- **V_kapasitor**: Tegangan kapasitor (V)
- **V_diode**: Tegangan dioda (V) 
- **I_rangkaian**: Arus rangkaian (A)
- **step_count**: Jumlah langkah simulasi
- **waktu_ms**: Waktu eksekusi (ms)

##  Visualisasi Hasil

### Analisis Rangkaian Individual
![Visualisasi Analisis Rangkaian RC-Dioda #1](https://i.imgur.com/Hc94oDk.png)

Grafik menunjukkan:
- **Respons Transient Tegangan Kapasitor**: Kurva pengisian/pengosongan
- **Tegangan Forward Dioda**: Karakteristik tegangan dioda
- **Arus Rangkaian**: Perubahan arus seiring waktu  
- **Karakteristik I-V Dioda**: Hubungan arus-tegangan dioda

### Perbandingan Multiple Rangkaian
![Perbandingan Respons Tegangan dan Arus](https://i.imgur.com/8qRoE8y.png)

Menampilkan perbandingan respons dari berbagai konfigurasi rangkaian dengan konstanta waktu (τ = R×C) yang berbeda.

##  Analisis yang Dihasilkan

1. **Konstanta Waktu (τ)**: Karakteristik waktu pengisian/pengosongan
2. **Efisiensi Pengisian**: Persentase tegangan final terhadap sumber
3. **Performa Algoritma**: Jumlah iterasi dan waktu komputasi
4. **Karakteristik Dioda**: Analisis tegangan forward dan kurva I-V

