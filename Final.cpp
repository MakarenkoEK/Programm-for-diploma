// Попытка1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>

void K(double* K_temp, int i, int g);
void U(double* U_new, double* U_old, double* Y_old, double h, double* K, double* X_temp, int i);
void Y(double* Y_new, double* Y_old, double h, double* U_new, double* U_old, double* X_temp, int i, double K_temp, int g);
double generateGaussian(double mean, double stdDev);
void Y_rung(double* Y_new, double* Y_old, double* X_temp, double K_temp, double* U_new, double* U_old, int M, int i, int g, double h);
double V_temp(double K_temp, double* U_new, double* U_old, double X_temp, double h, double Y_old, int i, int g, int M);
//void U_rung(double Y_old_for_rung, double* U_new, double* K_temp, double X_temp, double* U_old, double h, int i, int g, int M);
double U_temp(double Y_old_for_rung, double K_temp, double X_temp, double h, double U_old, int i, int g, int b);
void K_rung(double* K_temp, int i, int g, int b);
void matrix_metod(double* Y_new, double* Y_old, double* U_new, double* U_old, double K_temp, int N, int j);

int main() { //начинаем мэйн функцию
	srand(time(NULL));

	FILE* p_eul_posit_1 = fopen("1plot_Eul_Kpositive.txt", "w");
	FILE* p_eul_negat_1 = fopen("1plot_Eul_Knegative.txt", "w");
	FILE* p_eul_natur_1 = fopen("1plot_Eul_Knatural.txt", "w");

	FILE* p_rung_posit_2 = fopen("2plot_Rung_Kpositive.txt", "w");
	FILE* p_rung_negat_2 = fopen("2plot_Rung_Knegative.txt", "w");
	FILE* p_rung_natur_2 = fopen("2plot_Rung_Knatural.txt", "w");

	FILE* p_matrix_posit_1 = fopen("2plot_Matrix_Kpositive.txt", "w");
	FILE* p_matrix_negat_1 = fopen("2plot_Matrix_Knegative.txt", "w");
	FILE* p_matrix_natur_1 = fopen("2plot_Matrix_Knatural.txt", "w");

	FILE* f_eul_posit = fopen("output_Eul_Kpositive.txt", "w");
	FILE* f_eul_negat = fopen("output_Eul_Knegative.txt", "w");
	FILE* f_eul_natur = fopen("output_Eul_Knatural.txt", "w");

	FILE* f_rung_posit = fopen("output_Rung_Kpositive.txt", "w");
	FILE* f_rung_negat = fopen("output_Rung_Knegative.txt", "w");
	FILE* f_rung_natur = fopen("output_Rung_Knatural.txt", "w");

	FILE* f_matrix_posit = fopen("output_Matrix_Kpositive.txt", "w");
	FILE* f_matrix_negat = fopen("output_Matrix_Knegative.txt", "w");
	FILE* f_matrix_natur = fopen("output_Matrix_Knatural.txt", "w");

	FILE* p_eul_posit_norma_3 = fopen("3plot_Eul_Kpositive_Gauss.txt", "w");
	FILE* p_eul_negat_norma_3 = fopen("3plot_Eul_Knegative_Gauss.txt", "w");
	FILE* p_eul_natur_norma_3 = fopen("3plot_Eul_Knatural_Gauss.txt", "w");

	FILE* p_rung_posit_norma_4 = fopen("4plot_Rung_Kpositive_Gauss.txt", "w");
	FILE* p_rung_negat_norma_4 = fopen("4plot_Rung_Knegative_Gauss.txt", "w");
	FILE* p_rung_natur_norma_4 = fopen("4plot_Rung_Knatural_Gauss.txt", "w");

	FILE* f_eul_posit_norma = fopen("output_Eul_Kpositive_Gauss.txt", "w");
	FILE* f_eul_negat_norma = fopen("output_Eul_Knegative_Gauss.txt", "w");
	FILE* f_eul_natur_norma = fopen("output_Eul_Knatural_Gauss.txt", "w");

	FILE* f_rung_posit_norma = fopen("output_Rung_Kpositive_Gauss.txt", "w");
	FILE* f_rung_negat_norma = fopen("output_Rung_Knegative_Gauss.txt", "w");
	FILE* f_rung_natur_norma = fopen("output_Rung_Knatural_Gauss.txt", "w");


	FILE* plt_ravn_eul = fopen("k_vs_plt_ravn_eul.txt", "w");
	FILE* plt_norm_eul = fopen("k_vs_plt_norm_eul.txt", "w");

	FILE* plt_ravn_rung = fopen("k_vs_plt_ravn_rung.txt", "w");
	FILE* plt_norm_rung = fopen("k_vs_plt_norm_rung.txt", "w");


	FILE* plt_ravn_matrix = fopen("k_vs_plt_ravn_matrix.txt", "w");

	//обьявляем переменные
	int N = 1000, M = 1000000;      // N - количество делений шага по х M - количество итераций 100 1000000
	int i = 0, j =0;
	int r = 0;
	int R = M / N;
	double rr = 0;
	double X_max = 1000;
	double h = X_max / (M * 1.0);
	double w = X_max / (N * 1.0);
	double y_start = 0;
	double u_start = 1;
	int b = 0;

	//printf("%.4f", h);
	double K_sum = 0;                       // для орицательных К
	double K_avg = 0;

	double* X_temp = new double[M];         //пересчитываемый х
	double* Y_old = new double[M];
	double* U_old = new double[M];
	double* Y_new = new double[M];			//массив новых игреков
	double* U_new = new double[M];			//массив новых юшек
	double* K_temp = new double[M];			//массив чиселок слупа
	double* K_negative = new double[M];     // массив негативной штуки
	double* ln_y = new double[M];           //массив логарифмов

	//чтобы сохранить все лн
	double* ln_y_e_pos = new double[M];
	double* ln_y_e_nat = new double[M];
	double* ln_y_e_neg = new double[M];

	double* ln_y_r_pos = new double[M];
	double* ln_y_r_nat = new double[M];
	double* ln_y_r_neg = new double[M];

	double* ln_y_m_pos = new double[M];
	double* ln_y_m_nat = new double[M];
	double* ln_y_m_neg = new double[M];

	double* ln_y_e_posg = new double[M];
	double* ln_y_e_natg = new double[M];
	double* ln_y_e_negg = new double[M];

	double* ln_y_r_posg = new double[M];
	double* ln_y_r_natg = new double[M];
	double* ln_y_r_negg = new double[M];


	//чтобы сохранить все игреки
	double* y_e_pos = new double[M];
	double* y_e_nat = new double[M];
	double* y_e_neg = new double[M];

	double* y_m_pos = new double[M];
	double* y_m_nat = new double[M];
	double* y_m_neg = new double[M];

	double* y_r_pos = new double[M];
	double* y_r_nat = new double[M];
	double* y_r_neg = new double[M];

	double* y_e_posg = new double[M];
	double* y_e_natg = new double[M];
	double* y_e_negg = new double[M];

	double* y_r_posg = new double[M];
	double* y_r_natg = new double[M];
	double* y_r_negg = new double[M];

	
	////////////////////////////////////////////
	//                                        // 
	//                  СТАРТ                 //
	//                                        //
	////////////////////////////////////////////

	////////////////////////////////////////////
	//  ЭЙЛЕР + РАВНОМЕРНОЕ + ПОЛОЖИТЕЛЬНЫЕ К //
	//                   G=0                  //
	////////////////////////////////////////////
	int g = 0;

	if (g == 0) {
		for (i = 0; i < M; i++) {               //чтобы МАССИВЫ БЫЛИ ПУСТЫМИ
			X_temp[i] = 0.0;
			Y_new[i] = 0.0;
			Y_old[i] = 0.0;
			U_new[i] = 0.0;
			U_old[i] = 0.0;
			K_temp[i] = 0.0;
			ln_y[i] = 0.0;

		}
		for (j = 0; j < N; j++) {
			K_temp[i] = 0.0;
		}
		Y_old[0] = y_start;  //начальные условия
		U_old[0] = u_start;

		for (j = 0; j < N; j++) {
			K(K_temp, j, g);
			K_negative[j] = -K_temp[j];
		}
		for (j = 0; j < N; j++) {
			
			//printf(" % .4f\n", K_temp[j]);
		}

		for (i = 0; i < M; i++) {   // ПЕРВЫЙ РАЗ ВЫЗЫВАЕМ ФУНКЦИЮ ДЛЯ ЯВНОГО МЕТОДА ЭЙЛЕРА
			X_temp[i + 1] = X_temp[i] + h;
			rr = i * 1.0 / (R*1.0);
			//printf("%.4f", rr);
			r = (int)floor(rr);
			//printf("%.d %.4f\n",r, rr);
			Y(Y_new, Y_old, h, U_new, U_old, X_temp, i, K_temp[r], g);
			//printf("%.4f", Y_new);
		}


		fprintf(p_eul_posit_1, " /ln(y)\n"); //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
		for (i = 0; i < M; i++) {
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			if (X_temp[i] != 0) {
				ln_y[i] = log(abs(Y_new[i]));
				ln_y_e_pos[i] = ln_y[i];
				y_e_pos[i] = Y_new[i];
			}
			else ln_y[i] = 0;
			fprintf(p_eul_posit_1, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
			fprintf(f_eul_posit, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i, X_temp[i], i, i, ln_y[i], i, Y_new[i], i, K_temp[r]);
			
		}


		/// /////////////////////////////////////////////////////////////////////////////////

		for (i = 0; i < M; i++) {               //чтобы МАССИВЫ БЫЛИ ПУСТЫМИ
			X_temp[i] = 0.0;
			Y_new[i] = 0.0;
			Y_old[i] = 0.0;
			U_new[i] = 0.0;
			U_old[i] = 0.0;
			
			//K_temp[i] = 0.0;
			ln_y[i] = 0.0;
		}
		Y_old[0] = y_start;  //начальные условия
		U_old[0] = u_start;
		printf("1");
		for (i = 0; i < M; i++) {   // ПЕРВЫЙ РАЗ ВЫЗЫВАЕМ ФУНКЦИЮ ДЛЯ ЯВНОГО МЕТОДА ЭЙЛЕРА
			X_temp[i + 1] = X_temp[i] + h;
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			Y(Y_new, Y_old, h, U_new, U_old, X_temp, i, K_negative[r], g);
			//printf("%.4f", Y_new);
		}

		fprintf(p_eul_negat_1, " /ln(y)\n"); //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
		for (i = 0; i < M; i++) {
			if (X_temp[i] != 0) {
				ln_y[i] = log(abs(Y_new[i]));
				ln_y_e_neg[i] = ln_y[i];
				y_e_neg[i] = Y_new[i];
			}
			else ln_y[i] = 0;
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			fprintf(p_eul_negat_1, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
			fprintf(f_eul_negat, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i, X_temp[i], i, i, ln_y[i], i, Y_new[i], i, K_negative[r]);
			
		}
	}
	///////////////////////////////////////////////////
	//  ЭЙЛЕР + РАВНОМЕРНОЕ + В Т.Ч. ОТРИЦАТЕЛЬНЫЕ К //
	//                   G=1                         //
	///////////////////////////////////////////////////

	// функция расчета среднего К, для того, чтобы сделать выгрузку с отрыцательными К тоже
	//то есть, если из всех К вычесть все К средние, то появятся так же и отрицательные К
	g = 1;
	if (g == 1) {
		for (i = 0; i < M; i++) {         //ОПУСТОШАЕМ МАССИВЫ
			X_temp[i] = 0.0;
				Y_new[i] = 0.0;
				Y_old[i] = 0.0;
				U_new[i] = 0.0;
				U_old[i] = 0.0;
				ln_y[i] = 0.0;
				//K_temp[i] = 0.0;
		}
		Y_old[0] = y_start;
		U_old[0] = u_start;                             //НАЧАЛЬНЫЕ УСЛОВИЯ

			for (i = 0; i < N; i++) {
				K_sum = K_sum + K_temp[i];
			}
		K_avg = K_sum / N;
		for (i = 0; i < N; i++) {
			K_temp[i] = K_temp[i] - K_avg;
		}
		for (i = 0; i < M; i++) {                  //ВЫЗОВ ФУНКЦИИ С ЯВНЫМ МЕТОДОМ ЭЙЛЕРА ВТОРОЙ раз
			X_temp[i + 1] = X_temp[i] + h;
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			Y(Y_new, Y_old, h, U_new, U_old, X_temp, i, K_temp[r], g);
		}

		fprintf(p_eul_natur_1, " /ln(y)\n");
		for (i = 0; i < M; i++) {                 //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
			ln_y[i] = log(abs(Y_new[i]));
			ln_y_e_nat[i] = ln_y[i];
			y_e_nat[i] = Y_new[i];
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			fprintf(p_eul_natur_1, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
			fprintf(f_eul_natur, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i, X_temp[i], i, i, ln_y[i], i, Y_new[i], i, K_temp[r]);
			
		}
	}
	printf("2");

	////////////////////////////////////////////
	//  ЭЙЛЕР + НОРМАЛЬНОЕ  + ПОЛОЖИТЕЛЬНЫЕ К //
	//                   G=2                  //
	////////////////////////////////////////////

	g = 2;
	if (g == 2) {
		for (i = 0; i < M; i++) {                           //ОПУСТОШАЕМ МАССИВЫ
			X_temp[i] = 0.0;
			Y_new[i] = 0.0;
			Y_old[i] = 0.0;
			U_new[i] = 0.0;
			U_old[i] = 0.0;
			K_temp[i] = 0.0;
			ln_y[i] = 0.0;
		}
		Y_old[0] = y_start;                                   //НАЧАЛЬНЫЕ УСЛОВИЯ
		U_old[0] = u_start;

		for (i = 0; i < N; i++) {
			K_temp[i] = 0.0;
		}

		for (i = 0; i < N; i++) {
			K(K_temp, i, g);
			K_negative[i] = -K_temp[i];
		}

		for (i = 0; i < M; i++) {
			X_temp[i + 1] = X_temp[i] + h;             //ВЫЗОВ ФУНКЦИИ С ЯВНЫМ МЕТОДОМ ЭЙЛЕРА ТРЕТИЙ РАЗ
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			Y(Y_new, Y_old, h, U_new, U_old, X_temp, i, K_temp[r], g);
		}
		fprintf(p_eul_posit_norma_3, " /ln(y)\n");

		for (i = 0; i < M; i++) {                    //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
			ln_y[i] = log(abs(Y_new[i]));
			ln_y_e_posg[i] = ln_y[i];
			y_e_posg[i] = Y_new[i];
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			fprintf(p_eul_posit_norma_3, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
			fprintf(f_eul_posit_norma, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i, X_temp[i], i, i, ln_y[i], i, Y_new[i], i, K_temp[r]);
		
		}
		printf("3");
		/// /////////////////////////////////////


		for (i = 0; i < M; i++) {               //чтобы МАССИВЫ БЫЛИ ПУСТЫМИ
			X_temp[i] = 0.0;
			Y_new[i] = 0.0;
			Y_old[i] = 0.0;
			U_new[i] = 0.0;
			U_old[i] = 0.0;
			K_temp[i] = 0.0;
			ln_y[i] = 0.0;
		}
		Y_old[0] = y_start;  //начальные условия
		U_old[0] = u_start;

		for (i = 0; i < M; i++) {   // ПЕРВЫЙ РАЗ ВЫЗЫВАЕМ ФУНКЦИЮ ДЛЯ ЯВНОГО МЕТОДА ЭЙЛЕРА
			X_temp[i + 1] = X_temp[i] + h;
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			Y(Y_new, Y_old, h, U_new, U_old, X_temp, i, K_negative[r], g);
			//printf("%.4f", Y_new);
		}

		fprintf(p_eul_negat_norma_3, " /ln(y)\n"); //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
		for (i = 0; i < M; i++) {
			if (X_temp[i] != 0) {
				ln_y[i] = log(abs(Y_new[i]));
				ln_y_e_negg[i] = ln_y[i];
				y_e_negg[i] = Y_new[i];
			}
			else ln_y[i] = 0;
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			fprintf(p_eul_negat_norma_3, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
			fprintf(f_eul_negat_norma, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i, X_temp[i], i, i, ln_y[i], i, Y_new[i], i, K_negative[r]);
		
		}
	}
	printf("4");
	//////////////////////////////////////////////////
	//  ЭЙЛЕР + НОРМАЛЬНОЕ  + В Т.Ч ОТРИЦАТЕЛЬНЫЕ К //
	//                   G=3                        //
	//////////////////////////////////////////////////

	g = 3;
	if (g == 3) {
		for (i = 0; i < M; i++) { //ОПУСТОШАЕМ МАССИВЫ КРОМЕ К
			X_temp[i] = 0.0;
			Y_new[i] = 0.0;
			Y_old[i] = 0.0;
			U_new[i] = 0.0;
			U_old[i] = 0.0;
			K_temp[i] = 0.0;
			ln_y[i] = 0.0;
		}
		Y_old[0] = y_start;               //НАЧАЛЬНЫЕ УСЛОВИЯ
		U_old[0] = u_start;

		for (i = 0; i < N; i++) {
			K_temp[i] = 0.0;
		}
		//K_sum = 0;                  //ВЕДЕМ РАССЧЕТ СРЕДНЕГО К
		//K_avg = 0;
		//for (i = 0; i < M; i++) {
		//	K_sum = K_sum + K_temp[i];
		//}
		//K_avg = K_sum / M;
		//printf("% .2f / % .2f", K_avg, K_sum);
		//for (i = 0; i < M; i++) {
		//	K_temp[i] = K_temp[i] - K_avg;
			//printf("% .2f ", K_temp[i]);
		//}

		for (i = 0; i < N; i++) {
			K(K_temp, i, g);
		}

		for (i = 0; i < M; i++) {
			X_temp[i + 1] = X_temp[i] + h;
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			Y(Y_new, Y_old, h, U_new, U_old, X_temp, i, K_temp[r], g); //ЗАПУСКАЕМ ФУНКЦИЮ С ЭЛЕРОМ В ЧЕТВЕРТЫЙ РАЗ
		}

		fprintf(p_eul_natur_norma_3, " /ln(y)\n");
		for (i = 0; i < M; i++) {
			ln_y[i] = log(abs(Y_new[i])); //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
			ln_y_e_natg[i] = ln_y[i];
			y_e_natg[i] = Y_new[i];
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			fprintf(p_eul_natur_norma_3, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
			fprintf(f_eul_natur_norma, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i, X_temp[i], i, i, ln_y[i], i, Y_new[i], i, K_temp[r]);
		
		}
	}
	printf("5");

	////////////////////////////////////////////////
	// СТАРТ ИСПОЛЬЗОВАНИЮ ФУНКЦИЙ С РУНГЕ-КУТТОЙ //
	////////////////////////////////////////////////


	////////////////////////////////////////////
	//  РУНГЕ + РАВНОМЕРНОЕ + ПОЛОЖИТЕЛЬНЫЕ К //
	//                   G=0                  //
	////////////////////////////////////////////
	
	g = 0;
	if (g == 0) {
		for (i = 0; i < M; i++) {              //ОПУСТОШАЕМ МАССИВЫ
			X_temp[i] = 0.0;
			Y_new[i] = 0.0;
			Y_old[i] = 0.0;
			U_new[i] = 0.0;
			U_old[i] = 0.0;
			K_temp[i] = 0.0;
			//printf("obnulen\n");
			ln_y[i] = 0.0;
		}

		Y_old[0] = y_start;                     //начальные условия
		U_old[0] = u_start;

		int b = 0;

		for (i = 0; i < N; i++) {
			K_rung(K_temp, i, g, b);
			K_negative[i] = -K_temp[i];
			//printf("k napoln\n");

		}
		b = 1;

		for (i = 0; i < M; i++) {
			X_temp[i + 1] = X_temp[i] + h;
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			//printf("rr %.4f r %.4d\n", rr, i);
			Y_rung(Y_new, Y_old, X_temp, K_temp[r], U_new, U_old, M, i, g, h);  //ВЫЗЫВАЕМ Ф-Ю РУНГЕ ПЕРВЫЙ РАЗ
			//printf("y-rung\n");
		}

		fprintf(p_rung_posit_2, " /ln(y)\n");
		for (i = 0; i < M; i++) {    //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
			ln_y[i] = log(abs(Y_new[i]));
			ln_y_r_pos[i] = ln_y[i];
			y_r_pos[i] = Y_new[i];
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			//printf("in file");
			fprintf(p_rung_posit_2, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
			fprintf(f_rung_posit, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i, X_temp[i], i, i, ln_y[i], i, Y_new[i], i, K_temp[r]);
		
		}

		for (i = 0; i < M; i++) {              //ОПУСТОШАЕМ МАССИВЫ
			X_temp[i] = 0.0;
			Y_new[i] = 0.0;
			Y_old[i] = 0.0;
			U_new[i] = 0.0;
			U_old[i] = 0.0;
			ln_y[i] = 0.0;
		}
		Y_old[0] = y_start;
		U_old[0] = u_start;

		for (i = 0; i < M; i++) {
			X_temp[i + 1] = X_temp[i] + h;
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			Y_rung(Y_new, Y_old, X_temp, K_negative[r], U_new, U_old, M, i, g, h);  //ВЫЗЫВАЕМ Ф-Ю РУНГЕ ПЕРВЫЙ РАЗ
		}

		fprintf(p_rung_negat_2, " /ln(y)\n");
		for (i = 0; i < M; i++) {    //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
			ln_y[i] = log(abs(Y_new[i]));
			ln_y_r_neg[i] = ln_y[i];
			y_r_neg[i] = Y_new[i];
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			fprintf(p_rung_negat_2, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
			fprintf(f_rung_negat, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i, X_temp[i], i, i, ln_y[i], i, Y_new[i], i, K_negative[r]);
		
		}
	}

	printf("6");
	///////////////////////////////////////////////////
	//  РУНГЕ + РАВНОМЕРНОЕ + В Т.Ч. ОТРИЦАТЕЛЬНЫЕ К //
	//                   G=1                         //
	///////////////////////////////////////////////////

	g = 1;
	if (g == 1) {
		for (i = 0; i < M; i++) { //ОПУСТОШАЕМ МАССИВЫ КРОМЕ К
			X_temp[i] = 0.0;
			Y_new[i] = 0.0;
			Y_old[i] = 0.0;
			U_new[i] = 0.0;
			U_old[i] = 0.0;
			ln_y[i] = 0.0;
		}
		Y_old[0] = y_start;
		U_old[0] = u_start;             //НАЧАЛЬНЫЕ УСЛОВИЯ

		K_sum = 0;
		for (i = 0; i < N; i++) {
			printf("bylo %.4f\n", K_temp[i]);
			K_sum = K_sum + K_temp[i];

		}
		K_avg = K_sum / N; 

		printf("%.4f %.4f\n", K_avg, K_sum);
		

		for (i = 0; i < N; i++) {
			K_temp[i] = K_temp[i] - K_avg;
			printf(" stalo %.4f\n", K_temp[i]);
		}

		for (i = 0; i < M; i++) {
			X_temp[i + 1] = X_temp[i] + h;
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			Y_rung(Y_new, Y_old, X_temp, K_temp[r], U_new, U_old, M, i, g, h);  //ВЫЗЫВАЕМ Ф-Ю РУНГЕ ВТОРОЙ РАЗ
		}

		fprintf(p_rung_natur_2, " /ln(y)\n");
		for (i = 0; i < M; i++) {
			ln_y[i] = log(abs(Y_new[i]));    //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
			ln_y_r_nat[i] = ln_y[i];
			y_r_nat[i] = Y_new[i];
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			fprintf(p_rung_natur_2, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
			fprintf(f_rung_natur, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i, X_temp[i], i, i, ln_y[i], i, Y_new[i], i, K_temp[r]);
	
		}
	}

	////////////////////////////////////////////
	//  РУНГЕ + НОРМАЛЬНОЕ  + ПОЛОЖИТЕЛЬНЫЕ К //
	//                   G=2                  //
	////////////////////////////////////////////

	g = 2;
	if (g == 2) {
		for (i = 0; i < M; i++) { //ОПУСТОШАЕМ МАССИВЫ
			X_temp[i] = 0.0;
			Y_new[i] = 0.0;
			Y_old[i] = 0.0;
			U_new[i] = 0.0;
			U_old[i] = 0.0;
			K_temp[i] = 0.0;
			ln_y[i] = 0.0;
		}
		Y_old[0] = y_start;
		U_old[0] = u_start;                       //НАЧАЛЬНЫЕ УСЛОВИЯ
		for (i = 0; i < N; i++) {
			K_temp[i] = 0.0;
		}
		b = 10000;
		for (i = 0; i < N; i++) {
			int b = 0;
			K_rung(K_temp, i, g, b);
			K_negative[i] = -K_temp[i];
		}
		b = 10000;
		for (i = 0; i < M; i++) {
			X_temp[i + 1] = X_temp[i] + h;
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			Y_rung(Y_new, Y_old, X_temp, K_temp[r], U_new, U_old, M, i, g, h);  //ВЫЗЫВАЕМ Ф-Ю РУНГЕ 3 РАЗ
			//printf(" Y= % .4f, U = %.4f\n", Y_old[i], U_old[i]);
		}

		fprintf(p_rung_posit_norma_4, " /ln(y)\n");
		for (i = 0; i < M; i++) {
			ln_y[i] = log(abs(Y_new[i]));     //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
			ln_y_r_posg[i] = ln_y[i];
			y_r_posg[i] = Y_new[i];
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			fprintf(p_rung_posit_norma_4, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
			fprintf(f_rung_posit_norma, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i, X_temp[i], i, i, ln_y[i], i, Y_new[i], i, K_temp[r]);
		}


		for (i = 0; i < M; i++) {              //ОПУСТОШАЕМ МАССИВЫ
			X_temp[i] = 0.0;
			Y_new[i] = 0.0;
			Y_old[i] = 0.0;
			U_new[i] = 0.0;
			U_old[i] = 0.0;
			ln_y[i] = 0.0;
		}
		Y_old[0] = y_start;
		U_old[0] = u_start;

		for (i = 0; i < M; i++) {
			X_temp[i + 1] = X_temp[i] + h;
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			Y_rung(Y_new, Y_old, X_temp, K_negative[r], U_new, U_old, M, i, g, h);  //ВЫЗЫВАЕМ Ф-Ю РУНГЕ ПЕРВЫЙ РАЗ
		}

		fprintf(p_rung_negat_norma_4, " /ln(y)\n");
		for (i = 0; i < M; i++) {    //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
			ln_y[i] = log(abs(Y_new[i]));
			ln_y_r_negg[i] = ln_y[i];
			y_r_negg[i] = Y_new[i];
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			fprintf(p_rung_negat_norma_4, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
			fprintf(f_rung_negat_norma, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i, X_temp[i], i, i, ln_y[i], i, Y_new[i], i, K_negative[r]);
		
		}
	}
	printf("7");
	//////////////////////////////////////////////////
	//  РУНГЕ + НОРМАЛЬНОЕ  + В Т.Ч ОТРИЦАТЕЛЬНЫЕ К //
	//                   G=3                        //
	//////////////////////////////////////////////////

	g = 3;
	if (g == 3) {
		for (i = 0; i < M; i++) {  //ОПУСТОШАЕМ МАССИВЫ КРОМЕ К
			X_temp[i] = 0.0;
			Y_new[i] = 0.0;
			Y_old[i] = 0.0;
			U_new[i] = 0.0;
			U_old[i] = 0.0;
			ln_y[i] = 0.0;
			K_temp[i] = 0.0;
		}
		Y_old[0] = y_start;
		U_old[0] = u_start;
		for (i = 0; i < N; i++) {
			K_temp[i] = 0.0;
		}
		//K_sum = 0;
		//K_avg = 0;
		//for (i = 0; i < M; i++) {
		//	K_sum = K_sum + K_temp[i];
		//}
		//K_avg = K_sum / M;
		//printf("% .2f / % .2f", K_avg, K_sum);
		//for (i = 0; i < M; i++) {
		//	K_temp[i] = K_temp[i] - K_avg;        // СЧИТАЕМ СРЕДНЕЕ
			//printf("% .2f ", K_temp[i]);   
		//}

		for (i = 0; i < N; i++) {
			K_rung(K_temp, i, g, b);
		}
		b = 1;
		for (i = 0; i < M; i++) {
			X_temp[i + 1] = X_temp[i] + h;
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			Y_rung(Y_new, Y_old, X_temp, K_temp[r], U_new, U_old, M, i, g, h); //ВЫЗЫВАЕМ Ф=Ю РУНГЕ 4 РАЗ
		}


		fprintf(p_rung_natur_norma_4, " /ln(y)\n");
		for (i = 0; i < M; i++) {
			ln_y[i] = log(abs(Y_new[i]));            //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
			ln_y_r_natg[i] = ln_y[i];
			y_r_natg[i] = Y_new[i];
			rr = i * 1.0 / (R * 1.0);
			r = (int)floor(rr);
			fprintf(p_rung_natur_norma_4, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
			fprintf(f_rung_natur_norma, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i, X_temp[i], i, i, ln_y[i], i, Y_new[i], i, K_temp[r]);
			
		}
	}
	printf("8");
	
	printf("9");

	/////////////// matrix
	g = 0;

	if (g == 0) {
		X_temp[0] = 0.0;
		for (i = 0; i < M; i++) {               //чтобы МАССИВЫ БЫЛИ ПУСТЫМИ
			X_temp[i + 1] = X_temp[i] + h;
			Y_new[i] = 0.0;
			Y_old[i] = 0.0;
			U_new[i] = 0.0;
			U_old[i] = 0.0;
			K_temp[i] = 0.0;
			ln_y[i] = 0.0;

		}
		for (j = 0; j < N; j++) {
			K_temp[i] = 0.0;
		}
		Y_old[0] = y_start;  //начальные условия
		U_old[0] = u_start;

		for (j = 0; j < N; j++) {
			K(K_temp, j, g);
			K_negative[j] = -K_temp[j];
		}
		for (j = 0; j < N; j++) {
			matrix_metod(Y_new, Y_old, U_new, U_old, K_temp[j], w, j);
		}

		fprintf(p_matrix_posit_1, " /ln(y)\n"); //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
		for (i = 0; i < N; i++) {
			//rr = i * 1.0 / (R * 1.0);
			//r = (int)floor(rr);
			for (j = 0; j < R; j++) {
				if (X_temp[i] != 0) {
					ln_y[i] = log(abs(Y_new[i]));
					ln_y_m_pos[i] = ln_y[i];
					y_m_pos[i] = Y_new[i];
				}
				else ln_y[i] = 0;
				fprintf(p_matrix_posit_1, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
				fprintf(f_matrix_posit, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i*R+j, X_temp[i*R+j], i*R + j, i*R +j, ln_y[i], i*R+j, Y_new[i], i*R+j, K_temp[i]);
			}
		}
	}
	//negative
	if (g == 0) {
		X_temp[0] = 0.0;
		for (i = 0; i < M; i++) {               //чтобы МАССИВЫ БЫЛИ ПУСТЫМИ
			X_temp[i + 1] = X_temp[i] + h;
			Y_new[i] = 0.0;
			Y_old[i] = 0.0;
			U_new[i] = 0.0;
			U_old[i] = 0.0;
			//K_temp[i] = 0.0;
			ln_y[i] = 0.0;

		}
		for (j = 0; j < N; j++) {
			//K_temp[i] = 0.0;
		}
		Y_old[0] = y_start;  //начальные условия
		U_old[0] = u_start;

		for (j = 0; j < N; j++) {
			matrix_metod(Y_new, Y_old, U_new, U_old, K_negative[j], w, j);
		}

		fprintf(p_matrix_negat_1, " /ln(y)\n"); //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
		for (i = 0; i < N; i++) {
			//rr = i * 1.0 / (R * 1.0);
			//r = (int)floor(rr);
			for (j = 0; j < R; j++) {
				if (X_temp[i] != 0) {
					ln_y[i] = log(abs(Y_new[i]));
					ln_y_m_neg[i] = ln_y[i];
					y_m_neg[i] = Y_new[i];
				}
				else ln_y[i] = 0;
				fprintf(p_matrix_negat_1, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
				fprintf(f_matrix_negat, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i * R + j, X_temp[i * R + j], i * R + j, i * R + j, ln_y[i], i * R + j, Y_new[i], i * R + j, K_negative[i]);
			}
		}
	}

	//both 
	if (g == 0) {
		X_temp[0] = 0.0;
		for (i = 0; i < M; i++) {               //чтобы МАССИВЫ БЫЛИ ПУСТЫМИ
			X_temp[i + 1] = X_temp[i] + h;
			Y_new[i] = 0.0;
			Y_old[i] = 0.0;
			U_new[i] = 0.0;
			U_old[i] = 0.0;
			//K_temp[i] = 0.0;
			ln_y[i] = 0.0;

		}
		for (j = 0; j < N; j++) {
			//K_temp[i] = 0.0;
		}
		Y_old[0] = y_start;  //начальные условия
		U_old[0] = u_start;

		K_sum = 0;

		for (i = 0; i < N; i++) {
			K_sum = K_sum + K_temp[i];
		}
		K_avg = K_sum / N;

		for (i = 0; i < N; i++) {
			K_temp[i] = K_temp[i] - K_avg;
		}

		for (j = 0; j < N; j++) {
			matrix_metod(Y_new, Y_old, U_new, U_old, K_temp[j], w, j);
		}

		fprintf(p_matrix_negat_1, " /ln(y)\n"); //ФУНКЦИЯ ЛОГАРИФМОВ И ВЫГРУЗКА В ФАЙЛЫ
		for (i = 0; i < N; i++) {
			//rr = i * 1.0 / (R * 1.0);
			//r = (int)floor(rr);
			for (j = 0; j < R; j++) {
				if (X_temp[i] != 0) {
					ln_y[i] = log(abs(Y_new[i]));
					ln_y_m_nat[i] = ln_y[i];
					y_m_nat[i] = Y_new[i];
				}
				else ln_y[i] = 0;
				fprintf(p_matrix_natur_1, "%.4f/%.0f\n", X_temp[i], ln_y[i]);
				fprintf(f_matrix_natur, "X_%.0d = %.4f, ln(y(x_%.0d)/x_%.0d = %.6f, y(x_%.0d) = %.6f, K(x_%.0d) = %.6f;\n", i * R + j, X_temp[i * R + j], i * R + j, i * R + j, ln_y[i], i * R + j, Y_new[i], i * R + j, K_temp[i]);
			}
		}
	}
	
	for (i = 0; i < M; i++) {
		
		fprintf(plt_ravn_eul, "% .4f % / % .4f / % .4f / %.4f\n", X_temp[i], ln_y_e_pos[i], ln_y_e_nat[i], ln_y_e_neg[i]);

		fprintf(plt_norm_eul, "% .4f % / % .4f / % .4f/ %.4f\n", X_temp[i], ln_y_e_posg[i], ln_y_e_natg[i], ln_y_e_negg[i]);

		fprintf(plt_ravn_rung, "% .4f % / % .4f / % .4f/ %.4f\n", X_temp[i], ln_y_r_pos[i], ln_y_r_nat[i], ln_y_r_neg[i]);

		fprintf(plt_norm_rung, "% .4f % / % .4f / % .4f / %.4f\n", X_temp[i], ln_y_r_posg[i], ln_y_r_natg[i], ln_y_r_negg[i]);
		i = i + 999;
	}
	for (i = 0; i < N; i++) {
		//for (j = 0; j < R; j++) {
			fprintf(plt_ravn_matrix, "% .0d % / % .4f / % .4f / %.4f\n", i, ln_y_m_pos[i], ln_y_m_nat[i], ln_y_m_neg[i]);
		//}
	}
}

void K(double* K_temp, int i, int g) {
	if ((g == 0 || g == 1)) {
		//K_temp[i] = 1;
		K_temp[i] = (rand() % 100000) / (100000 * 1.0);
		//K_temp[i] = rand();
	}
	else if (g == 2 || g == 3) {
		double mean = 0;
		if (g == 2) {
			mean = 1;
		}
		else if (g == 3) {
			mean = 0;
		}
		double stdDev = 0.3;
		K_temp[i] = generateGaussian(mean, stdDev);
	}
}

void Y(double* Y_new, double* Y_old, double h, double* U_new, double* U_old, double* X_temp, int i, double K_temp, int g) {
	
	Y_old[i + 1] = Y_old[i] + U_old[i] * h; // явный метод эйлера 
	U_old[i + 1] = U_old[i] - K_temp * Y_old[i] * h;
	Y_old[i + 1] = Y_old[i] + (U_old[i + 1] + U_old[i]) * h / 2.0;
	U_old[i + 1] = U_old[i] - h * (K_temp * Y_old[i] + K_temp * Y_old[i + 1])/2.0;

	
	U_new[i] = U_old[i+1];
	Y_new[i] = Y_old[i+1];

}


void Y_rung(double* Y_new, double* Y_old, double* X_temp, double K_temp, double* U_new, double* U_old, int M, int i, int g, double h) {
	double K0 = 0.0;
	double K1 = 0.0;
	double K2 = 0.0;
	double K3 = 0.0;
	double Q0 = 0.0;
	double Q1 = 0.0;
	double Q2 = 0.0;
	double Q3 = 0.0; //массивы для рунге 
	double Y_old_for_rung = Y_old[i];
	int b = 0;

	Q0 = U_temp(Y_old_for_rung, K_temp, X_temp[i], h, Y_old[i], i, g, b);
	K0 = V_temp(K_temp, U_new, U_old, X_temp[i], h, U_old[i], i, g, M);

	Q1 = U_temp(Y_old_for_rung, K_temp, X_temp[i], h, Y_old[i] + h * K0 / 2, i, g, b);
	K1 = V_temp(K_temp, U_new, U_old, X_temp[i] + h / 2.0, h, U_old[i] + h * Q0 / 2.0, i, g, M);

	Q2 = U_temp(Y_old_for_rung, K_temp, X_temp[i], h, Y_old[i] + h * K1 / 2, i, g, b);
	K2 = V_temp(K_temp, U_new, U_old, X_temp[i] + h / 2, h, U_old[i] + h * Q1 / 2, i, g, M);

	Q3 = U_temp(Y_old_for_rung, K_temp, X_temp[i], h, Y_old[i] + h * K2, i, g, b);
	K3 = V_temp(K_temp, U_new, U_old, X_temp[i] + h, h, U_old[i] + h * Q2, i, g, M);

	U_old[i + 1] = U_old[i] + h * (Q0 + 2 * Q1 + 2 * Q2 + Q3) / 6;
	Y_old[i + 1] = Y_old[i] + h * (K0 + 2 * K1 + 2 * K2 + K3) / 6;

	Y_new[i] = Y_old[i];

}

double V_temp(double K_temp, double* U_new, double* U_old, double X_temp, double h, double Y_old, int i, int g, int M) {
	double K = Y_old;
	return K;
}

//FU(Y_old_for_rung, K_temp[i], X_temp[i], h, U_old[i] + Y4[i] * h, i, g, b);
double U_temp(double Y_old_for_rung, double K_temp, double X_temp, double h, double FFF, int i, int g, int b) {

	double KKK = K_temp;
	return -KKK * FFF;
}


void K_rung(double* K_temp, int i, int g, int b) {
	K(K_temp, i, g);
	//printf("%.4f",b);
}

double generateGaussian(double mean, double stdDev) {
	static double spare;
	static bool hasSpare = false;
	if (hasSpare) {
		hasSpare = false;
		return spare * stdDev + mean;
	}
	else {
		double u, v, s;
		do {
			u = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
			v = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
			s = u * u + v * v;
		} while (s >= 1.0 || s == 0.0);
		s = sqrt(-2.0 * log(s) / s);
		spare = v * s;
		hasSpare = true;
		return mean + stdDev * u * s;
	}
}

void matrix_metod(double* Y_new, double* Y_old, double* U_new, double* U_old, double K_temp, int w, int i) {
	//double N = 1.0 / w;
	int N = w;
	//printf("here");
	if (K_temp >= 0) {
		Y_old[i + 1] = Y_old[i] * cos(sqrt(K_temp)*N) + U_old[i] * (1 / (sqrt(K_temp))) * sin(sqrt(K_temp) * N) ;
		U_old[i + 1] = -Y_old[i] * (sqrt(K_temp)*N) * sin(sqrt(K_temp)*N) + U_old[i]* N * cos(sqrt(K_temp)* N) ;
		printf("%.4f %.4f\n", Y_old[i], U_old[i]);
		Y_new[i + 1] = Y_old[i + 1];
		U_new[i + 1] = U_old[i + 1];
	}
	else {
		Y_old[i+1] = Y_old[i] * cosh(sqrt(-K_temp)* N) + U_old[i] *(1/ (sqrt(-K_temp))) *sinh(sqrt(-K_temp)* N) ;
		U_old[i+1] = Y_old[i] * (sqrt(-K_temp)*N) * sinh(sqrt(-K_temp)*N) + U_old[i] * N * cosh(sqrt(-K_temp)* N) ;
		Y_new[i + 1] = Y_old[i + 1];
		U_new[i + 1] = U_old[i + 1];
	}
}

