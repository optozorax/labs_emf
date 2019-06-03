#pragma once

/** Определения типов. */
#include <functional>
#include <chrono>
#include <mutex>
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace placeholders;

/* Для плотных матриц и векторов используется библиотека Eigen. */
typedef Eigen::MatrixXd matrix_t; /// Плотная матрица
typedef Eigen::VectorXd vector_t; /// Вектор

/* Тип 1D, 2D, 3D функций. */
typedef function<double(double)> function_1d_t; 
typedef function<double(double, double)> function_2d_t;
typedef function<double(double, double, double)> function_3d_t;

/** Считает время выполнения функции f в микросекундах. */
inline double calc_time_microseconds(const function<void(void)>& f) {
	using namespace chrono;
	auto start = high_resolution_clock::now();
	f();
	auto end = high_resolution_clock::now();
	return duration_cast<microseconds>(end - start).count();;
}

/** Выводит на экран процент завершенной работы. Использует мьютексы для защиты cout при использовании несколькими потоками */
inline void write_percent(double percent) {
	static mutex m;
	lock_guard<mutex> g(m);
	cout << "\r" << setprecision(2) << fixed << setw(5) << percent * 100 << "%";
}