#pragma once

/** Определения типов. */
#include <functional>
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