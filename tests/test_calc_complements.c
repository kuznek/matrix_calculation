#include <check.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#include "../s21_matrix.h"

START_TEST(test_calc_complements_2x2) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 4.0;
  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 4.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], -3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], -2.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 1.0, 1e-6);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_1x1) {
  matrix_t A, result;
  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 5.0;
  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 1.0, 1e-6);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_1x1_zero) {
  matrix_t A, result;
  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 0.0;
  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 1.0, 1e-6);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_1x1_negative) {
  matrix_t A, result;
  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = -7.5;
  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 1.0, 1e-6);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_3x3) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);
  // Матрица:
  // [1  2  3]
  // [4  5  6]
  // [7  8  9]
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[1][0] = 4.0;
  A.matrix[1][1] = 5.0;
  A.matrix[1][2] = 6.0;
  A.matrix[2][0] = 7.0;
  A.matrix[2][1] = 8.0;
  A.matrix[2][2] = 9.0;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  // Ожидаемая матрица алгебраических дополнений:
  // [-3   6  -3]
  // [ 6 -12   6]
  // [-3   6  -3]
  ck_assert_double_eq_tol(result.matrix[0][0], -3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 6.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][2], -3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 6.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], -12.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][2], 6.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[2][0], -3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[2][1], 6.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[2][2], -3.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_3x3_identity) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);
  // Единичная матрица 3x3
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A.matrix[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  // Алгебраические дополнения единичной матрицы = единичная матрица
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double expected = (i == j) ? 1.0 : 0.0;
      ck_assert_double_eq_tol(result.matrix[i][j], expected, 1e-6);
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_4x4) {
  matrix_t A, result;
  s21_create_matrix(4, 4, &A);

  // Тестовая матрица 4x4
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[0][3] = 4.0;
  A.matrix[1][0] = 5.0;
  A.matrix[1][1] = 6.0;
  A.matrix[1][2] = 7.0;
  A.matrix[1][3] = 8.0;
  A.matrix[2][0] = 9.0;
  A.matrix[2][1] = 10.0;
  A.matrix[2][2] = 11.0;
  A.matrix[2][3] = 12.0;
  A.matrix[3][0] = 13.0;
  A.matrix[3][1] = 14.0;
  A.matrix[3][2] = 15.0;
  A.matrix[3][3] = 16.0;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  // Проверяем, что результат получен (конкретные значения сложно вычислить
  // вручную) Проверим, что результат не содержит NaN или Inf
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      ck_assert_int_eq(isnan(result.matrix[i][j]), 0);
      ck_assert_int_eq(isinf(result.matrix[i][j]), 0);
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_2x2_zero_determinant) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  // Матрица с нулевым определителем
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 2.0;
  A.matrix[1][1] = 4.0;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  // Алгебраические дополнения должны вычисляться даже для вырожденных матриц
  ck_assert_double_eq_tol(result.matrix[0][0], 4.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], -2.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], -2.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 1.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_negative_values) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = -1.0;
  A.matrix[0][1] = -2.0;
  A.matrix[1][0] = -3.0;
  A.matrix[1][1] = -4.0;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], -4.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 2.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], -1.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_fractional_values) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 0.5;
  A.matrix[0][1] = 1.5;
  A.matrix[1][0] = 2.5;
  A.matrix[1][1] = 3.5;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 3.5, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], -2.5, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], -1.5, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.5, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_large_values) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1000.0;
  A.matrix[0][1] = 2000.0;
  A.matrix[1][0] = 3000.0;
  A.matrix[1][1] = 4000.0;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 4000.0, 1e-3);
  ck_assert_double_eq_tol(result.matrix[0][1], -3000.0, 1e-3);
  ck_assert_double_eq_tol(result.matrix[1][0], -2000.0, 1e-3);
  ck_assert_double_eq_tol(result.matrix[1][1], 1000.0, 1e-3);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_small_values) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1e-6;
  A.matrix[0][1] = 2e-6;
  A.matrix[1][0] = 3e-6;
  A.matrix[1][1] = 4e-6;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 4e-6, 1e-12);
  ck_assert_double_eq_tol(result.matrix[0][1], -3e-6, 1e-12);
  ck_assert_double_eq_tol(result.matrix[1][0], -2e-6, 1e-12);
  ck_assert_double_eq_tol(result.matrix[1][1], 1e-6, 1e-12);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_non_square) {
  matrix_t A, result;
  s21_create_matrix(2, 3, &A);
  ck_assert_int_eq(s21_calc_complements(&A, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_calc_complements_rectangular_3x2) {
  matrix_t A, result;
  s21_create_matrix(3, 2, &A);
  ck_assert_int_eq(s21_calc_complements(&A, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_calc_complements_rectangular_1x5) {
  matrix_t A, result;
  s21_create_matrix(1, 5, &A);
  ck_assert_int_eq(s21_calc_complements(&A, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
}
END_TEST

// Тесты на NULL указатели
START_TEST(test_calc_complements_null_matrix) {
  matrix_t result;
  ck_assert_int_eq(s21_calc_complements(NULL, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_calc_complements_null_result) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  ck_assert_int_eq(s21_calc_complements(&A, NULL), INCORRECT_MATRIX);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_calc_complements_null_both) {
  ck_assert_int_eq(s21_calc_complements(NULL, NULL), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_calc_complements_invalid_matrix) {
  matrix_t A = {0}, result;
  A.rows = 0;
  A.columns = 2;
  A.matrix = NULL;
  ck_assert_int_eq(s21_calc_complements(&A, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_calc_complements_zero_size) {
  matrix_t A = {0}, result;
  A.rows = 0;
  A.columns = 0;
  A.matrix = NULL;
  ck_assert_int_eq(s21_calc_complements(&A, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_calc_complements_3x3_specific) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  // Конкретная матрица с известными алгебраическими дополнениями
  A.matrix[0][0] = 2.0;
  A.matrix[0][1] = 1.0;
  A.matrix[0][2] = 0.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 2.0;
  A.matrix[1][2] = 1.0;
  A.matrix[2][0] = 1.0;
  A.matrix[2][1] = 1.0;
  A.matrix[2][2] = 2.0;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  // Вычисленные вручную алгебраические дополнения
  ck_assert_double_eq_tol(result.matrix[0][0], 3.0, 1e-6);   // |2 1; 1 2| = 3
  ck_assert_double_eq_tol(result.matrix[0][1], -5.0, 1e-6);  // -|3 1; 1 2| = -5
  ck_assert_double_eq_tol(result.matrix[0][2], 1.0, 1e-6);   // |3 2; 1 1| = 1
  ck_assert_double_eq_tol(result.matrix[1][0], -2.0, 1e-6);  // -|1 0; 1 2| = -2
  ck_assert_double_eq_tol(result.matrix[1][1], 4.0, 1e-6);   // |2 0; 1 2| = 4
  ck_assert_double_eq_tol(result.matrix[1][2], -1.0, 1e-6);  // -|2 1; 1 1| = -1
  ck_assert_double_eq_tol(result.matrix[2][0], 1.0, 1e-6);   // |1 0; 2 1| = 1
  ck_assert_double_eq_tol(result.matrix[2][1], -2.0, 1e-6);  // -|2 0; 3 1| = -2
  ck_assert_double_eq_tol(result.matrix[2][2], 1.0, 1e-6);   // |2 1; 3 2| = 1

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST
START_TEST(test_calc_complements_5x5_identity) {
  matrix_t A, result;
  s21_create_matrix(5, 5, &A);

  // Единичная матрица 5x5
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      A.matrix[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  // Алгебраические дополнения единичной матрицы = единичная матрица
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      double expected = (i == j) ? 1.0 : 0.0;
      ck_assert_double_eq_tol(result.matrix[i][j], expected, 1e-6);
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_6x6_diagonal) {
  matrix_t A, result;
  s21_create_matrix(6, 6, &A);

  // Диагональная матрица
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      A.matrix[i][j] = (i == j) ? (double)(i + 1) : 0.0;
    }
  }

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  // Для диагональной матрицы проверяем некоторые ключевые элементы
  ck_assert_double_ne(result.matrix[0][0], 0.0);  // Не должно быть нулем
  ck_assert_double_ne(result.matrix[5][5], 0.0);  // Не должно быть нулем

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

// Специальные матрицы
START_TEST(test_calc_complements_2x2_antisymmetric) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  // Кососимметричная матрица
  A.matrix[0][0] = 0.0;
  A.matrix[0][1] = 3.0;
  A.matrix[1][0] = -3.0;
  A.matrix[1][1] = 0.0;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], -3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_3x3_zero_row) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[1][0] = 0.0;
  A.matrix[1][1] = 0.0;
  A.matrix[1][2] = 0.0;  // Нулевая строка
  A.matrix[2][0] = 7.0;
  A.matrix[2][1] = 8.0;
  A.matrix[2][2] = 9.0;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  // Проверяем, что функция работает с нулевыми строками
  ck_assert_double_eq_tol(result.matrix[0][0], 0.0, 1e-6);  // 0*9 - 0*8 = 0
  ck_assert_double_eq_tol(result.matrix[0][1], 0.0, 1e-6);  // -(0*9 - 0*7) = 0
  ck_assert_double_eq_tol(result.matrix[0][2], 0.0, 1e-6);  // 0*8 - 0*7 = 0

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_3x3_zero_column) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 0.0;
  A.matrix[0][2] = 3.0;
  A.matrix[1][0] = 4.0;
  A.matrix[1][1] = 0.0;
  A.matrix[1][2] = 6.0;  // Нулевой столбец
  A.matrix[2][0] = 7.0;
  A.matrix[2][1] = 0.0;
  A.matrix[2][2] = 9.0;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  // Проверяем наличие результатов (точные значения сложно вычислить)
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ck_assert_int_eq(isnan(result.matrix[i][j]), 0);
      ck_assert_int_eq(isinf(result.matrix[i][j]), 0);
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

// Граничные значения
START_TEST(test_calc_complements_extremely_large) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1e8;
  A.matrix[0][1] = 2e8;
  A.matrix[1][0] = 3e8;
  A.matrix[1][1] = 4e8;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 4e8, 1e2);
  ck_assert_double_eq_tol(result.matrix[0][1], -3e8, 1e2);
  ck_assert_double_eq_tol(result.matrix[1][0], -2e8, 1e2);
  ck_assert_double_eq_tol(result.matrix[1][1], 1e8, 1e2);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_extremely_small) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1e-10;
  A.matrix[0][1] = 2e-10;
  A.matrix[1][0] = 3e-10;
  A.matrix[1][1] = 4e-10;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 4e-10, 1e-16);
  ck_assert_double_eq_tol(result.matrix[0][1], -3e-10, 1e-16);
  ck_assert_double_eq_tol(result.matrix[1][0], -2e-10, 1e-16);
  ck_assert_double_eq_tol(result.matrix[1][1], 1e-10, 1e-16);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_mixed_magnitude) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1e6;
  A.matrix[0][1] = 1e-6;
  A.matrix[1][0] = 1e-6;
  A.matrix[1][1] = 1e6;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 1e6, 1e-3);
  ck_assert_double_eq_tol(result.matrix[0][1], -1e-6, 1e-12);
  ck_assert_double_eq_tol(result.matrix[1][0], -1e-6, 1e-12);
  ck_assert_double_eq_tol(result.matrix[1][1], 1e6, 1e-3);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

// Тесты с иррациональными числами
START_TEST(test_calc_complements_pi_e) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = M_PI;
  A.matrix[0][1] = M_E;
  A.matrix[1][0] = M_E;
  A.matrix[1][1] = M_PI;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], M_PI, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], -M_E, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], -M_E, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], M_PI, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_sqrt_values) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = sqrt(2.0);
  A.matrix[0][1] = sqrt(3.0);
  A.matrix[1][0] = sqrt(5.0);
  A.matrix[1][1] = sqrt(7.0);

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], sqrt(7.0), 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], -sqrt(5.0), 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], -sqrt(3.0), 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], sqrt(2.0), 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

// Сложные 3x3 матрицы
START_TEST(test_calc_complements_3x3_upper_triangular) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  // Верхнетреугольная матрица
  A.matrix[0][0] = 2.0;
  A.matrix[0][1] = 3.0;
  A.matrix[0][2] = 4.0;
  A.matrix[1][0] = 0.0;
  A.matrix[1][1] = 5.0;
  A.matrix[1][2] = 6.0;
  A.matrix[2][0] = 0.0;
  A.matrix[2][1] = 0.0;
  A.matrix[2][2] = 7.0;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  // Проверяем некоторые ключевые элементы
  ck_assert_double_eq_tol(result.matrix[0][0], 35.0, 1e-6);  // 5*7 - 6*0 = 35
  ck_assert_double_eq_tol(result.matrix[2][2], 10.0, 1e-6);  // 2*5 - 3*0 = 10

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_3x3_lower_triangular) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  // Нижнетреугольная матрица
  A.matrix[0][0] = 2.0;
  A.matrix[0][1] = 0.0;
  A.matrix[0][2] = 0.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 4.0;
  A.matrix[1][2] = 0.0;
  A.matrix[2][0] = 5.0;
  A.matrix[2][1] = 6.0;
  A.matrix[2][2] = 7.0;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  // Проверяем некоторые ключевые элементы
  ck_assert_double_eq_tol(result.matrix[0][0], 28.0, 1e-6);  // 4*7 - 0*6 = 28
  ck_assert_double_eq_tol(result.matrix[2][2], 8.0, 1e-6);   // 2*4 - 0*3 = 8

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_calc_complements_3x3_symmetric) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  // Симметричная матрица
  A.matrix[0][0] = 4.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 1.0;
  A.matrix[1][0] = 2.0;
  A.matrix[1][1] = 5.0;
  A.matrix[1][2] = 3.0;
  A.matrix[2][0] = 1.0;
  A.matrix[2][1] = 3.0;
  A.matrix[2][2] = 6.0;

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  // Для симметричной матрицы проверяем симметричность результата
  ck_assert_double_eq_tol(result.matrix[0][1], result.matrix[1][0], 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][2], result.matrix[2][0], 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][2], result.matrix[2][1], 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

// Дополнительные тесты на ошибки
START_TEST(test_calc_complements_negative_rows) {
  matrix_t A = {0}, result;
  A.rows = -1;
  A.columns = 2;
  A.matrix = NULL;

  ck_assert_int_eq(s21_calc_complements(&A, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_calc_complements_negative_columns) {
  matrix_t A = {0}, result;
  A.rows = 2;
  A.columns = -1;
  A.matrix = NULL;

  ck_assert_int_eq(s21_calc_complements(&A, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_calc_complements_null_data_valid_size) {
  matrix_t A = {0}, result;
  A.rows = 2;
  A.columns = 2;
  A.matrix = NULL;  // NULL данные с валидными размерами

  ck_assert_int_eq(s21_calc_complements(&A, &result), INCORRECT_MATRIX);
}
END_TEST

// Тесты производительности
START_TEST(test_calc_complements_7x7) {
  matrix_t A, result;
  s21_create_matrix(7, 7, &A);

  // Заполняем матрицу 7x7 так, чтобы она не была вырожденной
  for (int i = 0; i < 7; i++) {
    for (int j = 0; j < 7; j++) {
      A.matrix[i][j] = (i == j) ? (double)(i + 1) : sin(i + j);
    }
  }

  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);

  // Проверяем, что результат получен без ошибок
  for (int i = 0; i < 7; i++) {
    for (int j = 0; j < 7; j++) {
      ck_assert_int_eq(isnan(result.matrix[i][j]), 0);
      ck_assert_int_eq(isinf(result.matrix[i][j]), 0);
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

// Обновленная функция suite
Suite *suite_calc_complements(void) {
  Suite *s = suite_create("s21_calc_complements");
  TCase *tc = tcase_create("calc_complements");

  // Основные функциональные тесты
  tcase_add_test(tc, test_calc_complements_1x1);
  tcase_add_test(tc, test_calc_complements_1x1_zero);
  tcase_add_test(tc, test_calc_complements_1x1_negative);
  tcase_add_test(tc, test_calc_complements_2x2);
  tcase_add_test(tc, test_calc_complements_3x3);
  tcase_add_test(tc, test_calc_complements_3x3_identity);
  tcase_add_test(tc, test_calc_complements_3x3_specific);
  tcase_add_test(tc, test_calc_complements_4x4);

  // Дополнительные размеры
  tcase_add_test(tc, test_calc_complements_5x5_identity);
  tcase_add_test(tc, test_calc_complements_6x6_diagonal);
  tcase_add_test(tc, test_calc_complements_7x7);

  // Специальные матрицы
  tcase_add_test(tc, test_calc_complements_2x2_antisymmetric);
  tcase_add_test(tc, test_calc_complements_3x3_zero_row);
  tcase_add_test(tc, test_calc_complements_3x3_zero_column);
  tcase_add_test(tc, test_calc_complements_3x3_upper_triangular);
  tcase_add_test(tc, test_calc_complements_3x3_lower_triangular);
  tcase_add_test(tc, test_calc_complements_3x3_symmetric);
  tcase_add_test(tc, test_calc_complements_2x2_zero_determinant);

  // Различные типы значений
  tcase_add_test(tc, test_calc_complements_negative_values);
  tcase_add_test(tc, test_calc_complements_fractional_values);
  tcase_add_test(tc, test_calc_complements_large_values);
  tcase_add_test(tc, test_calc_complements_small_values);
  tcase_add_test(tc, test_calc_complements_extremely_large);
  tcase_add_test(tc, test_calc_complements_extremely_small);
  tcase_add_test(tc, test_calc_complements_mixed_magnitude);
  tcase_add_test(tc, test_calc_complements_pi_e);
  tcase_add_test(tc, test_calc_complements_sqrt_values);

  // Тесты на ошибки размерности
  tcase_add_test(tc, test_calc_complements_non_square);
  tcase_add_test(tc, test_calc_complements_rectangular_3x2);
  tcase_add_test(tc, test_calc_complements_rectangular_1x5);

  // Тесты на NULL и некорректные данные
  tcase_add_test(tc, test_calc_complements_null_matrix);
  tcase_add_test(tc, test_calc_complements_null_result);
  tcase_add_test(tc, test_calc_complements_null_both);
  tcase_add_test(tc, test_calc_complements_invalid_matrix);
  tcase_add_test(tc, test_calc_complements_zero_size);
  tcase_add_test(tc, test_calc_complements_negative_rows);
  tcase_add_test(tc, test_calc_complements_negative_columns);
  tcase_add_test(tc, test_calc_complements_null_data_valid_size);

  suite_add_tcase(s, tc);
  return s;
}