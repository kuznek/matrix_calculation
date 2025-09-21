#include <check.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "../s21_matrix.h"

START_TEST(test_transpose_basic) {
  matrix_t A, result;
  s21_create_matrix(2, 3, &A);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[1][0] = 4.0;
  A.matrix[1][1] = 5.0;
  A.matrix[1][2] = 6.0;

  ck_assert_int_eq(s21_transpose(&A, &result), OK);
  ck_assert_int_eq(result.rows, 3);
  ck_assert_int_eq(result.columns, 2);

  ck_assert_double_eq_tol(result.matrix[0][0], 1.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 4.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 2.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 5.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[2][0], 3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[2][1], 6.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_transpose_square) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 4.0;

  ck_assert_int_eq(s21_transpose(&A, &result), OK);
  ck_assert_int_eq(result.rows, 2);
  ck_assert_int_eq(result.columns, 2);

  ck_assert_double_eq_tol(result.matrix[0][0], 1.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 2.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 4.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_transpose_1x1) {
  matrix_t A, result;
  s21_create_matrix(1, 1, &A);

  A.matrix[0][0] = 42.0;

  ck_assert_int_eq(s21_transpose(&A, &result), OK);
  ck_assert_int_eq(result.rows, 1);
  ck_assert_int_eq(result.columns, 1);
  ck_assert_double_eq_tol(result.matrix[0][0], 42.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_transpose_1x5) {
  matrix_t A, result;
  s21_create_matrix(1, 5, &A);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[0][3] = 4.0;
  A.matrix[0][4] = 5.0;

  ck_assert_int_eq(s21_transpose(&A, &result), OK);
  ck_assert_int_eq(result.rows, 5);
  ck_assert_int_eq(result.columns, 1);

  ck_assert_double_eq_tol(result.matrix[0][0], 1.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 2.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[2][0], 3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[3][0], 4.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[4][0], 5.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_transpose_5x1) {
  matrix_t A, result;
  s21_create_matrix(5, 1, &A);

  A.matrix[0][0] = 10.0;
  A.matrix[1][0] = 20.0;
  A.matrix[2][0] = 30.0;
  A.matrix[3][0] = 40.0;
  A.matrix[4][0] = 50.0;

  ck_assert_int_eq(s21_transpose(&A, &result), OK);
  ck_assert_int_eq(result.rows, 1);
  ck_assert_int_eq(result.columns, 5);

  ck_assert_double_eq_tol(result.matrix[0][0], 10.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 20.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][2], 30.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][3], 40.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][4], 50.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_transpose_3x3_identity) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  // Единичная матрица
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A.matrix[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }

  ck_assert_int_eq(s21_transpose(&A, &result), OK);
  ck_assert_int_eq(result.rows, 3);
  ck_assert_int_eq(result.columns, 3);

  // Транспонированная единичная матрица = единичная матрица
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

START_TEST(test_transpose_3x3_zero) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  // Нулевая матрица
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A.matrix[i][j] = 0.0;
    }
  }

  ck_assert_int_eq(s21_transpose(&A, &result), OK);
  ck_assert_int_eq(result.rows, 3);
  ck_assert_int_eq(result.columns, 3);

  // Транспонированная нулевая матрица = нулевая матрица
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j], 0.0, 1e-6);
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_transpose_negative_values) {
  matrix_t A, result;
  s21_create_matrix(2, 3, &A);

  A.matrix[0][0] = -1.0;
  A.matrix[0][1] = -2.0;
  A.matrix[0][2] = -3.0;
  A.matrix[1][0] = -4.0;
  A.matrix[1][1] = -5.0;
  A.matrix[1][2] = -6.0;

  ck_assert_int_eq(s21_transpose(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], -1.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], -4.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], -2.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], -5.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[2][0], -3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[2][1], -6.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_transpose_fractional_values) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1.5;
  A.matrix[0][1] = 2.7;
  A.matrix[1][0] = 3.14;
  A.matrix[1][1] = 4.999;

  ck_assert_int_eq(s21_transpose(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 1.5, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 3.14, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 2.7, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 4.999, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_transpose_large_values) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1e6;
  A.matrix[0][1] = 2e6;
  A.matrix[1][0] = 3e6;
  A.matrix[1][1] = 4e6;

  ck_assert_int_eq(s21_transpose(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 1e6, 1e-3);
  ck_assert_double_eq_tol(result.matrix[0][1], 3e6, 1e-3);
  ck_assert_double_eq_tol(result.matrix[1][0], 2e6, 1e-3);
  ck_assert_double_eq_tol(result.matrix[1][1], 4e6, 1e-3);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_transpose_small_values) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1e-9;
  A.matrix[0][1] = 2e-9;
  A.matrix[1][0] = 3e-9;
  A.matrix[1][1] = 4e-9;

  ck_assert_int_eq(s21_transpose(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 1e-9, 1e-15);
  ck_assert_double_eq_tol(result.matrix[0][1], 3e-9, 1e-15);
  ck_assert_double_eq_tol(result.matrix[1][0], 2e-9, 1e-15);
  ck_assert_double_eq_tol(result.matrix[1][1], 4e-9, 1e-15);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_transpose_mixed_values) {
  matrix_t A, result;
  s21_create_matrix(3, 2, &A);

  A.matrix[0][0] = 0.0;
  A.matrix[0][1] = -1.5;
  A.matrix[1][0] = 1000.0;
  A.matrix[1][1] = 0.001;
  A.matrix[2][0] = -99.99;
  A.matrix[2][1] = 42.0;

  ck_assert_int_eq(s21_transpose(&A, &result), OK);
  ck_assert_int_eq(result.rows, 2);
  ck_assert_int_eq(result.columns, 3);

  ck_assert_double_eq_tol(result.matrix[0][0], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 1000.0, 1e-3);
  ck_assert_double_eq_tol(result.matrix[0][2], -99.99, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], -1.5, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.001, 1e-9);
  ck_assert_double_eq_tol(result.matrix[1][2], 42.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_transpose_symmetric_matrix) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  // Симметричная матрица
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[1][0] = 2.0;
  A.matrix[1][1] = 4.0;
  A.matrix[1][2] = 5.0;
  A.matrix[2][0] = 3.0;
  A.matrix[2][1] = 5.0;
  A.matrix[2][2] = 6.0;

  ck_assert_int_eq(s21_transpose(&A, &result), OK);

  // Транспонированная симметричная матрица равна исходной
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j], A.matrix[i][j], 1e-6);
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_transpose_double_transpose) {
  matrix_t A, result1, result2;
  s21_create_matrix(2, 3, &A);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[1][0] = 4.0;
  A.matrix[1][1] = 5.0;
  A.matrix[1][2] = 6.0;

  // Первое транспонирование
  ck_assert_int_eq(s21_transpose(&A, &result1), OK);

  // Второе транспонирование
  ck_assert_int_eq(s21_transpose(&result1, &result2), OK);

  // Двойное транспонирование должно вернуть исходную матрицу
  ck_assert_int_eq(result2.rows, A.rows);
  ck_assert_int_eq(result2.columns, A.columns);

  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < A.columns; j++) {
      ck_assert_double_eq_tol(result2.matrix[i][j], A.matrix[i][j], 1e-6);
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&result1);
  s21_remove_matrix(&result2);
}
END_TEST

START_TEST(test_transpose_large_matrix) {
  matrix_t A, result;
  s21_create_matrix(4, 5, &A);

  // Заполняем матрицу последовательными числами
  double value = 1.0;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 5; j++) {
      A.matrix[i][j] = value++;
    }
  }

  ck_assert_int_eq(s21_transpose(&A, &result), OK);
  ck_assert_int_eq(result.rows, 5);
  ck_assert_int_eq(result.columns, 4);

  // Проверяем правильность транспонирования
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 5; j++) {
      ck_assert_double_eq_tol(result.matrix[j][i], A.matrix[i][j], 1e-6);
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

// Тесты на ошибки
START_TEST(test_transpose_null_matrix) {
  matrix_t result;
  ck_assert_int_eq(s21_transpose(NULL, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_transpose_null_result) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  ck_assert_int_eq(s21_transpose(&A, NULL), INCORRECT_MATRIX);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_transpose_null_both) {
  ck_assert_int_eq(s21_transpose(NULL, NULL), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_transpose_invalid_matrix) {
  matrix_t A = {0}, result;
  A.rows = 0;
  A.columns = 2;
  A.matrix = NULL;
  ck_assert_int_eq(s21_transpose(&A, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_transpose_zero_size) {
  matrix_t A = {0}, result;
  A.rows = 0;
  A.columns = 0;
  A.matrix = NULL;
  ck_assert_int_eq(s21_transpose(&A, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_transpose_null_matrix_data) {
  matrix_t A = {0}, result;
  A.rows = 2;
  A.columns = 2;
  A.matrix = NULL;
  ck_assert_int_eq(s21_transpose(&A, &result), INCORRECT_MATRIX);
}
END_TEST

Suite *suite_transpose(void) {
  Suite *s = suite_create("s21_transpose");
  TCase *tc = tcase_create("transpose");

  // Основные функциональные тесты
  tcase_add_test(tc, test_transpose_basic);
  tcase_add_test(tc, test_transpose_square);
  tcase_add_test(tc, test_transpose_1x1);
  tcase_add_test(tc, test_transpose_1x5);
  tcase_add_test(tc, test_transpose_5x1);
  tcase_add_test(tc, test_transpose_large_matrix);

  // Специальные матрицы
  tcase_add_test(tc, test_transpose_3x3_identity);
  tcase_add_test(tc, test_transpose_3x3_zero);
  tcase_add_test(tc, test_transpose_symmetric_matrix);

  // Различные типы значений
  tcase_add_test(tc, test_transpose_negative_values);
  tcase_add_test(tc, test_transpose_fractional_values);
  tcase_add_test(tc, test_transpose_large_values);
  tcase_add_test(tc, test_transpose_small_values);
  tcase_add_test(tc, test_transpose_mixed_values);

  // Математические свойства
  tcase_add_test(tc, test_transpose_double_transpose);

  // Тесты на ошибки
  tcase_add_test(tc, test_transpose_null_matrix);
  tcase_add_test(tc, test_transpose_null_result);
  tcase_add_test(tc, test_transpose_null_both);
  tcase_add_test(tc, test_transpose_invalid_matrix);
  tcase_add_test(tc, test_transpose_zero_size);
  tcase_add_test(tc, test_transpose_null_matrix_data);

  suite_add_tcase(s, tc);
  return s;
}
