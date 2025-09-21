#include <check.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "../s21_matrix.h"

START_TEST(test_sum_matrix_basic) {
  matrix_t A, B, result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 4.0;

  B.matrix[0][0] = 5.0;
  B.matrix[0][1] = 6.0;
  B.matrix[1][0] = 7.0;
  B.matrix[1][1] = 8.0;

  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 6.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 8.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 10.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 12.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_sum_matrix_1x1) {
  matrix_t A, B, result;
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &B);

  A.matrix[0][0] = 3.5;
  B.matrix[0][0] = 2.5;

  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 6.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_sum_matrix_3x3) {
  matrix_t A, B, result;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &B);

  // Заполняем матрицы
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A.matrix[i][j] = i * 3 + j + 1;
      B.matrix[i][j] = (i * 3 + j + 1) * 10;
    }
  }

  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), OK);

  // Проверяем результат
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double expected = A.matrix[i][j] + B.matrix[i][j];
      ck_assert_double_eq_tol(result.matrix[i][j], expected, 1e-6);
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_sum_matrix_rectangular) {
  matrix_t A, B, result;
  s21_create_matrix(2, 4, &A);
  s21_create_matrix(2, 4, &B);

  // Матрица A
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[0][3] = 4.0;
  A.matrix[1][0] = 5.0;
  A.matrix[1][1] = 6.0;
  A.matrix[1][2] = 7.0;
  A.matrix[1][3] = 8.0;

  // Матрица B
  B.matrix[0][0] = 0.1;
  B.matrix[0][1] = 0.2;
  B.matrix[0][2] = 0.3;
  B.matrix[0][3] = 0.4;
  B.matrix[1][0] = 0.5;
  B.matrix[1][1] = 0.6;
  B.matrix[1][2] = 0.7;
  B.matrix[1][3] = 0.8;

  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 1.1, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 2.2, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][2], 3.3, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][3], 4.4, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 5.5, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 6.6, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][2], 7.7, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][3], 8.8, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_sum_matrix_negative_values) {
  matrix_t A, B, result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = -1.0;
  A.matrix[0][1] = -2.0;
  A.matrix[1][0] = -3.0;
  A.matrix[1][1] = -4.0;

  B.matrix[0][0] = 1.0;
  B.matrix[0][1] = 2.0;
  B.matrix[1][0] = 3.0;
  B.matrix[1][1] = 4.0;

  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_sum_matrix_zero_matrix) {
  matrix_t A, B, result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 5.0;
  A.matrix[0][1] = 10.0;
  A.matrix[1][0] = 15.0;
  A.matrix[1][1] = 20.0;

  // B - нулевая матрица (calloc автоматически инициализирует нулями)

  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 5.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 10.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 15.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 20.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_sum_matrix_large_values) {
  matrix_t A, B, result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1e6;
  A.matrix[0][1] = 2e6;
  A.matrix[1][0] = 3e6;
  A.matrix[1][1] = 4e6;

  B.matrix[0][0] = 5e6;
  B.matrix[0][1] = 6e6;
  B.matrix[1][0] = 7e6;
  B.matrix[1][1] = 8e6;

  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 6e6, 1e-3);
  ck_assert_double_eq_tol(result.matrix[0][1], 8e6, 1e-3);
  ck_assert_double_eq_tol(result.matrix[1][0], 10e6, 1e-3);
  ck_assert_double_eq_tol(result.matrix[1][1], 12e6, 1e-3);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_sum_matrix_small_values) {
  matrix_t A, B, result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1e-9;
  A.matrix[0][1] = 2e-9;
  A.matrix[1][0] = 3e-9;
  A.matrix[1][1] = 4e-9;

  B.matrix[0][0] = 5e-9;
  B.matrix[0][1] = 6e-9;
  B.matrix[1][0] = 7e-9;
  B.matrix[1][1] = 8e-9;

  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 6e-9, 1e-15);
  ck_assert_double_eq_tol(result.matrix[0][1], 8e-9, 1e-15);
  ck_assert_double_eq_tol(result.matrix[1][0], 10e-9, 1e-15);
  ck_assert_double_eq_tol(result.matrix[1][1], 12e-9, 1e-15);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_sum_matrix_mixed_values) {
  matrix_t A, B, result;
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(2, 3, &B);

  A.matrix[0][0] = 0.0;
  A.matrix[0][1] = -100.5;
  A.matrix[0][2] = 1e6;
  A.matrix[1][0] = 1e-9;
  A.matrix[1][1] = 42.0;
  A.matrix[1][2] = -3.14;

  B.matrix[0][0] = 1.0;
  B.matrix[0][1] = 100.5;
  B.matrix[0][2] = -1e6;
  B.matrix[1][0] = -1e-9;
  B.matrix[1][1] = -42.0;
  B.matrix[1][2] = 3.14;

  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 1.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][2], 0.0, 1e-3);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.0, 1e-15);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][2], 0.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_sum_matrix_commutative) {
  matrix_t A, B, result1, result2;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1.5;
  A.matrix[0][1] = 2.7;
  A.matrix[1][0] = 3.9;
  A.matrix[1][1] = 4.1;

  B.matrix[0][0] = 9.8;
  B.matrix[0][1] = 7.6;
  B.matrix[1][0] = 5.4;
  B.matrix[1][1] = 3.2;

  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result1), OK);
  ck_assert_int_eq(s21_sum_matrix(&B, &A, &result2), OK);

  // Проверяем коммутативность A + B = B + A
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      ck_assert_double_eq_tol(result1.matrix[i][j], result2.matrix[i][j], 1e-6);
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result1);
  s21_remove_matrix(&result2);
}
END_TEST

// Тесты на ошибки
START_TEST(test_sum_matrix_different_size) {
  matrix_t A, B, result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(3, 3, &B);
  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_sum_matrix_different_rows) {
  matrix_t A, B, result;
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(3, 3, &B);
  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_sum_matrix_different_columns) {
  matrix_t A, B, result;
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(3, 3, &B);
  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_sum_matrix_null_first) {
  matrix_t B, result;
  s21_create_matrix(2, 2, &B);
  ck_assert_int_eq(s21_sum_matrix(NULL, &B, &result), INCORRECT_MATRIX);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_sum_matrix_null_second) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);
  ck_assert_int_eq(s21_sum_matrix(&A, NULL, &result), INCORRECT_MATRIX);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_sum_matrix_null_result) {
  matrix_t A, B;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  ck_assert_int_eq(s21_sum_matrix(&A, &B, NULL), INCORRECT_MATRIX);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_sum_matrix_null_all) {
  ck_assert_int_eq(s21_sum_matrix(NULL, NULL, NULL), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_sum_matrix_invalid_first) {
  matrix_t A = {0}, B, result;
  A.rows = 0;
  A.columns = 2;
  A.matrix = NULL;
  s21_create_matrix(2, 2, &B);
  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), INCORRECT_MATRIX);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_sum_matrix_invalid_second) {
  matrix_t A, B = {0}, result;
  s21_create_matrix(2, 2, &A);
  B.rows = 2;
  B.columns = 0;
  B.matrix = NULL;
  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), INCORRECT_MATRIX);
  s21_remove_matrix(&A);
}
END_TEST

Suite *suite_sum_matrix(void) {
  Suite *s = suite_create("s21_sum_matrix");
  TCase *tc = tcase_create("sum_matrix");

  // Основные функциональные тесты
  tcase_add_test(tc, test_sum_matrix_basic);
  tcase_add_test(tc, test_sum_matrix_1x1);
  tcase_add_test(tc, test_sum_matrix_3x3);
  tcase_add_test(tc, test_sum_matrix_rectangular);

  // Специальные значения
  tcase_add_test(tc, test_sum_matrix_negative_values);
  tcase_add_test(tc, test_sum_matrix_zero_matrix);
  tcase_add_test(tc, test_sum_matrix_large_values);
  tcase_add_test(tc, test_sum_matrix_small_values);
  tcase_add_test(tc, test_sum_matrix_mixed_values);

  // Математические свойства
  tcase_add_test(tc, test_sum_matrix_commutative);

  // Тесты на ошибки
  tcase_add_test(tc, test_sum_matrix_different_size);
  tcase_add_test(tc, test_sum_matrix_different_rows);
  tcase_add_test(tc, test_sum_matrix_different_columns);
  tcase_add_test(tc, test_sum_matrix_null_first);
  tcase_add_test(tc, test_sum_matrix_null_second);
  tcase_add_test(tc, test_sum_matrix_null_result);
  tcase_add_test(tc, test_sum_matrix_null_all);
  tcase_add_test(tc, test_sum_matrix_invalid_first);
  tcase_add_test(tc, test_sum_matrix_invalid_second);

  suite_add_tcase(s, tc);
  return s;
}
