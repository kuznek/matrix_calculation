#include <check.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "../s21_matrix.h"

START_TEST(test_mult_number_basic) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 4.0;
  ck_assert_int_eq(s21_mult_number(&A, 2.5, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 2.5, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 5.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 7.5, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 10.0, 1e-6);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_number_zero) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 4.0;
  ck_assert_int_eq(s21_mult_number(&A, 0.0, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.0, 1e-6);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

// Дополнительные тесты
START_TEST(test_mult_number_negative) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = -2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = -4.0;
  ck_assert_int_eq(s21_mult_number(&A, -3.0, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], -3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 6.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], -9.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 12.0, 1e-6);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_number_very_small) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 4.0;
  double small_num = 1e-15;
  ck_assert_int_eq(s21_mult_number(&A, small_num, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], small_num, 1e-20);
  ck_assert_double_eq_tol(result.matrix[1][1], 4.0 * small_num, 1e-20);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_number_very_large) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 4.0;
  double large_num = 1e10;
  ck_assert_int_eq(s21_mult_number(&A, large_num, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], large_num, 1e4);
  ck_assert_double_eq_tol(result.matrix[1][1], 4.0 * large_num, 1e4);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_number_infinity) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 0.0;
  A.matrix[1][0] = -1.0;
  A.matrix[1][1] = 2.0;
  ck_assert_int_eq(s21_mult_number(&A, INFINITY, &result), OK);
  ck_assert_int_eq(isinf(result.matrix[0][0]), 1);
  ck_assert_int_eq(isnan(result.matrix[0][1]), 1);  // 0 * inf = nan
  ck_assert_int_eq(isinf(result.matrix[1][0]) && result.matrix[1][0] < 0, 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_number_nan) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 4.0;
  ck_assert_int_eq(s21_mult_number(&A, NAN, &result), OK);
  ck_assert_int_eq(isnan(result.matrix[0][0]), 1);
  ck_assert_int_eq(isnan(result.matrix[1][1]), 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_number_1x1) {
  matrix_t A, result;
  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 42.0;
  ck_assert_int_eq(s21_mult_number(&A, 0.5, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 21.0, 1e-6);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_number_null_matrix) {
  matrix_t result;
  ck_assert_int_eq(s21_mult_number(NULL, 2.0, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_mult_number_null_result) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  ck_assert_int_eq(s21_mult_number(&A, 2.0, NULL), INCORRECT_MATRIX);
  s21_remove_matrix(&A);
}
END_TEST

Suite *suite_mult_number(void) {
  Suite *s = suite_create("s21_mult_number");
  TCase *tc = tcase_create("mult_number");
  tcase_add_test(tc, test_mult_number_basic);
  tcase_add_test(tc, test_mult_number_zero);
  tcase_add_test(tc, test_mult_number_negative);
  tcase_add_test(tc, test_mult_number_very_small);
  tcase_add_test(tc, test_mult_number_very_large);
  tcase_add_test(tc, test_mult_number_infinity);
  tcase_add_test(tc, test_mult_number_nan);
  tcase_add_test(tc, test_mult_number_1x1);
  tcase_add_test(tc, test_mult_number_null_matrix);
  tcase_add_test(tc, test_mult_number_null_result);
  suite_add_tcase(s, tc);
  return s;
}
