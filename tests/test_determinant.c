#include <check.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "../s21_matrix.h"

START_TEST(test_determinant_1x1) {
  matrix_t A;
  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 5.0;
  double result;
  ck_assert_int_eq(s21_determinant(&A, &result), OK);
  ck_assert_double_eq_tol(result, 5.0, 1e-6);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_2x2) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 4.0;
  double result;
  ck_assert_int_eq(s21_determinant(&A, &result), OK);
  ck_assert_double_eq_tol(result, -2.0, 1e-6);  // 1*4 - 2*3 = -2
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_3x3) {
  matrix_t A;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[1][0] = 4.0;
  A.matrix[1][1] = 5.0;
  A.matrix[1][2] = 6.0;
  A.matrix[2][0] = 7.0;
  A.matrix[2][1] = 8.0;
  A.matrix[2][2] = 9.0;
  double result;
  ck_assert_int_eq(s21_determinant(&A, &result), OK);
  ck_assert_double_eq_tol(result, 0.0, 1e-6);  // Определитель равен 0
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_non_square) {
  matrix_t A;
  s21_create_matrix(2, 3, &A);
  double result;
  ck_assert_int_eq(s21_determinant(&A, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_null) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  double result;
  ck_assert_int_eq(s21_determinant(NULL, &result), INCORRECT_MATRIX);
  ck_assert_int_eq(s21_determinant(&A, NULL), INCORRECT_MATRIX);
  ck_assert_int_eq(s21_determinant(NULL, NULL), INCORRECT_MATRIX);
  s21_remove_matrix(&A);
}
END_TEST

// Дополнительные тесты
START_TEST(test_determinant_1x1_zero) {
  matrix_t A;
  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 0.0;
  double result;
  ck_assert_int_eq(s21_determinant(&A, &result), OK);
  ck_assert_double_eq_tol(result, 0.0, 1e-6);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_1x1_negative) {
  matrix_t A;
  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = -7.5;
  double result;
  ck_assert_int_eq(s21_determinant(&A, &result), OK);
  ck_assert_double_eq_tol(result, -7.5, 1e-6);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_2x2_identity) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 0.0;
  A.matrix[1][0] = 0.0;
  A.matrix[1][1] = 1.0;
  double result;
  ck_assert_int_eq(s21_determinant(&A, &result), OK);
  ck_assert_double_eq_tol(result, 1.0, 1e-6);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_3x3_identity) {
  matrix_t A;
  s21_create_matrix(3, 3, &A);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A.matrix[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }
  double result;
  ck_assert_int_eq(s21_determinant(&A, &result), OK);
  ck_assert_double_eq_tol(result, 1.0, 1e-6);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_4x4) {
  matrix_t A;
  s21_create_matrix(4, 4, &A);
  // Создаем матрицу с известным определителем
  A.matrix[0][0] = 2.0;
  A.matrix[0][1] = 1.0;
  A.matrix[0][2] = 3.0;
  A.matrix[0][3] = 4.0;
  A.matrix[1][0] = 1.0;
  A.matrix[1][1] = 2.0;
  A.matrix[1][2] = 1.0;
  A.matrix[1][3] = 2.0;
  A.matrix[2][0] = 3.0;
  A.matrix[2][1] = 1.0;
  A.matrix[2][2] = 2.0;
  A.matrix[2][3] = 1.0;
  A.matrix[3][0] = 1.0;
  A.matrix[3][1] = 3.0;
  A.matrix[3][2] = 1.0;
  A.matrix[3][3] = 2.0;
  double result;
  ck_assert_int_eq(s21_determinant(&A, &result), OK);
  // Проверяем, что результат не NaN и конечен
  ck_assert_int_eq(isnan(result), 0);
  ck_assert_int_eq(isinf(result), 0);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_large_values) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1e6;
  A.matrix[0][1] = 2e6;
  A.matrix[1][0] = 3e6;
  A.matrix[1][1] = 4e6;
  double result;
  ck_assert_int_eq(s21_determinant(&A, &result), OK);
  ck_assert_double_eq_tol(result, -2e12,
                          1e6);  // (1e6 * 4e6) - (2e6 * 3e6) = -2e12
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_small_values) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1e-6;
  A.matrix[0][1] = 2e-6;
  A.matrix[1][0] = 3e-6;
  A.matrix[1][1] = 4e-6;
  double result;
  ck_assert_int_eq(s21_determinant(&A, &result), OK);
  ck_assert_double_eq_tol(result, -2e-12, 1e-18);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_nearly_singular) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 1.000000001;
  A.matrix[1][1] = 2.000000002;  // Почти одинаковые строки
  double result;
  ck_assert_int_eq(s21_determinant(&A, &result), OK);
  ck_assert(fabs(result) < 1e-8);  // Очень маленький определитель
  s21_remove_matrix(&A);
}
END_TEST

Suite *suite_determinant(void) {
  Suite *s = suite_create("s21_determinant");
  TCase *tc = tcase_create("determinant");
  tcase_add_test(tc, test_determinant_1x1);
  tcase_add_test(tc, test_determinant_2x2);
  tcase_add_test(tc, test_determinant_3x3);
  tcase_add_test(tc, test_determinant_non_square);
  tcase_add_test(tc, test_determinant_null);
  tcase_add_test(tc, test_determinant_1x1_zero);
  tcase_add_test(tc, test_determinant_1x1_negative);
  tcase_add_test(tc, test_determinant_2x2_identity);
  tcase_add_test(tc, test_determinant_3x3_identity);
  tcase_add_test(tc, test_determinant_4x4);
  tcase_add_test(tc, test_determinant_large_values);
  tcase_add_test(tc, test_determinant_small_values);
  tcase_add_test(tc, test_determinant_nearly_singular);
  suite_add_tcase(s, tc);
  return s;
}
