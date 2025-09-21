#include <check.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "../s21_matrix.h"

START_TEST(test_eq_matrix_equal) {
  matrix_t A, B;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1.123456;
  A.matrix[0][1] = 2.654321;
  A.matrix[1][0] = 3.789012;
  A.matrix[1][1] = 4.345678;
  B.matrix[0][0] = 1.123456;
  B.matrix[0][1] = 2.654321;
  B.matrix[1][0] = 3.789012;
  B.matrix[1][1] = 4.345678;
  ck_assert_int_eq(s21_eq_matrix(&A, &B), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_eq_matrix_equal_within_precision) {
  matrix_t A, B;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  // Различие меньше 1e-6 - должны быть равны
  A.matrix[0][0] = 1.1234560;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 4.0;
  B.matrix[0][0] = 1.1234565;
  B.matrix[0][1] = 2.0;  // различие 5e-7 < 1e-6
  B.matrix[1][0] = 3.0;
  B.matrix[1][1] = 4.0;
  ck_assert_int_eq(s21_eq_matrix(&A, &B), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_eq_matrix_not_equal_precision_boundary) {
  matrix_t A, B;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  // Гарантированное различие больше 1e-6
  A.matrix[0][0] = 1.000000;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 4.0;
  B.matrix[0][0] = 1.000002;
  B.matrix[0][1] = 2.0;  // различие 2e-6 > 1e-6
  B.matrix[1][0] = 3.0;
  B.matrix[1][1] = 4.0;
  ck_assert_int_eq(s21_eq_matrix(&A, &B), FAILURE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_eq_matrix_not_equal) {
  matrix_t A, B;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 4.0;
  B.matrix[0][0] = 1.0;
  B.matrix[0][1] = 2.0;
  B.matrix[1][0] = 3.0;
  B.matrix[1][1] = 4.1;
  ck_assert_int_eq(s21_eq_matrix(&A, &B), FAILURE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_eq_matrix_different_size) {
  matrix_t A, B;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(3, 3, &B);
  ck_assert_int_eq(s21_eq_matrix(&A, &B), FAILURE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_eq_matrix_null) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  ck_assert_int_eq(s21_eq_matrix(NULL, &A), FAILURE);
  ck_assert_int_eq(s21_eq_matrix(&A, NULL), FAILURE);
  ck_assert_int_eq(s21_eq_matrix(NULL, NULL), FAILURE);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_eq_matrix_exactly_1e6_difference) {
  matrix_t A, B;
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &B);
  A.matrix[0][0] = 1.0;
  B.matrix[0][0] = 1.000010;  // Гарантированно больше 1e-6
  ck_assert_int_eq(s21_eq_matrix(&A, &B), FAILURE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_eq_matrix_less_than_1e6_difference) {
  matrix_t A, B;
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &B);
  A.matrix[0][0] = 1.0;
  B.matrix[0][0] = 1.0000005;  // 5e-7 < 1e-6
  ck_assert_int_eq(s21_eq_matrix(&A, &B), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_eq_matrix_negative_values_precision) {
  matrix_t A, B;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = -1.123456;
  A.matrix[0][1] = -2.789012;
  A.matrix[1][0] = -3.456789;
  A.matrix[1][1] = -4.012345;
  B.matrix[0][0] = -1.123456;
  B.matrix[0][1] = -2.789012;
  B.matrix[1][0] = -3.456789;
  B.matrix[1][1] = -4.012345;
  ck_assert_int_eq(s21_eq_matrix(&A, &B), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_eq_matrix_zero_values) {
  matrix_t A, B;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 0.0;
  A.matrix[0][1] = 0.000010;  // 1e-5 > 1e-6
  A.matrix[1][0] = -0.000010;
  A.matrix[1][1] = 0.0;
  B.matrix[0][0] = 0.0;
  B.matrix[0][1] = 0.0;
  B.matrix[1][0] = 0.0;
  B.matrix[1][1] = 0.0;
  ck_assert_int_eq(s21_eq_matrix(&A, &B), FAILURE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_eq_matrix_very_small_numbers) {
  matrix_t A, B;
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &B);
  A.matrix[0][0] = 1e-10;
  B.matrix[0][0] = 1.0000001e-10;  // Очень маленькое различие
  ck_assert_int_eq(s21_eq_matrix(&A, &B), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_eq_matrix_large_numbers_precision) {
  matrix_t A, B;
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &B);
  A.matrix[0][0] = 1000000.123456;
  B.matrix[0][0] = 1000000.123466;  // Различие 1e-5 > 1e-6
  ck_assert_int_eq(s21_eq_matrix(&A, &B), FAILURE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

Suite *suite_eq_matrix(void) {
  Suite *s = suite_create("s21_eq_matrix");
  TCase *tc = tcase_create("eq_matrix");
  tcase_add_test(tc, test_eq_matrix_equal);
  tcase_add_test(tc, test_eq_matrix_equal_within_precision);
  tcase_add_test(tc, test_eq_matrix_not_equal_precision_boundary);
  tcase_add_test(tc, test_eq_matrix_not_equal);
  tcase_add_test(tc, test_eq_matrix_different_size);
  tcase_add_test(tc, test_eq_matrix_null);
  tcase_add_test(tc, test_eq_matrix_exactly_1e6_difference);
  tcase_add_test(tc, test_eq_matrix_less_than_1e6_difference);
  tcase_add_test(tc, test_eq_matrix_negative_values_precision);
  tcase_add_test(tc, test_eq_matrix_zero_values);
  tcase_add_test(tc, test_eq_matrix_very_small_numbers);
  tcase_add_test(tc, test_eq_matrix_large_numbers_precision);
  suite_add_tcase(s, tc);
  return s;
}
