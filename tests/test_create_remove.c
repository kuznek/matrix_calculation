#include <check.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>

#include "../s21_matrix.h"

// Базовые тесты
START_TEST(test_create_matrix_basic) {
  matrix_t result;
  ck_assert_int_eq(s21_create_matrix(3, 3, &result), OK);
  ck_assert_int_eq(result.rows, 3);
  ck_assert_int_eq(result.columns, 3);
  ck_assert_ptr_nonnull(result.matrix);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_create_matrix_invalid) {
  matrix_t result;
  ck_assert_int_eq(s21_create_matrix(0, 3, &result), INCORRECT_MATRIX);
  ck_assert_int_eq(s21_create_matrix(3, 0, &result), INCORRECT_MATRIX);
  ck_assert_int_eq(s21_create_matrix(-1, 3, &result), INCORRECT_MATRIX);
  ck_assert_int_eq(s21_create_matrix(3, -1, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_create_matrix_null) {
  ck_assert_int_eq(s21_create_matrix(3, 3, NULL), INCORRECT_MATRIX);
}
END_TEST

// Дополнительные тесты на граничные случаи
START_TEST(test_create_matrix_1x1) {
  matrix_t result;
  ck_assert_int_eq(s21_create_matrix(1, 1, &result), OK);
  ck_assert_int_eq(result.rows, 1);
  ck_assert_int_eq(result.columns, 1);
  ck_assert_ptr_nonnull(result.matrix);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_create_matrix_large) {
  matrix_t result;
  // Тест на большую матрицу (но не слишком большую для памяти)
  ck_assert_int_eq(s21_create_matrix(100, 100, &result), OK);
  ck_assert_int_eq(result.rows, 100);
  ck_assert_int_eq(result.columns, 100);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_create_matrix_extreme_negative) {
  matrix_t result;
  ck_assert_int_eq(s21_create_matrix(INT_MIN, 3, &result), INCORRECT_MATRIX);
  ck_assert_int_eq(s21_create_matrix(3, INT_MIN, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_create_matrix_rectangle_extreme) {
  matrix_t result;
  ck_assert_int_eq(s21_create_matrix(1, 1000, &result), OK);
  s21_remove_matrix(&result);
  ck_assert_int_eq(s21_create_matrix(1000, 1, &result), OK);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_remove_matrix_basic) {
  matrix_t matrix;
  s21_create_matrix(3, 3, &matrix);
  s21_remove_matrix(&matrix);
  ck_assert_ptr_null(matrix.matrix);
  ck_assert_int_eq(matrix.rows, 0);
  ck_assert_int_eq(matrix.columns, 0);
}
END_TEST

START_TEST(test_remove_matrix_null) {
  s21_remove_matrix(NULL);  // Не должно вызывать ошибку
  ck_assert_int_eq(1, 1);  // Тест пройден, если не было crash
}
END_TEST

START_TEST(test_remove_matrix_double_free) {
  matrix_t matrix;
  s21_create_matrix(3, 3, &matrix);
  s21_remove_matrix(&matrix);
  // Повторный вызов не должен вызывать ошибку
  s21_remove_matrix(&matrix);
  ck_assert_int_eq(1, 1);
}
END_TEST

START_TEST(test_remove_matrix_null_inside) {
  matrix_t matrix = {0};
  matrix.matrix = NULL;
  matrix.rows = 0;
  matrix.columns = 0;
  s21_remove_matrix(&matrix);  // Не должно вызывать ошибку
  ck_assert_int_eq(1, 1);
}
END_TEST

Suite *suite_create_remove(void) {
  Suite *s = suite_create("s21_create_remove");
  TCase *tc = tcase_create("create_remove");
  tcase_add_test(tc, test_create_matrix_basic);
  tcase_add_test(tc, test_create_matrix_invalid);
  tcase_add_test(tc, test_create_matrix_null);
  tcase_add_test(tc, test_create_matrix_1x1);
  tcase_add_test(tc, test_create_matrix_large);
  tcase_add_test(tc, test_create_matrix_extreme_negative);
  tcase_add_test(tc, test_create_matrix_rectangle_extreme);
  tcase_add_test(tc, test_remove_matrix_basic);
  tcase_add_test(tc, test_remove_matrix_null);
  tcase_add_test(tc, test_remove_matrix_double_free);
  tcase_add_test(tc, test_remove_matrix_null_inside);
  suite_add_tcase(s, tc);
  return s;
}
