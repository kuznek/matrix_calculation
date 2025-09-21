#include <check.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "../s21_matrix.h"

START_TEST(test_mult_matrix_basic) {
  matrix_t A, B, result;
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(3, 2, &B);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[1][0] = 4.0;
  A.matrix[1][1] = 5.0;
  A.matrix[1][2] = 6.0;

  B.matrix[0][0] = 7.0;
  B.matrix[0][1] = 8.0;
  B.matrix[1][0] = 9.0;
  B.matrix[1][1] = 10.0;
  B.matrix[2][0] = 11.0;
  B.matrix[2][1] = 12.0;

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 58.0, 1e-6);  // 1*7 + 2*9 + 3*11
  ck_assert_double_eq_tol(result.matrix[0][1], 64.0,
                          1e-6);  // 1*8 + 2*10 + 3*12
  ck_assert_double_eq_tol(result.matrix[1][0], 139.0,
                          1e-6);  // 4*7 + 5*9 + 6*11
  ck_assert_double_eq_tol(result.matrix[1][1], 154.0,
                          1e-6);  // 4*8 + 5*10 + 6*12

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_matrix_1x1) {
  matrix_t A, B, result;
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &B);

  A.matrix[0][0] = 5.0;
  B.matrix[0][0] = 3.0;

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 15.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_matrix_1x1_zero) {
  matrix_t A, B, result;
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &B);

  A.matrix[0][0] = 7.0;
  B.matrix[0][0] = 0.0;

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 0.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_matrix_square_2x2) {
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

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 19.0, 1e-6);  // 1*5 + 2*7 = 19
  ck_assert_double_eq_tol(result.matrix[0][1], 22.0, 1e-6);  // 1*6 + 2*8 = 22
  ck_assert_double_eq_tol(result.matrix[1][0], 43.0, 1e-6);  // 3*5 + 4*7 = 43
  ck_assert_double_eq_tol(result.matrix[1][1], 50.0, 1e-6);  // 3*6 + 4*8 = 50

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_matrix_square_3x3) {
  matrix_t A, B, result;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &B);

  // Инициализация A
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[1][0] = 4.0;
  A.matrix[1][1] = 5.0;
  A.matrix[1][2] = 6.0;
  A.matrix[2][0] = 7.0;
  A.matrix[2][1] = 8.0;
  A.matrix[2][2] = 9.0;

  // Инициализация B (единичная матрица)
  B.matrix[0][0] = 1.0;
  B.matrix[0][1] = 0.0;
  B.matrix[0][2] = 0.0;
  B.matrix[1][0] = 0.0;
  B.matrix[1][1] = 1.0;
  B.matrix[1][2] = 0.0;
  B.matrix[2][0] = 0.0;
  B.matrix[2][1] = 0.0;
  B.matrix[2][2] = 1.0;

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);

  // Умножение на единичную матрицу должно дать исходную матрицу
  ck_assert_double_eq_tol(result.matrix[0][0], 1.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 2.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][2], 3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 4.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 5.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][2], 6.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[2][0], 7.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[2][1], 8.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[2][2], 9.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_matrix_vector_multiplication) {
  matrix_t A, B, result;
  s21_create_matrix(3, 1, &A);  // Столбец-вектор
  s21_create_matrix(1, 3, &B);  // Строка-вектор

  A.matrix[0][0] = 1.0;
  A.matrix[1][0] = 2.0;
  A.matrix[2][0] = 3.0;

  B.matrix[0][0] = 4.0;
  B.matrix[0][1] = 5.0;
  B.matrix[0][2] = 6.0;

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_int_eq(result.rows, 3);
  ck_assert_int_eq(result.columns, 3);

  // Результат: внешнее произведение векторов
  ck_assert_double_eq_tol(result.matrix[0][0], 4.0, 1e-6);   // 1*4
  ck_assert_double_eq_tol(result.matrix[0][1], 5.0, 1e-6);   // 1*5
  ck_assert_double_eq_tol(result.matrix[0][2], 6.0, 1e-6);   // 1*6
  ck_assert_double_eq_tol(result.matrix[1][0], 8.0, 1e-6);   // 2*4
  ck_assert_double_eq_tol(result.matrix[1][1], 10.0, 1e-6);  // 2*5
  ck_assert_double_eq_tol(result.matrix[1][2], 12.0, 1e-6);  // 2*6
  ck_assert_double_eq_tol(result.matrix[2][0], 12.0, 1e-6);  // 3*4
  ck_assert_double_eq_tol(result.matrix[2][1], 15.0, 1e-6);  // 3*5
  ck_assert_double_eq_tol(result.matrix[2][2], 18.0, 1e-6);  // 3*6

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_matrix_row_vector_dot_column_vector) {
  matrix_t A, B, result;
  s21_create_matrix(1, 3, &A);  // Строка-вектор
  s21_create_matrix(3, 1, &B);  // Столбец-вектор

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;

  B.matrix[0][0] = 4.0;
  B.matrix[1][0] = 5.0;
  B.matrix[2][0] = 6.0;

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_int_eq(result.rows, 1);
  ck_assert_int_eq(result.columns, 1);

  // Результат: скалярное произведение = 1*4 + 2*5 + 3*6 = 32
  ck_assert_double_eq_tol(result.matrix[0][0], 32.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_matrix_negative_values) {
  matrix_t A, B, result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = -1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = -3.0;
  A.matrix[1][1] = 4.0;

  B.matrix[0][0] = 5.0;
  B.matrix[0][1] = -6.0;
  B.matrix[1][0] = -7.0;
  B.matrix[1][1] = 8.0;

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], -19.0,
                          1e-6);  // -1*5 + 2*(-7) = -19
  ck_assert_double_eq_tol(result.matrix[0][1], 22.0,
                          1e-6);  // -1*(-6) + 2*8 = 22
  ck_assert_double_eq_tol(result.matrix[1][0], -43.0,
                          1e-6);  // -3*5 + 4*(-7) = -43
  ck_assert_double_eq_tol(result.matrix[1][1], 50.0,
                          1e-6);  // -3*(-6) + 4*8 = 50

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_matrix_fractional_values) {
  matrix_t A, B, result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 0.5;
  A.matrix[0][1] = 1.5;
  A.matrix[1][0] = 2.5;
  A.matrix[1][1] = 3.5;

  B.matrix[0][0] = 0.2;
  B.matrix[0][1] = 0.4;
  B.matrix[1][0] = 0.6;
  B.matrix[1][1] = 0.8;

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 1.0,
                          1e-6);  // 0.5*0.2 + 1.5*0.6 = 1.0
  ck_assert_double_eq_tol(result.matrix[0][1], 1.4,
                          1e-6);  // 0.5*0.4 + 1.5*0.8 = 1.4
  ck_assert_double_eq_tol(result.matrix[1][0], 2.6,
                          1e-6);  // 2.5*0.2 + 3.5*0.6 = 2.6
  ck_assert_double_eq_tol(result.matrix[1][1], 3.8,
                          1e-6);  // 2.5*0.4 + 3.5*0.8 = 3.8

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_matrix_large_values) {
  matrix_t A, B, result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1000.0;
  A.matrix[0][1] = 2000.0;
  A.matrix[1][0] = 3000.0;
  A.matrix[1][1] = 4000.0;

  B.matrix[0][0] = 1.0;
  B.matrix[0][1] = 2.0;
  B.matrix[1][0] = 3.0;
  B.matrix[1][1] = 4.0;

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 7000.0,
                          1e-3);  // 1000*1 + 2000*3
  ck_assert_double_eq_tol(result.matrix[0][1], 10000.0,
                          1e-3);  // 1000*2 + 2000*4
  ck_assert_double_eq_tol(result.matrix[1][0], 15000.0,
                          1e-3);  // 3000*1 + 4000*3
  ck_assert_double_eq_tol(result.matrix[1][1], 22000.0,
                          1e-3);  // 3000*2 + 4000*4

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_matrix_small_values) {
  matrix_t A, B, result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1e-6;
  A.matrix[0][1] = 2e-6;
  A.matrix[1][0] = 3e-6;
  A.matrix[1][1] = 4e-6;

  B.matrix[0][0] = 1e6;
  B.matrix[0][1] = 2e6;
  B.matrix[1][0] = 3e6;
  B.matrix[1][1] = 4e6;

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 7.0,
                          1e-6);  // 1e-6*1e6 + 2e-6*3e6 = 7
  ck_assert_double_eq_tol(result.matrix[0][1], 10.0,
                          1e-6);  // 1e-6*2e6 + 2e-6*4e6 = 10
  ck_assert_double_eq_tol(result.matrix[1][0], 15.0,
                          1e-6);  // 3e-6*1e6 + 4e-6*3e6 = 15
  ck_assert_double_eq_tol(result.matrix[1][1], 22.0,
                          1e-6);  // 3e-6*2e6 + 4e-6*4e6 = 22

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_matrix_zero_matrix) {
  matrix_t A, B, result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  // A - обычная матрица, B - нулевая
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = 4.0;

  B.matrix[0][0] = 0.0;
  B.matrix[0][1] = 0.0;
  B.matrix[1][0] = 0.0;
  B.matrix[1][1] = 0.0;

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_matrix_rectangular_long) {
  matrix_t A, B, result;
  s21_create_matrix(1, 4, &A);  // Длинная строка
  s21_create_matrix(4, 1, &B);  // Длинный столбец

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[0][3] = 4.0;

  B.matrix[0][0] = 2.0;
  B.matrix[1][0] = 3.0;
  B.matrix[2][0] = 4.0;
  B.matrix[3][0] = 5.0;

  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_int_eq(result.rows, 1);
  ck_assert_int_eq(result.columns, 1);

  // Результат: 1*2 + 2*3 + 3*4 + 4*5 = 2 + 6 + 12 + 20 = 40
  ck_assert_double_eq_tol(result.matrix[0][0], 40.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

// Тесты на ошибки
START_TEST(test_mult_matrix_incompatible) {
  matrix_t A, B, result;
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(2, 2, &B);  // A.columns != B.rows (3 != 2)
  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_mult_matrix_incompatible_1x1) {
  matrix_t A, B, result;
  s21_create_matrix(1, 2, &A);
  s21_create_matrix(1, 1, &B);  // A.columns != B.rows (2 != 1)
  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_mult_matrix_incompatible_large) {
  matrix_t A, B, result;
  s21_create_matrix(4, 5, &A);
  s21_create_matrix(3, 2, &B);  // A.columns != B.rows (5 != 3)
  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_mult_matrix_null_first) {
  matrix_t B, result;
  s21_create_matrix(2, 2, &B);
  ck_assert_int_eq(s21_mult_matrix(NULL, &B, &result), INCORRECT_MATRIX);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_mult_matrix_null_second) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);
  ck_assert_int_eq(s21_mult_matrix(&A, NULL, &result), INCORRECT_MATRIX);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_mult_matrix_null_result) {
  matrix_t A, B;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  ck_assert_int_eq(s21_mult_matrix(&A, &B, NULL), INCORRECT_MATRIX);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_mult_matrix_null_all) {
  ck_assert_int_eq(s21_mult_matrix(NULL, NULL, NULL), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_mult_matrix_invalid_first) {
  matrix_t A = {0}, B, result;
  A.rows = 0;
  A.columns = 2;
  A.matrix = NULL;
  s21_create_matrix(2, 2, &B);
  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), INCORRECT_MATRIX);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_mult_matrix_invalid_second) {
  matrix_t A, B = {0}, result;
  s21_create_matrix(2, 2, &A);
  B.rows = 2;
  B.columns = 0;
  B.matrix = NULL;
  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), INCORRECT_MATRIX);
  s21_remove_matrix(&A);
}
END_TEST

Suite *suite_mult_matrix(void) {
  Suite *s = suite_create("s21_mult_matrix");
  TCase *tc = tcase_create("mult_matrix");

  // Основные тесты функциональности
  tcase_add_test(tc, test_mult_matrix_basic);
  tcase_add_test(tc, test_mult_matrix_1x1);
  tcase_add_test(tc, test_mult_matrix_1x1_zero);
  tcase_add_test(tc, test_mult_matrix_square_2x2);
  tcase_add_test(tc, test_mult_matrix_square_3x3);

  // Тесты векторных операций
  tcase_add_test(tc, test_mult_matrix_vector_multiplication);
  tcase_add_test(tc, test_mult_matrix_row_vector_dot_column_vector);
  tcase_add_test(tc, test_mult_matrix_rectangular_long);

  // Тесты на различные типы значений
  tcase_add_test(tc, test_mult_matrix_negative_values);
  tcase_add_test(tc, test_mult_matrix_fractional_values);
  tcase_add_test(tc, test_mult_matrix_large_values);
  tcase_add_test(tc, test_mult_matrix_small_values);
  tcase_add_test(tc, test_mult_matrix_zero_matrix);

  // Тесты на ошибки размерности
  tcase_add_test(tc, test_mult_matrix_incompatible);
  tcase_add_test(tc, test_mult_matrix_incompatible_1x1);
  tcase_add_test(tc, test_mult_matrix_incompatible_large);

  // Тесты на NULL указатели и некорректные данные
  tcase_add_test(tc, test_mult_matrix_null_first);
  tcase_add_test(tc, test_mult_matrix_null_second);
  tcase_add_test(tc, test_mult_matrix_null_result);
  tcase_add_test(tc, test_mult_matrix_null_all);
  tcase_add_test(tc, test_mult_matrix_invalid_first);
  tcase_add_test(tc, test_mult_matrix_invalid_second);

  suite_add_tcase(s, tc);
  return s;
}
