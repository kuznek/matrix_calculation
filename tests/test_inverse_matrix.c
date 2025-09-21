#include <check.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "../s21_matrix.h"

// Основные функциональные тесты
START_TEST(test_inverse_matrix_1x1) {
  matrix_t A, result;
  s21_create_matrix(1, 1, &A);

  A.matrix[0][0] = 4.0;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 0.25, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_1x1_negative) {
  matrix_t A, result;
  s21_create_matrix(1, 1, &A);

  A.matrix[0][0] = -8.0;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], -0.125, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_1x1_fraction) {
  matrix_t A, result;
  s21_create_matrix(1, 1, &A);

  A.matrix[0][0] = 0.25;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 4.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_1x1_large) {
  matrix_t A, result;
  s21_create_matrix(1, 1, &A);

  A.matrix[0][0] = 1e6;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 1e-6, 1e-12);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_1x1_small) {
  matrix_t A, result;
  s21_create_matrix(1, 1, &A);

  A.matrix[0][0] = 1e-6;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);
  ck_assert_double_eq_tol(result.matrix[0][0], 1e6, 1e-3);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_2x2_basic) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 2.0;
  A.matrix[0][1] = 5.0;
  A.matrix[1][0] = 6.0;
  A.matrix[1][1] = 3.0;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  // Проверяем результат: det = 2*3 - 5*6 = -24
  // A^-1 = (1/det) * [3 -5; -6 2] = [-1/8, 5/24; 1/4, -1/12]
  ck_assert_double_eq_tol(result.matrix[0][0], -0.125, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 0.208333, 1e-5);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.25, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], -0.083333, 1e-5);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_2x2_identity) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 0.0;
  A.matrix[1][0] = 0.0;
  A.matrix[1][1] = 1.0;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 1.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 1.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_2x2_diagonal) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 3.0;
  A.matrix[0][1] = 0.0;
  A.matrix[1][0] = 0.0;
  A.matrix[1][1] = 4.0;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 1.0 / 3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.25, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_2x2_negative_elements) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = -1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 3.0;
  A.matrix[1][1] = -4.0;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  // det = (-1)*(-4) - 2*3 = 4 - 6 = -2
  // A^-1 = (1/-2) * [-4 -2; -3 -1] = [2, 1; 1.5, 0.5]
  ck_assert_double_eq_tol(result.matrix[0][0], 2.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], 1.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 1.5, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.5, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_2x2_fractional) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 0.5;
  A.matrix[0][1] = 1.5;
  A.matrix[1][0] = 2.5;
  A.matrix[1][1] = 0.25;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  // det = 0.5*0.25 - 1.5*2.5 = 0.125 - 3.75 = -3.625
  double det = -3.625;
  ck_assert_double_eq_tol(result.matrix[0][0], 0.25 / det, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], -1.5 / det, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], -2.5 / det, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.5 / det, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_3x3_identity) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A.matrix[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  // Обратная от единичной матрицы - единичная матрица
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

START_TEST(test_inverse_matrix_3x3_diagonal) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 2.0;
  A.matrix[0][1] = 0.0;
  A.matrix[0][2] = 0.0;
  A.matrix[1][0] = 0.0;
  A.matrix[1][1] = 3.0;
  A.matrix[1][2] = 0.0;
  A.matrix[2][0] = 0.0;
  A.matrix[2][1] = 0.0;
  A.matrix[2][2] = 4.0;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 0.5, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 1.0 / 3.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[2][2], 0.25, 1e-6);

  // Проверяем, что остальные элементы равны нулю
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i != j) {
        ck_assert_double_eq_tol(result.matrix[i][j], 0.0, 1e-6);
      }
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_3x3_complex) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 0.0;
  A.matrix[1][0] = 0.0;
  A.matrix[1][1] = 1.0;
  A.matrix[1][2] = 1.0;
  A.matrix[2][0] = 1.0;
  A.matrix[2][1] = 0.0;
  A.matrix[2][2] = 2.0;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  // (det = 4):
  ck_assert_double_eq_tol(result.matrix[0][0], 0.5, 1e-6);    // 2/4 = 0.5 ✓
  ck_assert_double_eq_tol(result.matrix[0][1], -1.0, 1e-6);   // -4/4 = -1.0 ✓
  ck_assert_double_eq_tol(result.matrix[0][2], 0.5, 1e-6);    // 2/4 = 0.5 ✓
  ck_assert_double_eq_tol(result.matrix[1][0], 0.25, 1e-6);   // 1/4 = 0.25
  ck_assert_double_eq_tol(result.matrix[1][1], 0.5, 1e-6);    // 2/4 = 0.5
  ck_assert_double_eq_tol(result.matrix[1][2], -0.25, 1e-6);  // -1/4 = -0.25
  ck_assert_double_eq_tol(result.matrix[2][0], -0.25, 1e-6);  // -1/4 = -0.25
  ck_assert_double_eq_tol(result.matrix[2][1], 0.5, 1e-6);    // 2/4 = 0.5
  ck_assert_double_eq_tol(result.matrix[2][2], 0.25, 1e-6);   // 1/4 = 0.25

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_4x4_identity) {
  matrix_t A, result;
  s21_create_matrix(4, 4, &A);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      A.matrix[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  // Проверяем, что результат - единичная матрица
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      double expected = (i == j) ? 1.0 : 0.0;
      ck_assert_double_eq_tol(result.matrix[i][j], expected, 1e-6);
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_4x4_diagonal) {
  matrix_t A, result;
  s21_create_matrix(4, 4, &A);

  // Заполняем нулями
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      A.matrix[i][j] = 0.0;
    }
  }

  A.matrix[0][0] = 2.0;
  A.matrix[1][1] = 0.5;
  A.matrix[2][2] = 4.0;
  A.matrix[3][3] = 0.25;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 0.5, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 2.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[2][2], 0.25, 1e-6);
  ck_assert_double_eq_tol(result.matrix[3][3], 4.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_5x5_identity) {
  matrix_t A, result;
  s21_create_matrix(5, 5, &A);

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      A.matrix[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  // Проверяем некоторые ключевые элементы
  for (int i = 0; i < 5; i++) {
    ck_assert_double_eq_tol(result.matrix[i][i], 1.0, 1e-6);
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

// Тесты с граничными значениями
START_TEST(test_inverse_matrix_large_values) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1e6;
  A.matrix[0][1] = 0.0;
  A.matrix[1][0] = 0.0;
  A.matrix[1][1] = 2e6;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 1e-6, 1e-12);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.5e-6, 1e-12);
  ck_assert_double_eq_tol(result.matrix[0][1], 0.0, 1e-12);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.0, 1e-12);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_small_values) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1e-3;
  A.matrix[0][1] = 0.0;
  A.matrix[1][0] = 0.0;
  A.matrix[1][1] = 2e-3;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 1e3, 1e-3);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.5e3, 1e-3);
  ck_assert_double_eq_tol(result.matrix[0][1], 0.0, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.0, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_inverse_matrix_too_small_values) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1e-15;
  A.matrix[0][1] = 0.0;
  A.matrix[1][0] = 0.0;
  A.matrix[1][1] = 2e-15;

  // Ожидаем CALCULATION_ERROR из-за потери точности
  ck_assert_int_eq(s21_inverse_matrix(&A, &result), CALCULATION_ERROR);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_inverse_matrix_mixed_scale) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1e6;
  A.matrix[0][1] = 0.0;
  A.matrix[1][0] = 0.0;
  A.matrix[1][1] = 1e-6;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  ck_assert_double_eq_tol(result.matrix[0][0], 1e-6, 1e-12);
  ck_assert_double_eq_tol(result.matrix[1][1], 1e6, 1e-3);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

// Математические свойства
START_TEST(test_inverse_matrix_double_inverse) {
  matrix_t A, inverse, double_inverse;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 3.0;
  A.matrix[0][1] = 1.0;
  A.matrix[1][0] = 2.0;
  A.matrix[1][1] = 4.0;

  ck_assert_int_eq(s21_inverse_matrix(&A, &inverse), OK);
  ck_assert_int_eq(s21_inverse_matrix(&inverse, &double_inverse), OK);

  // (A^-1)^-1 = A
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      ck_assert_double_eq_tol(double_inverse.matrix[i][j], A.matrix[i][j],
                              1e-6);
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&inverse);
  s21_remove_matrix(&double_inverse);
}
END_TEST

START_TEST(test_inverse_matrix_multiply_check) {
  matrix_t A, inverse, product;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 2.0;
  A.matrix[0][1] = 1.0;
  A.matrix[0][2] = 0.0;
  A.matrix[1][0] = 1.0;
  A.matrix[1][1] = 2.0;
  A.matrix[1][2] = 1.0;
  A.matrix[2][0] = 0.0;
  A.matrix[2][1] = 1.0;
  A.matrix[2][2] = 2.0;

  ck_assert_int_eq(s21_inverse_matrix(&A, &inverse), OK);
  ck_assert_int_eq(s21_mult_matrix(&A, &inverse, &product), OK);

  // A * A^-1 = I
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double expected = (i == j) ? 1.0 : 0.0;
      ck_assert_double_eq_tol(product.matrix[i][j], expected, 1e-6);
    }
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&inverse);
  s21_remove_matrix(&product);
}
END_TEST

START_TEST(test_inverse_matrix_orthogonal) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  // Матрица поворота на 45 градусов (ортогональная)
  double cos45 = 1.0 / sqrt(2.0);
  double sin45 = 1.0 / sqrt(2.0);

  A.matrix[0][0] = cos45;
  A.matrix[0][1] = -sin45;
  A.matrix[1][0] = sin45;
  A.matrix[1][1] = cos45;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  // Для ортогональной матрицы A^-1 = A^T
  ck_assert_double_eq_tol(result.matrix[0][0], cos45, 1e-6);
  ck_assert_double_eq_tol(result.matrix[0][1], sin45, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][0], -sin45, 1e-6);
  ck_assert_double_eq_tol(result.matrix[1][1], cos45, 1e-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

// Тесты на вырожденные матрицы (singular)
START_TEST(test_inverse_matrix_1x1_zero) {
  matrix_t A, result;
  s21_create_matrix(1, 1, &A);

  A.matrix[0][0] = 0.0;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), CALCULATION_ERROR);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_inverse_matrix_2x2_singular) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 2.0;
  A.matrix[1][1] = 4.0;  // det = 0

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), CALCULATION_ERROR);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_inverse_matrix_3x3_singular_proportional_rows) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[1][0] = 2.0;
  A.matrix[1][1] = 4.0;
  A.matrix[1][2] = 6.0;  // 2 * первая строка
  A.matrix[2][0] = 5.0;
  A.matrix[2][1] = 6.0;
  A.matrix[2][2] = 7.0;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), CALCULATION_ERROR);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_inverse_matrix_3x3_singular_zero_row) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[1][0] = 0.0;
  A.matrix[1][1] = 0.0;
  A.matrix[1][2] = 0.0;  // нулевая строка
  A.matrix[2][0] = 7.0;
  A.matrix[2][1] = 8.0;
  A.matrix[2][2] = 9.0;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), CALCULATION_ERROR);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_inverse_matrix_3x3_singular_arithmetic_progression) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[1][0] = 4.0;
  A.matrix[1][1] = 5.0;
  A.matrix[1][2] = 6.0;
  A.matrix[2][0] = 7.0;
  A.matrix[2][1] = 8.0;
  A.matrix[2][2] = 9.0;  // det = 0

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), CALCULATION_ERROR);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_inverse_matrix_4x4_singular) {
  matrix_t A, result;
  s21_create_matrix(4, 4, &A);

  // Заполняем так, чтобы определитель был 0
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      A.matrix[i][j] = i + j;
    }
  }

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), CALCULATION_ERROR);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_inverse_matrix_nearly_singular) {
  matrix_t A, result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[1][0] = 1.0000001;
  A.matrix[1][1] = 2.0000002;  // Очень маленький det

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), CALCULATION_ERROR);

  s21_remove_matrix(&A);
}
END_TEST

// Тесты на неквадратные матрицы
START_TEST(test_inverse_matrix_rectangular_2x3) {
  matrix_t A, result;
  s21_create_matrix(2, 3, &A);

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), CALCULATION_ERROR);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_inverse_matrix_rectangular_3x2) {
  matrix_t A, result;
  s21_create_matrix(3, 2, &A);

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), CALCULATION_ERROR);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_inverse_matrix_rectangular_1x5) {
  matrix_t A, result;
  s21_create_matrix(1, 5, &A);

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), CALCULATION_ERROR);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_inverse_matrix_rectangular_5x1) {
  matrix_t A, result;
  s21_create_matrix(5, 1, &A);

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), CALCULATION_ERROR);

  s21_remove_matrix(&A);
}
END_TEST

// Тесты на NULL указатели
START_TEST(test_inverse_matrix_null_input) {
  matrix_t result;

  ck_assert_int_eq(s21_inverse_matrix(NULL, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_inverse_matrix_null_result) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);

  ck_assert_int_eq(s21_inverse_matrix(&A, NULL), INCORRECT_MATRIX);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_inverse_matrix_null_both) {
  ck_assert_int_eq(s21_inverse_matrix(NULL, NULL), INCORRECT_MATRIX);
}
END_TEST

// Тесты на невалидные матрицы
START_TEST(test_inverse_matrix_invalid_rows) {
  matrix_t A = {0}, result;
  A.rows = 0;
  A.columns = 2;
  A.matrix = NULL;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_inverse_matrix_invalid_columns) {
  matrix_t A = {0}, result;
  A.rows = 2;
  A.columns = 0;
  A.matrix = NULL;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_inverse_matrix_invalid_negative_size) {
  matrix_t A = {0}, result;
  A.rows = -1;
  A.columns = -1;
  A.matrix = NULL;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_inverse_matrix_null_data) {
  matrix_t A = {0}, result;
  A.rows = 2;
  A.columns = 2;
  A.matrix = NULL;

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), INCORRECT_MATRIX);
}
END_TEST

// Тест производительности для больших матриц
START_TEST(test_inverse_matrix_6x6_identity) {
  matrix_t A, result;
  s21_create_matrix(6, 6, &A);

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      A.matrix[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }

  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);

  // Проверяем диагональные элементы
  for (int i = 0; i < 6; i++) {
    ck_assert_double_eq_tol(result.matrix[i][i], 1.0, 1e-6);
  }

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

Suite *suite_inverse_matrix(void) {
  Suite *s = suite_create("s21_inverse_matrix");
  TCase *tc = tcase_create("inverse_matrix");

  // Основные функциональные тесты
  tcase_add_test(tc, test_inverse_matrix_1x1);
  tcase_add_test(tc, test_inverse_matrix_1x1_negative);
  tcase_add_test(tc, test_inverse_matrix_1x1_fraction);
  tcase_add_test(tc, test_inverse_matrix_1x1_large);
  tcase_add_test(tc, test_inverse_matrix_1x1_small);

  tcase_add_test(tc, test_inverse_matrix_2x2_basic);
  tcase_add_test(tc, test_inverse_matrix_2x2_identity);
  tcase_add_test(tc, test_inverse_matrix_2x2_diagonal);
  tcase_add_test(tc, test_inverse_matrix_2x2_negative_elements);
  tcase_add_test(tc, test_inverse_matrix_2x2_fractional);

  tcase_add_test(tc, test_inverse_matrix_3x3_identity);
  tcase_add_test(tc, test_inverse_matrix_3x3_diagonal);
  tcase_add_test(tc, test_inverse_matrix_3x3_complex);

  tcase_add_test(tc, test_inverse_matrix_4x4_identity);
  tcase_add_test(tc, test_inverse_matrix_4x4_diagonal);
  tcase_add_test(tc, test_inverse_matrix_5x5_identity);
  tcase_add_test(tc, test_inverse_matrix_6x6_identity);

  // Граничные значения
  tcase_add_test(tc, test_inverse_matrix_large_values);
  tcase_add_test(tc, test_inverse_matrix_small_values);
  tcase_add_test(tc, test_inverse_matrix_too_small_values);
  tcase_add_test(tc, test_inverse_matrix_mixed_scale);

  // Математические свойства
  tcase_add_test(tc, test_inverse_matrix_double_inverse);
  tcase_add_test(tc, test_inverse_matrix_multiply_check);
  tcase_add_test(tc, test_inverse_matrix_orthogonal);

  // Вырожденные матрицы
  tcase_add_test(tc, test_inverse_matrix_1x1_zero);
  tcase_add_test(tc, test_inverse_matrix_2x2_singular);
  tcase_add_test(tc, test_inverse_matrix_3x3_singular_proportional_rows);
  tcase_add_test(tc, test_inverse_matrix_3x3_singular_zero_row);
  tcase_add_test(tc, test_inverse_matrix_3x3_singular_arithmetic_progression);
  tcase_add_test(tc, test_inverse_matrix_4x4_singular);
  tcase_add_test(tc, test_inverse_matrix_nearly_singular);

  // Неквадратные матрицы
  tcase_add_test(tc, test_inverse_matrix_rectangular_2x3);
  tcase_add_test(tc, test_inverse_matrix_rectangular_3x2);
  tcase_add_test(tc, test_inverse_matrix_rectangular_1x5);
  tcase_add_test(tc, test_inverse_matrix_rectangular_5x1);

  // Тесты на ошибки
  tcase_add_test(tc, test_inverse_matrix_null_input);
  tcase_add_test(tc, test_inverse_matrix_null_result);
  tcase_add_test(tc, test_inverse_matrix_null_both);
  tcase_add_test(tc, test_inverse_matrix_invalid_rows);
  tcase_add_test(tc, test_inverse_matrix_invalid_columns);
  tcase_add_test(tc, test_inverse_matrix_invalid_negative_size);
  tcase_add_test(tc, test_inverse_matrix_null_data);

  suite_add_tcase(s, tc);
  return s;
}
