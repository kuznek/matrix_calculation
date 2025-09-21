#include "s21_matrix.h"

// создание матрицы rows*columns с динамическим выделением памяти
int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int status = OK;
  if (rows <= 0 || columns <= 0 || result == NULL) {
    status = INCORRECT_MATRIX;
  }
  if (status == OK) {
    // массив указателей на строки
    result->matrix = (double **)calloc(rows, sizeof(double *));
    if (result->matrix == NULL) {
      status = INCORRECT_MATRIX;
    }
  }
  int allocation_success = status == OK ? 1 : 0;
  for (int i = 0; i < rows && status == OK; i++) {
    // строки
    result->matrix[i] = (double *)calloc(columns, sizeof(double));
    if (result->matrix[i] == NULL) {
      for (int j = 0; j < i; j++) {
        free(result->matrix[j]);
      }
      free(result->matrix);
      result->matrix = NULL;
      allocation_success = 0;
      status = INCORRECT_MATRIX;
    }
  }

  if (allocation_success) {
    result->rows = rows;
    result->columns = columns;
  }
  return status;
}
// освобождение памяти
void s21_remove_matrix(matrix_t *A) {
  if (A != NULL && A->matrix != NULL) {
    for (int i = 0; i < A->rows; i++) {
      if (A->matrix[i] != NULL) {
        free(A->matrix[i]);
        A->matrix[i] = NULL;
      }
    }
    free(A->matrix);
    A->matrix = NULL;
    A->rows = 0;
    A->columns = 0;
  }
}

// A(i,j) == B(i,j) точность до 6 знака (включительно) после запятой
int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int status = SUCCESS;
  if (!A || !B || !A->matrix || !B->matrix) {
    status = FAILURE;
  } else if (A->rows != B->rows || A->columns != B->columns) {
    status = FAILURE;
  } else {
    for (int i = 0; i < A->rows && status != FAILURE; i++) {
      for (int j = 0; j < A->columns && status != FAILURE; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-6) {
          status = FAILURE;
        }
      }
    }
  }

  return status;
}

// C(i,j) = A(i,j) + B(i,j) одинаковые размеры
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  check_matrices(&status, 2, A, B);
  check_matrices_same_size(&status, A, B);
  if (!result) {
    status = INCORRECT_MATRIX;
  }
  if (status == OK) {
    status = s21_create_matrix(A->rows, B->columns, result);
  }
  if (status == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
    }
  }
  return status;
}
// C(i,j) = A(i,j) - B(i,j) одинаковые размеры
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  check_matrices(&status, 2, A, B);
  check_matrices_same_size(&status, A, B);
  if (!result) {
    status = INCORRECT_MATRIX;
  }
  if (status == OK) {
    status = s21_create_matrix(A->rows, B->columns, result);
  }
  if (status == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  }
  return status;
}
// B = λ × A(i,j)
int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int status = OK;
  check_matrices(&status, 1, A);
  if (!result) {
    status = INCORRECT_MATRIX;
  }
  if (status == OK) {
    status = s21_create_matrix(A->rows, A->columns, result);
  }
  if (status == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = number * A->matrix[i][j];
      }
    }
  }
  return status;
}
// C(i,j) = A(i,k) × B(k,j) строка * столбец - скалярное произведение
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  check_matrices(&status, 2, A, B);
  if (!result) {
    status = INCORRECT_MATRIX;
  }
  if (status == OK) {
    if (A->columns != B->rows) {
      status = CALCULATION_ERROR;
    }
  }
  if (status == OK) {
    // итоговая матрица i × j
    status = s21_create_matrix(A->rows, B->columns, result);
  }
  if (status == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        double sum = 0.0;
        // скалярное произведение векторов
        for (int k = 0; k < A->columns; k++) {
          sum += A->matrix[i][k] * B->matrix[k][j];
        }
        result->matrix[i][j] = sum;
      }
    }
  }
  return status;
}

// B[j][i] = A[i][j]
int s21_transpose(matrix_t *A, matrix_t *result) {
  int status = OK;
  check_matrices(&status, 1, A);
  if (!result) {
    status = INCORRECT_MATRIX;
  }
  if (status == OK) {
    // итоговая матрица j × i
    status = s21_create_matrix(A->columns, A->rows, result);
  }
  if (status == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }
  return status;
}
// вычисление определителя путём разложения по строке при размерах больше 2х2
// (только квадратные матрицы)
int s21_determinant(matrix_t *A, double *result) {
  int status = OK;
  check_matrices(&status, 1, A);
  if (!result) {
    status = INCORRECT_MATRIX;
  }
  if (status == OK) {
    if (A->rows != A->columns) {
      status = CALCULATION_ERROR;
    }
  }
  if (status == OK) {
    *result = 0.0;
    if (A->rows == 1) {
      *result = A->matrix[0][0];
    } else if (A->rows == 2) {
      // произведение элементов главной диагонали - произведение побочной
      *result =
          A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
    } else {
      // разложение по первой строке
      for (int j = 0; j < A->columns; j++) {
        double minor_det = calculate_minor(A, 0, j);
        // алгебраическое дополнение: (-1)^(i+j) * minor
        double cofactor = ((0 + j) % 2 == 0 ? 1 : -1) * minor_det;
        *result += A->matrix[0][j] * cofactor;
      }
    }
  }
  return status;
}
// вычисление матрицы алгебраических дополнений (только квадратные матрицы)
int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int status = OK;
  check_matrices(&status, 1, A);
  if (!result) {
    status = INCORRECT_MATRIX;
  }
  if (status == OK) {
    if (A->rows != A->columns) {
      status = CALCULATION_ERROR;
    }
  }
  if (status == OK) {
    status = s21_create_matrix(A->rows, A->columns, result);
  }
  if (status == OK) {
    if (A->rows == 1) {
      // матрица 1x1 алгебраическое дополнение: 1
      result->matrix[0][0] = 1.0;
    } else {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          double minor_det = calculate_minor(A, i, j);
          // алгебраическое дополнение: (-1)^(i+j) * minor
          result->matrix[i][j] = ((i + j) % 2 == 0 ? 1.0 : -1.0) * minor_det;
        }
      }
    }
  }

  return status;
}

// обратная матрица A^(-1) = (1/|A|) * M^T
int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int status = OK;
  check_matrices(&status, 1, A);
  if (!result) {
    status = INCORRECT_MATRIX;
  }
  if (status == OK) {
    if (A->rows != A->columns) {
      status = CALCULATION_ERROR;
    }
  }
  double det;
  if (status == OK) {
    status = s21_determinant(A, &det);
    if (fabs(det) < EPS) status = CALCULATION_ERROR;
  }
  if (status == OK) {
    matrix_t complements, transpose;
    status = s21_calc_complements(A, &complements);
    if (status == OK) status = s21_transpose(&complements, &transpose);
    if (status == OK) status = s21_mult_number(&transpose, 1.0 / det, result);
    s21_remove_matrix(&complements);
    s21_remove_matrix(&transpose);
  }
  return status;
}
