#include "s21_matrix.h"

// void print_matrix(matrix_t* result) {
//     for (int i = 0; i < result->rows; i++) {
//         for (int j = 0; j < result->columns; j++) {
//             printf("%.2f ", result->matrix[i][j]);
//         }
//         putchar('\n');
//     }
// }

// проверка переменного количества матриц на валидность (вспомогательная)
void check_matrices(int *status, int n, ...) {
  va_list matrices;
  va_start(matrices, n);
  for (int i = 0; i < n && *status != INCORRECT_MATRIX; i++) {
    matrix_t *current = va_arg(matrices, matrix_t *);
    if (current == NULL || !current->matrix) {
      *status = INCORRECT_MATRIX;
    }
  }
  va_end(matrices);
}

// проверка матриц на одинаковые размеры (вспомогательная)
void check_matrices_same_size(int *status, matrix_t *A, matrix_t *B) {
  if (*status == OK) {
    if (A->rows != B->rows || A->columns != B->columns) {
      *status = CALCULATION_ERROR;
    }
  }
}

// cоздание матрицы-минора (исключение строки и столбца)
int create_minor_matrix(matrix_t *A, int exclude_row, int exclude_col,
                        matrix_t *minor) {
  int status = OK;

  if (!A || !A->matrix || exclude_row < 0 || exclude_row >= A->rows ||
      exclude_col < 0 || exclude_col >= A->columns || A->rows != A->columns ||
      A->rows <= 1) {
    status = INCORRECT_MATRIX;
  } else {
    // итоговая размером (n-1)x(n-1)
    status = s21_create_matrix(A->rows - 1, A->columns - 1, minor);

    int minor_row = 0;
    for (int i = 0; i < A->rows && status == OK; i++) {
      if (i == exclude_row) continue;

      int minor_col = 0;
      for (int j = 0; j < A->columns; j++) {
        if (j == exclude_col) continue;

        minor->matrix[minor_row][minor_col] = A->matrix[i][j];
        minor_col++;
      }
      minor_row++;
    }
  }

  return status;
}

// вычисление минора элемента A[row][col]
double calculate_minor(matrix_t *A, int row, int col) {
  matrix_t minor;
  double minor_det = 0.0;

  if (create_minor_matrix(A, row, col, &minor) == OK) {
    s21_determinant(&minor, &minor_det);
    s21_remove_matrix(&minor);
  }

  return minor_det;
}