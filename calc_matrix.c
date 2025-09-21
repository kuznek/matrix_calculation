#include <stdio.h>

#include "s21_matrix.c"
#include "s21_matrix_auxiliary.c"

int fill_matrix(int* rows, int* cols, matrix_t* A) {
  if (*rows == 0 || *cols == 0) {
    printf("Размеры матрицы\nСтроки: ");
    scanf("%d", rows);
    printf("Столбцы: ");
    scanf("%d", cols);
  }
  if (s21_create_matrix(*rows, *cols, A) != 0) return 1;
  printf("Матрица: \n");
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      scanf("%lf", &A->matrix[i][j]);
    }
  }
  return 0;
}

void print_matrix(matrix_t* result) {
  for (int i = 0; i < result->rows; i++) {
    for (int j = 0; j < result->columns; j++) {
      printf("%lf ", result->matrix[i][j]);
    }
    putchar('\n');
  }
}

int main() {
  printf(
      "\
1) Вычисление определителя                         |A|\n\
2) Сумма матриц\n\
3) Разность матриц\n\
4) Cравнение матриц (точность 1е-6)    \n\
5) Умножение матрицы на число                       B(i,j) = λ × A(i,j)\n\
6) Умножение матриц                                 C(i,j) = A(i,k) × B(k,j)\n\
7) Транспонирование матрицы                         B(j,i) = A(i,j)\n\
8) Вычисление матрицы алгебраических дополнений     (-1)^(i+j) * minor\n\
9) Вычисление обратной матрицы                      A^(-1) = (1/|A|) * M^T\n");

  int oper = -1;
  printf("Номер операции: ");
  scanf("%d", &oper);
  if (oper < 0 || oper > 9) {
    printf("Значение должно быть от 0 до 9\n");
    return 1;
  }
  matrix_t A;
  matrix_t B;
  matrix_t C;
  double result;
  int rows = 0;
  int cols = 0;
  switch (oper) {
    case 1:
      if (fill_matrix(&rows, &cols, &A) != 0) return 1;
      if (s21_determinant(&A, &result)) return 1;
      printf("%f\n", result);
      s21_remove_matrix(&A);
      break;
    case 2:
      if (fill_matrix(&rows, &cols, &A) != 0) return 1;
      if (fill_matrix(&rows, &cols, &B) != 0) return 1;
      if (s21_sum_matrix(&A, &B, &C)) return 1;
      printf("Результат: \n");
      print_matrix(&C);
      s21_remove_matrix(&A);
      s21_remove_matrix(&B);
      break;
    case 3:
      if (fill_matrix(&rows, &cols, &A) != 0) return 1;
      if (fill_matrix(&rows, &cols, &B) != 0) return 1;
      if (s21_sub_matrix(&A, &B, &C)) return 1;
      printf("Результат: \n");
      print_matrix(&C);
      s21_remove_matrix(&A);
      s21_remove_matrix(&B);
      break;
    case 4:
      if (fill_matrix(&rows, &cols, &A) != 0) return 1;
      if (fill_matrix(&rows, &cols, &B) != 0) return 1;
      int status = s21_eq_matrix(&A, &B);
      status == 1 ? printf("Матрицы равны\n") : printf("Матрицы не равны\n");
      s21_remove_matrix(&A);
      s21_remove_matrix(&B);
      break;
    case 5:
      if (fill_matrix(&rows, &cols, &A) != 0) return 1;
      printf("Число: ");
      double num = 0;
      scanf("%lf", &num);
      if (s21_mult_number(&A, num, &B) != 0) return 1;
      printf("Результат: \n");
      print_matrix(&B);
      s21_remove_matrix(&A);
      s21_remove_matrix(&B);
      break;
    case 6:
      if (fill_matrix(&rows, &cols, &A) != 0) return 1;
      rows = 0;
      cols = 0;
      if (fill_matrix(&rows, &cols, &B) != 0) return 1;
      if (s21_mult_matrix(&A, &B, &C) != 0) return 1;
      printf("Результат: \n");
      print_matrix(&C);
      s21_remove_matrix(&A);
      s21_remove_matrix(&B);
      s21_remove_matrix(&C);
      break;
    case 7:
      if (fill_matrix(&rows, &cols, &A) != 0) return 1;
      if (s21_transpose(&A, &B) != 0) return 1;
      printf("Результат: \n");
      print_matrix(&B);
      s21_remove_matrix(&A);
      s21_remove_matrix(&B);
      break;
    case 8:
      if (fill_matrix(&rows, &cols, &A) != 0) return 1;
      if (s21_calc_complements(&A, &B) != 0) return 1;
      printf("Результат: \n");
      print_matrix(&B);
      s21_remove_matrix(&A);
      s21_remove_matrix(&B);
      break;
    case 9:
      if (fill_matrix(&rows, &cols, &A) != 0) return 1;
      if (s21_inverse_matrix(&A, &B) != 0) return 1;
      printf("Результат: \n");
      print_matrix(&B);
      s21_remove_matrix(&A);
      s21_remove_matrix(&B);
      break;
    default:
      break;
  }

  return 0;
}