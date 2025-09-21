#ifndef S21_MATRIX_H
#define S21_MATRIX_H

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#define EPS 1e-12

#define OK 0
#define INCORRECT_MATRIX 1
#define CALCULATION_ERROR 2

#define SUCCESS 1
#define FAILURE 0

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

// вспомогательные
// void print_matrix(matrix_t *result);
void check_matrices(int *status, int n, ...);
void check_matrices_same_size(int *status, matrix_t *A, matrix_t *B);
int create_minor_matrix(matrix_t *A, int exclude_row, int exclude_col,
                        matrix_t *minor);
double calculate_minor(matrix_t *A, int row, int col);

// основные
int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);

int s21_eq_matrix(matrix_t *A, matrix_t *B);

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);

int s21_transpose(matrix_t *A, matrix_t *result);

int s21_determinant(matrix_t *A, double *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

#endif
