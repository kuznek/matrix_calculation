#include <check.h>
#include <stdlib.h>

#include "test_calc_complements.c"
#include "test_create_remove.c"
#include "test_determinant.c"
#include "test_eq_matrix.c"
#include "test_inverse_matrix.c"
#include "test_mult_matrix.c"
#include "test_mult_number.c"
#include "test_sub_matrix.c"
#include "test_sum_matrix.c"
#include "test_transpose.c"

int main(void) {
  int total_tests = 0;
  int total_failures = 0;
  int suite_count = 0;

  Suite *s;
  SRunner *runner;

  typedef struct {
    const char *name;
    int tests;
    int failures;
  } suite_result_t;

  suite_result_t results[10];

  // Массив функций для создания suite
  struct {
    const char *name;
    Suite *(*create_func)(void);
  } suites[] = {{"s21_create_remove", suite_create_remove},
                {"s21_eq_matrix", suite_eq_matrix},
                {"s21_sum_matrix", suite_sum_matrix},
                {"s21_sub_matrix", suite_sub_matrix},
                {"s21_mult_number", suite_mult_number},
                {"s21_mult_matrix", suite_mult_matrix},
                {"s21_transpose", suite_transpose},
                {"s21_determinant", suite_determinant},
                {"s21_calc_complements", suite_calc_complements},
                {"s21_inverse_matrix", suite_inverse_matrix}};

  // Запуск всех suite
  for (size_t i = 0; i < sizeof(suites) / sizeof(suites[0]); i++) {
    s = suites[i].create_func();
    runner = srunner_create(s);
    srunner_run_all(runner, CK_NORMAL);

    results[suite_count].name = suites[i].name;
    results[suite_count].tests = srunner_ntests_run(runner);
    results[suite_count].failures = srunner_ntests_failed(runner);

    total_tests += results[suite_count].tests;
    total_failures += results[suite_count].failures;

    srunner_free(runner);
    if (results[suite_count].failures > 0) {
      printf("❌ SUITE %s HAS FAILURES!\n", suites[i].name);
    }
    suite_count++;
  }

  printf("\n");
  printf("========================================\n");
  printf("           ОБЩАЯ СТАТИСТИКА             \n");
  printf("========================================\n");
  printf("Всего test suites: %d\n", suite_count);
  printf("Всего тестов:      %d\n", total_tests);
  printf("Успешных:          %d\n", total_tests - total_failures);
  printf("Неуспешных:        %d\n", total_failures);

  double success_rate =
      total_tests > 0
          ? (double)(total_tests - total_failures) / total_tests * 100.0
          : 0.0;
  printf("\n%.1f%% ", success_rate);
  if (total_failures == 0) {
    printf("✅!\n");
  } else {
    printf("❌ И не было печали...\n");
  }

  printf("========================================\n");
  printf("            РАЗБОР ПОЛЁТОВ              \n");
  printf("========================================\n");

  for (int i = 0; i < suite_count; i++) {
    printf("%-25s: %2d тестов", results[i].name, results[i].tests);
    if (results[i].failures == 0) {
      printf(" ✅\n");
    } else {
      printf(" ❌ (%d ошибок)\n", results[i].failures);
    }
  }
  printf("========================================\n");

  return EXIT_SUCCESS;
}
