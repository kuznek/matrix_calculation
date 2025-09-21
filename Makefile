CC = gcc
CFLAGS = -Wall -Wextra -Werror -std=c11
LDFLAGS = -lcheck -lm -lpthread -lsubunit
SRC = s21_matrix.c s21_matrix_auxiliary.c
OBJ = $(SRC:.c=.o)
TEST_C = tests/test_matrix_main.c
EXEC = test_s21_matrix
MANUAL = calc_matrix.c
MANUAL_EXEC = calc_matrix
LIB = s21_matrix.a
GCOV_DIR = gcov_report

# имя ядра
UNAME_S := $(shell uname -s)

# macOS без -lsubunit
ifeq ($(UNAME_S), Darwin)
    LDFLAGS = -lcheck -lm -lpthread
endif

all: $(LIB)

$(LIB): $(OBJ)
	ar rcs $@ $^
	ranlib $@

build_test: $(LIB)
	$(CC) $(CFLAGS) $(TEST_C) $< -o $(EXEC) $(LDFLAGS)

test: build_test
	./$(EXEC)

manual:
	$(CC) $(CFLAGS) $(MANUAL) -o $(MANUAL_EXEC)

%.o: %.c clean
	$(CC) $(CFLAGS) -c $< -o $@

# компиляция + линковка с --coverage
gcov_report: CFLAGS += --coverage
gcov_report: LDFLAGS += --coverage
gcov_report: build_test
	mkdir -p $(GCOV_DIR)
	./$(EXEC)  # запуск тестов для генерации .gcda
	lcov --capture --directory . --output-file $(GCOV_DIR)/coverage.info
	lcov --remove $(GCOV_DIR)/coverage.info "*/tests/*" -o $(GCOV_DIR)/coverage_filtered.info
	genhtml $(GCOV_DIR)/coverage_filtered.info --output-directory $(GCOV_DIR)
	echo "Отчёт сгенерирован в папке $(GCOV_DIR)/"

valgrind_report: build_test
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./$(EXEC)

leaks_report: build_test
	leaks -atExit -- ./$(EXEC)

style_check:
	clang-format -n -style=Google $(SRC) $(MANUAL) *.h tests/*.c
style_format:
	clang-format -i -style=Google $(SRC) $(MANUAL) *.h tests/*.c

cppcheck:
	cppcheck --std=c11 --enable=warning --enable=performance --enable=portability --enable=unusedFunction $(SRC) $(MANUAL) *.h tests/*.c

clean:
	@rm -rf *.o *.a $(EXEC) *.gcno *.gcda *.info $(GCOV_DIR) .DS_Store