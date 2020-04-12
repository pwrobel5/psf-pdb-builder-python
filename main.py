import sys

CORRECT_ARGC_VALUE = 2


def main():
    if len(sys.argv) != CORRECT_ARGC_VALUE:
        input_file_name = input("Enter input file name: ")
    else:
        input_file_name = sys.argv[1]


if __name__ == "__main__":
    main()