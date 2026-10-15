"""Check docstring diagnostics recorded by the shared library build."""
import argparse

from build_state import prepared_build


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', default='build', help='Prepared library directory')
    args = parser.parse_args()
    with prepared_build(args.build) as snapshot:
        problems = snapshot.manifest['docstrings']
        for problem in problems:
            print(problem)
        print(f'{len(problems)} missing or partial docstrings')
        return 1 if problems else 0


if __name__ == '__main__':
    raise SystemExit(main())
