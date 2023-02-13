#!/usr/bin/env python

# Pull the script part out of a YAML file and put it into a script
# of a similar name, so that it can be used for testing purposes.
# if the YAML is called MyDataFxn.yaml, the script file will be
# called MyDataFxn_script.py.  The test script can then import
# whatever it wants from the Python file and you can be sure that
# you're testing what's in the YAML file.

import argparse
import re
import shutil


from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Extract script from YAML'
                                                 ' file.')
    parser.add_argument('-I', '--input-dir', dest='input_dir',
                        required=True,
                        help='Name of directory containing the YAML files.')
    parser.add_argument('-O', '--output-dir', dest='output_dir',
                        required=True,
                        help='Name of directory for output scripts')
    parser.add_argument('-Y', '--yaml-file', dest='yaml_files',
                        required=True, action='append',
                        help='Name of YAML file to process.  Multiple'
                             ' instances allowed.')
    args = parser.parse_args()
    return args


def extract_script(filename: str) -> list[str]:
    """
    Extract the lines in the YAML, if any between the 'script: |' line
    and the first line which starts at position 0.  The script lines
    are all indented by 2 spaces at least.  Removes those 2 spaces and
    returns a list of lines, still with the \n at the end, ready to be
    dumped unchanged into the new script file.
    """
    script_lines = []
    blank_start = re.compile(r'\s')
    two_space = re.compile(r'  ')
    with open(filename, 'r') as f:
        reading_script = False
        for nextline in f:
            if not blank_start.match(nextline):
                reading_script = False
            if reading_script:
                if two_space.match(nextline):
                    script_lines.append(nextline[2:])
                else:
                    script_lines.append(nextline)

            if nextline.startswith('script: |'):
                reading_script = True

    return script_lines


def write_script_file(script_lines: list[str], filename: str,
                      outdir: str) -> None:
    outfile = Path(outdir) / f'{Path(filename).stem}_script.py'
    print(outfile)
    with open(outfile, 'w') as f:
        f.write(''.join(script_lines))


def main() -> None:
    args = parse_args()
    yaml_dir = Path(args.input_dir)
    for yf in args.yaml_files:
        f = Path(args.input_dir) / yf
        print(yf)
        if f.is_file() and f.suffix == '.yaml':
            print(f'Doing file {f}')
            script_lines = extract_script(f)
            if not script_lines:
                print(f'{f} had no script')
            else:
                write_script_file(script_lines, f, args.output_dir)


if __name__ == '__main__':
    main()
