import subprocess

# List of Python files to execute in sequence
python_files = ["aP_H.py", "aP_H2.py", "aP_U.py", "aP_UH.py", "aP_V.py", "aP_V2.py", "aP_VH.py"]

# Execute each file in sequence
for file in python_files:
    with open(file, 'r') as f:
        code = compile(f.read(), file, 'exec')
        exec(code)