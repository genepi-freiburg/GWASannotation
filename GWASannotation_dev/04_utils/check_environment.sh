#!/bin/bash

# Function to check if a Python module is installed and the Python version is correct
function check_python_env() {
    # Check Python version
    if python --version 2>&1 | grep -q "Python 3"; then
            echo "Python 3 is installed"
    else
        echo "Python 3 is not installed"
        exit 1
    fi
    
    # Check if the module is installed
    if python -c "import importlib.util; exit(0 if importlib.util.find_spec('$1') else 1)"; then
        echo "$1 is installed"
    else
        echo "$1 is not installed"
        exit 1
    fi
}

if check_python_env "sklearn.linear_model"; then
    echo "Environment check passed"
else
    echo "Environment check failed"
fi
