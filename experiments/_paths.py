import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ROOT, "models"))
INPUTS = os.path.join(ROOT, "data", "inputs")
OUTPUTS = os.path.join(ROOT, "data", "outputs")
os.makedirs(OUTPUTS, exist_ok=True)
