import os


def ensure_directories_exist(paths):
    for path in paths:
        os.makedirs(path, exist_ok=True)