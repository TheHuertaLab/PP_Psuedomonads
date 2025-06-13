import os

def add_basename_to_files(root_dir):
    for subdir, _, files in os.walk(root_dir):
        # Skip the root directory
        if subdir == root_dir:
            continue

        print(f"Processing directory: {subdir}")
        basename = os.path.basename(subdir)
        print(f"Basename = {basename}")
        first_part = basename.split('.')[0]
        print(f"First part = {first_part}")

        for file_name in files:
            old_file_path = os.path.join(subdir, file_name)
            new_file_name = f"{first_part}_{file_name}"
            new_file_path = os.path.join(subdir, new_file_name)

            try:
                os.rename(old_file_path, new_file_path)
                print(f"Renamed: {old_file_path} to {new_file_path}")
            except OSError as e:
                print(f"Error renaming {old_file_path} to {new_file_path}: {e}")
if __name__ == "__main__":
    root_directory = "/Users/pjoglekar/work/reads/trycycler/subsets/"
    if os.path.exists(root_directory):
        add_basename_to_files(root_directory)
    else:
        print(f"Directory does not exist: {root_directory}")       