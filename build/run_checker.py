import os

def check_slurm_outputs(directory, job_id, max_id):
    missing_files = []
    incomplete_files = []

    for i in range(0, max_id + 1):
        filename = f"slurm-{job_id}_{i}.out"
        filepath = os.path.join(directory, filename)

        # Check if file exists
        if not os.path.exists(filepath):
            missing_files.append(i)
            continue

        # Check last line for "Finished successfully"
        try:
            with open(filepath, "r") as f:
                lines = f.read().strip().splitlines()
                if not lines or not any("Finished successfully" in line for line in lines):
                    incomplete_files.append(i)
        except Exception as e:
            print(f"Error reading {filepath}: {e}")
            incomplete_files.append(i)
    # Print results
    print("Missing files:", ",".join(map(str, sorted(missing_files))))
    print("Incomplete files:", ",".join(map(str, sorted(incomplete_files))))


if __name__ == "__main__":
    directory = f"/mnt/home/suzuekar/Avida/build/"
    check_slurm_outputs(directory, 10814483, 319)
