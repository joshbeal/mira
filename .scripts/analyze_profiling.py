import re
from collections import defaultdict


def extract_time(time_string):
    match = re.search(r"[\.\s]*([\d.]+)(s?)$", time_string)
    if match:
        return float(match.group(1))
    else:
        print(f"Warning: Could not extract time from '{time_string}'")
        return 0.0


def extract_memory_usage(memory_file_path):
    """Extract memory usage from the memory log file."""
    try:
        with open(memory_file_path, "r") as file:
            content = file.read()

        # Extract the "Maximum resident set size (kbytes)" line
        memory_match = re.search(r"Maximum resident set size \(kbytes\):\s+(\d+)", content)
        if memory_match:
            memory_kb = int(memory_match.group(1))
            # Convert memory usage from KB to GB
            memory_gb = memory_kb / (1024**2)  # 1 GB = 1024^2 KB
            return memory_gb
        else:
            print(f"Warning: Could not extract memory usage from '{memory_file_path}'.")
            return 0.0
    except FileNotFoundError:
        print(f"Memory file '{memory_file_path}' not found.")
        return 0.0


def parse_log_file(file_path):
    with open(file_path, "r") as file:
        content = file.read()

    steps = re.findall(
        r"Start: ivc_(?:new|fold_step) step=(\d+)(.*?)End: ivc_(?:new|fold_step) step=\1", content, re.DOTALL
    )

    results = []
    for step, step_content in steps:
        step_number = int(step)

        if step_number == 0:
            continue

        step_data = {"step": step_number}

        for circuit in ["primary", "secondary"]:
            circuit_content = re.search(rf"Start: {circuit}(.*?)End: {circuit}", step_content, re.DOTALL)
            if circuit_content:
                circuit_data = circuit_content.group(1)

                # Witness generation time
                witness_gen = re.search(
                    r"Start: circuit_collect_witness.*?End: circuit_collect_witness.*?([\.\s\d.]+)",
                    circuit_data,
                    re.DOTALL,
                )
                witness_gen_time = extract_time(witness_gen.group(1)) if witness_gen else 0

                # Proving time components (witness_commit + commit_cross_terms + marginal prove work)
                proving_time = 0
                commit_cross_terms_total = 0  # Track total commit_cross_terms contribution
                witness_commit_total = 0  # Track total witness_commit contribution

                # Extract all witness_commit times
                witness_commits = re.findall(
                    r"Start: witness_commit.*?End: witness_commit.*?([\.\s\d.]+)", circuit_data, re.DOTALL
                )
                for commit in witness_commits:
                    witness_commit_total += extract_time(commit)

                # Extract total prove time and account for marginal work
                prove_blocks = re.findall(r"Start: prove(.*?)End: prove.*?([\.\s\d.]+)", circuit_data, re.DOTALL)
                for prove_block, prove_time in prove_blocks:
                    total_prove_time = extract_time(prove_time)

                    # Extract all commit_cross_terms within this prove block
                    commit_cross_terms = re.findall(
                        r"Start: commit_cross_terms.*?End: commit_cross_terms.*?([\.\s\d.]+)", prove_block, re.DOTALL
                    )
                    for cross_terms in commit_cross_terms:
                        commit_cross_terms_total += extract_time(cross_terms)

                    # Marginal work is total prove time minus the known commit_cross_terms and witness_commit
                    marginal_work = total_prove_time - (commit_cross_terms_total + witness_commit_total)

                    # If there's marginal work, add it
                    if marginal_work > 0:
                        proving_time += marginal_work

                # Add witness generation and proving time components to step data
                proving_time += commit_cross_terms_total + witness_commit_total

                step_data[f"{circuit}_witness_gen"] = witness_gen_time
                step_data[f"{circuit}_proving"] = proving_time  # Total proving time
                step_data[f"{circuit}_commit_cross_terms"] = commit_cross_terms_total
                step_data[f"{circuit}_witness_commit"] = witness_commit_total

        results.append(step_data)

    # Now search for ivc_verify at the end of the log
    ivc_verify = re.search(r"Start: ivc_verify.*?End: ivc_verify.*?([\.\s\d.]+)s", content, re.DOTALL)
    ivc_verify_time = extract_time(ivc_verify.group(1)) if ivc_verify else 0

    return results, ivc_verify_time


def calculate_averages(results):
    if not results:
        return {}

    total_steps = len(results)
    averages = defaultdict(float)

    for step_data in results:
        for key, value in step_data.items():
            if key != "step":
                averages[key] += value

    for key in averages:
        averages[key] /= total_steps

    averages["total_witness_gen"] = averages["primary_witness_gen"] + averages["secondary_witness_gen"]
    averages["total_proving"] = averages["primary_proving"] + averages["secondary_proving"]
    averages["total_commit_cross_terms"] = (
        averages["primary_commit_cross_terms"] + averages["secondary_commit_cross_terms"]
    )
    averages["total_witness_commit"] = averages["primary_witness_commit"] + averages["secondary_witness_commit"]

    return dict(averages)


def main(log_file_path, memory_file_path):
    try:
        # Parse the log file for timing information
        results, ivc_verify_time = parse_log_file(log_file_path)
        print(f"Parsed {len(results)} steps from the log file, excluding step 0.")

        if not results:
            print("No data could be extracted from the log file.")
            return

        print("\nRaw results:")
        for result in results:
            print(result)

        averages = calculate_averages(results)

        print("\n===== Average Times per Step =====")
        print("\nPrimary Circuit:")
        print(f"  Proving Time: {averages['primary_proving']:.2f}s")
        print(f"    - Commit Cross Terms Time: {averages['primary_commit_cross_terms']:.2f}s")
        print(f"    - Witness Commit Time: {averages['primary_witness_commit']:.2f}s")

        print("\nSecondary Circuit:")
        print(f"  Proving Time: {averages['secondary_proving']:.2f}s")

        print("\nIVC Verify Time:")
        print(f"  {ivc_verify_time:.2f}s")

        print("\n===== Total Times =====")
        print(f"\nTotal Proving Time: {averages['total_proving']:.2f}s")
        print(f"  - Total Commit Cross Terms Time: {averages['total_commit_cross_terms']:.2f}s")
        print(f"  - Total Witness Commit Time: {averages['total_witness_commit']:.2f}s")

        # Extract and display the memory usage from the memory file
        memory_usage_gb = extract_memory_usage(memory_file_path)
        print(f"\nMaximum Memory Usage: {memory_usage_gb:.2f} GB")

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        import traceback

        print(traceback.format_exc())


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage: python analyze_profiling.py <log_file_path> <memory_file_path>")
        sys.exit(1)

    log_file = sys.argv[1]
    memory_file = sys.argv[2]

    main(log_file, memory_file)
