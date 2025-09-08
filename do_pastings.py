from paste import *

cases = [13]

for case in cases:
    input_file = f"./j4k5/rk4j5_18.txt"
    graphs = read_graphs_from_file(input_file)
    print(f"Read {len(graphs)} graphs from {input_file}.")
    results = get_all_pastes(graphs)

    for size, pastes in results.items():
        output_file = f"pastes/case{case}/d{size}_case{case}_pastes.g6"
        write_graphs_to_file(pastes, output_file)
        print(f"Processed {input_file} for d = {size} and saved pastes to {output_file}.")
